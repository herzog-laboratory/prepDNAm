# Simulate raw type II betas for a sample with a known rho under the model
# beta_obs = (rho * beta_signal + beta_e) / (1 + rho).
simulate_beta <- function(rho, n = 20000, beta_e = 0.3309135,
                          beta_s1 = 0.019570, beta_s2 = 0.970650, sd = 0.02) {
  signal <- sample(c(beta_s1, beta_s2), n, replace = TRUE)
  obs <- (rho * signal + beta_e) / (1 + rho) + stats::rnorm(n, sd = sd)
  pmin(pmax(obs, 0), 1)
}

test_that("recovers a known rho from simulated betas", {
  set.seed(42)
  true_rho <- c(a = 2, b = 5, c = 10)
  beta <- sapply(true_rho, simulate_beta)
  # add some type I probes that must be ignored
  beta <- rbind(beta, matrix(stats::runif(300), ncol = 3, dimnames = list(NULL, names(true_rho))))
  probe_type <- rep(c("II", "I"), c(20000, 100))

  res <- rhoEstimator(beta, probe_type = probe_type)

  expect_s3_class(res, "data.frame")
  expect_equal(rownames(res), names(true_rho))
  expect_named(res, c("mode0", "mode1", "rho"))
  expect_equal(res$rho, unname(true_rho), tolerance = 0.1)
  expect_true(all(res$mode0 < res$mode1))
})

test_that("samples with too many NAs or a single mode get NA", {
  set.seed(1)
  ok <- simulate_beta(5)
  mostly_na <- ok; mostly_na[1:15000] <- NA
  unimodal <- rep(0.5, 20000)
  beta <- cbind(ok = ok, mostly_na = mostly_na, unimodal = unimodal)

  res <- rhoEstimator(beta, probe_type = rep("II", nrow(beta)))
  expect_false(is.na(res["ok", "rho"]))
  expect_true(all(is.na(unlist(res["mostly_na", ]))))
  expect_true(all(is.na(unlist(res["unimodal", ]))))
})

test_that("works with a single sample", {
  set.seed(2)
  beta <- matrix(simulate_beta(3), ncol = 1, dimnames = list(NULL, "s1"))
  res <- rhoEstimator(beta, probe_type = rep("II", nrow(beta)))
  expect_equal(nrow(res), 1L)
  expect_equal(res$rho, 3, tolerance = 0.1)
})

test_that("input validation", {
  beta <- matrix(stats::runif(20), ncol = 2)
  expect_error(rhoEstimator(beta), "probe_type")
  expect_error(rhoEstimator(beta, probe_type = rep("II", 3)), "one entry per row")
  expect_error(rhoEstimator(beta, probe_type = rep("III", 10)), "\"I\" and \"II\"")
  expect_error(rhoEstimator(as.data.frame(beta), probe_type = rep("II", 10)), "MethylSet")
  expect_error(rhoEstimator(beta, probe_type = rep("I", 10)), "no type II probes")
})

test_that("MethylSet input reproduces eutopsQC 1.2.2 values on minfiData::MsetEx", {
  skip_if_not_installed("minfiData")
  skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
  MsetEx <- NULL
  utils::data("MsetEx", package = "minfiData", envir = environment())

  res <- rhoEstimator(MsetEx)

  # reference values computed with eutopsQC 1.2.2 rhoEstimator(MsetEx, anno = "450k")
  legacy <- data.frame(
    mode0 = c(0.05675146771, 0.06849315068, 0.05675146771, 0.08023483366, 0.08023483366, 0.06457925636),
    mode1 = c(0.8512720157, 0.8590998043, 0.8806262231, 0.8669275930, 0.8571428571, 0.8825831703),
    rho   = c(4.76, 4.85, 6.32, 4.94, 4.53, 6.20),
    row.names = c("5723646052_R02C02", "5723646052_R04C01", "5723646052_R05C02",
                  "5723646053_R04C02", "5723646053_R05C02", "5723646053_R06C02"))
  expect_equal(res, legacy, tolerance = 1e-8)
})

test_that("EPIC v2: probe types come from the manifest, not the shorter annotation", {
  skip_if_not_installed("IlluminaHumanMethylationEPICv2manifest")
  skip_if_not_installed("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  suppressWarnings(suppressMessages(withr::local_package("IlluminaHumanMethylationEPICv2anno.20a1.hg38")))
  types <- manifest_types_for("IlluminaHumanMethylationEPICv2manifest")
  # probes in the manifest but not in the 20a1 annotation (as in real EPIC v2 data)
  ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38, what = "Manifest")
  unannotated <- setdiff(names(types), rownames(ann))
  expect_gt(length(unannotated), 0)
  set.seed(1)
  probes <- c(sample(names(types), 5000), unannotated[1:50])
  n <- length(probes)
  meth <- matrix(rexp(n * 3, 1 / 2000), n, 3, dimnames = list(probes, paste0("s", 1:3)))
  unmeth <- matrix(rexp(n * 3, 1 / 2000), n, 3, dimnames = dimnames(meth))
  mset <- minfi::MethylSet(Meth = meth, Unmeth = unmeth,
                           annotation = c(array = "IlluminaHumanMethylationEPICv2", annotation = "20a1.hg38"))
  expect_lt(length(minfi::getProbeType(mset)), n)          # the minfi behaviour behind the bug
  expect_identical(mset_probe_types(mset), unname(types[probes]))
  expect_no_error(r <- rhoEstimator(mset))
  expect_identical(nrow(r), 3L)
  skip_if_not_installed("ChAMP")
  expect_no_error(b <- suppressWarnings(suppressMessages(normalise_bmiq(mset))))  # random data: degenerate fits warn
  expect_identical(dim(b), dim(meth))
})
