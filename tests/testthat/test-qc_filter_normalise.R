test_that("no samples are removed with default thresholds", {
  run_dir <- copy_run_450k()
  expect_identical(unname(qc_filter_samples(run_dir)), c(TRUE, TRUE))
  for (b in c("5723646052", "5723646053")) {
    rm <- read_samples_removed(run_dir, b)
    expect_identical(nrow(rm), 0L)
    expect_named(rm, c("basename", "batch", "reason", "value", "threshold"))
  }
  expect_identical(unname(qc_filter_samples(run_dir)), c(FALSE, FALSE))
  expect_match(readLines(file.path(run_dir, "log.txt")), "0 of 3 samples removed", all = FALSE)
})

test_that("intensity, failed-probe and exclude reasons are recorded with values", {
  b <- fixture_basenames_450k()
  run_dir <- copy_run_450k(intensity_threshold = 11.8, failed_probe_threshold = 0.0015,
                           exclude = b[1], pheno = fixture_pheno_df())
  qc_filter_samples(run_dir)
  rm <- rbind(read_samples_removed(run_dir, "5723646052"), read_samples_removed(run_dir, "5723646053"))
  expect_true(all(rm$reason %in% c("excluded", "low_intensity", "failed_probes")))

  ex <- rm[rm$reason == "excluded", ]
  expect_identical(ex$basename, b[1])
  expect_true(is.na(ex$value) && is.na(ex$threshold))

  qc <- rbind(readRDS(file.path(batch_dir(run_dir, "5723646052"), "qc.rds")),
              readRDS(file.path(batch_dir(run_dir, "5723646053"), "qc.rds")))
  expected_low <- qc$basename[pmin(qc$mMed, qc$uMed) < 11.8]
  low <- rm[rm$reason == "low_intensity", ]
  expect_setequal(low$basename, expected_low)
  expect_true(length(expected_low) > 0 && length(expected_low) < 6)
  expect_equal(low$value, pmin(qc$mMed, qc$uMed)[match(low$basename, qc$basename)])
  expect_true(all(low$threshold == 11.8))

  detP <- cbind(readRDS(file.path(batch_dir(run_dir, "5723646052"), "detP.rds")),
                readRDS(file.path(batch_dir(run_dir, "5723646053"), "detP.rds")))
  frac <- colMeans(detP > 0.01)
  fp <- rm[rm$reason == "failed_probes", ]
  expect_setequal(fp$basename, names(frac)[frac > 0.0015])
  expect_true(nrow(fp) > 0 && nrow(fp) < 6)
  expect_equal(fp$value, unname(frac[fp$basename]))
  expect_match(readLines(file.path(run_dir, "log.txt")), "excluded: 1", all = FALSE)
})

test_that("beta input only applies the exclude list", {
  skip_if_no_fixture_data()
  run_dir <- file.path(withr::local_tempdir(), "run")
  b <- fixture_basenames_450k()
  qc_init(suppressMessages(qc_config(input = fixture_beta()$rds, run_dir = run_dir, array = "450k",
                                     batch = "none", exclude = b[2:3], intensity_threshold = 100)))
  qc_read(run_dir)
  qc_filter_samples(run_dir)
  rm <- read_samples_removed(run_dir, "all")
  expect_identical(rm$basename, sort(b[2:3]))
  expect_identical(unique(rm$reason), "excluded")
  expect_match(readLines(file.path(run_dir, "log.txt")), "intensity filter skipped", all = FALSE)
})

test_that("qc_normalise writes BMIQ betas for the retained samples and adds bmiq densities", {
  skip_if_not_installed("ChAMP")
  b <- fixture_basenames_450k()
  run_dir <- copy_run_450k(exclude = b[1], pheno = fixture_pheno_df())
  qc_filter_samples(run_dir)
  expect_identical(unname(qc_normalise(run_dir)), c(TRUE, TRUE))
  expect_identical(unname(qc_normalise(run_dir)), c(FALSE, FALSE))

  bdir <- batch_dir(run_dir, "5723646052")
  beta <- readRDS(file.path(bdir, "beta_norm.rds"))
  expect_identical(dim(beta), c(485512L, 2L))
  expect_identical(colnames(beta), b[2:3])
  expect_true(all(beta >= 0 & beta <= 1))
  # same as calling normalise_bmiq() on the Noob mset directly
  mset <- readRDS(file.path(bdir, "mset.rds"))
  ref <- normalise_bmiq(mset[, b[2:3]])
  expect_equal(beta, ref, tolerance = 1e-12)
  d <- readRDS(file.path(bdir, "densities.rds"))
  expect_named(d, c("x", "raw", "noob", "bmiq"))
  expect_identical(dim(d$bmiq), c(512L, 2L))
  expect_identical(colnames(d$bmiq), b[2:3])
  # the other batch keeps all three samples
  expect_identical(ncol(readRDS(file.path(batch_dir(run_dir, "5723646053"), "beta_norm.rds"))), 3L)
  expect_error(qc_normalise(copy_run_450k()), "Run qc_filter_samples\\(\\) first")
})

test_that("qc_normalise on beta input: forced on/off and the auto heuristic", {
  skip_if_no_fixture_data()
  skip_if_not_installed("ChAMP")
  b <- fixture_basenames_450k()
  beta_noob <- fixture_beta()$beta[, b[1:2]]

  run_for <- function(beta, normalise, env = parent.frame()) {
    f <- withr::local_tempfile(fileext = ".rds", .local_envir = env); saveRDS(beta, f)
    run_dir <- file.path(withr::local_tempdir(.local_envir = env), "run")
    qc_init(suppressMessages(qc_config(input = f, run_dir = run_dir, array = "450k", batch = "none",
                                       normalise = normalise)))
    qc_read(run_dir); qc_filter_samples(run_dir); qc_normalise(run_dir)
    run_dir
  }
  # off: output is the input
  r_off <- run_for(beta_noob, FALSE)
  expect_equal(readRDS(file.path(batch_dir(r_off, "all"), "beta_norm.rds")), beta_noob)
  expect_null(readRDS(file.path(batch_dir(r_off, "all"), "densities.rds"))$bmiq)
  expect_match(readLines(file.path(r_off, "log.txt")), "BMIQ off \\(`normalise = FALSE`\\)", all = FALSE)

  # auto on Noob betas: unnormalised -> BMIQ runs
  r_auto <- run_for(beta_noob, "auto")
  log <- readLines(file.path(r_auto, "log.txt"))
  expect_match(log, "BMIQ on \\(input looks unnormalised\\)", all = FALSE)
  beta_auto <- readRDS(file.path(batch_dir(r_auto, "all"), "beta_norm.rds"))
  expect_false(isTRUE(all.equal(beta_auto, beta_noob)))
  types <- manifest_probe_types(array_info("450k"))[rownames(beta_noob)]
  expect_equal(beta_auto, normalise_bmiq(beta_noob, probe_type = unname(types)), tolerance = 1e-12)

  # auto on already-BMIQ'd betas: normalised -> BMIQ skipped
  r_auto2 <- run_for(beta_auto, "auto")
  expect_match(readLines(file.path(r_auto2, "log.txt")), "BMIQ off \\(input looks normalised\\)", all = FALSE)
  expect_equal(readRDS(file.path(batch_dir(r_auto2, "all"), "beta_norm.rds")), beta_auto)
})

test_that("type_ii_peak_shift measures the shift", {
  set.seed(3)
  n1 <- 5000; n2 <- 15000
  b_I <- c(stats::rbeta(n1 / 2, 2, 40), stats::rbeta(n1 / 2, 40, 2))
  b_II <- c(stats::rbeta(n2 / 2, 4, 20), stats::rbeta(n2 / 2, 20, 4))
  beta <- cbind(s1 = c(b_I, b_II), s2 = c(b_I, b_I[1:n2 %% n1 + 1]))   # s2: type II drawn like type I
  pt <- rep(c("I", "II"), c(n1, n2))
  expect_gt(type_ii_peak_shift(beta[, "s1", drop = FALSE], pt), 0.08)
  expect_lt(type_ii_peak_shift(beta[, "s2", drop = FALSE], pt), 0.02)
})
