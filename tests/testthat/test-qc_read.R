# Initialised run for the 450k IDAT fixture (2 chips -> 2 batches).
init_450k <- function(..., layout = "nested", batch = "dir") {
  run_dir <- file.path(withr::local_tempdir(.local_envir = parent.frame()), "run")
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(layout), run_dir = run_dir,
                                    array = "450k", batch = batch, ...))
  qc_init(cfg)
  run_dir
}

test_that("IDAT batches produce the full bundle", {
  skip_if_no_fixture_data()
  run_dir <- init_450k(save_rs = TRUE)
  ran <- qc_read(run_dir)
  expect_identical(ran, c(`5723646052` = TRUE, `5723646053` = TRUE))

  bdir <- batch_dir(run_dir, "5723646052")
  files <- c("mset.rds", "detP.rds", "qc.rds", "rho.rds", "sex_raw.rds", "rs.rds", "ctrls.rds", "densities.rds")
  expect_true(all(file.exists(file.path(bdir, files))))
  expect_false(file.exists(file.path(bdir, "shiny.rds")))
  b <- read_samples(run_dir); b <- b$basename[b$batch == "5723646052"]

  mset <- readRDS(file.path(bdir, "mset.rds"))
  expect_s4_class(mset, "MethylSet")
  expect_identical(dim(mset), c(485512L, 3L))
  expect_identical(colnames(mset), b)
  expect_match(minfi::preprocessMethod(mset)[1], "Noob")
  expect_identical(unname(minfi::annotation(mset)), c("IlluminaHumanMethylation450k", "ilmn12.hg19"))

  detP <- readRDS(file.path(bdir, "detP.rds"))
  expect_identical(dim(detP), c(485512L, 3L))
  expect_identical(colnames(detP), b)
  expect_true(all(detP >= 0 & detP <= 1))

  qc <- readRDS(file.path(bdir, "qc.rds"))
  expect_named(qc, c("basename", "mMed", "uMed", "snp_outlier"))
  expect_identical(qc$basename, b)
  expect_true(all(qc$mMed > 9 & qc$uMed > 9))
  skip_if_not_installed("ewastools")
  expect_false(anyNA(qc$snp_outlier))

  rho <- readRDS(file.path(bdir, "rho.rds"))
  expect_identical(rownames(rho), b)
  expect_named(rho, c("mode0", "mode1", "rho"))

  sex <- readRDS(file.path(bdir, "sex_raw.rds"))
  expect_named(sex, c("basename", "xMed", "yMed"))
  expect_identical(sex$basename, b)

  rs <- readRDS(file.path(bdir, "rs.rds"))
  expect_identical(nrow(rs), 65L)
  expect_true(all(startsWith(rownames(rs), "rs")))

  ctrls <- readRDS(file.path(bdir, "ctrls.rds"))
  expect_named(ctrls, c("info", "red", "green"))
  expect_identical(nrow(ctrls$info), nrow(ctrls$red))
  expect_identical(colnames(ctrls$green), b)

  d <- readRDS(file.path(bdir, "densities.rds"))
  expect_length(d$x, 512)
  expect_identical(dim(d$raw), c(512L, 3L))
  expect_identical(dim(d$noob), c(512L, 3L))
  expect_identical(colnames(d$noob), b)

  info <- read_run_info(run_dir)
  expect_identical(info$stages[["qc_read:5723646052"]]$status, "ok")
  expect_identical(info$stages[["qc_read:5723646053"]]$status, "ok")
})

test_that("qc_read matches minfi run directly, is idempotent and honours force/batch", {
  skip_if_no_fixture_data()
  run_dir <- init_450k(save_sex = FALSE)
  expect_identical(unname(qc_read(run_dir, batch = "5723646053")), TRUE)
  expect_false(file.exists(file.path(batch_dir(run_dir, "5723646052"), "mset.rds")))
  expect_identical(unname(qc_read(run_dir, batch = "5723646053")), FALSE)
  expect_identical(unname(qc_read(run_dir)), c(TRUE, FALSE))
  expect_identical(unname(qc_read(run_dir, batch = "5723646053", force = TRUE)), TRUE)
  expect_error(qc_read(run_dir, batch = "nope"), "Unknown batch")
  expect_false(file.exists(file.path(batch_dir(run_dir, "5723646053"), "sex_raw.rds")))

  # reference: the same computation with minfi
  b <- read_samples(run_dir); b <- b[b$batch == "5723646053", ]
  rg <- minfi::read.metharray(b$path, force = TRUE)
  ref_qc <- minfi::getQC(minfi::preprocessRaw(rg))
  qc <- readRDS(file.path(batch_dir(run_dir, "5723646053"), "qc.rds"))
  expect_equal(qc$mMed, unname(ref_qc$mMed))
  expect_equal(qc$uMed, unname(ref_qc$uMed))
  ref_beta <- minfi::getBeta(minfi::preprocessNoob(rg, dyeCorr = TRUE, dyeMethod = "single"))
  beta <- minfi::getBeta(readRDS(file.path(batch_dir(run_dir, "5723646053"), "mset.rds")))
  expect_equal(unname(beta), unname(ref_beta))
  ref_detP <- minfi::detectionP(rg, type = "m+u")
  expect_equal(unname(readRDS(file.path(batch_dir(run_dir, "5723646053"), "detP.rds"))), unname(ref_detP))
  ref_sex <- minfi::getSex(minfi::mapToGenome(minfi::preprocessRaw(minfi::read.metharray(read_samples(run_dir)$path, force = TRUE))))
  run_dir2 <- init_450k(batch = "none")
  qc_read(run_dir2)
  sex <- readRDS(file.path(batch_dir(run_dir2, "all"), "sex_raw.rds"))
  expect_equal(sex$xMed, unname(ref_sex$xMed))
  expect_equal(sex$yMed, unname(ref_sex$yMed))
})

test_that("gzipped and flat layouts read identically", {
  skip_if_no_fixture_data()
  run_gz <- init_450k(layout = "gz", batch = "chip", save_sex = FALSE)
  run_flat <- init_450k(layout = "flat", batch = "chip", save_sex = FALSE)
  qc_read(run_gz, batch = "5723646052")
  qc_read(run_flat, batch = "5723646052")
  expect_identical(readRDS(file.path(batch_dir(run_gz, "5723646052"), "qc.rds")),
                   readRDS(file.path(batch_dir(run_flat, "5723646052"), "qc.rds")))
})

test_that("a wrong `array` for the IDATs is caught", {
  skip_if_no_fixture_data()
  skip_if_not_installed("IlluminaHumanMethylationEPICmanifest")
  run_dir <- file.path(withr::local_tempdir(), "run")
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = run_dir, array = "EPIC",
                                    batch = "none", save_sex = FALSE))
  qc_init(cfg)
  expect_error(qc_read(run_dir), "look like IlluminaHumanMethylation450k but `array` is EPIC")
  expect_true(file.exists(file.path(run_dir, "FAILED_qc_read:all.txt")))
})

test_that("signal-matrix input reproduces the IDAT raw metrics", {
  skip_if_no_fixture_data()
  run_dir <- file.path(withr::local_tempdir(), "run")
  cfg <- suppressMessages(qc_config(input = fixture_sigmatrix(), run_dir = run_dir, array = "450k",
                                    batch = "chip", save_rs = TRUE))
  qc_init(cfg)
  expect_identical(unname(qc_read(run_dir)), c(TRUE, TRUE))
  bdir <- batch_dir(run_dir, "5723646052")
  expect_true(all(file.exists(file.path(bdir, c("mset.rds", "detP.rds", "qc.rds", "rho.rds", "sex_raw.rds", "densities.rds")))))
  expect_false(file.exists(file.path(bdir, "ctrls.rds")))
  expect_false(file.exists(file.path(bdir, "rs.rds")))   # fixture has no rs rows
  mset <- readRDS(file.path(bdir, "mset.rds"))
  expect_identical(dim(mset), c(485512L, 3L))
  expect_identical(unname(minfi::preprocessMethod(mset)[1]), "upstream (signal matrix)")
  d <- readRDS(file.path(bdir, "densities.rds"))
  expect_null(d$noob)

  # same raw M/U as the IDAT path -> identical QC, rho and detP
  run_idat <- init_450k(batch = "chip")
  qc_read(run_idat, batch = "5723646052")
  idir <- batch_dir(run_idat, "5723646052")
  qc_sig <- readRDS(file.path(bdir, "qc.rds")); qc_idat <- readRDS(file.path(idir, "qc.rds"))
  expect_equal(qc_sig[, c("mMed", "uMed")], qc_idat[, c("mMed", "uMed")], tolerance = 1e-10)
  expect_true(all(is.na(qc_sig$snp_outlier)))
  expect_equal(readRDS(file.path(bdir, "rho.rds")), readRDS(file.path(idir, "rho.rds")), tolerance = 1e-10)
  expect_equal(readRDS(file.path(bdir, "detP.rds")), readRDS(file.path(idir, "detP.rds")), tolerance = 1e-10)
  expect_equal(readRDS(file.path(bdir, "sex_raw.rds")), readRDS(file.path(idir, "sex_raw.rds")), tolerance = 1e-10)
})

test_that("signal matrix without p-values logs a warning and skips detP", {
  skip_if_no_fixture_data()
  hdr <- strsplit(readLines(fixture_sigmatrix(), n = 1), "\t")[[1]]
  keep <- !grepl("Detection Pval", hdr)
  f <- withr::local_tempfile(fileext = ".tsv")
  dt <- data.table::fread(fixture_sigmatrix(), select = which(keep))
  dt <- dt[seq(1, nrow(dt), by = 25), ]   # every 25th row: both probe types, 19421 rows
  data.table::fwrite(dt, f, sep = "\t")
  run_dir <- file.path(withr::local_tempdir(), "run")
  cfg <- suppressMessages(qc_config(input = f, run_dir = run_dir, array = "450k", batch = "none",
                                    save_sex = FALSE))
  qc_init(cfg)
  qc_read(run_dir)
  expect_false(file.exists(file.path(batch_dir(run_dir, "all"), "detP.rds")))
  expect_identical(nrow(readRDS(file.path(batch_dir(run_dir, "all"), "mset.rds"))), 19421L)
  log <- readLines(file.path(run_dir, "log.txt"))
  expect_match(log, "no detection p-value columns", all = FALSE)
  expect_match(log, "only 19421 of 485512", all = FALSE)
})

test_that("beta input is stored per batch, restricted to manifest probes", {
  skip_if_no_fixture_data()
  for (path in c(fixture_beta()$rds, fixture_beta()$csv)) {
    run_dir <- file.path(withr::local_tempdir(), "run")
    cfg <- suppressMessages(qc_config(input = path, run_dir = run_dir, array = "450k",
                                      pheno = fixture_pheno_df(), batch = "plate"))
    qc_init(cfg)
    expect_identical(unname(qc_read(run_dir)), c(TRUE, TRUE))
    bdir <- batch_dir(run_dir, "5723646052")
    expect_identical(list.files(bdir), c("beta_input.rds", "densities.rds"))
    beta <- readRDS(file.path(bdir, "beta_input.rds"))
    expect_identical(dim(beta), c(485512L, 3L))
    expect_true(anyNA(beta))
    expect_equal(beta, fixture_beta()$beta[, colnames(beta)], tolerance = 1e-12)
  }
})

test_that("beta input outside [0, 1] and unknown probes are handled", {
  skip_if_no_fixture_data()
  beta <- fixture_beta()$beta[1:1000, ]
  rownames(beta)[1:10] <- paste0("bogus", 1:10)
  f <- withr::local_tempfile(fileext = ".rds"); saveRDS(beta, f)
  run_dir <- file.path(withr::local_tempdir(), "run")
  qc_init(suppressMessages(qc_config(input = f, run_dir = run_dir, array = "450k", batch = "none")))
  qc_read(run_dir)
  expect_identical(nrow(readRDS(file.path(batch_dir(run_dir, "all"), "beta_input.rds"))), 990L)
  expect_match(readLines(file.path(run_dir, "log.txt")), "10 rows .* dropped", all = FALSE)

  mvals <- log2(beta / (1 - beta))
  f2 <- withr::local_tempfile(fileext = ".rds"); saveRDS(mvals, f2)
  run_dir2 <- file.path(withr::local_tempdir(), "run")
  qc_init(suppressMessages(qc_config(input = f2, run_dir = run_dir2, array = "450k", batch = "none")))
  expect_error(qc_read(run_dir2), "betas must lie in \\[0, 1\\]")
})

test_that("beta_densities handles NA-only samples and the grid", {
  set.seed(1)
  beta <- cbind(a = stats::runif(500), b = NA_real_)
  d <- beta_densities(beta)
  expect_identical(dim(d), c(512L, 2L))
  expect_false(anyNA(d[, "a"]))
  expect_true(all(is.na(d[, "b"])))
  expect_equal(range(density_grid()), c(-0.05, 1.05))
})
