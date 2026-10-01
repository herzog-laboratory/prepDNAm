cfg_450k <- function(run = "run", ...) {
  suppressMessages(qc_config(
    input = fixture_idat_450k(), run_dir = file.path(withr::local_tempdir(.local_envir = parent.frame()), run),
    array = "450k", ...))
}

test_that("qc_init writes the run directory layout", {
  skip_if_no_fixture_data()
  cfg <- cfg_450k(batch = "dir", pheno = fixture_pheno("csv"))
  expect_true(qc_init(cfg))
  expect_true(dir.exists(cfg$run_dir))
  expect_true(all(file.exists(file.path(cfg$run_dir, c("config.yml", "samples.csv", "run_info.json", "log.txt")))))
  expect_false(file.exists(file.path(cfg$run_dir, "FAILED_qc_init.txt")))

  # config round-trips
  on_disk <- read_config(cfg$run_dir)
  expect_s3_class(on_disk, "prepDNAm_config")
  expect_length(config_diff(cfg, on_disk), 0)
  expect_identical(on_disk$cores, 1L)
  expect_identical(on_disk$seed, 1L)
  expect_null(on_disk$exclude)
  expect_identical(on_disk$detp_threshold, 0.01)

  # sample sheet
  s <- read_samples(cfg$run_dir)
  expect_identical(s$basename, fixture_basenames_450k())
  expect_type(s$exclude, "logical")
  expect_type(s$batch, "character")
  expect_identical(s$batch, as.character(s$plate))   # plate == chip == dir here
  expect_true(all(c("sex", "age", "group") %in% names(s)))

  # run info
  info <- read_run_info(cfg$run_dir)
  expect_identical(info$prepDNAm_version, cfg$prepDNAm_version)
  expect_identical(info$n_samples, 6L)
  expect_identical(info$batches, list(`5723646052` = 3L, `5723646053` = 3L))
  expect_identical(info$stages$qc_init$status, "ok")
  expect_true(is.character(info$stages$qc_init$started))
  expect_true(is.character(info$stages$qc_init$finished))
  expect_identical(info$md5$samples, unname(tools::md5sum(file.path(cfg$run_dir, "samples.csv"))))
  expect_identical(info$md5$pheno, unname(tools::md5sum(fixture_pheno("csv"))))
  expect_true("prepDNAm" %in% names(info$packages))
  expect_true("minfi" %in% names(info$packages))

  # log
  log <- readLines(file.path(cfg$run_dir, "log.txt"))
  expect_match(log, "qc_init: start", all = FALSE)
  expect_match(log, "6 samples in 2 batch", all = FALSE)
  expect_match(log, "qc_init: done", all = FALSE)
})

test_that("a pheno data frame is materialised as pheno.csv", {
  skip_if_no_fixture_data()
  cfg <- cfg_450k(batch = "plate", pheno = fixture_pheno_df())
  qc_init(cfg)
  expect_true(file.exists(file.path(cfg$run_dir, "pheno.csv")))
  expect_identical(read_config(cfg$run_dir)$pheno, file.path(cfg$run_dir, "pheno.csv"))
  expect_identical(read_samples(cfg$run_dir)$age, fixture_pheno_df()$age)
  # re-running with the same data frame is a no-op, a changed one an error
  expect_false(qc_init(cfg))
  cfg2 <- cfg; cfg2$pheno$age[1] <- 99
  expect_error(qc_init(cfg2), "Differing fields: pheno")
})

test_that("re-running is idempotent and a different config is refused", {
  skip_if_no_fixture_data()
  cfg <- cfg_450k(batch = "chip")
  expect_true(qc_init(cfg))
  mtime <- file.mtime(file.path(cfg$run_dir, "samples.csv"))
  expect_false(qc_init(cfg))
  expect_identical(file.mtime(file.path(cfg$run_dir, "samples.csv")), mtime)
  expect_match(readLines(file.path(cfg$run_dir, "log.txt")), "skipping", all = FALSE)

  cfg2 <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = cfg$run_dir,
                                     array = "450k", batch = "none", detp_threshold = 0.05))
  expect_error(qc_init(cfg2), "different configuration")
  expect_error(qc_init(cfg2), "Differing fields: batch, detp_threshold")
  # force overwrites
  expect_true(qc_init(cfg2, force = TRUE))
  expect_identical(read_config(cfg$run_dir)$batch, "none")
  expect_identical(unique(read_samples(cfg$run_dir)$batch), "all")
})

test_that("a stage failure leaves a FAILED marker and is recorded", {
  skip_if_no_fixture_data()
  # pheno lists a sample that is not on disk -> discover_samples() errors inside the stage
  cfg <- cfg_450k(batch = "none", pheno = fixture_pheno_df("missing"))
  expect_error(qc_init(cfg), class = "prepDNAm_stage_error")
  expect_error(qc_init(cfg), "Stage qc_init failed: .*not found in `input`")
  failed <- file.path(cfg$run_dir, "FAILED_qc_init.txt")
  expect_true(file.exists(failed))
  expect_match(readLines(failed), "5723646053_R01C01", all = FALSE)
  expect_match(readLines(failed), "Traceback", all = FALSE)
  expect_identical(read_run_info(cfg$run_dir)$stages$qc_init$status, "failed")
  expect_false(file.exists(file.path(cfg$run_dir, "samples.csv")))
  expect_match(readLines(file.path(cfg$run_dir, "log.txt")), "FAILED", all = FALSE)

  # fixing the problem and re-running resumes and clears the marker
  cfg_ok <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = cfg$run_dir,
                                       array = "450k", batch = "none", pheno = fixture_pheno_df()))
  expect_true(qc_init(cfg_ok))
  expect_false(file.exists(failed))
  expect_identical(read_run_info(cfg$run_dir)$stages$qc_init$status, "ok")
})

test_that("read_config / read_samples need an initialised run dir", {
  d <- withr::local_tempdir()
  expect_error(read_config(d), "not a run directory")
  expect_error(read_samples(d), "no samples.csv")
})

test_that("config_diff tolerates yaml round-trip type changes", {
  a <- structure(list(x = 1L, y = NULL, z = "a", created = "t1"), class = "prepDNAm_config")
  b <- structure(list(x = 1, z = "a", created = "t2"), class = "prepDNAm_config")
  expect_length(config_diff(a, b), 0)
  b$z <- "b"
  expect_identical(config_diff(a, b), "z")
})

test_that("run_stage skips, runs, seeds and records", {
  d <- withr::local_tempdir()
  out <- file.path(d, "out.txt")
  ran <- 0
  body <- function() { ran <<- ran + 1; writeLines(as.character(stats::runif(1)), out) }
  expect_true(run_stage(d, "demo", out, force = FALSE, seed = 7, body))
  first <- readLines(out)
  expect_false(run_stage(d, "demo", out, force = FALSE, seed = 7, body))
  expect_identical(ran, 1)
  expect_true(run_stage(d, "demo", out, force = TRUE, seed = 7, body))
  expect_identical(readLines(out), first)   # same seed, same draw
  info <- read_run_info(d)
  expect_identical(info$stages$demo$status, "ok")
  expect_true(is.numeric(info$stages$demo$seconds))
})

test_that("qc_init validates its arguments", {
  expect_error(qc_init(list()), "must be a `prepDNAm_config`")
  skip_if_no_fixture_data()
  expect_error(qc_init(cfg_450k(batch = "none"), force = "yes"), "must be TRUE or FALSE")
})
