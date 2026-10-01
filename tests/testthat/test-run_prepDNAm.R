test_that("run_prepDNAm runs end to end, resumes, and summarises", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- file.path(withr::local_tempdir(), "run")
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = run_dir, array = "450k",
                                    batch = "chip", pheno = fixture_pheno_df()))
  res <- run_prepDNAm(cfg)
  expect_s3_class(res, "prepDNAm_run")
  expect_true(res$complete)
  expect_identical(res$n_samples, 6L)
  expect_identical(res$n_batches, 2L)
  expect_identical(res$n_samples_retained, 6L)
  expect_gt(res$n_probes_retained, 400000)
  expect_true(all(res$stages$status == "ok"))
  expect_setequal(res$stages$stage,
                  c("qc_init", "qc_read:5723646052", "qc_read:5723646053",
                    "qc_filter_samples:5723646052", "qc_filter_samples:5723646053",
                    "qc_normalise:5723646052", "qc_normalise:5723646053", "qc_merge", "qc_estimates"))
  expect_output(print(res), "status:   complete")
  expect_output(print(res), "6 in 2 batch\\(es\\), 6 retained")
  expect_true(file.exists(file.path(run_dir, "beta_merged.rds")))
  log <- readLines(file.path(run_dir, "log.txt"))
  expect_match(log, "run_prepDNAm: finished", all = FALSE)

  # resume from the run directory: every stage skips
  n_log <- length(log)
  res2 <- run_prepDNAm(run_dir)
  expect_true(res2$complete)
  new_lines <- readLines(file.path(run_dir, "log.txt"))[-seq_len(n_log)]
  expect_identical(sum(grepl("skipping", new_lines)), 9L)
  expect_false(any(grepl(": start$", new_lines)))
})

test_that("a failed stage stops the run and the summary says so", {
  skip_if_no_fixture_data()
  run_dir <- file.path(withr::local_tempdir(), "run")
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = run_dir, array = "450k",
                                    batch = "none", pheno = fixture_pheno_df("missing")))
  expect_error(run_prepDNAm(cfg), class = "prepDNAm_stage_error")
  res <- run_summary(run_dir)
  expect_false(res$complete)
  expect_identical(res$failed, "FAILED_qc_init.txt")
  expect_output(print(res), "FAILED - FAILED_qc_init.txt")
  expect_null(res$n_probes_retained)
  expect_true(is.na(res$n_samples))
  expect_identical(read_config(run_dir)$batch, "none")
})

test_that("run_prepDNAm validates its input and dependencies", {
  expect_error(run_prepDNAm(42), "must be a `prepDNAm_config`")
  expect_error(run_prepDNAm(withr::local_tempdir()), "not a run directory")
  skip_if_no_fixture_data()
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = file.path(withr::local_tempdir(), "r"),
                                    array = "450k", batch = "none", impute = FALSE, save_sex = FALSE))
  expect_setequal(check_run_dependencies(cfg), c("IlluminaHumanMethylation450kmanifest", "ChAMP"))
  cfg$normalise <- FALSE; cfg$input_type <- "beta"
  expect_identical(check_run_dependencies(cfg), "IlluminaHumanMethylation450kmanifest")
})

test_that("the CLI script parses its arguments", {
  script <- system.file("scripts", "run_prepDNAm.R", package = "prepDNAm")
  if (!nzchar(script)) script <- file.path("..", "..", "inst", "scripts", "run_prepDNAm.R")
  expect_true(file.exists(script))
  out <- suppressWarnings(system2(file.path(R.home("bin"), "Rscript"), script, stdout = TRUE, stderr = TRUE))
  expect_match(out, "Usage:", all = FALSE)
  expect_identical(attr(out, "status"), 2L)
})
