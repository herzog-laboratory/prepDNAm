test_that("qc_estimates adds estimates to pheno_postQC and estimates.csv", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  before <- readRDS(file.path(run_dir, "pheno_postQC.rds"))
  expect_true(qc_estimates(run_dir))
  est <- utils::read.csv(file.path(run_dir, "estimates.csv"), check.names = FALSE)
  pheno <- readRDS(file.path(run_dir, "pheno_postQC.rds"))
  expect_identical(est$basename, before$basename)
  expect_identical(names(pheno), c(names(before), setdiff(names(est), "basename")))
  expect_identical(pheno[names(before)], before)
  expect_true(all(c("cg05575921", "cg16867657") %in% names(est)))
  if (requireNamespace("EpiDISH", quietly = TRUE)) {
    expect_true(all(c("Epi", "Fib", "IC", "CD4Tnv", "Mono") %in% names(pheno)))
    immune <- c("CD4Tnv", "Baso", "CD4Tmem", "Bmem", "Bnv", "Treg", "CD8Tmem", "CD8Tnv", "Eos", "NK", "Neu", "Mono")
    expect_equal(pheno$IC, unname(rowSums(pheno[immune])), tolerance = 1e-6)
    expect_equal(pheno$Epi + pheno$Fib + pheno$IC, rep(1, nrow(pheno)), tolerance = 1e-6)
  }
  if (requireNamespace("WIDclocks", quietly = TRUE)) expect_true(is.numeric(pheno$wid_general_age))
  # idempotent: skip, and a forced re-run does not duplicate columns
  expect_false(qc_estimates(run_dir))
  qc_estimates(run_dir, force = TRUE)
  expect_identical(names(readRDS(file.path(run_dir, "pheno_postQC.rds"))), names(pheno))
  # re-merging drops the estimates and lets the stage run again
  qc_merge(run_dir, force = TRUE)
  expect_false(file.exists(file.path(run_dir, "estimates.csv")))
  expect_true(qc_estimates(run_dir))
  # export metadata carries the estimates
  meta <- utils::read.csv(export_beta_csv(run_dir, out_dir = withr::local_tempdir())[2], check.names = FALSE)
  expect_true(all(setdiff(names(est), "basename") %in% names(meta)))
})

test_that("qc_estimates refuses to overwrite sample sheet columns", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  pheno <- fixture_pheno_df(); pheno$cg05575921 <- 1
  run_dir <- copy_run_450k(pheno = pheno, failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  expect_error(qc_estimates(run_dir), "cg05575921")
  expect_error(qc_estimates(copy_run_450k()), "qc_merge")
})
