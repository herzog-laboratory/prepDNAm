# Extracted from test-qc_estimates.R:16

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "prepDNAm", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
skip_if_no_fixture_data()
skip_if_not_installed("ChAMP")
skip_if_not_installed("impute")
run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
qc_filter_samples(run_dir)
qc_normalise(run_dir)
qc_merge(run_dir)
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
    expect_equal(pheno$IC, rowSums(pheno[immune]), tolerance = 1e-6)
    expect_equal(pheno$Epi + pheno$Fib + pheno$IC, rep(1, nrow(pheno)), tolerance = 1e-6)
  }
