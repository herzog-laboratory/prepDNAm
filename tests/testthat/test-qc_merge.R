test_that("qc_merge writes the final outputs with consistent dimensions", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- fixture_merged_450k()
  b <- fixture_basenames_450k()
  expect_true(all(file.exists(file.path(run_dir, c(
    "beta_merged.rds", "probes_removed.csv", "imputed.csv", "pheno_postQC.rds", "sample_qc.csv", "rs_merged.rds")))))
  expect_false(file.exists(file.path(run_dir, "beta_merged_compatible.rds")))
  expect_false(file.exists(file.path(run_dir, "imputed_mask.csv")))

  beta <- readRDS(file.path(run_dir, "beta_merged.rds"))
  expect_identical(colnames(beta), b)
  expect_false(anyNA(beta))
  expect_true(all(beta >= 0 & beta <= 1))
  removed <- utils::read.csv(file.path(run_dir, "probes_removed.csv"), stringsAsFactors = FALSE)
  expect_identical(nrow(beta) + length(unique(removed$probe)), 485512L)
  expect_false(any(removed$probe %in% rownames(beta)))
  expect_true(all(removed$reason %in% c("chrY", "snp", "non_cpg", "zhou", "detP", "na_fraction")))
  # list-based removals are exactly the intersection with the array's lists
  for (r in c("chrY", "snp", "non_cpg", "zhou")) {
    lst <- probe_list(array_info("450k")$probe_lists[[r]])
    expect_setequal(removed$probe[removed$reason == r], intersect(rownames(fixture_beta()$beta), lst))
  }
  expect_gt(sum(removed$reason == "detP"), 0)

  imp <- utils::read.csv(file.path(run_dir, "imputed.csv"), stringsAsFactors = FALSE)
  expect_identical(imp$basename, b)
  detP <- cbind(readRDS(file.path(batch_dir(run_dir, "5723646052"), "detP.rds")),
                readRDS(file.path(batch_dir(run_dir, "5723646053"), "detP.rds")))
  expect_equal(imp$n_masked, unname(colSums(detP[, b] > 0.01)))
  expect_true(all(imp$n_missing > 0 & imp$n_missing <= imp$n_masked))
  expect_true(all(imp$imputed))
  # a probe is removed for detP only if it fails in > 50% of samples
  fails <- rowSums(detP[, b] > 0.01)
  expect_setequal(removed$probe[removed$reason == "detP"], names(fails)[fails / 6 > 0.5])

  pheno <- readRDS(file.path(run_dir, "pheno_postQC.rds"))
  expect_identical(pheno$basename, b)
  expect_true(all(c("batch", "mMed", "uMed", "rho", "snp_outlier", "frac_failed_probes",
                    "n_missing", "xMed", "yMed", "predicted_sex") %in% names(pheno)))
  expect_false(any(c("removed", "reasons") %in% names(pheno)))
  qc <- utils::read.csv(file.path(run_dir, "sample_qc.csv"), stringsAsFactors = FALSE)
  expect_identical(nrow(qc), 6L)
  expect_false(any(qc$removed))
  expect_equal(qc$mMed, pheno$mMed)

  # sex as minfi calls it on all samples at once
  rg <- minfi::read.metharray(read_samples(run_dir)$path, force = TRUE)
  ref <- suppressWarnings(minfi::getSex(minfi::mapToGenome(minfi::preprocessRaw(rg))))
  expect_identical(pheno$predicted_sex, unname(ref$predictedSex))
  expect_equal(pheno$xMed, unname(ref$xMed))
  # ... and matches the sample sheet's sex for these samples
  expect_identical(pheno$predicted_sex, fixture_pheno_df()$sex)

  rs <- readRDS(file.path(run_dir, "rs_merged.rds"))
  expect_identical(dim(rs), c(65L, 6L))
  info <- read_run_info(run_dir)
  expect_identical(info$n_samples_retained, 6L)
  expect_identical(info$n_probes_retained, nrow(beta))
  expect_identical(info$stages$qc_merge$status, "ok")
  expect_false(qc_merge(run_dir))
})

test_that("the result does not depend on batching", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  by_dir <- fixture_merged_450k()
  run_dir <- file.path(withr::local_tempdir(), "run")
  qc_init(suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = run_dir, array = "450k",
                                     batch = "none", save_rs = TRUE, failed_sample_threshold = 0.5)))
  qc_read(run_dir); qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  expect_true(all(utils::read.csv(file.path(run_dir, "imputed.csv"))$imputed))
  expect_identical(unique(read_samples(run_dir)$batch), "all")

  expect_equal(readRDS(file.path(run_dir, "beta_merged.rds")), readRDS(file.path(by_dir, "beta_merged.rds")))
  rm1 <- utils::read.csv(file.path(run_dir, "probes_removed.csv")); rm2 <- utils::read.csv(file.path(by_dir, "probes_removed.csv"))
  expect_identical(rm1[order(rm1$probe, rm1$reason), ], rm2[order(rm2$probe, rm2$reason), ], ignore_attr = TRUE)
  expect_identical(utils::read.csv(file.path(run_dir, "imputed.csv")), utils::read.csv(file.path(by_dir, "imputed.csv")))
  p1 <- readRDS(file.path(run_dir, "pheno_postQC.rds")); p2 <- readRDS(file.path(by_dir, "pheno_postQC.rds"))
  expect_equal(p1[, setdiff(names(p1), "batch")], p2[, setdiff(names(p2), "batch")])
})

test_that("removed samples are absent from beta/pheno but present in sample_qc.csv", {
  skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  b <- fixture_basenames_450k()
  run_dir <- copy_run_450k(exclude = b[1], intensity_threshold = 11.8, pheno = fixture_pheno_df())
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  removed <- unique(rbind(read_samples_removed(run_dir, "5723646052"),
                          read_samples_removed(run_dir, "5723646053"))$basename)
  expect_true(length(removed) >= 2)
  beta <- readRDS(file.path(run_dir, "beta_merged.rds"))
  expect_identical(colnames(beta), setdiff(b, removed))
  pheno <- readRDS(file.path(run_dir, "pheno_postQC.rds"))
  expect_identical(pheno$basename, setdiff(b, removed))
  expect_true(all(c("plate", "sex", "age", "group", "person") %in% names(pheno)))
  qc <- utils::read.csv(file.path(run_dir, "sample_qc.csv"), stringsAsFactors = FALSE)
  expect_identical(qc$basename[qc$removed], removed)
  expect_match(qc$reasons[qc$basename == b[1]], "excluded")
  expect_true(all(nzchar(qc$reasons[qc$removed])))
  expect_true(all(qc$reasons[!qc$removed] == ""))
  expect_identical(colnames(readRDS(file.path(run_dir, "rs_merged.rds"))), setdiff(b, removed))
})

test_that("impute = FALSE leaves NAs; save_na_mask records them", {
  skip_if_not_installed("ChAMP")
  run_dir <- copy_run_450k(impute = FALSE, save_na_mask = TRUE, failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  beta <- readRDS(file.path(run_dir, "beta_merged.rds"))
  expect_true(anyNA(beta))
  imp <- utils::read.csv(file.path(run_dir, "imputed.csv"))
  expect_false(any(imp$imputed))
  expect_equal(imp$n_missing, unname(colSums(is.na(beta))))
  mask <- utils::read.csv(file.path(run_dir, "imputed_mask.csv"), stringsAsFactors = FALSE)
  expect_identical(nrow(mask), sum(is.na(beta)))
  expect_true(all(is.na(beta[cbind(match(mask$probe, rownames(beta)), match(mask$basename, colnames(beta)))])))
  expect_match(readLines(file.path(run_dir, "log.txt")), "left in place", all = FALSE)
})

test_that("probe lists can be switched off", {
  skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- copy_run_450k(remove_chrY = FALSE, remove_zhou = FALSE)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  removed <- utils::read.csv(file.path(run_dir, "probes_removed.csv"), stringsAsFactors = FALSE)
  expect_false(any(removed$reason %in% c("chrY", "zhou")))
  expect_true(all(c("snp", "non_cpg") %in% removed$reason))
})

test_that("call_sex follows minfi's cutoff and warns on disagreement", {
  expect_identical(call_sex(c(12, 12, 12), c(8, 11.5, 12)), c("F", "M", "M"))
  expect_identical(call_sex(c(12, NA), c(8, 8)), c("F", NA))
  expect_identical(call_sex(12, 8), "F")
})

test_that("qc_merge needs normalised batches", {
  run_dir <- copy_run_450k()
  expect_error(qc_merge(run_dir), "Run qc_filter_samples\\(\\) first|Run qc_normalise\\(\\) first")
})
