merged_run_with_pheno <- function(env = parent.frame()) {
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5, env = env)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  run_dir
}

test_that("export_beta_csv writes the V-2410 matrix and metadata", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- merged_run_with_pheno()
  # one missing value to check the NaN coding
  b <- readRDS(file.path(run_dir, "beta_merged.rds")); b[1, 1] <- NA
  saveRDS(b, file.path(run_dir, "beta_merged.rds"))
  files <- export_beta_csv(run_dir, prefix = "study")
  expect_identical(basename(files), c("study_methylation.csv", "study_metadata.csv"))
  expect_true(all(startsWith(files, file.path(run_dir, "export"))))

  lines <- readLines(files[1], n = 3)
  expect_identical(strsplit(lines[1], ",")[[1]], c("cpgSite", colnames(b)))
  expect_match(lines[2], "^cg[0-9]+,NaN,")
  m <- utils::read.csv(files[1], check.names = FALSE)
  expect_identical(dim(m), c(nrow(b), ncol(b) + 1L))
  v <- as.matrix(m[, -1])
  expect_true(all(v >= 0 & v <= 1, na.rm = TRUE))
  expect_true(all(abs(v * 1000 - round(v * 1000)) < 1e-6, na.rm = TRUE))
  expect_equal(v[2, ], round(b[2, ], 3), ignore_attr = TRUE)

  meta <- utils::read.csv(files[2], colClasses = "character")
  expect_identical(names(meta)[1:4], c("SampleID", "Sex", "Age", "PredictedSex"))
  pheno <- fixture_pheno_df()
  pheno <- pheno[match(meta$SampleID, pheno$basename), ]
  expect_identical(meta$Sex, ifelse(normalise_sex(pheno$sex) == "M", "1", "0"))
  expect_equal(as.numeric(meta$Age), pheno$age)
  expect_identical(meta$PredictedSex, meta$Sex)
  expect_true(all(c("person", "batch") %in% names(meta)))
  expect_false(any(c("path", "exclude", "sex", "age") %in% names(meta)))
})

test_that("export_beta_csv splits files, checks sample IDs", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- merged_run_with_pheno()
  out <- withr::local_tempdir()
  files <- export_beta_csv(run_dir, out_dir = out, prefix = "x", max_samples = 4)
  expect_identical(basename(files), c("x_methylation_part1.csv", "x_methylation_part2.csv", "x_metadata.csv"))
  expect_identical(length(strsplit(readLines(files[1], n = 1), ",")[[1]]), 5L)
  expect_identical(length(strsplit(readLines(files[2], n = 1), ",")[[1]]), 3L)
  expect_error(export_beta_csv(run_dir, out_dir = out, sample_id = "person"), "unique")
  expect_error(export_beta_csv(run_dir, out_dir = out, sample_id = "nope"), "No column")
  expect_error(export_beta_csv(copy_run_450k()), "qc_merge")
})


test_that("compatible = TRUE exports beta_merged_compatible", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- merged_run_with_pheno()
  expect_error(export_beta_csv(run_dir, compatible = TRUE), "compatible_subset")
  # stand-in for an EPIC v2 run: a compatible matrix with v1 IDs
  b <- readRDS(file.path(run_dir, "beta_merged.rds"))[1:50, ]
  saveRDS(b, file.path(run_dir, "beta_merged_compatible.rds"))
  files <- export_beta_csv(run_dir, compatible = TRUE)
  expect_match(basename(files), "_compatible_(methylation|metadata)\\.csv$")
  m <- utils::read.csv(files[1], check.names = FALSE)
  expect_identical(m$cpgSite, rownames(b))
})

test_that("export_csv = TRUE exports at the end of run_prepDNAm, without re-running", {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
  run_dir <- merged_run_with_pheno()
  expect_false(dir.exists(file.path(run_dir, "export")))
  # switching the export on for a finished run is not a configuration change
  cfg <- read_config(run_dir); cfg$export_csv <- TRUE
  expect_no_error(res <- run_prepDNAm(cfg))
  expect_true(res$complete)
  expect_true(file.exists(file.path(run_dir, "export", sprintf("%s_metadata.csv", cfg$run_name))))
  expect_identical(res$stages$status[res$stages$stage == "export_csv"], "ok")
  expect_false("export_csv:compatible" %in% res$stages$stage)
  # resume: export skips
  n_log <- length(readLines(file.path(run_dir, "log.txt")))
  run_prepDNAm(cfg)
  new_lines <- readLines(file.path(run_dir, "log.txt"))[-seq_len(n_log)]
  expect_match(new_lines, "export_csv: outputs exist, skipping", all = FALSE)
})
