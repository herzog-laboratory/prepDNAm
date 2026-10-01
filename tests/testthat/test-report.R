skip_if_no_report_run <- function() {
  skip_if_no_fixture_data(); skip_if_not_installed("ChAMP"); skip_if_not_installed("impute")
}

test_that("qc_data loads a merged run and reports what is present", {
  skip_if_no_report_run()
  rd <- qc_data(fixture_merged_450k())
  expect_s3_class(rd, "prepDNAm_data")
  expect_identical(rd$batches, c("5723646052", "5723646053"))
  expect_true(all(rd$has[c("ctrls", "densities", "rho", "qc", "detP", "sex", "rs", "beta_merged")]))
  expect_false(rd$has[["beta_compatible"]])
  expect_identical(nrow(rd$sample_qc), 6L)
  expect_identical(qc_data(rd), rd)
  # before merge: sample_qc/pheno absent
  rd0 <- qc_data(copy_run_450k())
  expect_null(rd0$sample_qc); expect_null(rd0$pheno)
  expect_false(rd0$has[["beta_merged"]])
})

test_that("tables: batches, samples, probes, thresholds", {
  skip_if_no_report_run()
  rd <- qc_data(fixture_merged_450k())
  b <- qc_table_batches(rd)
  expect_named(b, c("batch", "n_samples", "n_removed", "n_retained", "pass_rate", "excluded", "low_intensity", "failed_probes"))
  expect_identical(b$n_samples, c(3L, 3L)); expect_identical(b$pass_rate, c(1, 1))
  s <- qc_table_samples(rd)
  expect_identical(s$basename, fixture_basenames_450k())
  expect_true(all(c("mMed", "rho", "predicted_sex", "removed", "reasons") %in% names(s)))
  expect_false(any(s$removed))
  p <- qc_table_probes(rd)
  expect_identical(p$reason[nrow(p)], "unique")
  expect_true(all(c("chrY", "snp", "non_cpg", "zhou", "detP") %in% p$reason))
  th <- qc_thresholds(rd)
  expect_identical(th$value[th$setting == "intensity_threshold"], "9.5")
  # before merge, samples table still works from per-batch removals
  run0 <- copy_run_450k(exclude = fixture_basenames_450k()[1])
  qc_filter_samples(run0)
  s0 <- qc_table_samples(run0)
  expect_identical(s0$reasons[1], "excluded"); expect_true(s0$removed[1])
  expect_identical(qc_table_batches(run0)$excluded, c(1L, 0L))
})

test_that("plots return ggplots (or NULL where data is absent)", {
  skip_if_no_report_run()
  rd <- qc_data(fixture_merged_450k())
  for (f in list(qc_plot_intensities, qc_plot_failed_probes, qc_plot_rho, qc_plot_sex)) {
    expect_s3_class(f(rd), "ggplot")
  }
  expect_setequal(qc_density_stages(rd), c("raw", "noob", "bmiq"))
  expect_s3_class(qc_plot_densities(rd, "raw"), "ggplot")
  expect_s3_class(qc_plot_densities(rd, "bmiq", batch = "5723646052"), "ggplot")
  thin <- qc_plot_densities(rd, "raw", n_points = 100)$data
  expect_identical(nrow(thin), 100L * 6L)
  expect_true(all(startsWith(thin$hover, as.character(thin$basename))))
  expect_null(qc_plot_densities(rd, "nope"))
  types <- qc_control_types(rd)
  expect_true(all(c("BISULFITE CONVERSION I", "NEGATIVE", "NORM_A") %in% types))
  expect_s3_class(qc_plot_controls(rd, "BISULFITE CONVERSION I"), "ggplot")
  expect_s3_class(qc_plot_controls(rd, "NORM_A"), "ggplot")
  expect_null(qc_plot_controls(rd, "NOT A TYPE"))
  # beta input: no intensities, controls or densities beyond raw
  run_b <- file.path(withr::local_tempdir(), "run")
  qc_init(suppressMessages(qc_config(input = fixture_beta()$rds, run_dir = run_b, array = "450k", batch = "none")))
  qc_read(run_b)
  rdb <- qc_data(run_b)
  expect_null(qc_plot_intensities(rdb)); expect_null(qc_plot_rho(rdb)); expect_null(qc_plot_sex(rdb))
  expect_identical(qc_control_types(rdb), character())
  expect_identical(qc_density_stages(rdb), "raw")
})

test_that("sex table compares against a reported sex column", {
  skip_if_no_report_run()
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  t <- qc_table_sex(run_dir)
  expect_identical(t$reported_sex, fixture_pheno_df()$sex)
  expect_false(any(t$sex_mismatch))
  pheno <- fixture_pheno_df(); pheno$sex[2] <- "male"   # sample 2 is female
  run2 <- copy_run_450k(pheno = pheno, failed_sample_threshold = 0.5)
  qc_filter_samples(run2); qc_normalise(run2); qc_merge(run2)
  t2 <- qc_table_sex(run2)
  expect_identical(which(t2$sex_mismatch), 2L)
  expect_error(qc_table_sex(run2, sex_col = "nope"), "No column")
  expect_identical(normalise_sex(c("F", "female", "M", "Male", "x", NA)), c("F", "F", "M", "M", NA, NA))
})

test_that("genotype agreement recovers the minfiData donors", {
  skip_if_no_report_run(); skip_if_not_installed("ewastools")
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  expect_message(g <- qc_snp_agreement(run_dir), "donor_id")
  expect_identical(g$donor_col, "person")
  expect_identical(nrow(g$conflicts), 0L)
  d <- g$donors
  # samples of the same person share an inferred donor
  per_person <- tapply(d$inferred_donor, d$reported_donor, function(x) length(unique(x)))
  expect_identical(as.vector(per_person), c(1L, 1L, 1L))
  expect_identical(length(unique(d$inferred_donor)), 3L)
  expect_true(all(d$n_samples_of_donor == 2))
  # a mislabelled sample is a conflict
  pheno <- fixture_pheno_df(); pheno$person[1] <- "id1"
  run2 <- copy_run_450k(pheno = pheno, failed_sample_threshold = 0.5)
  qc_filter_samples(run2); qc_normalise(run2); qc_merge(run2)
  expect_gt(nrow(suppressMessages(qc_snp_agreement(run2))$conflicts), 0)
  # without rs
  expect_message(expect_null(qc_snp_agreement(copy_run_450k())), "rs_merged")
})

test_that("biology helpers degrade gracefully", {
  skip_if_no_report_run()
  rd <- qc_data(fixture_merged_450k())
  mk <- qc_index_cpgs(rd)
  expect_true(is.null(mk) || all(c("basename") %in% names(mk)))
  if (!requireNamespace("EpiDISH", quietly = TRUE)) expect_message(expect_null(qc_estimate_ic(rd)), "EpiDISH")
  if (!requireNamespace("WIDclocks", quietly = TRUE)) expect_message(expect_null(qc_estimate_age(rd)), "WIDclocks")
  # with the packages installed, called from the package namespace (not attached)
  if (requireNamespace("EpiDISH", quietly = TRUE)) {
    ic <- qc_estimate_ic(rd)
    expect_named(ic, c("basename", "Epi", "Fib", "IC"))
    expect_equal(rowSums(ic[, -1]), rep(1, 6), tolerance = 1e-6)
    h <- qc_estimate_ic(rd, subtypes = TRUE)
    expect_identical(names(h)[1:3], c("basename", "Epi", "Fib"))
    expect_identical(ncol(h), 15L)                       # Epi, Fib + 12 immune types
    expect_equal(rowSums(h[, -1]), rep(1, 6), tolerance = 1e-6)
    expect_equal(h$Epi, ic$Epi, tolerance = 1e-6)
    expect_s3_class(qc_plot_ic(rd, subtypes = TRUE), "ggplot")
  }
  if (requireNamespace("WIDclocks", quietly = TRUE)) {
    age <- suppressMessages(qc_estimate_age(rd))
    expect_identical(age$basename, colnames(readRDS(file.path(rd$run_dir, "beta_merged.rds"))))
    expect_true(all(is.finite(age$wid_general_age)))
  }
  for (f in list(qc_plot_ic, qc_plot_age, qc_plot_index_cpgs)) {
    p <- suppressMessages(f(rd))
    expect_true(is.null(p) || inherits(p, "ggplot"))
  }
  if (!is.null(mk)) expect_s3_class(qc_plot_index_cpgs(rd), "ggplot")
  expect_error(qc_plot_age(rd, age_col = "nope"), "No column")
})

test_that("SNP scatter plots per donor", {
  skip_if_no_report_run(); skip_if_not_installed("ewastools")
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  pl <- suppressMessages(qc_plot_snp_pairs(run_dir))
  expect_length(pl, 3L)                       # three persons, two samples each
  expect_true(all(vapply(pl, inherits, logical(1), "ggplot")))
  expect_identical(length(unique(pl[[1]]$data$pair)), 1L)
  expect_match(as.character(pl[[1]]$data$pair[1]), "r = ")
  expect_length(suppressMessages(qc_plot_snp_pairs(run_dir, donors = names(pl)[1])), 1L)
  # without a donor column, inferred donors are used
  pheno <- fixture_pheno_df(); pheno$person <- NULL
  run2 <- copy_run_450k(pheno = pheno, failed_sample_threshold = 0.5)
  qc_filter_samples(run2); qc_normalise(run2); qc_merge(run2)
  pl2 <- qc_plot_snp_pairs(run2)
  expect_length(pl2, 3L)
  expect_match(pl2[[1]]$labels$title, "inferred donor")
  expect_message(expect_null(qc_plot_snp_pairs(copy_run_450k())), "rs_merged")
  f <- suppressMessages(qc_snp_pairs_pdf(run_dir, file = withr::local_tempfile(fileext = ".pdf")))
  expect_true(file.exists(f) && file.size(f) > 0)
})

test_that("donor column: donor_id preferred, aliases with a message", {
  s <- data.frame(basename = "a", person = "p1", donor_id = "d1")
  expect_identical(find_donor_col(s, NULL), "donor_id")
  expect_message(expect_identical(find_donor_col(s[, 1:2], NULL), "person"), "donor_id")
  expect_null(find_donor_col(s[, 1, drop = FALSE], NULL))
  expect_identical(find_donor_col(s, "person"), "person")
  expect_error(find_donor_col(s, "nope"), "No column")
})

test_that("dimensionality reduction and association screen", {
  skip_if_no_report_run()
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  dr <- qc_dimred(run_dir, n_probes = 5000)
  expect_s3_class(dr, "prepDNAm_dimred")
  expect_identical(dr$n_probes, 5000L)
  expect_identical(nrow(dr$coords), 6L)
  expect_true(all(c("PC1", "PC5") %in% names(dr$coords)))
  expect_false("PC6" %in% names(dr$coords))
  expect_named(dr$var_explained, paste0("PC", 1:5))
  expect_true(all(c("batch", "sex", "age", "group", "person", "sentrix_id", "sentrix_position") %in% dr$variables))
  expect_false(any(c("basename", "path", "removed", "reasons", "mMed", "uMed", "rho", "xMed",
                     "frac_failed_probes", "predicted_sex") %in% dr$variables))
  # technical first, then biological
  expect_identical(unname(dr$groups[c("batch", "plate", "sentrix_id", "sentrix_position", "sex", "age", "person")]),
                   rep(c("technical", "biological"), c(4, 3)))
  expect_identical(dr$variables, names(dr$groups))
  expect_true(all(diff(match(dr$groups, c("technical", "biological"))) >= 0))
  expect_identical(unname(qc_dimred(run_dir, n_probes = 5000, technical = "group")$groups["group"]), "technical")
  expect_identical(add_sentrix_columns(data.frame(basename = "a_b", Sentrix_ID = 1))$basename, "a_b")
  s2 <- add_sentrix_columns(data.frame(basename = c("206467010001_R01C01", "206467010001_R02C01")))
  expect_identical(s2$sentrix_position, c("R01C01", "R02C01"))
  expect_identical(dim(dr$pvals), c(length(dr$variables), 5L))
  expect_true(all(dr$pvals >= 0 & dr$pvals <= 1, na.rm = TRUE))
  expect_s3_class(qc_plot_dimred(dr, "pca"), "ggplot")
  expect_s3_class(qc_plot_dimred(dr, "pca", colour = "sex"), "ggplot")
  # categorical variables with many values: no legend, count in the title
  many <- qc_plot_dimred(dr, "pca", colour = "sex", max_legend = 1)
  expect_match(many$labels$title, "2 values")
  expect_identical(many$scales$get_scales("colour")$guide, "none")
  expect_identical(qc_plot_dimred(dr, "pca", colour = "sex")$scales$get_scales("colour")$guide, "legend")
  expect_identical(truncate_label(c("short", strrep("x", 30)), 25), c("short", paste0(strrep("x", 24), "\u2026")))
  pa <- qc_plot_pc_association(dr)
  expect_s3_class(pa, "ggplot")
  expect_setequal(levels(pa$data$group), c("technical", "biological"))
  expect_setequal(levels(pa$data$variable), dr$variables)
  expect_s3_class(qc_plot_pc_association(dr, cluster = FALSE), "ggplot")
  if (dr$has_umap) expect_s3_class(qc_plot_dimred(dr, "umap"), "ggplot") else expect_null(qc_plot_dimred(dr, "umap"))
  expect_error(qc_dimred(copy_run_450k()), "beta_merged")
})

test_that("qc_report renders the site", {
  skip_if_no_report_run()
  skip_if_not_installed("rmarkdown"); skip_if_not_installed("plotly"); skip_if_not_installed("DT")
  skip_if_not(rmarkdown::pandoc_available("2.0"), "pandoc not available")
  run_dir <- copy_run_450k(pheno = fixture_pheno_df(), failed_sample_threshold = 0.5)
  qc_filter_samples(run_dir); qc_normalise(run_dir); qc_merge(run_dir)
  expect_true(qc_report(run_dir))
  pages <- c("index", "1-samples", "2-controls", "3-intensities", "4-beta-densities", "5-identity",
             "6-dimred", "7-run")
  expect_true(all(file.exists(file.path(run_dir, "report", paste0(pages, ".html")))))
  idx <- readLines(file.path(run_dir, "report", "index.html"))
  expect_match(idx, "6 / 6", all = FALSE)
  expect_match(idx, "probes retained", all = FALSE)
  expect_false(any(grepl("## Error", readLines(file.path(run_dir, "report", "5-identity.html")))))
  expect_false(qc_report(run_dir))
  expect_identical(read_run_info(run_dir)$stages$qc_report$status, "ok")
  # self-contained by default: no shared library folder, and every page with
  # plotly widgets embeds plotly.js (widgets printed in loops used to lose it)
  expect_false(dir.exists(file.path(run_dir, "report", "site_libs")))
  for (pg in pages) {
    html <- readLines(file.path(run_dir, "report", paste0(pg, ".html")), warn = FALSE)
    if (any(grepl("plotly html-widget", html, fixed = TRUE))) {
      expect_true(any(grepl("plotly.js v", html, fixed = TRUE)), label = paste(pg, "embeds plotly.js"))
    }
  }
  expect_false(any(grepl("<h1", grep("samples retained", idx, value = TRUE))))
  expect_match(idx, "Input type", all = FALSE)
  # options: shared libraries, SNP scatter plots
  expect_true(qc_report(run_dir, force = TRUE, self_contained = FALSE, donor_col = "person"))
  expect_true(dir.exists(file.path(run_dir, "report", "site_libs")))
  expect_match(readLines(file.path(run_dir, "report", "5-identity.html")), "column <code>person</code>", all = FALSE)
})

test_that("qc_report explains a missing pandoc", {
  skip_if_no_fixture_data()
  run_dir <- copy_run_450k()
  withr::local_envvar(RSTUDIO_PANDOC = "/nonexistent", PATH = "/nonexistent")
  skip_if(rmarkdown::pandoc_available("2.0"), "pandoc still found")
  expect_error(qc_report(run_dir), "RSTUDIO_PANDOC")
})
