#' Export the merged beta matrix and metadata as CSV
#'
#' Writes the final beta matrix and sample metadata of a run as CSV files in
#' the DNA Methylation Array Data Standard V-2410, the standard format for
#' DNAm array data read by `biolearn`
#' (<https://bio-learn.github.io/methylation-standard.html>):
#'
#' * methylation matrix: CpG sites in rows (first column `cpgSite`), one
#'   column per sample with the sample ID as header, betas rounded to three
#'   decimals as the standard requires (`digits`), missing values as `NaN`; at most `max_samples` samples per
#'   file (`<prefix>_methylation_part<k>.csv` if split, else
#'   `<prefix>_methylation.csv`).
#' * metadata (`<prefix>_metadata.csv`): one row per sample with `SampleID`,
#'   `Sex` (0 = female, 1 = male, `NaN` unknown, from the reported sex
#'   column), `Age` (years, `NaN` unknown), `PredictedSex` (same coding, from
#'   chrX/chrY intensities), the other sample sheet columns and, if
#'   [qc_estimates()] has run, the cell-type, age and marker CpG estimates.
#'
#' Probe IDs are written as in `beta_merged.rds`: for EPIC v2 and mouse with
#' their design suffix (`cg00000029_TC21`), so replicate probes of a CpG stay
#' separate rows. For one row per CpG with EPIC v1 IDs, export the compatible
#' matrix.
#'
#' For EPIC v2, `compatible = TRUE` exports `beta_merged_compatible` instead
#' (written by [qc_merge()] with `compatible_subset = TRUE`, the default for
#' EPIC v2): only probes with identical sequence on EPIC v1 and v2, with EPIC
#' v1 IDs (see [subset_versionshared_CpGs()]), for combining with EPIC v1 data
#' or using tools trained on EPIC v1. File names then get `_compatible`.
#'
#' `run_prepDNAm()` calls this at the end of a run if the configuration has
#' `export_csv = TRUE` (for EPIC v2 both the full and the compatible matrix).
#'
#' @param x Run directory or [qc_data()] object of a merged run.
#' @param out_dir Output directory; defaults to `export/` in the run
#'   directory.
#' @param prefix File name prefix; defaults to the run name.
#' @param sample_id Sample sheet column with the sample IDs to use as column
#'   headers (e.g. GEO accessions, which the standard requires for GEO data);
#'   defaults to `basename`. Must be unique.
#' @param sex_col,age_col Sample sheet columns with reported sex and age;
#'   `NULL` to auto-detect (`sex`/`gender`, `age`).
#' @param max_samples Maximum number of samples per methylation file (the
#'   standard says 1000).
#' @param compatible Export the EPIC v1/v2 compatible matrix (EPIC v2 runs
#'   only).
#' @param digits Decimal places to round betas to; the standard asks for at
#'   most 3. `NULL` writes full precision (15 significant digits), so values
#'   match `beta_merged.rds` (or `beta_merged_compatible.rds`); the file then
#'   no longer follows the standard strictly.
#' @return The paths of the files written, invisibly.
#' @examples
#' \dontrun{
#' export_beta_csv("runs/study1", sample_id = "geo_accession")
#' export_beta_csv("runs/study1", out_dir = "shared/study1", prefix = "study1")
#' export_beta_csv("runs/epicv2", compatible = TRUE)   # EPIC v1 IDs
#' }
#' @export
export_beta_csv <- function(x, out_dir = NULL, prefix = NULL, sample_id = "basename",
                            sex_col = NULL, age_col = NULL, max_samples = 1000,
                            compatible = FALSE, digits = 3) {
  rd <- qc_data(x)
  if (!rd$has[["beta_merged"]]) rlang::abort("No beta_merged.rds; run qc_merge() first.")
  setup_logging(rd$run_dir)
  check_string(sample_id, "sample_id")
  check_number(max_samples, "max_samples", min = 1)
  check_flag(compatible, "compatible")
  if (!is.null(digits)) check_number(digits, "digits", min = 0)
  if (compatible && !rd$has[["beta_compatible"]]) {
    rlang::abort(c("No beta_merged_compatible.rds in this run.",
                   i = "It is written for EPIC v2 runs with `compatible_subset = TRUE` (the default)."))
  }
  if (is.null(out_dir)) out_dir <- run_file(rd$run_dir, "export")
  if (is.null(prefix)) prefix <- default_export_prefix(rd$config, compatible)
  check_string(prefix, "prefix")

  beta <- report_beta(rd, compatible = compatible)
  s <- read_samples(rd$run_dir)
  s <- s[match(colnames(beta), s$basename), , drop = FALSE]
  if (!sample_id %in% names(s)) rlang::abort(sprintf("No column \"%s\" in the sample sheet.", sample_id))
  ids <- as.character(s[[sample_id]])
  if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) {
    rlang::abort(sprintf("Sample IDs in column \"%s\" must be unique and non-missing.", sample_id))
  }

  # ---- methylation matrix ----
  beta <- pmin(pmax(beta, 0), 1)
  if (!is.null(digits)) beta <- round(beta, digits)
  colnames(beta) <- ids

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  chunks <- split(seq_len(ncol(beta)), ceiling(seq_len(ncol(beta)) / max_samples))
  files <- if (length(chunks) == 1) {
    file.path(out_dir, sprintf("%s_methylation.csv", prefix))
  } else {
    file.path(out_dir, sprintf("%s_methylation_part%d.csv", prefix, seq_along(chunks)))
  }
  for (k in seq_along(chunks)) {
    dt <- data.table::data.table(cpgSite = rownames(beta), beta[, chunks[[k]], drop = FALSE])
    data.table::fwrite(dt, files[k], na = "NaN", quote = FALSE)
  }

  # ---- metadata ----
  sex_code <- function(v) ifelse(is.na(v), "NaN", ifelse(v == "M", "1", "0"))
  col_sex <- find_sex_col(s, sex_col)
  col_age <- find_age_col(s, age_col)
  age <- if (is.null(col_age)) rep(NA_real_, nrow(s)) else suppressWarnings(as.numeric(s[[col_age]]))
  qc <- rd$sample_qc
  predicted <- if (!is.null(qc) && "predicted_sex" %in% names(qc)) {
    qc$predicted_sex[match(s$basename, qc$basename)]
  } else rep(NA_character_, nrow(s))
  meta <- data.frame(
    SampleID = ids,
    Sex = sex_code(if (is.null(col_sex)) rep(NA_character_, nrow(s)) else normalise_sex(s[[col_sex]])),
    Age = ifelse(is.na(age), "NaN", format(age, trim = TRUE, drop0trailing = TRUE)),
    PredictedSex = sex_code(predicted),
    stringsAsFactors = FALSE)
  other <- setdiff(names(s), c(sample_id, col_sex, col_age, "path", "exclude", "pheno_path", "pheno_dir"))
  if (sample_id != "basename") other <- union("basename", other)
  meta <- cbind(meta, s[, other, drop = FALSE])
  est_file <- run_file(rd$run_dir, "estimates.csv")
  if (file.exists(est_file)) {
    est <- utils::read.csv(est_file, stringsAsFactors = FALSE, check.names = FALSE,
                           colClasses = c(basename = "character"))
    meta <- cbind(meta, est[match(s$basename, est$basename), setdiff(names(est), "basename"), drop = FALSE])
  }
  meta_file <- file.path(out_dir, sprintf("%s_metadata.csv", prefix))
  data.table::fwrite(meta, meta_file, na = "")

  log_info(sprintf("export_beta_csv: %d CpGs x %d samples written to %s", nrow(beta), ncol(beta), out_dir))
  invisible(c(files, meta_file))
}

default_export_prefix <- function(cfg, compatible = FALSE) {
  paste0(gsub("[^A-Za-z0-9._-]+", "_", cfg$run_name), if (compatible) "_compatible")
}

# Export stage run by run_prepDNAm() when the config has export_csv = TRUE:
# the full matrix, plus the compatible one where it exists. The metadata file
# is written last, so it marks a finished export.
export_stage <- function(run_dir, cfg, force = FALSE) {
  compat <- c(FALSE, if (isTRUE(cfg$compatible_subset)) TRUE)
  for (cmp in compat) {
    stage <- if (cmp) "export_csv:compatible" else "export_csv"
    marker <- run_file(run_dir, "export", sprintf("%s_metadata.csv", default_export_prefix(cfg, cmp)))
    run_stage(run_dir, stage, marker, force, seed = NULL,
              body = function() export_beta_csv(run_dir, compatible = cmp,
                                                digits = export_digits_value(cfg$export_digits %||% 3)))
  }
  invisible(NULL)
}

# export_digits of the config: a number, or "full" for no rounding (NULL).
export_digits_value <- function(x) {
  if (is.character(x) && tolower(x) %in% c("full", "none")) NULL else x
}
