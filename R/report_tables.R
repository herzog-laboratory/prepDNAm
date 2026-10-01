#' Report tables
#'
#' Data frames behind the report pages; usable directly for custom reports.
#' All take a run directory or a [qc_data()] object.
#'
#' * `qc_table_batches()`: one row per batch with sample counts, removals by
#'   reason and pass rate.
#' * `qc_table_samples()`: one row per sample: sample sheet columns, QC
#'   metrics (`mMed`, `uMed`, `frac_failed_probes`, `rho`, `snp_outlier`,
#'   `n_missing`, `predicted_sex`), `removed` and `reasons`. Before
#'   [qc_merge()] only the sample sheet and per-batch removals are available.
#' * `qc_table_probes()`: probe removal counts by reason (a probe can have
#'   several reasons; `unique` is the number of distinct probes removed).
#' * `qc_thresholds()`: the thresholds and switches of the run.
#'
#' @param x Run directory or [qc_data()] object.
#' @return A data frame.
#' @examples
#' \dontrun{
#' qc_table_batches("runs/study1")
#' }
#' @name qc_tables
NULL

#' @rdname qc_tables
#' @export
qc_table_batches <- function(x) {
  rd <- qc_data(x)
  removed <- batch_removals(rd)
  reasons <- c("excluded", "low_intensity", "failed_probes")
  out <- do.call(rbind, lapply(rd$batches, function(b) {
    n <- sum(rd$samples$batch == b)
    r <- removed[removed$batch == b, , drop = FALSE]
    n_rm <- length(unique(r$basename))
    counts <- vapply(reasons, function(z) sum(r$reason == z), integer(1))
    data.frame(batch = b, n_samples = n, n_removed = n_rm, n_retained = n - n_rm,
               pass_rate = if (n > 0) (n - n_rm) / n else NA_real_,
               t(counts), stringsAsFactors = FALSE)
  }))
  rownames(out) <- NULL
  out
}

#' @rdname qc_tables
#' @export
qc_table_samples <- function(x) {
  rd <- qc_data(x)
  s <- rd$samples
  if (!is.null(rd$sample_qc)) {
    m <- rd$sample_qc[, setdiff(names(rd$sample_qc), "batch"), drop = FALSE]
    s <- merge(s, m, by = "basename", sort = FALSE)
    s <- s[match(rd$samples$basename, s$basename), , drop = FALSE]
  } else {
    removed <- batch_removals(rd)
    s$removed <- s$basename %in% removed$basename
    s$reasons <- vapply(s$basename, function(b) {
      paste(sort(unique(removed$reason[removed$basename == b])), collapse = ";")
    }, character(1))
  }
  rownames(s) <- NULL
  s
}

#' @rdname qc_tables
#' @export
qc_table_probes <- function(x) {
  rd <- qc_data(x)
  if (is.null(rd$probes_removed)) return(NULL)
  order <- c("chrY", "snp", "non_cpg", "zhou", "detP", "na_fraction", "not_in_all_batches")
  counts <- table(factor(rd$probes_removed$reason, levels = order))
  out <- data.frame(reason = names(counts), n_probes = as.integer(counts), stringsAsFactors = FALSE)
  out <- out[out$n_probes > 0 | out$reason %in% c("detP", "na_fraction"), , drop = FALSE]
  rbind(out, data.frame(reason = "unique", n_probes = length(unique(rd$probes_removed$probe))))
}

#' @rdname qc_tables
#' @export
qc_thresholds <- function(x) {
  rd <- qc_data(x)
  cfg <- rd$config
  data.frame(
    setting = c("intensity_threshold", "detp_threshold", "failed_probe_threshold",
                "failed_sample_threshold", "probe_na_threshold", "remove_chrY", "remove_snp",
                "remove_non_cpg", "remove_zhou", "impute", "normalise", "compatible_subset"),
    value = c(cfg$intensity_threshold, cfg$detp_threshold, cfg$failed_probe_threshold,
              cfg$failed_sample_threshold, cfg$probe_na_threshold, cfg$remove_chrY, cfg$remove_snp,
              cfg$remove_non_cpg, cfg$remove_zhou, cfg$impute,
              if (cfg$input_type == "beta") as.character(cfg$normalise) else "BMIQ",
              cfg$compatible_subset),
    description = c("min median log2 intensity (methylated and unmethylated) per sample",
                    "detection p-value above which a probe fails in a sample",
                    "max fraction of failed probes per sample",
                    "max fraction of samples in which a probe may fail",
                    "max fraction of samples in which a probe may be missing (before imputation)",
                    "remove chrY probes", "remove probes overlapping SNPs", "remove non-CpG probes",
                    "remove cross-reactive/unreliable probes (Zhou et al. 2017)",
                    "kNN imputation of missing values", "probe-type normalisation",
                    "EPIC v1/v2 compatible subset written"),
    stringsAsFactors = FALSE)
}

# samples_removed.csv of every batch, stacked (empty data frame if none yet).
batch_removals <- function(rd) {
  parts <- lapply(rd$batches, function(b) {
    p <- file.path(batch_dir(rd$run_dir, b), "samples_removed.csv")
    if (file.exists(p)) read_samples_removed(rd$run_dir, b) else NULL
  })
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (length(parts) == 0) {
    return(data.frame(basename = character(), batch = character(), reason = character(),
                      value = numeric(), threshold = numeric(), stringsAsFactors = FALSE))
  }
  do.call(rbind, parts)
}
