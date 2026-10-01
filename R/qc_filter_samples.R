#' Flag samples to remove from a batch
#'
#' Third stage. For each batch, decides which samples are dropped before
#' normalisation and writes `batches/<batch>/samples_removed.csv` with one row
#' per sample and reason:
#'
#' | `reason` | condition | `value` | `threshold` |
#' |---|---|---|---|
#' | `excluded` | listed in `exclude` / flagged in the sample sheet | `NA` | `NA` |
#' | `low_intensity` | median methylated or unmethylated log2 intensity below `intensity_threshold` | `min(mMed, uMed)` | `intensity_threshold` |
#' | `failed_probes` | fraction of probes with detection p-value above `detp_threshold` exceeds `failed_probe_threshold` | that fraction | `failed_probe_threshold` |
#'
#' Intensity and detection-p-value criteria need `qc.rds` / `detP.rds` from
#' [qc_read()] and are skipped (with a log message) when absent, e.g. for
#' beta input. The file is written even when no sample is removed.
#'
#' @inheritParams qc_read
#' @return Named logical vector, one element per batch: `TRUE` if run,
#'   `FALSE` if skipped; invisibly.
#' @examples
#' \dontrun{
#' qc_filter_samples("runs/study1")
#' read.csv("runs/study1/batches/plate1/samples_removed.csv")
#' }
#' @export
qc_filter_samples <- function(run_dir, batch = NULL, force = FALSE) {
  cfg <- read_config(run_dir)
  samples <- read_samples(run_dir)
  check_flag(force, "force")
  batch <- check_batches(batch, samples)

  ran <- vapply(batch, function(b) {
    bdir <- batch_dir(run_dir, b)
    bs <- samples[samples$batch == b, , drop = FALSE]
    out <- file.path(bdir, "samples_removed.csv")
    run_stage(run_dir, paste0("qc_filter_samples:", b), out, force, seed = NULL, body = function() {
      removed <- filter_batch_samples(cfg, bs, bdir, b)
      utils::write.csv(removed, out, row.names = FALSE)
      n <- length(unique(removed$basename))
      log_info(sprintf("qc_filter_samples:%s: %d of %d samples removed%s", b, n, nrow(bs),
                       if (n > 0) paste0(" (", paste(sprintf("%s: %d", names(table(removed$reason)),
                                                              table(removed$reason)), collapse = ", "), ")")
                       else ""))
    })
  }, logical(1))
  invisible(ran)
}

# One row per (sample, reason).
filter_batch_samples <- function(cfg, bs, bdir, b) {
  rows <- list()
  add <- function(basename, reason, value, threshold) {
    if (length(basename) == 0) return()
    rows[[length(rows) + 1]] <<- data.frame(
      basename = basename, batch = b, reason = reason,
      value = value, threshold = threshold, stringsAsFactors = FALSE)
  }
  add(bs$basename[bs$exclude], "excluded", NA_real_, NA_real_)

  qc_file <- file.path(bdir, "qc.rds")
  if (file.exists(qc_file)) {
    qc <- readRDS(qc_file)
    low <- pmin(qc$mMed, qc$uMed)
    bad <- !is.na(low) & low < cfg$intensity_threshold
    add(qc$basename[bad], "low_intensity", low[bad], cfg$intensity_threshold)
  } else {
    log_info(sprintf("qc_filter_samples:%s: no qc.rds, intensity filter skipped", b))
  }

  detp_file <- file.path(bdir, "detP.rds")
  if (file.exists(detp_file)) {
    detP <- readRDS(detp_file)
    frac <- colMeans(detP > cfg$detp_threshold, na.rm = TRUE)
    bad <- frac > cfg$failed_probe_threshold
    add(colnames(detP)[bad], "failed_probes", unname(frac[bad]), cfg$failed_probe_threshold)
  } else {
    log_info(sprintf("qc_filter_samples:%s: no detP.rds, failed-probe filter skipped", b))
  }

  if (length(rows) == 0) {
    return(data.frame(basename = character(), batch = character(), reason = character(),
                      value = numeric(), threshold = numeric(), stringsAsFactors = FALSE))
  }
  removed <- do.call(rbind, rows)
  removed[order(match(removed$basename, bs$basename), removed$reason), , drop = FALSE]
}

# Read samples_removed.csv of a batch.
read_samples_removed <- function(run_dir, batch) {
  path <- file.path(batch_dir(run_dir, batch), "samples_removed.csv")
  if (!file.exists(path)) {
    rlang::abort(sprintf("Batch %s has no samples_removed.csv. Run qc_filter_samples() first.", batch))
  }
  utils::read.csv(path, stringsAsFactors = FALSE,
                  colClasses = c(basename = "character", batch = "character", reason = "character",
                                 value = "numeric", threshold = "numeric"))
}

# Validate / default the `batch` argument of the per-batch stages.
check_batches <- function(batch, samples) {
  batches <- unique(samples$batch)
  if (is.null(batch)) return(batches)
  unknown <- setdiff(batch, batches)
  if (length(unknown) > 0) {
    rlang::abort(c(sprintf("Unknown batch(es): %s", paste(unknown, collapse = ", ")),
                   i = sprintf("Batches in samples.csv: %s", paste(batches, collapse = ", "))))
  }
  batch
}
