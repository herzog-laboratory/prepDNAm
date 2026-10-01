#' Load what the report pages need from a run directory
#'
#' Reads the configuration, sample sheet, per-sample QC metrics, probe and
#' sample removal records and per-batch bundles of a run once, so that the
#' report helper functions (`qc_table_*()`, `qc_plot_*()`) can be called
#' repeatedly without re-reading. Large matrices (`beta_merged`) are read on
#' demand by the functions that need them.
#'
#' @param run_dir Run directory.
#' @return A list of class `prepDNAm_data` with elements `run_dir`, `config`,
#'   `info` (run_info.json), `array` (see [array_info()]), `samples`,
#'   `sample_qc` (NULL before [qc_merge()]), `probes_removed`, `imputed`,
#'   `pheno` (post-QC), `batches` (character), and `has`, a named logical
#'   vector saying which optional items exist (`ctrls`, `densities`, `rho`,
#'   `sex`, `rs`, `beta_merged`, `beta_compatible`, `detP`).
#' @examples
#' \dontrun{
#' rd <- qc_data("runs/study1")
#' rd$has
#' }
#' @export
qc_data <- function(run_dir) {
  if (inherits(run_dir, "prepDNAm_data")) return(run_dir)
  cfg <- read_config(run_dir)
  samples <- read_samples(run_dir)
  batches <- unique(samples$batch)
  csv <- function(f, ...) {
    p <- run_file(run_dir, f)
    if (file.exists(p)) utils::read.csv(p, stringsAsFactors = FALSE, ...) else NULL
  }
  rds <- function(f) { p <- run_file(run_dir, f); if (file.exists(p)) readRDS(p) else NULL }
  in_batches <- function(f) all(file.exists(file.path(batch_dir(run_dir, batches), f)))
  in_any_batch <- function(f) any(file.exists(file.path(batch_dir(run_dir, batches), f)))

  sample_qc <- csv("sample_qc.csv", colClasses = c(basename = "character", batch = "character"))
  if (!is.null(sample_qc)) sample_qc$reasons[is.na(sample_qc$reasons)] <- ""
  structure(list(
    run_dir = run_dir,
    config = cfg,
    info = read_run_info(run_dir),
    array = array_info(cfg$array),
    samples = samples,
    sample_qc = sample_qc,
    probes_removed = csv("probes_removed.csv"),
    imputed = csv("imputed.csv", colClasses = c(basename = "character")),
    pheno = rds("pheno_postQC.rds"),
    batches = batches,
    has = c(
      ctrls = in_batches("ctrls.rds"),
      densities = in_batches("densities.rds"),
      rho = in_batches("rho.rds"),
      qc = in_batches("qc.rds"),
      detP = in_any_batch("detP.rds"),
      sex = in_batches("sex_raw.rds"),
      rs = file.exists(run_file(run_dir, "rs_merged.rds")),
      beta_merged = file.exists(run_file(run_dir, "beta_merged.rds")),
      beta_compatible = file.exists(run_file(run_dir, "beta_merged_compatible.rds")))
  ), class = "prepDNAm_data")
}

# Read one per-batch item across batches as a named list (NULL where absent).
batch_items <- function(rd, file) {
  stats::setNames(lapply(rd$batches, function(b) {
    p <- file.path(batch_dir(rd$run_dir, b), file)
    if (file.exists(p)) readRDS(p) else NULL
  }), rd$batches)
}

# The beta matrix used for indices trained on EPIC v1 (compatible subset for
# EPIC v2), or the merged matrix.
report_beta <- function(rd, compatible = FALSE) {
  f <- if (compatible && rd$has[["beta_compatible"]]) "beta_merged_compatible.rds" else "beta_merged.rds"
  readRDS(run_file(rd$run_dir, f))
}

# Sample sheet columns worth colouring/testing by: 2..(n-1) unique non-NA
# values, not identifiers, paths, dates or bookkeeping columns.
pheno_variables <- function(pheno) {
  skip <- c("basename", "path", "exclude", "removed", "reasons", "pheno_path", "pheno_dir")
  ok <- vapply(names(pheno), function(v) {
    if (v %in% skip) return(FALSE)
    x <- pheno[[v]]
    if (inherits(x, c("Date", "POSIXt", "difftime"))) return(FALSE)
    u <- length(unique(x[!is.na(x)]))
    if (u < 2) return(FALSE)
    if (is.character(x) || is.factor(x) || is.logical(x)) return(u < nrow(pheno))
    is.numeric(x)
  }, logical(1))
  names(pheno)[ok]
}
