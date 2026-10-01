#' Run the whole pipeline
#'
#' Runs [qc_init()], [qc_read()], [qc_filter_samples()], [qc_normalise()] and
#' [qc_merge()] in order, then [qc_estimates()] (unless `estimates = FALSE`)
#' and [export_beta_csv()] if the configuration has `export_csv = TRUE`. Every stage skips work whose outputs already exist,
#' so calling this again on the same run directory resumes after a crash or
#' continues an interrupted run; `force = TRUE` recomputes everything.
#'
#' Required packages (array manifest / annotation, `ChAMP`, `impute`) are
#' checked before anything runs.
#'
#' @param x A [qc_config()] object, or the path of an existing run directory
#'   to resume (its `config.yml` is used).
#' @param force Recompute all stages even if their outputs exist.
#' @return A `prepDNAm_run` summary (see [run_summary()]), invisibly.
#' @examples
#' \dontrun{
#' cfg <- qc_config(input = "idats/", run_dir = "runs/study1", array = "EPIC",
#'                  pheno = "pheno.csv", cores = 8)
#' run_prepDNAm(cfg)
#' run_prepDNAm("runs/study1")   # resume
#' }
#' @export
run_prepDNAm <- function(x, force = FALSE) {
  cfg <- if (inherits(x, "prepDNAm_config")) {
    x
  } else if (is.character(x) && length(x) == 1 && dir.exists(x)) {
    read_config(x)
  } else {
    rlang::abort("`x` must be a `prepDNAm_config` or the path of an initialised run directory.")
  }
  check_flag(force, "force")
  check_run_dependencies(cfg)
  run_dir <- cfg$run_dir

  t0 <- Sys.time()
  qc_init(cfg, force = force)
  setup_logging(run_dir)
  log_info(sprintf("run_prepDNAm: %s (prepDNAm %s)", cfg$run_name, cfg$prepDNAm_version))
  qc_read(run_dir, force = force)
  qc_filter_samples(run_dir, force = force)
  qc_normalise(run_dir, force = force)
  qc_merge(run_dir, force = force)
  # runs without the setting (created before it existed) get estimates too
  if (!isFALSE(cfg$estimates)) qc_estimates(run_dir, force = force)
  if (isTRUE(cfg$export_csv)) export_stage(run_dir, cfg, force = force)
  log_info(sprintf("run_prepDNAm: finished in %.1f min", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  invisible(run_summary(run_dir))
}

# Abort up front if a package a run will need is missing.
check_run_dependencies <- function(cfg) {
  info <- array_info(cfg$array)
  needed <- c(info$manifest_pkg,
              if (cfg$save_sex) info$anno_pkg,
              if (cfg$input_type != "beta" || !isFALSE(cfg$normalise)) "ChAMP",
              if (cfg$impute) "impute",
              if (isTRUE(cfg$save_shiny)) "shinyMethyl")
  suppressWarnings(for (p in needed) requireNamespace(p, quietly = TRUE))
  rlang::check_installed(needed, reason = "to run this configuration")
  if (info$has_snp_probes && cfg$input_type != "beta" &&
      !requireNamespace("ewastools", quietly = TRUE)) {
    rlang::inform("Package ewastools is not installed; the SNP outlier metric will be NA. Install with remotes::install_github(\"hhhh5/ewastools\").")
  }
  invisible(needed)
}

#' Summarise a run directory
#'
#' Collects what `run_info.json`, `samples.csv` and the stage outputs say about
#' a run: sample and probe counts, stage status and timings, and any failed
#' stage.
#'
#' @param run_dir Run directory.
#' @return A list of class `prepDNAm_run` with a print method.
#' @examples
#' \dontrun{
#' run_summary("runs/study1")
#' }
#' @export
run_summary <- function(run_dir) {
  cfg <- read_config(run_dir)
  info <- read_run_info(run_dir)
  samples <- if (file.exists(run_file(run_dir, "samples.csv"))) read_samples(run_dir) else NULL
  stages <- info$stages %||% list()
  stage_df <- data.frame(
    stage = names(stages),
    status = vapply(stages, function(s) s$status %||% NA_character_, character(1)),
    seconds = vapply(stages, function(s) as.numeric(s$seconds %||% NA_real_), numeric(1)),
    stringsAsFactors = FALSE)
  rownames(stage_df) <- NULL
  failed <- list.files(run_dir, pattern = "^FAILED_.*\\.txt$")
  structure(list(
    run_dir = run_dir,
    run_name = cfg$run_name,
    array = cfg$array,
    input_type = cfg$input_type,
    prepDNAm_version = info$prepDNAm_version %||% cfg$prepDNAm_version,
    n_samples = if (is.null(samples)) NA_integer_ else nrow(samples),
    n_batches = if (is.null(samples)) NA_integer_ else length(unique(samples$batch)),
    n_samples_retained = info$n_samples_retained,
    n_probes_retained = info$n_probes_retained,
    n_probes_removed = info$n_probes_removed,
    stages = stage_df,
    failed = failed,
    complete = file.exists(run_file(run_dir, "beta_merged.rds")) && length(failed) == 0
  ), class = "prepDNAm_run")
}

#' @export
print.prepDNAm_run <- function(x, ...) {
  cat(sprintf("<prepDNAm run> %s  [%s, %s, prepDNAm %s]\n", x$run_name, x$array, x$input_type,
              x$prepDNAm_version))
  cat(sprintf("  run_dir:  %s\n", x$run_dir))
  cat(sprintf("  samples:  %d in %d batch(es)%s\n", x$n_samples, x$n_batches,
              if (!is.null(x$n_samples_retained)) sprintf(", %d retained", x$n_samples_retained) else ""))
  if (!is.null(x$n_probes_retained)) {
    cat(sprintf("  probes:   %d retained, %d removed\n", x$n_probes_retained, x$n_probes_removed))
  }
  cat(sprintf("  status:   %s\n", if (x$complete) "complete" else if (length(x$failed) > 0) {
    paste("FAILED -", paste(x$failed, collapse = ", "))
  } else "incomplete"))
  if (nrow(x$stages) > 0) {
    cat("  stages:\n")
    cat(sprintf("    %-32s %-8s %8s\n", x$stages$stage, x$stages$status,
                ifelse(is.na(x$stages$seconds), "", sprintf("%.1f s", x$stages$seconds))), sep = "")
  }
  invisible(x)
}

`%||%` <- function(a, b) if (is.null(a)) b else a
