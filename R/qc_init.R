#' Initialise a run directory
#'
#' First stage of a run. Creates `run_dir` and writes
#' \describe{
#'   \item{`config.yml`}{the validated configuration ([qc_config()]); a
#'     `pheno` data frame is written to `pheno.csv` and referenced from there}
#'   \item{`samples.csv`}{the sample sheet ([discover_samples()])}
#'   \item{`run_info.json`}{package and R version, platform, package versions,
#'     seed, checksums of the sample sheet and pheno file, and per-stage
#'     status/timings maintained by every later stage}
#'   \item{`log.txt`}{the run log (appended to by all stages)}
#' }
#'
#' Like every stage, `qc_init()` is idempotent: if the files exist and the
#' configuration matches the one on disk, nothing is rewritten and the run can
#' be resumed. A different configuration for an existing run directory is an
#' error unless `force = TRUE`, which rewrites the files (later stages then
#' need `force = TRUE` as well to recompute).
#'
#' @param config A [qc_config()] object.
#' @param force Rewrite the files even if they exist.
#' @return `TRUE` if the stage ran, `FALSE` if it was skipped because its
#'   outputs exist; invisibly.
#' @examples
#' \dontrun{
#' cfg <- qc_config(input = "idats/", run_dir = "runs/study1", array = "EPIC")
#' qc_init(cfg)
#' read_samples("runs/study1")
#' }
#' @export
qc_init <- function(config, force = FALSE) {
  if (!inherits(config, "prepDNAm_config")) {
    rlang::abort("`config` must be a `prepDNAm_config` object from `qc_config()`.")
  }
  check_flag(force, "force")
  run_dir <- config$run_dir
  dir.create(run_dir, showWarnings = FALSE, recursive = FALSE)
  if (!dir.exists(run_dir)) {
    rlang::abort(sprintf("Could not create run directory %s.", run_dir))
  }

  config_file <- run_file(run_dir, "config.yml")
  outputs <- run_file(run_dir, c("config.yml", "samples.csv", "run_info.json"))

  # A pheno data frame is materialised as pheno.csv so the run directory is
  # self-contained; the config on disk references that file.
  cfg <- config
  pheno_df <- NULL
  if (is.data.frame(cfg$pheno)) {
    pheno_df <- cfg$pheno
    cfg$pheno <- run_file(run_dir, "pheno.csv")
  }

  # A run directory belongs to one configuration - once qc_init has completed.
  # After a failed init (config.yml written, samples.csv not) any config may
  # be used, e.g. a corrected pheno file.
  if (file.exists(config_file) && file.exists(run_file(run_dir, "samples.csv"))) {
    on_disk <- read_config(run_dir)
    diff <- config_diff(cfg, on_disk)
    if (!is.null(pheno_df) && file.exists(cfg$pheno) && !"pheno" %in% diff) {
      tmp <- withr::local_tempfile(fileext = ".csv")
      utils::write.csv(pheno_df, tmp, row.names = FALSE)
      if (unname(tools::md5sum(tmp)) != unname(tools::md5sum(cfg$pheno))) diff <- c(diff, "pheno")
    }
    if (length(diff) > 0 && !force) {
      rlang::abort(c(
        sprintf("%s was initialised with a different configuration.", run_dir),
        x = sprintf("Differing fields: %s", paste(diff, collapse = ", ")),
        i = "Use a new `run_dir`, or `force = TRUE` to overwrite (and re-run later stages with `force = TRUE`)."))
    }
    if (!identical(on_disk$prepDNAm_version, config$prepDNAm_version)) {
      setup_logging(run_dir)
      log_warn(sprintf("qc_init: run directory was created with prepDNAm %s, now running %s",
                       on_disk$prepDNAm_version, config$prepDNAm_version))
    }
  }

  ran <- run_stage(run_dir, "qc_init", outputs, force, seed = NULL, body = function() {
    # written first so that a failed run still records what was attempted
    write_config(cfg, config_file)
    if (!is.null(pheno_df)) {
      utils::write.csv(pheno_df, cfg$pheno, row.names = FALSE)
      log_info("qc_init: pheno data frame written to pheno.csv")
    }
    samples <- discover_samples(cfg)
    log_info(sprintf("qc_init: %d samples in %d batch(es) (%s), %d excluded, input type %s, array %s",
                     nrow(samples), length(unique(samples$batch)), cfg$batch,
                     sum(samples$exclude), cfg$input_type, cfg$array))
    utils::write.csv(samples, run_file(run_dir, "samples.csv"), row.names = FALSE)

    info <- read_run_info(run_dir)
    info$run_name <- cfg$run_name
    info$prepDNAm_version <- cfg$prepDNAm_version
    info$r_version <- R.version.string
    info$platform <- R.version$platform
    info$host <- unname(Sys.info()[["nodename"]])
    info$user <- unname(Sys.info()[["user"]])
    info$seed <- cfg$seed
    info$cores <- cfg$cores
    info$input_type <- cfg$input_type
    info$array <- cfg$array
    info$n_samples <- nrow(samples)
    info$batches <- as.list(table(samples$batch))
    info$md5 <- list(
      samples = unname(tools::md5sum(run_file(run_dir, "samples.csv"))),
      pheno = if (is.null(cfg$pheno)) NULL else unname(tools::md5sum(cfg$pheno)))
    info$packages <- package_versions()
    write_run_info(run_dir, info)
  })
  invisible(ran)
}

# Versions of the packages loaded in this session (prepDNAm and everything
# it pulled in).
package_versions <- function() {
  pk <- tryCatch(sessioninfo::package_info("loaded"), error = function(e) NULL)
  versions <- if (is.null(pk)) list() else stats::setNames(as.list(pk$loadedversion), pk$package)
  versions$prepDNAm <- as.character(utils::packageVersion("prepDNAm"))
  versions[order(names(versions))]
}
