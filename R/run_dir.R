# Run-directory infrastructure shared by all stages (DESIGN.md §2, §7):
# config.yml / samples.csv / run_info.json / log.txt / FAILED_<stage>.txt,
# and run_stage(), which wraps a stage body with idempotency, logging,
# timing and failure markers.

LOG_NS <- "prepDNAm"

# ---- paths -------------------------------------------------------------------

run_file <- function(run_dir, ...) file.path(run_dir, ...)

#' Path of a batch directory inside a run directory
#'
#' @param run_dir Run directory.
#' @param batch Batch label as in `samples.csv`.
#' @return `run_dir/batches/<batch>`.
#' @keywords internal
#' @export
batch_dir <- function(run_dir, batch) {
  file.path(run_dir, "batches", batch)
}

# ---- config.yml --------------------------------------------------------------

write_config <- function(config, path) {
  x <- unclass(config)
  yaml::write_yaml(x, path)
  invisible(path)
}

#' Read the configuration and sample sheet of a run
#'
#' @param run_dir Run directory written by [qc_init()].
#' @return `read_config()`: the `prepDNAm_config`; `read_samples()`: the
#'   sample sheet data frame (`basename`, `batch`, `exclude`, `path`, pheno
#'   columns).
#' @export
read_config <- function(run_dir) {
  path <- run_file(run_dir, "config.yml")
  if (!file.exists(path)) {
    rlang::abort(sprintf("%s is not a run directory (no config.yml). Run qc_init() first.", run_dir))
  }
  x <- yaml::read_yaml(path)
  # yaml turns integers into integers; keep the types qc_config() produces
  x$cores <- as.integer(x$cores)
  x$seed <- as.integer(x$seed)
  structure(x, class = "prepDNAm_config")
}

#' @rdname read_config
#' @export
read_samples <- function(run_dir) {
  path <- run_file(run_dir, "samples.csv")
  if (!file.exists(path)) {
    rlang::abort(sprintf("%s has no samples.csv. Run qc_init() first.", run_dir))
  }
  s <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE,
                       colClasses = c(basename = "character", batch = "character",
                                      exclude = "logical", path = "character"))
  s
}

# Fields that may legitimately differ between two configs of the same run.
# Not compared when re-using a run directory: bookkeeping, and options that
# only add outputs after the merge (switching them on needs no re-run).
config_volatile <- c("created", "prepDNAm_version", "estimates", "export_csv", "export_digits")

# Names of fields whose values differ (ignoring volatile ones).
config_diff <- function(a, b) {
  a <- unclass(a); b <- unclass(b)
  fields <- setdiff(union(names(a), names(b)), config_volatile)
  fields[!vapply(fields, function(f) {
    x <- a[[f]]; y <- b[[f]]
    # yaml round-trips: NULL <-> missing, integer <-> double
    if (is.null(x) && is.null(y)) return(TRUE)
    if (is.numeric(x) && is.numeric(y)) return(isTRUE(all.equal(as.numeric(x), as.numeric(y))))
    identical(x, y)
  }, logical(1))]
}

# ---- run_info.json -----------------------------------------------------------

read_run_info <- function(run_dir) {
  path <- run_file(run_dir, "run_info.json")
  if (!file.exists(path)) return(list())
  jsonlite::read_json(path)
}

write_run_info <- function(run_dir, info) {
  jsonlite::write_json(info, run_file(run_dir, "run_info.json"),
                       auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null")
  invisible(info)
}

# Merge `...` into the record of one stage in run_info.json.
update_stage_info <- function(run_dir, stage, ...) {
  info <- read_run_info(run_dir)
  if (is.null(info$stages)) info$stages <- list()
  cur <- info$stages[[stage]]
  if (is.null(cur)) cur <- list()
  new <- list(...)
  cur[names(new)] <- new
  info$stages[[stage]] <- cur
  write_run_info(run_dir, info)
}

now_iso <- function() format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")

# ---- logging -----------------------------------------------------------------

# Route this package's logger namespace to run_dir/log.txt, and to the console
# unless option `prepDNAm.quiet` is TRUE.
setup_logging <- function(run_dir) {
  logfile <- run_file(run_dir, "log.txt")
  appender <- if (isTRUE(getOption("prepDNAm.quiet", FALSE))) {
    logger::appender_file(logfile)
  } else {
    logger::appender_tee(logfile)
  }
  logger::log_threshold(logger::INFO, namespace = LOG_NS)
  logger::log_formatter(logger::formatter_paste, namespace = LOG_NS)
  logger::log_layout(logger::layout_glue_generator(
    format = "{format(time, \"%Y-%m-%d %H:%M:%S\")} [{level}] {msg}"), namespace = LOG_NS)
  logger::log_appender(appender, namespace = LOG_NS)
  invisible(logfile)
}

log_info <- function(...) logger::log_info(..., namespace = LOG_NS)
log_warn <- function(...) logger::log_warn(..., namespace = LOG_NS)
log_error <- function(...) logger::log_error(..., namespace = LOG_NS)

# ---- stage runner ------------------------------------------------------------

# Run one stage body with the shared bookkeeping:
#  * skip if every path in `outputs` exists (unless `force`)
#  * remove a stale FAILED_<stage>.txt, log start/end and record timings
#    in run_info.json
#  * set the seed
#  * on error: write FAILED_<stage>.txt (message + traceback), mark the stage
#    failed and re-throw with the stage name
# Returns TRUE if the body ran, FALSE if skipped (invisibly).
run_stage <- function(run_dir, stage, outputs, force, seed, body) {
  setup_logging(run_dir)
  failed_file <- run_file(run_dir, sprintf("FAILED_%s.txt", stage))
  if (!force && length(outputs) > 0 && all(file.exists(outputs))) {
    log_info(sprintf("%s: outputs exist, skipping (force = TRUE to redo)", stage))
    return(invisible(FALSE))
  }
  if (file.exists(failed_file)) file.remove(failed_file)
  log_info(sprintf("%s: start", stage))
  update_stage_info(run_dir, stage, status = "running", started = now_iso(), finished = NULL)
  if (!is.null(seed)) set.seed(seed)

  t0 <- Sys.time()
  trace <- NULL
  tryCatch(
    withCallingHandlers(body(), error = function(e) {
      trace <<- rlang::trace_back()
    }),
    error = function(e) {
      msg <- conditionMessage(e)
      writeLines(c(
        sprintf("Stage %s failed at %s", stage, now_iso()),
        "", msg, "", "Traceback:",
        if (!is.null(trace)) format(trace) else "<not available>"),
        failed_file)
      update_stage_info(run_dir, stage, status = "failed", finished = now_iso(), error = msg)
      log_error(sprintf("%s: FAILED - %s", stage, msg))
      rlang::abort(c(
        sprintf("Stage %s failed: %s", stage, msg),
        i = sprintf("Details in %s. Fix the problem and re-run to resume.", failed_file)),
        parent = e, class = "prepDNAm_stage_error")
    })
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  update_stage_info(run_dir, stage, status = "ok", finished = now_iso(),
                    seconds = round(secs, 1))
  log_info(sprintf("%s: done in %.1f s", stage, secs))
  invisible(TRUE)
}
