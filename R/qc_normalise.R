#' Normalise the retained samples of a batch
#'
#' Fourth stage. For each batch, takes the samples not listed in
#' `samples_removed.csv`, applies BMIQ ([normalise_bmiq()]) and writes
#' `batches/<batch>/beta_norm.rds` (probes x retained samples). The BMIQ
#' density curves are added to `densities.rds` as `bmiq`.
#'
#' For IDAT and signal-matrix input BMIQ always runs on the Noob / supplied
#' intensities. For beta input `normalise` decides: `TRUE`/`FALSE` force it,
#' `"auto"` runs BMIQ only if the type II peaks are visibly shifted relative
#' to the type I peaks (median absolute shift of the unmethylated or
#' methylated peak across samples above 0.05), i.e. the matrix does not look
#' normalised already. The decision is logged. When BMIQ is not run,
#' `beta_norm.rds` is the input matrix restricted to the retained samples.
#'
#' @inheritParams qc_read
#' @return Named logical vector, one element per batch: `TRUE` if run,
#'   `FALSE` if skipped; invisibly.
#' @examples
#' \dontrun{
#' qc_normalise("runs/study1")
#' beta <- readRDS("runs/study1/batches/plate1/beta_norm.rds")
#' }
#' @export
qc_normalise <- function(run_dir, batch = NULL, force = FALSE) {
  cfg <- read_config(run_dir)
  samples <- read_samples(run_dir)
  check_flag(force, "force")
  batch <- check_batches(batch, samples)
  info <- array_info(cfg$array)

  ran <- vapply(batch, function(b) {
    bdir <- batch_dir(run_dir, b)
    out <- file.path(bdir, "beta_norm.rds")
    run_stage(run_dir, paste0("qc_normalise:", b), out, force, seed = cfg$seed, body = function() {
      removed <- unique(read_samples_removed(run_dir, b)$basename)
      inp <- batch_betas(cfg, info, bdir, b)
      keep <- setdiff(colnames(inp$beta), removed)
      beta <- inp$beta[, keep, drop = FALSE]
      log_info(sprintf("qc_normalise:%s: %d of %d samples retained", b, length(keep), ncol(inp$beta)))

      run_bmiq <- if (cfg$input_type == "beta") decide_bmiq(cfg$normalise, beta, inp$probe_type, b) else TRUE
      if (run_bmiq && ncol(beta) > 0) {
        log_info(sprintf("qc_normalise:%s: BMIQ on %d samples (%d cores)", b, ncol(beta), cfg$cores))
        beta <- normalise_bmiq(beta, probe_type = inp$probe_type, cores = cfg$cores)
      } else if (run_bmiq) {
        log_warn(sprintf("qc_normalise:%s: no samples left to normalise", b))
      }
      saveRDS(beta, out)

      dens_file <- file.path(bdir, "densities.rds")
      dens <- if (file.exists(dens_file)) readRDS(dens_file) else list(x = density_grid(), raw = NULL, noob = NULL)
      dens$bmiq <- if (run_bmiq && ncol(beta) > 0) beta_densities(beta) else NULL
      saveRDS(dens, dens_file)
    })
  }, logical(1))
  invisible(ran)
}

# Beta matrix and probe types of a batch, from mset.rds or beta_input.rds.
batch_betas <- function(cfg, info, bdir, b) {
  if (cfg$input_type == "beta") {
    path <- file.path(bdir, "beta_input.rds")
    if (!file.exists(path)) rlang::abort(sprintf("Batch %s has no beta_input.rds. Run qc_read() first.", b))
    beta <- readRDS(path)
    types <- manifest_probe_types(info)
    list(beta = beta, probe_type = unname(types[rownames(beta)]))
  } else {
    path <- file.path(bdir, "mset.rds")
    if (!file.exists(path)) rlang::abort(sprintf("Batch %s has no mset.rds. Run qc_read() first.", b))
    mset <- readRDS(path)
    list(beta = minfi::getBeta(mset), probe_type = mset_probe_types(mset))
  }
}

# Named character vector probe -> "I"/"II" for an array's manifest.
manifest_probe_types <- function(info) manifest_types_for(info$manifest_pkg)

manifest_types_for <- function(manifest_pkg) {
  cached(paste0("manifest_types:", manifest_pkg), function() {
    rlang::check_installed(manifest_pkg)
    e <- new.env()
    utils::data(list = manifest_pkg, package = manifest_pkg, envir = e)
    m <- get(manifest_pkg, e)
    i <- minfi::getProbeInfo(m, "I")$Name
    ii <- minfi::getProbeInfo(m, "II")$Name
    stats::setNames(rep(c("I", "II"), c(length(i), length(ii))), c(i, ii))
  })
}

# Probe types ("I"/"II") of the rows of a MethylSet, from its array's
# manifest, matched by probe name. Not minfi::getProbeType(): that reads the
# annotation package, which for EPIC v2 (20a1) lacks 6,915 manifest probes and
# so returns fewer entries than the object has rows.
mset_probe_types <- function(mset) {
  types <- manifest_types_for(paste0(minfi::annotation(mset)[["array"]], "manifest"))
  out <- unname(types[rownames(mset)])
  if (anyNA(out)) {
    rlang::abort(sprintf("%d probe(s) of the data are not in the %s manifest.",
                         sum(is.na(out)), minfi::annotation(mset)[["array"]]))
  }
  out
}

# Whether to run BMIQ on beta input: TRUE/FALSE as configured, or for "auto"
# the peak-shift heuristic.
decide_bmiq <- function(normalise, beta, probe_type, b) {
  if (isTRUE(normalise) || isFALSE(normalise)) {
    log_info(sprintf("qc_normalise:%s: BMIQ %s (`normalise = %s`)", b,
                     if (normalise) "on" else "off", normalise))
    return(normalise)
  }
  if (ncol(beta) == 0) return(FALSE)
  shift <- type_ii_peak_shift(beta, probe_type)
  run <- is.finite(shift) && shift > 0.05
  log_info(sprintf(
    "qc_normalise:%s: `normalise = \"auto\"`: median type II peak shift %.3f -> BMIQ %s", b,
    shift, if (run) "on (input looks unnormalised)" else "off (input looks normalised)"))
  run
}

# Median over samples of the larger absolute difference between the type II
# and type I unmethylated / methylated peak locations.
type_ii_peak_shift <- function(beta, probe_type) {
  grid <- 512L
  dx <- seq(0, 1, length.out = grid)
  peaks <- function(x) {
    x <- x[!is.na(x)]
    if (length(x) < 100) return(c(NA_real_, NA_real_))
    two_main_modes(stats::density(x, bw = 0.01, n = grid, from = 0, to = 1)$y, dx)
  }
  shifts <- vapply(seq_len(ncol(beta)), function(j) {
    p1 <- peaks(beta[probe_type == "I", j])
    p2 <- peaks(beta[probe_type == "II", j])
    max(abs(p2 - p1))
  }, numeric(1))
  stats::median(shifts, na.rm = TRUE)
}
