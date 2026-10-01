#' Read one batch of raw data and extract QC metrics
#'
#' Second stage. For each batch, reads the raw input and writes a bundle of
#' fixed-name files into `run_dir/batches/<batch>/`. Items that an input type
#' cannot provide are simply absent; later stages and report pages check for
#' presence.
#'
#' | File | Content | idat | sigmatrix | beta |
#' |---|---|---|---|---|
#' | `mset.rds` | `minfi` `MethylSet`; Noob background/dye corrected for IDATs, as supplied otherwise | yes | yes | -- |
#' | `beta_input.rds` | beta matrix as supplied, restricted to manifest probes | -- | -- | yes |
#' | `detP.rds` | detection p-values (probes x samples) | yes | if present | -- |
#' | `qc.rds` | per-sample metrics: `mMed`, `uMed` (median log2 intensities, raw), `snp_outlier` (`ewastools`) | yes | yes | -- |
#' | `rho.rds` | signal-to-noise estimate ([rhoEstimator()], raw betas) | yes | yes | -- |
#' | `sex_raw.rds` | median log2 copy number of chrX / chrY probes per sample (`xMed`, `yMed`); sex is called on all samples at merge | if `save_sex` | if `save_sex` | -- |
#' | `rs.rds` | SNP probe betas | if `save_rs` | if `save_rs` and rs rows present | -- |
#' | `ctrls.rds` | control probe intensities (`info`, `red`, `green`) | yes | -- | -- |
#' | `densities.rds` | per-sample beta densities on a fixed grid (`x`, `raw`, `noob`) for the report | yes | raw only | raw only |
#' | `shiny.rds` | `shinyMethyl::shinySummarize()` of the batch | if `save_shiny` | -- | -- |
#'
#' Samples flagged `exclude` in the sample sheet are read and measured like
#' the others so that the report can show why they were dropped; they are
#' removed by `qc_filter_samples()`.
#'
#' @param run_dir Run directory initialised by [qc_init()].
#' @param batch Batch label(s) to read; `NULL` for all batches in turn.
#' @param force Recompute even if the batch bundle exists.
#' @return Named logical vector, one element per batch: `TRUE` if read,
#'   `FALSE` if skipped; invisibly.
#' @examples
#' \dontrun{
#' qc_read("runs/study1")                 # all batches
#' qc_read("runs/study1", batch = "plate2", force = TRUE)
#' }
#' @export
qc_read <- function(run_dir, batch = NULL, force = FALSE) {
  cfg <- read_config(run_dir)
  samples <- read_samples(run_dir)
  check_flag(force, "force")
  batch <- check_batches(batch, samples)
  info <- array_info(cfg$array)
  rlang::check_installed(info$manifest_pkg, reason = sprintf("to read %s data", info$label))
  if (cfg$save_sex) {
    rlang::check_installed(info$anno_pkg, reason = "for sex estimation (`save_sex = TRUE`)")
  }

  ran <- vapply(batch, function(b) {
    bdir <- batch_dir(run_dir, b)
    bs <- samples[samples$batch == b, , drop = FALSE]
    outputs <- file.path(bdir, batch_outputs(cfg))
    run_stage(run_dir, paste0("qc_read:", b), outputs, force, seed = cfg$seed, body = function() {
      dir.create(bdir, recursive = TRUE, showWarnings = FALSE)
      log_info(sprintf("qc_read:%s: %d samples, input type %s", b, nrow(bs), cfg$input_type))
      switch(cfg$input_type,
        idat = read_batch_idat(cfg, info, bs, bdir, b),
        sigmatrix = read_batch_sigmatrix(cfg, info, bs, bdir, b),
        beta = read_batch_beta(cfg, info, bs, bdir, b))
    })
  }, logical(1))
  invisible(ran)
}

# Files a batch bundle must contain for the stage to count as done.
batch_outputs <- function(cfg) {
  switch(cfg$input_type,
    idat = c("mset.rds", "detP.rds", "qc.rds", "rho.rds", "ctrls.rds", "densities.rds",
             if (cfg$save_sex) "sex_raw.rds", if (cfg$save_rs) "rs.rds",
             if (isTRUE(cfg$save_shiny)) "shiny.rds"),
    sigmatrix = c("mset.rds", "qc.rds", "rho.rds", "densities.rds",
                  if (cfg$save_sex) "sex_raw.rds"),
    beta = c("beta_input.rds", "densities.rds"))
}

# ---- idat --------------------------------------------------------------------

read_batch_idat <- function(cfg, info, bs, bdir, b) {
  out <- function(f) file.path(bdir, f)
  log_info(sprintf("qc_read:%s: reading IDATs", b))
  rgset <- minfi::read.metharray(bs$path, force = TRUE)
  colnames(rgset) <- bs$basename
  guessed <- unname(minfi::annotation(rgset)["array"])
  if (guessed %in% arrays$minfi_array && guessed != info$minfi_array) {
    rlang::abort(c(
      sprintf("The IDATs look like %s but `array` is %s.", guessed, cfg$array),
      i = "Check the `array` argument of qc_config()."))
  }
  rgset@annotation <- c(array = info$minfi_array, annotation = info$minfi_annotation)

  # control probes
  ctrl_info <- minfi::getProbeInfo(rgset, type = "Control")
  ctrl_info <- ctrl_info[ctrl_info$Address %in% rownames(rgset) & !is.na(ctrl_info$Color), ]
  saveRDS(list(
    info = as.data.frame(ctrl_info),
    red = minfi::getRed(rgset)[ctrl_info$Address, , drop = FALSE],
    green = minfi::getGreen(rgset)[ctrl_info$Address, , drop = FALSE]), out("ctrls.rds"))

  # detection p-values
  log_info(sprintf("qc_read:%s: detection p-values", b))
  detP <- minfi::detectionP(rgset, type = "m+u")
  saveRDS(detP, out("detP.rds"))

  # raw intensities: median intensity QC, rho, sex
  mset_raw <- minfi::preprocessRaw(rgset)
  qc <- as.data.frame(minfi::getQC(mset_raw))
  qc <- data.frame(basename = rownames(qc), mMed = qc$mMed, uMed = qc$uMed,
                   stringsAsFactors = FALSE)
  qc$snp_outlier <- if (info$has_snp_probes) {
    snp_beta <- minfi::getSnpBeta(rgset)
    if (cfg$save_rs) saveRDS(snp_beta, out("rs.rds"))
    snp_outlier_metric(snp_beta, b)
  } else NA_real_
  saveRDS(qc, out("qc.rds"))

  log_info(sprintf("qc_read:%s: signal-to-noise (rho)", b))
  saveRDS(rhoEstimator(mset_raw), out("rho.rds"))

  if (cfg$save_sex) {
    saveRDS(sex_medians(mset_raw), out("sex_raw.rds"))
  }

  if (isTRUE(cfg$save_shiny)) {
    rlang::check_installed("shinyMethyl", reason = "for `save_shiny = TRUE`")
    saveRDS(shinyMethyl::shinySummarize(rgset), out("shiny.rds"))
  }

  # Noob background + dye-bias correction (single-sample: independent of the
  # other samples in the batch, so failed samples can be dropped afterwards)
  log_info(sprintf("qc_read:%s: Noob", b))
  mset <- minfi::preprocessNoob(rgset, dyeCorr = TRUE, dyeMethod = "single")
  saveRDS(mset, out("mset.rds"))

  saveRDS(list(
    x = density_grid(),
    raw = beta_densities(minfi::getBeta(mset_raw)),
    noob = beta_densities(minfi::getBeta(mset))), out("densities.rds"))
  invisible(NULL)
}

# ---- sigmatrix ---------------------------------------------------------------

read_batch_sigmatrix <- function(cfg, info, bs, bdir, b) {
  out <- function(f) file.path(bdir, f)
  log_info(sprintf("qc_read:%s: reading signal matrix", b))
  sig <- read_sigmatrix(cfg, info)
  keep <- match(bs$basename, colnames(sig$meth))
  if (anyNA(keep)) {
    rlang::abort(sprintf("Samples not found in the signal matrix: %s",
                         paste(bs$basename[is.na(keep)], collapse = ", ")))
  }
  mset <- minfi::MethylSet(
    Meth = sig$meth[, keep, drop = FALSE], Unmeth = sig$unmeth[, keep, drop = FALSE],
    annotation = c(array = info$minfi_array, annotation = info$minfi_annotation),
    preprocessMethod = c(rg.norm = "upstream (signal matrix)",
                         minfi = as.character(utils::packageVersion("minfi")),
                         manifest = as.character(utils::packageVersion(info$manifest_pkg))))
  saveRDS(mset, out("mset.rds"))

  if (!is.null(sig$detP)) {
    saveRDS(sig$detP[, keep, drop = FALSE], out("detP.rds"))
  } else {
    log_warn(sprintf("qc_read:%s: no detection p-value columns; probe/sample failure filters will be skipped", b))
  }

  qc <- as.data.frame(minfi::getQC(mset))
  qc <- data.frame(basename = rownames(qc), mMed = qc$mMed, uMed = qc$uMed,
                   stringsAsFactors = FALSE)
  qc$snp_outlier <- NA_real_
  if (!is.null(sig$rs)) {
    rs <- sig$rs[, keep, drop = FALSE]
    if (cfg$save_rs) saveRDS(rs, out("rs.rds"))
    qc$snp_outlier <- snp_outlier_metric(rs, b)
  }
  saveRDS(qc, out("qc.rds"))

  log_info(sprintf("qc_read:%s: signal-to-noise (rho)", b))
  saveRDS(rhoEstimator(mset), out("rho.rds"))

  if (cfg$save_sex) {
    saveRDS(sex_medians(mset), out("sex_raw.rds"))
  }

  saveRDS(list(x = density_grid(), raw = beta_densities(minfi::getBeta(mset)), noob = NULL),
          out("densities.rds"))
  invisible(NULL)
}

# Read the whole signal-intensity file once per session (memoised by path and
# mtime): list(meth, unmeth, detP or NULL, rs or NULL) with manifest CpG
# probes in rows and parsed sample names in columns.
read_sigmatrix <- function(cfg, info) {
  cached(cfg$input, function() {
    header <- read_header(cfg$input)
    col_for <- function(pattern) {
      cols <- header[grepl(pattern, header, fixed = TRUE)]
      stats::setNames(cols, header_samples(header, pattern))
    }
    u_cols <- col_for(cfg$unmeth_pattern)
    m_cols <- col_for(cfg$meth_pattern)
    p_cols <- if (!is.null(cfg$detp_pattern)) col_for(cfg$detp_pattern) else character()
    samples <- names(u_cols)

    dt <- data.table::fread(cfg$input, select = unname(c(header[1], u_cols, m_cols, p_cols)),
                            data.table = FALSE, showProgress = FALSE)
    ids <- as.character(dt[[1]])
    to_matrix <- function(cols) {
      m <- as.matrix(dt[, cols[samples], drop = FALSE])
      storage.mode(m) <- "double"
      dimnames(m) <- list(ids, samples)
      m
    }
    unmeth <- to_matrix(u_cols)
    meth <- to_matrix(m_cols)
    detP <- if (length(p_cols) > 0) to_matrix(p_cols) else NULL

    # rs probes (if supplied) are kept apart; everything else is restricted to
    # the manifest's CpG probes
    is_rs <- startsWith(ids, "rs")
    rs <- if (any(is_rs)) meth[is_rs, , drop = FALSE] / (meth[is_rs, , drop = FALSE] + unmeth[is_rs, , drop = FALSE] + 100) else NULL
    cpg <- manifest_probes(info)
    keep <- ids %in% cpg
    n_drop <- sum(!keep & !is_rs)
    if (n_drop > 0) {
      log_warn(sprintf("qc_read: %d rows of the signal matrix are not %s CpG probes and were dropped", n_drop, info$label))
    }
    if (sum(keep) < 0.5 * length(cpg)) {
      log_warn(sprintf("qc_read: only %d of %d %s probes present in the signal matrix", sum(keep), length(cpg), info$label))
    }
    list(meth = meth[keep, , drop = FALSE], unmeth = unmeth[keep, , drop = FALSE],
         detP = if (!is.null(detP)) detP[keep, , drop = FALSE] else NULL, rs = rs)
  })
}

# ---- beta --------------------------------------------------------------------

read_batch_beta <- function(cfg, info, bs, bdir, b) {
  out <- function(f) file.path(bdir, f)
  log_info(sprintf("qc_read:%s: reading beta matrix", b))
  beta <- read_beta_matrix(cfg, info)
  keep <- match(bs$basename, colnames(beta))
  if (anyNA(keep)) {
    rlang::abort(sprintf("Samples not found in the beta matrix: %s",
                         paste(bs$basename[is.na(keep)], collapse = ", ")))
  }
  beta <- beta[, keep, drop = FALSE]
  saveRDS(beta, out("beta_input.rds"))
  saveRDS(list(x = density_grid(), raw = beta_densities(beta), noob = NULL), out("densities.rds"))
  invisible(NULL)
}

# Read a beta matrix file once per session: numeric matrix restricted to
# manifest CpG probes, samples in columns.
read_beta_matrix <- function(cfg, info) {
  cached(cfg$input, function() {
    ext <- tolower(tools::file_ext(sub("\\.gz$", "", cfg$input)))
    beta <- if (ext == "rds") {
      readRDS(cfg$input)
    } else if (ext %in% c("rdata", "rda")) {
      load_single_object(cfg$input)
    } else {
      dt <- data.table::fread(cfg$input, data.table = FALSE, showProgress = FALSE)
      m <- as.matrix(dt[, -1, drop = FALSE])
      rownames(m) <- as.character(dt[[1]])
      m
    }
    if (is.data.frame(beta)) beta <- as.matrix(beta)
    if (!is.matrix(beta) || is.null(rownames(beta)) || is.null(colnames(beta))) {
      rlang::abort("The beta matrix must have probe IDs as rownames (or first column) and samples as column names.")
    }
    storage.mode(beta) <- "double"
    rng <- range(beta, na.rm = TRUE)
    if (rng[1] < 0 || rng[2] > 1) {
      rlang::abort(c(
        sprintf("Values in the beta matrix range from %.3g to %.3g; betas must lie in [0, 1].", rng[1], rng[2]),
        i = "M-values are not supported; convert with 2^M / (2^M + 1) first."))
    }
    cpg <- manifest_probes(info)
    keep <- rownames(beta) %in% cpg
    if (sum(!keep) > 0) {
      log_warn(sprintf("qc_read: %d rows of the beta matrix are not %s CpG probes and were dropped", sum(!keep), info$label))
    }
    if (sum(keep) < 0.5 * length(cpg)) {
      log_warn(sprintf("qc_read: only %d of %d %s probes present in the beta matrix", sum(keep), length(cpg), info$label))
    }
    beta[keep, , drop = FALSE]
  })
}

# ---- shared helpers ----------------------------------------------------------

# CpG probe names (type I + II) of an array's manifest.
manifest_probes <- function(info) {
  cached(paste0("manifest:", info$array), function() {
    e <- new.env()
    utils::data(list = info$manifest_pkg, package = info$manifest_pkg, envir = e)
    m <- get(info$manifest_pkg, e)
    c(minfi::getProbeInfo(m, "I")$Name, minfi::getProbeInfo(m, "II")$Name)
  })
}

# Median log2 copy number of chrX and chrY probes per sample (what
# minfi::getSex() computes before clustering).
sex_medians <- function(mset) {
  gm <- minfi::mapToGenome(mset)
  cn <- minfi::getCN(gm)
  chr <- as.character(GenomicRanges::seqnames(gm))
  data.frame(
    basename = colnames(cn),
    xMed = matrixStats::colMedians(cn[chr == "chrX", , drop = FALSE], na.rm = TRUE),
    yMed = matrixStats::colMedians(cn[chr == "chrY", , drop = FALSE], na.rm = TRUE),
    stringsAsFactors = FALSE)
}

# ewastools SNP outlier metric per sample; NA if ewastools is not installed.
snp_outlier_metric <- function(snp_beta, b) {
  if (!requireNamespace("ewastools", quietly = TRUE)) {
    log_warn(sprintf("qc_read:%s: ewastools not installed; snp_outlier not computed", b))
    return(rep(NA_real_, ncol(snp_beta)))
  }
  genotypes <- ewastools::call_genotypes(snp_beta, learn = FALSE)
  unname(ewastools::snp_outliers(genotypes))
}

density_grid <- function(n = 512L) seq(-0.05, 1.05, length.out = n)

# Per-sample beta densities on density_grid(): matrix n_grid x n_samples.
beta_densities <- function(beta, bw = 0.02) {
  grid <- density_grid()
  d <- apply(beta, 2, function(x) {
    if (sum(!is.na(x)) < 2) return(rep(NA_real_, length(grid)))
    stats::density(x, bw = bw, from = grid[1], to = grid[length(grid)], n = length(grid),
                   na.rm = TRUE)$y
  })
  d <- matrix(d, nrow = length(grid), dimnames = list(NULL, colnames(beta)))
  d
}

# Session-level memoisation keyed by a path (invalidated when the file changes)
# or an arbitrary key.
.cache <- new.env(parent = emptyenv())
cached <- function(key, build) {
  stamp <- if (file.exists(key)) paste(normalizePath(key), file.mtime(key)) else key
  hit <- .cache[[key]]
  if (!is.null(hit) && identical(hit$stamp, stamp)) return(hit$value)
  value <- build()
  .cache[[key]] <- list(stamp = stamp, value = value)
  value
}
