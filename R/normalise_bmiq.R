#' BMIQ normalisation of type I / type II probe bias
#'
#' Beta-mixture quantile normalisation (BMIQ, Teschendorff et al. 2013) adjusts
#' the beta distribution of Infinium type II probes to match that of type I
#' probes, sample by sample. This is a thin, array-agnostic wrapper around the
#' BMIQ implementation in `ChAMP`: the probe design vector is taken from the
#' array manifest (matched by probe name; or supplied explicitly), so it
#' works for 450k, EPIC, EPIC v2 and mouse arrays alike, whereas
#' `ChAMP::champ.norm()` only knows 450k and EPIC v1.
#'
#' Zeros and ones are replaced by `1e-6` and `1 - 1e-6` before fitting, as BMIQ
#' requires values in the open interval (0, 1). Missing values are excluded
#' from the fit for that sample and returned as `NA`.
#'
#' `champ.BMIQ()` subsamples probes for the mixture fit using an internal fixed
#' seed. The RNG kind is pinned to R's default (Mersenne-Twister) while it
#' runs, so results are deterministic and identical for any value of `cores`
#' and to those of `ChAMP::champ.norm()` run interactively; the caller's RNG
#' state is restored afterwards.
#'
#' @param x A `minfi` `MethylSet` (or `GenomicMethylSet`), typically after
#'   [minfi::preprocessNoob()], or a numeric matrix of beta values (probes in
#'   rows, samples in columns). For a matrix, `probe_type` is required.
#' @param probe_type Character vector, one entry per row of `x`, giving the
#'   Infinium design type (`"I"` or `"II"`). Ignored (taken from the manifest)
#'   when `x` is a `MethylSet`.
#' @param cores Number of samples to normalise in parallel. Capped at
#'   [parallelly::availableCores()].
#' @param nfit Number of type I and type II probes subsampled to fit the beta
#'   mixture models. `ChAMP`'s default.
#' @return A numeric matrix of normalised beta values with the same dimensions
#'   and dimnames as the input betas.
#' @references
#' The method and original implementation: Teschendorff AE et al. (2013). A
#' beta-mixture quantile normalization method for correcting probe design bias
#' in Illumina Infinium 450 k DNA methylation data. *Bioinformatics* 29(2),
#' 189-196. \doi{10.1093/bioinformatics/bts680}
#'
#' The implementation called here (`ChAMP`): Morris TJ et al. (2014). ChAMP:
#' 450k Chip Analysis Methylation Pipeline. *Bioinformatics* 30(3), 428-430.
#' \doi{10.1093/bioinformatics/btt684}; Tian Y et al. (2017). ChAMP: updated
#' methylation analysis pipeline for Illumina BeadChips. *Bioinformatics*
#' 33(24), 3982-3984. \doi{10.1093/bioinformatics/btx513}
#' @examples
#' \dontrun{
#' mset <- minfi::preprocessNoob(minfiData::RGsetEx)
#' beta_norm <- normalise_bmiq(mset, cores = 2)
#' }
#' @export
normalise_bmiq <- function(x, probe_type = NULL, cores = 1L, nfit = 10000L) {
  # ChAMP's namespace load emits import-replacement warnings of its own
  suppressWarnings(requireNamespace("ChAMP", quietly = TRUE))
  rlang::check_installed("ChAMP", reason = "for BMIQ normalisation")

  if (methods::is(x, "MethylSet") || methods::is(x, "GenomicMethylSet")) {
    probe_type <- mset_probe_types(x)
    x <- minfi::getBeta(x)
  } else if (!is.matrix(x) || !is.numeric(x)) {
    rlang::abort("`x` must be a MethylSet or a numeric matrix of beta values.")
  }
  if (is.null(probe_type)) {
    rlang::abort("`probe_type` is required when `x` is a matrix.")
  }
  if (length(probe_type) != nrow(x)) {
    rlang::abort("`probe_type` must have one entry per row of `x`.")
  }
  if (!all(probe_type %in% c("I", "II"))) {
    rlang::abort("`probe_type` must contain only \"I\" and \"II\".")
  }
  if (is.null(colnames(x))) {
    colnames(x) <- seq_len(ncol(x))
  }
  design <- ifelse(probe_type == "I", 1, 2)

  n_zero <- sum(x == 0, na.rm = TRUE)
  n_one <- sum(x == 1, na.rm = TRUE)
  if (n_zero + n_one > 0) {
    x[x == 0] <- 1e-6
    x[x == 1] <- 1 - 1e-6
    rlang::inform(sprintf(
      "BMIQ: replaced %d zeros and %d ones with 1e-6 / 1 - 1e-6.", n_zero, n_one))
  }

  cores <- min(cores, parallelly::availableCores())
  bp <- if (cores > 1 && .Platform$OS.type != "windows") {
    BiocParallel::MulticoreParam(cores)
  } else if (cores > 1) {
    BiocParallel::SnowParam(cores)
  } else {
    BiocParallel::SerialParam()
  }

  cols <- BiocParallel::bplapply(colnames(x), function(id) {
    bmiq_one(x[, id], design, id = id, nfit = nfit)
  }, BPPARAM = bp)

  out <- do.call(cbind, cols)
  dimnames(out) <- dimnames(x)
  out
}

# BMIQ for one sample; NAs are left out of the fit and returned as NA.
# champ.BMIQ() print()s progress for every sample, which is swallowed here, and
# emits bare warnings from max()/min() when a mixture component ends up empty;
# those are re-issued once, naming the sample.
bmiq_one <- function(beta, design, id, nfit) {
  # champ.BMIQ() calls set.seed() before subsampling; the result still depends
  # on the RNG kind, which BiocParallel workers change to L'Ecuyer-CMRG.
  withr::local_preserve_seed()
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  ok <- !is.na(beta)
  # champ.BMIQ() samples `nfit` probes of each type without replacement
  nfit <- min(nfit, sum(design[ok] == 1), sum(design[ok] == 2))
  fit_warnings <- character()
  res <- tryCatch(
    withCallingHandlers(
      {
        utils::capture.output(
          res <- ChAMP:::champ.BMIQ(beta[ok], design[ok], nfit = nfit,
                                    plots = FALSE, sampleID = id))
        res
      },
      warning = function(w) {
        fit_warnings <<- c(fit_warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) {
      rlang::abort(sprintf(
        "BMIQ failed for sample %s (poor quality samples with no bimodal beta distribution cannot be normalised).",
        id), parent = e)
    })
  if (length(fit_warnings) > 0) {
    rlang::warn(c(
      sprintf("BMIQ: degenerate mixture fit for sample %s; check its beta distribution.", id),
      i = paste(unique(fit_warnings), collapse = "; ")))
  }
  out <- rep(NA_real_, length(beta))
  out[ok] <- res$nbeta
  out
}
