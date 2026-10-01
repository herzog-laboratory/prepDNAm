#' Estimate the signal-to-noise ratio (rho) of each sample
#'
#' Estimates a per-sample signal-to-noise ratio from the raw beta distribution
#' of type II probes. The two main modes of the distribution (unmethylated and
#' methylated peaks) are located and compared with the modes expected from the
#' model `beta_obs = (rho * beta_signal + beta_e) / (1 + rho)`, where
#' `beta_signal` is `beta_s1` for the unmethylated and `beta_s2` for the
#' methylated peak and `beta_e` is the beta of pure noise. `rho` is the value
#' in `[0, 100]` (step 0.01) that minimises the squared distance between
#' expected and observed modes.
#'
#' @param x A `minfi` `MethylSet` (or `GenomicMethylSet`) of raw, un-normalised
#'   data, or a numeric matrix of raw beta values (probes in rows, samples in
#'   columns). For a matrix, `probe_type` is required.
#' @param probe_type Character vector, one entry per row of `x`, giving the
#'   Infinium design type (`"I"` or `"II"`). Ignored (taken from the manifest)
#'   when `x` is a `MethylSet`.
#' @param beta_e Beta value of pure noise.
#' @param beta_s1,beta_s2 Beta values of the unmethylated and methylated peaks
#'   of a noise-free sample.
#' @param max_na Samples with a larger fraction of missing type II betas are
#'   skipped and get `NA` modes and `rho`.
#' @return A data frame with one row per sample (rownames are sample names) and
#'   columns `mode0`, `mode1` (location of the two main modes of the type II
#'   beta density) and `rho`. `rho` is `NA` if fewer than two modes are found.
#' @examples
#' \dontrun{
#' rhoEstimator(minfiData::MsetEx)
#' }
#' @export
rhoEstimator <- function(x,
                         probe_type = NULL,
                         beta_e = 0.3309135,
                         beta_s1 = 0.019570,
                         beta_s2 = 0.970650,
                         max_na = 0.1) {
  if (methods::is(x, "MethylSet") || methods::is(x, "GenomicMethylSet")) {
    probe_type <- mset_probe_types(x)
    x <- minfi::getBeta(x)
  } else if (!is.matrix(x) || !is.numeric(x)) {
    rlang::abort("`x` must be a MethylSet or a numeric matrix of raw beta values.")
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

  beta_II <- x[probe_type == "II", , drop = FALSE]
  if (nrow(beta_II) == 0) {
    rlang::abort("`x` has no type II probes; rho cannot be estimated.")
  }

  # type II beta density on a fixed grid, then the two highest local maxima
  grid <- 512L
  dens <- matrix(NA_real_, nrow = grid, ncol = ncol(beta_II))
  for (j in seq_len(ncol(beta_II))) {
    if (mean(is.na(beta_II[, j])) <= max_na) {
      dens[, j] <- stats::density(beta_II[, j], na.rm = TRUE, bw = 0.01,
                                  n = grid, from = 0, to = 1)$y
    }
  }
  dx <- seq(0, 1, length.out = grid)
  modes <- t(apply(dens, 2, two_main_modes, dx = dx))
  dimnames(modes) <- list(colnames(x), c("mode0", "mode1"))

  # expected mode locations as a function of rho; pick the closest rho
  rho_grid <- seq(0, 100, by = 0.01)
  b0 <- (rho_grid * beta_s1 + beta_e) / (1 + rho_grid)
  b1 <- (rho_grid * beta_s2 + beta_e) / (1 + rho_grid)
  rho <- vapply(seq_len(nrow(modes)), function(i) {
    if (anyNA(modes[i, ])) return(NA_real_)
    L <- (b0 - modes[i, "mode0"])^2 + (b1 - modes[i, "mode1"])^2
    rho_grid[which.min(L)]
  }, numeric(1))

  data.frame(modes, rho = rho)
}

# Location of the two highest local maxima of a density evaluated on `dx`;
# c(NA, NA) if fewer than two are found. Peaks below 1% of the maximum are
# ignored: the FFT-based density oscillates around zero in the tails, which
# would otherwise count as modes.
two_main_modes <- function(y, dx) {
  if (anyNA(y)) return(c(NA_real_, NA_real_))
  dy <- diff(y)
  d2y <- diff(dy)
  turning <- which(diff(sign(dy)) != 0)
  peaks <- turning[d2y[turning] < 0]
  peaks <- peaks[y[peaks] > 0.01 * max(y)]
  if (length(peaks) < 2) return(c(NA_real_, NA_real_))
  top2 <- peaks[order(y[peaks], decreasing = TRUE)[1:2]]
  sort(dx[top2])
}
