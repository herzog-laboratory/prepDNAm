#' Merge batches into the final beta matrix
#'
#' Fifth stage. Combines the normalised batches, removes probes, masks and
#' imputes failed values, calls sex, and writes into `run_dir`:
#' \describe{
#'   \item{`beta_merged.rds`}{probes x retained samples (sample-sheet order)}
#'   \item{`beta_merged_compatible.rds`}{EPIC v2 with `compatible_subset`:
#'     the same restricted to EPIC v1/v2 compatible probes with v1 IDs
#'     ([subset_versionshared_CpGs()])}
#'   \item{`probes_removed.csv`}{`probe`, `reason`: `chrY`, `snp`, `non_cpg`,
#'     `zhou` (array probe lists), `detP` (failed in more than
#'     `failed_sample_threshold` of samples), `na_fraction` (missing in more
#'     than `probe_na_threshold` of samples), `not_in_all_batches`}
#'   \item{`imputed.csv`}{per sample: `n_masked` (values set to `NA` for
#'     detection p-value above `detp_threshold`), `n_missing` (missing before
#'     imputation, masked values included) and `frac_missing`; `imputed_mask.csv`
#'     (`probe`, `basename`) with `save_na_mask = TRUE`}
#'   \item{`pheno_postQC.rds`}{sample sheet rows of the retained samples plus
#'     per-sample metrics (`mMed`, `uMed`, `rho`, `snp_outlier`,
#'     `frac_failed_probes`, `n_missing`, `xMed`, `yMed`, `predicted_sex`)}
#'   \item{`sample_qc.csv`}{the same metrics for all samples, with `removed`
#'     and `reasons`}
#'   \item{`rs_merged.rds`}{SNP probe betas of the retained samples, if saved}
#' }
#'
#' Batches are only a unit of work: probe failure fractions, the missing-value
#' threshold, imputation (`impute::impute.knn()`, k = 10) and the sex call are
#' all computed on the pooled retained samples, so the result does not depend
#' on how the samples were batched. Sex is called as in [minfi::getSex()]:
#' `yMed - xMed < -2` is female.
#'
#' @param run_dir Run directory with normalised batches ([qc_normalise()]).
#' @param force Recompute even if the outputs exist.
#' @return `TRUE` if the stage ran, `FALSE` if skipped; invisibly.
#' @examples
#' \dontrun{
#' qc_merge("runs/study1")
#' beta <- readRDS("runs/study1/beta_merged.rds")
#' pheno <- readRDS("runs/study1/pheno_postQC.rds")
#' }
#' @export
qc_merge <- function(run_dir, force = FALSE) {
  cfg <- read_config(run_dir)
  samples <- read_samples(run_dir)
  check_flag(force, "force")
  info <- array_info(cfg$array)
  outputs <- run_file(run_dir, c("beta_merged.rds", "probes_removed.csv", "imputed.csv",
                                 "pheno_postQC.rds", "sample_qc.csv",
                                 if (cfg$compatible_subset) "beta_merged_compatible.rds"))

  run_stage(run_dir, "qc_merge", outputs, force, seed = cfg$seed, body = function() {
    # pheno_postQC.rds is rewritten without estimates: let qc_estimates() redo them
    unlink(run_file(run_dir, "estimates.csv"))
    batches <- unique(samples$batch)
    b <- lapply(batches, function(x) read_batch_bundle(run_dir, x))
    names(b) <- batches

    # ---- retained samples, sample-sheet order ----
    removed <- do.call(rbind, lapply(b, `[[`, "removed"))
    retained <- samples$basename[!samples$basename %in% removed$basename]
    if (length(retained) == 0) rlang::abort("No samples retained after sample filtering.")
    log_info(sprintf("qc_merge: %d of %d samples retained across %d batch(es)",
                     length(retained), nrow(samples), length(batches)))

    # ---- inner join of the probe sets, cbind ----
    probe_sets <- lapply(b, function(x) rownames(x$beta))
    common <- Reduce(intersect, probe_sets)
    dropped <- setdiff(unique(unlist(probe_sets)), common)
    probes_removed <- data.frame(probe = dropped,
                                 reason = rep("not_in_all_batches", length(dropped)),
                                 stringsAsFactors = FALSE)
    if (length(dropped) > 0) {
      log_warn(sprintf("qc_merge: %d probes not present in every batch were dropped", length(dropped)))
    }
    beta <- do.call(cbind, lapply(b, function(x) x$beta[common, , drop = FALSE]))
    beta <- beta[, retained, drop = FALSE]

    # ---- mask failed values (detP), pooled failure fraction per probe ----
    has_detp <- !vapply(b, function(x) is.null(x$detP), logical(1))
    n_masked <- stats::setNames(rep(NA_integer_, length(retained)), retained)
    fail_count <- rep(0, length(common)); fail_n <- 0
    for (x in b[has_detp]) {
      cols <- intersect(colnames(x$detP), retained)
      if (length(cols) == 0) next
      failed <- x$detP[common, cols, drop = FALSE] > cfg$detp_threshold
      beta[, cols][failed] <- NA
      n_masked[cols] <- colSums(failed)
      fail_count <- fail_count + rowSums(failed)
      fail_n <- fail_n + length(cols)
    }
    if (fail_n > 0) {
      bad <- common[fail_count / fail_n > cfg$failed_sample_threshold]
      probes_removed <- rbind(probes_removed, data.frame(
        probe = bad, reason = rep("detP", length(bad)), stringsAsFactors = FALSE))
      log_info(sprintf("qc_merge: %d values masked (detP > %s); %d probes failed in > %s of samples",
                       sum(n_masked, na.rm = TRUE), cfg$detp_threshold, length(bad), cfg$failed_sample_threshold))
    } else {
      log_info("qc_merge: no detection p-values; detP masking and probe failure filter skipped")
    }

    # ---- array probe lists ----
    enabled <- c(chrY = cfg$remove_chrY, snp = cfg$remove_snp,
                 non_cpg = cfg$remove_non_cpg, zhou = cfg$remove_zhou)
    for (reason in names(enabled)[enabled]) {
      if (!reason %in% names(info$probe_lists)) next
      hits <- intersect(common, probe_list(info$probe_lists[[reason]]))
      probes_removed <- rbind(probes_removed, data.frame(
        probe = hits, reason = rep(reason, length(hits)), stringsAsFactors = FALSE))
      log_info(sprintf("qc_merge: %d %s probes removed", length(hits), reason))
    }
    beta <- beta[!rownames(beta) %in% probes_removed$probe, , drop = FALSE]

    # ---- probes with too many missing values ----
    na_frac <- rowMeans(is.na(beta))
    bad <- rownames(beta)[na_frac > cfg$probe_na_threshold]
    probes_removed <- rbind(probes_removed, data.frame(
      probe = bad, reason = rep("na_fraction", length(bad)), stringsAsFactors = FALSE))
    beta <- beta[!rownames(beta) %in% bad, , drop = FALSE]
    log_info(sprintf("qc_merge: %d probes missing in > %s of samples removed; %d probes retained",
                     length(bad), cfg$probe_na_threshold, nrow(beta)))
    utils::write.csv(probes_removed, run_file(run_dir, "probes_removed.csv"), row.names = FALSE)

    # ---- imputation ----
    missing <- is.na(beta)
    n_missing <- colSums(missing)
    if (cfg$save_na_mask) {
      idx <- which(missing, arr.ind = TRUE)
      utils::write.csv(data.frame(probe = rownames(beta)[idx[, 1]], basename = colnames(beta)[idx[, 2]],
                                  stringsAsFactors = FALSE),
                       run_file(run_dir, "imputed_mask.csv"), row.names = FALSE)
    }
    imputed <- FALSE
    if (cfg$impute && any(missing)) {
      rlang::check_installed("impute", reason = "for kNN imputation (`impute = TRUE`)")
      log_info(sprintf("qc_merge: imputing %d missing values (impute.knn, k = 10)", sum(missing)))
      utils::capture.output(beta <- impute::impute.knn(beta, k = 10, rowmax = 0.8)$data)
      imputed <- TRUE
    } else if (any(missing)) {
      log_info(sprintf("qc_merge: %d missing values left in place (`impute = FALSE`)", sum(missing)))
    }
    utils::write.csv(data.frame(
      basename = retained, n_masked = unname(n_masked[retained]),
      n_missing = unname(n_missing[retained]),
      frac_missing = unname(n_missing[retained]) / nrow(beta),
      imputed = imputed, stringsAsFactors = FALSE),
      run_file(run_dir, "imputed.csv"), row.names = FALSE)
    saveRDS(beta, run_file(run_dir, "beta_merged.rds"))

    if (cfg$compatible_subset) {
      compat <- subset_versionshared_CpGs(beta, cfg$array)
      saveRDS(compat, run_file(run_dir, "beta_merged_compatible.rds"))
      log_info(sprintf("qc_merge: %d EPIC v1/v2 compatible probes in beta_merged_compatible", nrow(compat)))
    }

    # ---- per-sample metrics, sex, pheno ----
    metrics <- sample_metrics(b, samples, cfg)
    metrics$n_missing <- unname(n_missing[metrics$basename])
    metrics$removed <- metrics$basename %in% removed$basename
    metrics$reasons <- vapply(metrics$basename, function(s) {
      paste(sort(unique(removed$reason[removed$basename == s])), collapse = ";")
    }, character(1))
    utils::write.csv(metrics, run_file(run_dir, "sample_qc.csv"), row.names = FALSE)

    pheno <- merge(samples, metrics[, setdiff(names(metrics), c("batch", "removed", "reasons"))],
                   by = "basename", sort = FALSE)
    pheno <- pheno[match(retained, pheno$basename), , drop = FALSE]
    rownames(pheno) <- NULL
    saveRDS(pheno, run_file(run_dir, "pheno_postQC.rds"))

    rs <- lapply(b, `[[`, "rs")
    if (all(!vapply(rs, is.null, logical(1)))) {
      rs <- do.call(cbind, rs)
      saveRDS(rs[, intersect(retained, colnames(rs)), drop = FALSE], run_file(run_dir, "rs_merged.rds"))
    }

    info_json <- read_run_info(run_dir)
    info_json$n_samples_retained <- length(retained)
    info_json$n_probes_retained <- nrow(beta)
    info_json$n_probes_removed <- length(unique(probes_removed$probe))
    write_run_info(run_dir, info_json)
  })
}

# Everything qc_merge needs from one batch directory.
read_batch_bundle <- function(run_dir, batch) {
  bdir <- batch_dir(run_dir, batch)
  rds <- function(f) { p <- file.path(bdir, f); if (file.exists(p)) readRDS(p) else NULL }
  beta <- rds("beta_norm.rds")
  if (is.null(beta)) rlang::abort(sprintf("Batch %s has no beta_norm.rds. Run qc_normalise() first.", batch))
  list(batch = batch, beta = beta, detP = rds("detP.rds"), qc = rds("qc.rds"), rho = rds("rho.rds"),
       sex_raw = rds("sex_raw.rds"), rs = rds("rs.rds"), removed = read_samples_removed(run_dir, batch))
}

# Per-sample metrics for every sample in the sample sheet (NA where a batch
# bundle lacks the item), with the pooled sex call.
sample_metrics <- function(b, samples, cfg) {
  m <- data.frame(basename = samples$basename, batch = samples$batch, stringsAsFactors = FALSE)
  pull <- function(item, col, key = "basename") {
    vals <- rep(NA_real_, nrow(m))
    for (x in b) {
      if (is.null(x[[item]])) next
      d <- x[[item]]
      keys <- if (key == "rownames") rownames(d) else d[[key]]
      i <- match(keys, m$basename)
      vals[i[!is.na(i)]] <- d[[col]][!is.na(i)]
    }
    vals
  }
  m$mMed <- pull("qc", "mMed")
  m$uMed <- pull("qc", "uMed")
  m$snp_outlier <- pull("qc", "snp_outlier")
  m$rho <- pull("rho", "rho", key = "rownames")
  m$frac_failed_probes <- NA_real_
  for (x in b) {
    if (is.null(x$detP)) next
    frac <- colMeans(x$detP > cfg$detp_threshold, na.rm = TRUE)
    m$frac_failed_probes[match(names(frac), m$basename)] <- unname(frac)
  }
  m$xMed <- pull("sex_raw", "xMed")
  m$yMed <- pull("sex_raw", "yMed")
  m$predicted_sex <- call_sex(m$xMed, m$yMed)
  m
}

# minfi::getSex() rule on pooled medians; NA where medians are missing.
call_sex <- function(xMed, yMed, cutoff = -2) {
  dd <- yMed - xMed
  sex <- ifelse(dd < cutoff, "F", "M")
  ok <- !is.na(dd)
  if (sum(ok) >= 2 && length(unique(dd[ok])) >= 2) {
    k <- stats::kmeans(dd[ok], centers = c(min(dd[ok]), max(dd[ok])))
    sex_k <- ifelse(k$cluster == which.min(k$centers), "F", "M")
    if (!identical(sex_k, sex[ok])) {
      log_warn("qc_merge: sex calls by cutoff and by clustering disagree; possibly only one sex present. Check xMed/yMed in sample_qc.csv.")
    }
  }
  sex
}

# A probe list shipped in data/, by object name.
probe_list <- function(name) {
  e <- new.env()
  utils::data(list = name, package = "prepDNAm", envir = e)
  as.character(get(name, e))
}
