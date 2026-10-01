#' Add cell-type, age and marker CpG estimates to the post-QC sample sheet
#'
#' Stage after [qc_merge()]. Computes, for the retained samples of a human
#' array run:
#' \describe{
#'   \item{`Epi`, `Fib`, `IC`}{epithelial, fibroblast and immune cell fractions
#'     (`EpiDISH`, `centEpiFibIC.m`)}
#'   \item{`CD4Tnv`, `Baso`, `CD4Tmem`, `Bmem`, `Bnv`, `Treg`, `CD8Tmem`,
#'     `CD8Tnv`, `Eos`, `NK`, `Neu`, `Mono`}{the immune fraction split into 12
#'     immune cell types (hierarchical EpiDISH, `cent12CT.m`; they sum to `IC`)}
#'   \item{`wid_general_age`}{WID general clock age (`WIDclocks`)}
#'   \item{`cg05575921`, `cg16867657`}{betas of the AHRR (smoking) and ELOVL2
#'     (age) marker CpGs}
#' }
#' and writes them to `estimates.csv` and as extra columns of
#' `pheno_postQC.rds`. Estimates whose package is not installed, or that are
#' not available for the array, are skipped (logged). EPIC v2 runs use the
#' EPIC v1/v2 compatible matrix, as these methods were trained on EPIC v1.
#'
#' [run_prepDNAm()] runs this stage when the configuration has `estimates =
#' TRUE` (the default). For a run finished without it, call `qc_estimates()`
#' directly; the merge is not repeated.
#'
#' @param run_dir Run directory of a merged run.
#' @param force Recompute even if `estimates.csv` exists.
#' @return `TRUE` if computed, `FALSE` if skipped; invisibly.
#' @seealso [qc_estimate_ic()], [qc_estimate_age()], [qc_index_cpgs()] for
#'   the individual estimates.
#' @examples
#' \dontrun{
#' qc_estimates("runs/study1")
#' pheno <- readRDS("runs/study1/pheno_postQC.rds")
#' }
#' @export
qc_estimates <- function(run_dir, force = FALSE) {
  cfg <- read_config(run_dir)
  check_flag(force, "force")
  if (!file.exists(run_file(run_dir, "beta_merged.rds"))) {
    rlang::abort("No beta_merged.rds; run qc_merge() first.")
  }
  out <- run_file(run_dir, "estimates.csv")
  ran <- run_stage(run_dir, "qc_estimates", out, force, seed = cfg$seed, body = function() {
    rd <- qc_data(run_dir)
    # the estimate functions explain skips with messages; send them to the log
    quietly <- function(expr) withCallingHandlers(expr, message = function(m) {
      log_info(sprintf("qc_estimates: %s", trimws(conditionMessage(m))))
      invokeRestart("muffleMessage")
    })
    pheno <- rd$pheno
    est <- data.frame(basename = pheno$basename, stringsAsFactors = FALSE)
    add <- function(df) {
      if (is.null(df)) return(invisible())
      df <- df[match(est$basename, df$basename), setdiff(names(df), "basename"), drop = FALSE]
      est <<- cbind(est, df)
    }

    ic <- quietly(qc_estimate_ic(rd, subtypes = TRUE))
    if (!is.null(ic)) {
      immune <- setdiff(names(ic), c("basename", "Epi", "Fib"))
      ic <- cbind(ic[, c("basename", "Epi", "Fib")], IC = rowSums(ic[, immune, drop = FALSE]),
                  ic[, immune, drop = FALSE])
    }
    add(ic)
    add(quietly(qc_estimate_age(rd)))
    mk <- quietly(qc_index_cpgs(rd))
    if (!is.null(mk)) attr(mk, "labels") <- NULL
    add(mk)

    # replace estimates from an earlier run of this stage; never overwrite
    # sample sheet columns
    new_cols <- setdiff(names(est), "basename")
    previous <- if (file.exists(out)) setdiff(names(utils::read.csv(out, nrows = 1, check.names = FALSE)), "basename")
    sheet_cols <- names(read_samples(run_dir))
    clash <- intersect(new_cols, sheet_cols)
    if (length(clash) > 0) {
      rlang::abort(sprintf("The sample sheet already has column(s) %s; rename them to add the estimates.",
                           paste(clash, collapse = ", ")))
    }
    pheno <- pheno[, setdiff(names(pheno), union(previous, new_cols)), drop = FALSE]
    pheno <- cbind(pheno, est[, new_cols, drop = FALSE])

    utils::write.csv(est, out, row.names = FALSE)
    saveRDS(pheno, run_file(run_dir, "pheno_postQC.rds"))
    log_info(sprintf("qc_estimates: %d estimate column(s) for %d samples added to pheno_postQC.rds",
                     length(new_cols), nrow(est)))
  })
  invisible(ran)
}
