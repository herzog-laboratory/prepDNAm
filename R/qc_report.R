#' Render the HTML QC report
#'
#' Final stage. Renders the report pages shipped in `inst/rmd/` into
#' `run_dir/report/` (open `report/index.html`). It only reads the run
#' directory, so it can be re-run at any time, also on a partially finished
#' run (pages then show what is missing).
#'
#' Needs `rmarkdown` (with pandoc), `plotly`, `DT` and `htmltools`. Optional
#' content depends on `ewastools` (genotype checks), `EpiDISH` (immune cell
#' fractions), `WIDclocks` (epigenetic age) and `uwot` (UMAP).
#'
#' The logic behind every page lives in package functions
#' (`qc_table_*()`, `qc_plot_*()`, [qc_dimred()], ...), so the same figures
#' can be produced in custom reports.
#'
#' @param run_dir Run directory.
#' @param force Re-render even if the report exists.
#' @param open Open `index.html` in the browser afterwards.
#' @param self_contained `TRUE` (default): every page is a single HTML file
#'   with its JavaScript embedded, so pages can be shared individually (each
#'   is a few MB). `FALSE`: pages share one `report/site_libs/` folder, which
#'   must be kept next to them; smaller in total.
#' @param sex_col,donor_col Sample sheet columns with reported sex and donor
#'   identifiers for the identity checks (see [qc_table_sex()],
#'   [qc_snp_agreement()]); `NULL` to auto-detect.
#' @return `TRUE` if rendered, `FALSE` if skipped; invisibly.
#' @examples
#' \dontrun{
#' qc_report("runs/study1", open = TRUE)
#' qc_report("runs/study1", force = TRUE, donor_col = "patient_id")
#' }
#' @export
qc_report <- function(run_dir, force = FALSE, open = FALSE, self_contained = TRUE,
                      sex_col = NULL, donor_col = NULL) {
  cfg <- read_config(run_dir)
  check_flag(force, "force")
  check_flag(open, "open")
  check_flag(self_contained, "self_contained")
  if (!is.null(sex_col)) check_string(sex_col, "sex_col")
  if (!is.null(donor_col)) check_string(donor_col, "donor_col")
  rlang::check_installed(c("rmarkdown", "knitr", "plotly", "DT", "htmltools"), reason = "to render the report")
  if (!rmarkdown::pandoc_available("2.0")) {
    rlang::abort(c(
      "pandoc (>= 2.0) is needed to render the report but was not found.",
      i = "RStudio, Positron and Quarto ship one: set the RSTUDIO_PANDOC environment variable to its directory, e.g. Sys.setenv(RSTUDIO_PANDOC = \"/Applications/quarto/bin/tools/aarch64\").",
      i = "Or install pandoc from https://pandoc.org and put it on the PATH."))
  }
  report_dir <- run_file(run_dir, "report")
  ran <- run_stage(run_dir, "qc_report", file.path(report_dir, "index.html"), force,
                   seed = cfg$seed, body = function() {
    src <- withr::local_tempdir("prepDNAm-report")
    file.copy(list.files(system.file("rmd", package = "prepDNAm"), full.names = TRUE), src)
    # read by _setup.R on every page
    yaml::write_yaml(list(run_dir = normalizePath(run_dir), sex_col = sex_col, donor_col = donor_col),
                     file.path(src, "_report_options.yml"))
    out <- file.path(src, "_site")
    if (self_contained) {
      pages <- list.files(src, pattern = "\\.Rmd$")
      for (p in pages) {
        log_info(sprintf("qc_report: rendering %s", p))
        rmarkdown::render(file.path(src, p), output_dir = out, envir = new.env(parent = globalenv()),
                          quiet = TRUE)
      }
    } else {
      log_info("qc_report: rendering pages")
      rmarkdown::render_site(input = src, envir = new.env(parent = globalenv()), quiet = TRUE)
    }
    # replace run_dir/report with the rendered pages
    if (dir.exists(report_dir)) unlink(report_dir, recursive = TRUE)
    dir.create(report_dir)
    file.copy(list.files(out, full.names = TRUE, all.files = TRUE, no.. = TRUE),
              report_dir, recursive = TRUE)
    log_info(sprintf("qc_report: written to %s", report_dir))
  })
  if (open && file.exists(file.path(report_dir, "index.html"))) {
    utils::browseURL(file.path(report_dir, "index.html"))
  }
  invisible(ran)
}
