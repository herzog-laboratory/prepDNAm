#' Build and validate a run configuration
#'
#' All options of a prepDNAm run live in one configuration object, validated
#' here once and later written to `run_dir/config.yml` by [qc_init()]. Every
#' stage reads from it; nothing is passed between stages through the R
#' session.
#'
#' @section Input types:
#' `input_type` is inferred from `input` when `NULL`: a directory is scanned
#' for IDAT files (`"idat"`); a file is a signal-intensity matrix
#' (`"sigmatrix"`) if its header contains `unmeth_pattern`, otherwise a beta
#' matrix (`"beta"`). `.rds`/`.Rdata` files are always beta matrices.
#'
#' @section Batches:
#' Samples are processed in batches, the unit of work for reading,
#' normalisation and imputation. `batch` is one of
#' \describe{
#'   \item{`"none"`}{one batch}
#'   \item{`"dir"`}{IDATs only: the immediate subdirectory of `input` a sample sits in}
#'   \item{`"chip"`}{the Sentrix ID parsed from the basename (`<chip>_<position>`)}
#'   \item{`"<column>"`}{any column of `pheno`}
#' }
#' When `NULL`: the `pheno` column `batch` if there is one, else `"dir"` if
#' `input` is an IDAT directory with subdirectories, else `"none"`.
#'
#' @param input Directory of IDAT files (scanned recursively; `.idat.gz`
#'   accepted) or a signal-intensity / beta matrix file.
#' @param run_dir Directory to write all outputs into. Created by [qc_init()];
#'   its parent must exist.
#' @param array Array name, see [supported_arrays()].
#' @param input_type `"idat"`, `"sigmatrix"`, `"beta"`, or `NULL` to infer.
#' @param pheno Optional sample annotation: a data frame or a path to a
#'   `.csv`, `.tsv`/`.txt`, `.xlsx`, `.rds` or `.Rdata` file. Must have a
#'   `basename` column matching the IDAT basenames / matrix column names
#'   exactly. Samples in `pheno` but not in `input` are an error; samples in
#'   `input` but not in `pheno` are dropped with a message. Standard columns
#'   `sex`, `age` and `donor_id` switch on checks in the report; see
#'   [sample_sheet] for their format.
#' @param exclude Optional samples to exclude from the start: a character
#'   vector of basenames, or a path to a text file with one basename per line
#'   (or a `.csv` with a `basename` column). A logical `exclude` column in
#'   `pheno` is honoured as well.
#' @param batch Batching mode, see Batches.
#' @param run_name Label for the run; defaults to `basename(input)`.
#' @param impute Impute missing betas (k-nearest neighbours) after probe
#'   filtering.
#' @param save_na_mask Record every imputed `(probe, sample)` pair in
#'   `imputed.csv` instead of only per-sample counts. Can be large.
#' @param normalise Beta input only: run BMIQ on the input betas. `"auto"`
#'   runs it only if the type I / type II distributions look unnormalised;
#'   `TRUE`/`FALSE` force it. IDAT and signal-matrix input are always
#'   normalised.
#' @param compatible_subset EPIC v2 only: also write `beta_merged_compatible`
#'   restricted to EPIC v1/v2 compatible probes with v1 IDs (see
#'   [subset_versionshared_CpGs()]). Defaults to `TRUE` for EPIC v2.
#' @param save_rs Save the SNP (`rs`) probe betas, used for the genotype
#'   checks in the report ([qc_snp_agreement()]). `NULL` (default): save them
#'   whenever the array and input type have SNP probes. `TRUE` is ignored,
#'   with a message, if they have none.
#' @param save_sex Estimate sex from X/Y intensities. Ignored, with a message,
#'   if not available for the array or input type.
#' @param save_shiny IDAT input only: also save a `shinyMethyl` summary of
#'   each batch (`batches/<batch>/shiny.rds`) for interactive exploration with
#'   `shinyMethyl::runShinyMethyl()`. Needs the `shinyMethyl` package.
#' @param estimates Add immune cell fractions (EpiDISH), epigenetic age
#'   (WIDclocks) and marker CpG betas of the retained samples to
#'   `pheno_postQC.rds` and `estimates.csv` after the merge; see
#'   [qc_estimates()]. Human arrays; skipped where a package is missing.
#' @param export_digits Decimal places of the exported betas (the standard
#'   asks for at most 3), or `"full"` for full precision, matching
#'   `beta_merged.rds`; see [export_beta_csv()].
#' @param export_csv Also export the final beta matrix and sample metadata as
#'   CSV (DNA Methylation Array Data Standard, the format read by
#'   `biolearn`) into `run_dir/export/` at the end of [run_prepDNAm()]; see
#'   [export_beta_csv()]. For EPIC v2 with `compatible_subset`, the EPIC
#'   v1/v2 compatible matrix is exported as well.
#' @param intensity_threshold Minimum median methylated and unmethylated
#'   log2 intensity for a sample to pass.
#' @param detp_threshold Detection p-value above which a probe is considered
#'   failed in a sample.
#' @param failed_probe_threshold Maximum fraction of failed probes per sample.
#' @param failed_sample_threshold Maximum fraction of samples in which a probe
#'   may fail before it is removed.
#' @param probe_na_threshold Probes missing in more than this fraction of
#'   samples are removed before imputation.
#' @param remove_chrY,remove_snp,remove_non_cpg,remove_zhou Remove the
#'   corresponding probe sets (chrY; probes overlapping SNPs; non-CpG probes;
#'   cross-reactive/unreliable probes of Zhou et al. 2017).
#' @param unmeth_pattern,meth_pattern,detp_pattern Signal-matrix input only:
#'   the column-name suffixes identifying unmethylated, methylated and
#'   detection p-value columns. The sample name is the column name with the
#'   suffix removed. Set `detp_pattern = NULL` if there are no p-values.
#' @param cores Number of parallel workers; capped at
#'   [parallelly::availableCores()].
#' @param seed Random seed, set once per stage.
#' @return A list of class `prepDNAm_config`.
#' @examples
#' \dontrun{
#' cfg <- qc_config(
#'   input = "path/to/idats", run_dir = "path/to/run", array = "EPICv2",
#'   pheno = "pheno.csv", batch = "chip", cores = 4)
#' cfg
#' }
#' @export
qc_config <- function(input,
                      run_dir,
                      array,
                      input_type = NULL,
                      pheno = NULL,
                      exclude = NULL,
                      batch = NULL,
                      run_name = NULL,
                      impute = TRUE,
                      save_na_mask = FALSE,
                      normalise = "auto",
                      compatible_subset = NULL,
                      save_rs = NULL,
                      save_sex = TRUE,
                      save_shiny = FALSE,
                      estimates = TRUE,
                      export_csv = FALSE,
                      export_digits = 3,
                      intensity_threshold = 9.5,
                      detp_threshold = 0.01,
                      failed_probe_threshold = 0.1,
                      failed_sample_threshold = 0.1,
                      probe_na_threshold = 0.8,
                      remove_chrY = TRUE,
                      remove_snp = TRUE,
                      remove_non_cpg = TRUE,
                      remove_zhou = TRUE,
                      unmeth_pattern = "Unmethylated signal",
                      meth_pattern = "Methylated signal",
                      detp_pattern = "Detection Pval",
                      cores = 1L,
                      seed = 1L) {
  # ---- input / run_dir ----
  check_string(input, "input")
  if (!file.exists(input)) {
    rlang::abort(sprintf("`input` does not exist: %s", input))
  }
  input <- normalizePath(input, mustWork = TRUE)

  check_string(run_dir, "run_dir")
  parent <- dirname(normalizePath(run_dir, mustWork = FALSE))
  if (!dir.exists(parent)) {
    rlang::abort(sprintf("Parent directory of `run_dir` does not exist: %s", parent))
  }
  run_dir <- file.path(normalizePath(parent), basename(run_dir))

  # ---- array ----
  info <- array_info(array)
  array <- info$array

  # ---- input type ----
  check_string(unmeth_pattern, "unmeth_pattern")
  check_string(meth_pattern, "meth_pattern")
  if (!is.null(detp_pattern)) check_string(detp_pattern, "detp_pattern")
  if (is.null(input_type)) {
    input_type <- detect_input_type(input, unmeth_pattern)
  } else {
    input_type <- rlang::arg_match(input_type, c("idat", "sigmatrix", "beta"))
    if (input_type == "idat" && !dir.exists(input)) {
      rlang::abort("`input_type = \"idat\"` needs `input` to be a directory.")
    }
    if (input_type != "idat" && dir.exists(input)) {
      rlang::abort(sprintf("`input_type = \"%s\"` needs `input` to be a file.", input_type))
    }
  }

  # ---- run name ----
  if (is.null(run_name)) {
    run_name <- if (input_type == "idat") basename(input) else tools::file_path_sans_ext(basename(input))
  }
  check_string(run_name, "run_name")

  # ---- pheno / exclude ----
  pheno_cols <- character()
  if (!is.null(pheno)) {
    if (is.data.frame(pheno)) {
      pheno <- check_pheno_df(pheno, "pheno")
      pheno_cols <- names(pheno)
    } else {
      check_string(pheno, "pheno")
      if (!file.exists(pheno)) {
        rlang::abort(sprintf("`pheno` file does not exist: %s", pheno))
      }
      pheno <- normalizePath(pheno, mustWork = TRUE)
      pheno_cols <- names(read_pheno(pheno))
    }
  }
  if (!is.null(exclude)) {
    if (is.character(exclude) && length(exclude) == 1 && file.exists(exclude)) {
      exclude <- normalizePath(exclude, mustWork = TRUE)
    } else if (!is.character(exclude) || anyNA(exclude)) {
      rlang::abort("`exclude` must be a character vector of basenames or a path to a file listing them.")
    }
  }

  # ---- batch ----
  if (is.null(batch)) {
    batch <- if ("batch" %in% pheno_cols) {
      "batch"
    } else if (input_type == "idat" && has_subdirs(input)) {
      "dir"
    } else {
      "none"
    }
    rlang::inform(sprintf("`batch` not set; using \"%s\".", batch))
  }
  check_string(batch, "batch")
  if (batch == "dir" && input_type != "idat") {
    rlang::abort("`batch = \"dir\"` is only available for IDAT input.")
  }
  if (!batch %in% c("none", "dir", "chip")) {
    if (is.null(pheno)) {
      rlang::abort(c(
        sprintf("`batch = \"%s\"` is not a batching mode and no `pheno` was given.", batch),
        i = "Use \"none\", \"dir\", \"chip\" or the name of a `pheno` column."))
    }
    if (!batch %in% pheno_cols) {
      rlang::abort(c(
        sprintf("`batch = \"%s\"` is neither a batching mode nor a `pheno` column.", batch),
        i = sprintf("`pheno` columns: %s.", paste(pheno_cols, collapse = ", "))))
    }
  }

  # ---- feature flags ----
  check_flag(impute, "impute")
  check_flag(save_na_mask, "save_na_mask")
  if (!(identical(normalise, "auto") || is_flag(normalise))) {
    rlang::abort("`normalise` must be \"auto\", TRUE or FALSE.")
  }
  if (is.null(compatible_subset)) {
    compatible_subset <- info$needs_compatible_subset
  }
  check_flag(compatible_subset, "compatible_subset")
  if (compatible_subset && !info$needs_compatible_subset) {
    rlang::abort(sprintf("`compatible_subset = TRUE` only applies to EPIC v2, not %s.", array))
  }
  rs_available <- info$has_snp_probes && input_type != "beta"
  if (is.null(save_rs)) save_rs <- rs_available
  check_flag(save_rs, "save_rs")
  check_flag(save_sex, "save_sex")
  check_flag(save_shiny, "save_shiny")
  check_flag(estimates, "estimates")
  check_flag(export_csv, "export_csv")
  if (!(identical(export_digits, "full") || (is.numeric(export_digits) && length(export_digits) == 1 &&
                                              !is.na(export_digits) && export_digits >= 0))) {
    rlang::abort("`export_digits` must be a number of decimal places or \"full\".")
  }
  if (save_shiny && input_type != "idat") {
    rlang::inform("`save_shiny` ignored: only available for IDAT input.")
    save_shiny <- FALSE
  }
  if (save_rs && !rs_available) {
    rlang::inform(sprintf("`save_rs` ignored: no SNP probes for %s / %s input.", array, input_type))
    save_rs <- FALSE
  }
  if (save_sex && (!info$has_sex_estimation || input_type == "beta")) {
    rlang::inform(sprintf("`save_sex` ignored: sex estimation not available for %s / %s input.", array, input_type))
    save_sex <- FALSE
  }

  # ---- thresholds ----
  check_number(intensity_threshold, "intensity_threshold", min = 0)
  check_number(detp_threshold, "detp_threshold", min = 0, max = 1)
  check_number(failed_probe_threshold, "failed_probe_threshold", min = 0, max = 1)
  check_number(failed_sample_threshold, "failed_sample_threshold", min = 0, max = 1)
  check_number(probe_na_threshold, "probe_na_threshold", min = 0, max = 1)
  check_flag(remove_chrY, "remove_chrY")
  check_flag(remove_snp, "remove_snp")
  check_flag(remove_non_cpg, "remove_non_cpg")
  check_flag(remove_zhou, "remove_zhou")
  if (length(info$probe_lists) == 0 &&
      any(remove_chrY, remove_snp, remove_non_cpg, remove_zhou)) {
    rlang::inform(sprintf("No probe-removal lists are available for %s; `remove_*` ignored.", array))
    remove_chrY <- remove_snp <- remove_non_cpg <- remove_zhou <- FALSE
  }

  # ---- compute ----
  check_number(cores, "cores", min = 1)
  cores <- as.integer(cores)
  available <- parallelly::availableCores()
  if (cores > available) {
    rlang::inform(sprintf("`cores = %d` capped at %d available.", cores, available))
    cores <- as.integer(available)
  }
  check_number(seed, "seed")
  seed <- as.integer(seed)

  structure(list(
    input = input,
    input_type = input_type,
    run_dir = run_dir,
    run_name = run_name,
    array = array,
    pheno = pheno,
    exclude = exclude,
    batch = batch,
    impute = impute,
    save_na_mask = save_na_mask,
    normalise = normalise,
    compatible_subset = compatible_subset,
    save_rs = save_rs,
    save_sex = save_sex,
    save_shiny = save_shiny,
    estimates = estimates,
    export_csv = export_csv,
    export_digits = export_digits,
    intensity_threshold = intensity_threshold,
    detp_threshold = detp_threshold,
    failed_probe_threshold = failed_probe_threshold,
    failed_sample_threshold = failed_sample_threshold,
    probe_na_threshold = probe_na_threshold,
    remove_chrY = remove_chrY,
    remove_snp = remove_snp,
    remove_non_cpg = remove_non_cpg,
    remove_zhou = remove_zhou,
    unmeth_pattern = unmeth_pattern,
    meth_pattern = meth_pattern,
    detp_pattern = detp_pattern,
    cores = cores,
    seed = seed,
    prepDNAm_version = as.character(utils::packageVersion("prepDNAm")),
    created = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  ), class = "prepDNAm_config")
}

#' @export
print.prepDNAm_config <- function(x, ...) {
  cat("<prepDNAm_config>\n")
  pheno <- if (is.null(x$pheno)) "none" else if (is.data.frame(x$pheno)) {
    sprintf("data frame (%d rows)", nrow(x$pheno))
  } else x$pheno
  exclude <- if (is.null(x$exclude)) "none" else if (length(x$exclude) == 1 && file.exists(x$exclude)) {
    x$exclude
  } else sprintf("%d samples", length(x$exclude))
  fields <- c(
    run_name = x$run_name, input = x$input, input_type = x$input_type,
    run_dir = x$run_dir, array = x$array, pheno = pheno, exclude = exclude,
    batch = x$batch,
    normalise = if (x$input_type == "beta") as.character(x$normalise) else "BMIQ",
    impute = as.character(x$impute),
    compatible_subset = as.character(x$compatible_subset),
    save_rs = as.character(x$save_rs), save_sex = as.character(x$save_sex),
    save_shiny = as.character(x$save_shiny),
    estimates = as.character(!isFALSE(x$estimates)),
    export_csv = as.character(isTRUE(x$export_csv)),
    export_digits = as.character(if (is.null(x$export_digits)) 3 else x$export_digits),
    thresholds = sprintf(
      "intensity %s, detP %s, failed probes %s, failed samples %s, probe NA %s",
      x$intensity_threshold, x$detp_threshold, x$failed_probe_threshold,
      x$failed_sample_threshold, x$probe_na_threshold),
    remove = paste(c("chrY", "snp", "non_cpg", "zhou")[
      c(x$remove_chrY, x$remove_snp, x$remove_non_cpg, x$remove_zhou)], collapse = ", "),
    cores = as.character(x$cores), seed = as.character(x$seed),
    version = x$prepDNAm_version, created = x$created)
  cat(sprintf("  %-18s %s\n", names(fields), fields), sep = "")
  invisible(x)
}

# ---- helpers -----------------------------------------------------------------

has_subdirs <- function(dir) {
  length(list.dirs(dir, recursive = FALSE)) > 0
}

# Directory -> idat; .rds/.Rdata -> beta; other files: sigmatrix if the header
# contains `unmeth_pattern`, else beta.
detect_input_type <- function(input, unmeth_pattern) {
  if (dir.exists(input)) return("idat")
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", input)))
  if (ext %in% c("rds", "rdata", "rda")) return("beta")
  header <- read_header(input)
  if (any(grepl(unmeth_pattern, header, fixed = TRUE))) "sigmatrix" else "beta"
}

is_flag <- function(x) is.logical(x) && length(x) == 1 && !is.na(x)

check_flag <- function(x, name) {
  if (!is_flag(x)) rlang::abort(sprintf("`%s` must be TRUE or FALSE.", name))
}

check_string <- function(x, name) {
  if (!is.character(x) || length(x) != 1 || is.na(x) || !nzchar(x)) {
    rlang::abort(sprintf("`%s` must be a single non-empty string.", name))
  }
}

check_number <- function(x, name, min = -Inf, max = Inf) {
  if (!is.numeric(x) || length(x) != 1 || is.na(x) || x < min || x > max) {
    range <- if (is.finite(min) && is.finite(max)) sprintf(" in [%s, %s]", min, max) else
      if (is.finite(min)) sprintf(" >= %s", min) else ""
    rlang::abort(sprintf("`%s` must be a single number%s.", name, range))
  }
}
