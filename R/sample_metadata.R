#' The sample sheet (`pheno`)
#'
#' The optional `pheno` table passed to [qc_config()] annotates the samples
#' of a run. Only `basename` is required; a few columns with standard names
#' switch on checks in the report and are carried into the exports. Any other
#' columns are kept, shown in the report and tested against the principal
#' components (page 6).
#'
#' @section Standard columns:
#' \describe{
#'   \item{`basename`}{Required. The sample identifier exactly as in the data:
#'     the IDAT file name without `_Grn.idat` / `_Red.idat` (e.g.
#'     `206467010001_R01C01`), or the column name of a signal / beta matrix.
#'     Unique, no missing values.}
#'   \item{`sex`}{Reported sex: `F` or `M` (`female` / `male` and any
#'     capitalisation are accepted). Leave empty, `NA` or `unknown` if not
#'     known. Numeric codes (`0`/`1`, `1`/`2`) are not accepted because
#'     conventions differ. Compared with the sex predicted from chrX/chrY
#'     (report page 5) and exported as `Sex` 0/1 by [export_beta_csv()].}
#'   \item{`age`}{Chronological age in years, as a number (decimals allowed,
#'     e.g. `0.5` for six months; no units such as `45y`, no ranges). Empty or
#'     `NA` if not known. Plotted against epigenetic age (page 5) and
#'     exported as `Age`.}
#'   \item{`donor_id`}{Identifier of the person a sample comes from; samples
#'     of the same donor (replicates, longitudinal samples) share it. Used for
#'     the genotype agreement check ([qc_snp_agreement()]) and
#'     [qc_plot_snp_pairs()].}
#'   \item{`batch`}{Optional batch label, used as the batching mode when
#'     `qc_config(batch = NULL)`.}
#'   \item{`exclude`}{Optional; `TRUE` / `yes` / `1` marks samples to exclude
#'     from the start (they are still read and shown in the report).}
#' }
#'
#' @section Checks:
#' When the run is initialised ([qc_init()]), the `sex` and `age` columns are
#' validated: unrecognised sex values, and ages that are not numbers or are
#' negative, are an error listing the offending samples; ages above 120 give a
#' warning. A text column of numbers (e.g. from Excel) is converted to
#' numeric.
#'
#' Common alternative column names are recognised, with a message suggesting
#' the standard name: `gender`, `reported_sex` for `sex`;
#' `age_at_sampling`, `chron_age`, `chronological_age` for `age`; `donor`,
#' `person`, `person_id`, `individual`, `subject`, `patient`, `participant`
#' for `donor_id`. Any other column can be named explicitly with the
#' `sex_col`, `age_col` and `donor_col` arguments of the report and export
#' functions.
#'
#' @section Example:
#' \preformatted{
#' basename,sex,age,donor_id,group
#' 206467010001_R01C01,F,34.5,D001,control
#' 206467010001_R02C01,M,51,D002,case
#' 206467010001_R03C01,F,,D001,control
#' 206467010001_R04C01,unknown,62,D003,case
#' }
#' @name sample_sheet
NULL

SEX_ALIASES <- c("gender", "reported_sex")
AGE_ALIASES <- c("age_at_sampling", "chron_age", "chronological_age")
DONOR_ALIASES <- c("donor", "person", "person_id", "individual", "subject", "patient", "participant")

# Reported sex as "F"/"M"/NA.
normalise_sex <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("f", "female", "w", "woman", "women")] <- "F"
  out[x %in% c("m", "male", "man", "men")] <- "M"
  out
}

# Values that mean "not known".
is_missing_value <- function(x) {
  is.na(x) | tolower(trimws(as.character(x))) %in% c("", "na", "nan", "unknown", "unk", "not known")
}

# Find a sample sheet column: an explicit name (must exist), else the
# standard name, else the first alias (with a message).
find_col <- function(s, col, standard, aliases, what) {
  if (!is.null(col)) {
    if (!col %in% names(s)) rlang::abort(sprintf("No column \"%s\" in the sample sheet.", col))
    return(col)
  }
  if (standard %in% names(s)) return(standard)
  hit <- names(s)[tolower(names(s)) %in% c(standard, aliases)]
  if (length(hit) == 0) return(NULL)
  rlang::inform(sprintf("Using sample sheet column \"%s\" as %s (name it `%s`, or pass `%s_col`, to be explicit).",
                        hit[1], what, standard, sub("_id$", "", standard)))
  hit[1]
}

find_sex_col <- function(s, sex_col = NULL) find_col(s, sex_col, "sex", SEX_ALIASES, "reported sex")
find_age_col <- function(s, age_col = NULL) find_col(s, age_col, "age", AGE_ALIASES, "age")
find_donor_col <- function(s, donor_col = NULL) find_col(s, donor_col, "donor_id", DONOR_ALIASES, "donor identifier")

# Validate the standard sex / age columns of a pheno table (see ?sample_sheet);
# returns it with a numeric-text age column converted to numeric.
check_sample_metadata <- function(df) {
  show <- function(i, values) {
    i <- utils::head(i, 10)
    stats::setNames(sprintf("%s: \"%s\"", df$basename[i], values[i]), rep("x", length(i)))
  }
  sex_col <- find_sex_col(df)
  if (!is.null(sex_col)) {
    v <- df[[sex_col]]
    bad <- which(!is_missing_value(v) & is.na(normalise_sex(v)))
    if (length(bad) > 0) {
      rlang::abort(c(
        sprintf("%d sample(s) have an unrecognised value in the sex column `%s`:", length(bad), sex_col),
        show(bad, as.character(v)),
        i = "Use F or M (female/male accepted); leave empty or write NA / unknown if not known.",
        i = if (is.numeric(v) || all(grepl("^[0-9]+$", trimws(v[!is_missing_value(v)])))) {
          "Numeric codes are not accepted because conventions differ (0/1, 1/2); recode to F/M."
        },
        i = "See ?sample_sheet."))
    }
  }
  age_col <- find_age_col(df)
  if (!is.null(age_col)) {
    v <- df[[age_col]]
    missing <- is_missing_value(v)
    num <- suppressWarnings(as.numeric(trimws(as.character(v))))
    bad <- which(!missing & is.na(num))
    if (length(bad) > 0) {
      rlang::abort(c(
        sprintf("%d sample(s) have a non-numeric value in the age column `%s`:", length(bad), age_col),
        show(bad, as.character(v)),
        i = "Give age in years as a number (decimals allowed, no units or ranges); leave empty or NA if not known.",
        i = "See ?sample_sheet."))
    }
    num[missing] <- NA_real_
    neg <- which(!is.na(num) & num < 0)
    if (length(neg) > 0) {
      rlang::abort(c(sprintf("%d sample(s) have a negative age in `%s`:", length(neg), age_col),
                     show(neg, as.character(v))))
    }
    old <- which(!is.na(num) & num > 120)
    if (length(old) > 0) {
      rlang::warn(sprintf("%d sample(s) have an age above 120 in `%s` (%s); ages should be in years.",
                          length(old), age_col, paste(utils::head(df$basename[old], 5), collapse = ", ")))
    }
    if (!is.numeric(v)) df[[age_col]] <- num
  }
  df
}
