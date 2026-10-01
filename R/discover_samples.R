#' Build the sample sheet for a run
#'
#' Finds the samples in `config$input`, joins the optional `pheno` table by
#' exact `basename`, flags excluded samples and assigns each sample to a batch
#' (see the Batches section of [qc_config()]). [qc_init()] writes the result
#' to `run_dir/samples.csv`; every later stage works from that file.
#'
#' For IDAT input the directory is scanned recursively for `*_Grn.idat` /
#' `*_Grn.idat.gz`; the basename is the file name with that suffix removed and
#' a matching `_Red` file must exist. For signal-matrix input the sample names
#' are parsed from the column headers (`<sample><unmeth_pattern>` etc.); for
#' beta input they are the column names.
#'
#' Samples listed in `pheno` that are not found in `input` are an error;
#' samples found in `input` but not listed in `pheno` are dropped with a
#' message.
#'
#' @param config A [qc_config()] object.
#' @return A data frame with one row per sample and columns `basename`,
#'   `batch`, `exclude` (logical), `path` (for IDATs the `minfi` basename,
#'   i.e. the file path without `_Grn.idat`; otherwise the input file), and
#'   any further `pheno` columns. Ordered as `pheno` if given, otherwise by
#'   basename.
#' @examples
#' \dontrun{
#' samples <- discover_samples(cfg)
#' table(samples$batch, samples$exclude)
#' }
#' @export
discover_samples <- function(config) {
  if (!inherits(config, "prepDNAm_config")) {
    rlang::abort("`config` must be a `prepDNAm_config` object from `qc_config()`.")
  }
  found <- switch(config$input_type,
    idat = scan_idats(config$input),
    sigmatrix = parse_sigmatrix_samples(config$input, config$unmeth_pattern,
                                        config$meth_pattern, config$detp_pattern),
    beta = parse_beta_samples(config$input))
  if (nrow(found) == 0) {
    rlang::abort(sprintf("No samples found in %s.", config$input))
  }

  # ---- pheno join (exact basename) ----
  if (!is.null(config$pheno)) {
    pheno <- if (is.data.frame(config$pheno)) {
      check_pheno_df(config$pheno, "pheno")
    } else {
      read_pheno(config$pheno)
    }
    pheno <- check_sample_metadata(pheno)
    missing <- setdiff(pheno$basename, found$basename)
    if (length(missing) > 0) {
      rlang::abort(c(
        sprintf("%d sample(s) in `pheno` were not found in `input`:", length(missing)),
        stats::setNames(utils::head(missing, 10), rep("x", min(10, length(missing)))),
        if (length(missing) > 10) sprintf("... and %d more.", length(missing) - 10)))
    }
    extra <- setdiff(found$basename, pheno$basename)
    if (length(extra) > 0) {
      rlang::inform(sprintf(
        "%d sample(s) in `input` are not in `pheno` and will be dropped: %s",
        length(extra), paste(utils::head(extra, 10), collapse = ", ")))
    }
    if ("batch" %in% names(pheno) && config$batch != "batch") {
      rlang::inform("`pheno` has a `batch` column but another batching mode was chosen; keeping it as `pheno_batch`.")
      names(pheno)[names(pheno) == "batch"] <- "pheno_batch"
    }
    samples <- merge(pheno, found, by = "basename", sort = FALSE)
    samples <- samples[match(pheno$basename, samples$basename), , drop = FALSE]
  } else {
    samples <- found[order(found$basename), , drop = FALSE]
  }
  rownames(samples) <- NULL

  # ---- exclude: file / vector, OR'd with a pheno `exclude` column ----
  exclude <- read_exclude(config$exclude)
  unknown <- setdiff(exclude, samples$basename)
  if (length(unknown) > 0) {
    rlang::warn(sprintf("%d basename(s) in `exclude` are not among the samples: %s",
                        length(unknown), paste(utils::head(unknown, 10), collapse = ", ")))
  }
  from_pheno <- if ("exclude" %in% names(samples)) samples$exclude else FALSE
  samples$exclude <- samples$basename %in% exclude | from_pheno

  # ---- batch ----
  samples$batch <- resolve_batch(samples, config)
  samples$dir <- NULL

  first <- c("basename", "batch", "exclude", "path")
  samples[, c(first, setdiff(names(samples), first)), drop = FALSE]
}

# ---- input scanning ----------------------------------------------------------

# One row per *_Grn.idat(.gz): basename, path (minfi basename), dir.
scan_idats <- function(input) {
  grn <- list.files(input, pattern = "_Grn\\.idat(\\.gz)?$", recursive = TRUE,
                    full.names = TRUE)
  basename <- sub("_Grn\\.idat(\\.gz)?$", "", basename(grn))
  red <- sub("_Grn\\.idat", "_Red.idat", grn)
  no_red <- !file.exists(red)
  if (any(no_red)) {
    rlang::abort(c(
      sprintf("%d sample(s) have a _Grn but no _Red IDAT:", sum(no_red)),
      stats::setNames(utils::head(basename[no_red], 10), rep("x", min(10, sum(no_red))))))
  }
  dup <- unique(basename[duplicated(basename)])
  if (length(dup) > 0) {
    rlang::abort(c(
      sprintf("%d basename(s) occur more than once under `input`:", length(dup)),
      stats::setNames(utils::head(dup, 10), rep("x", min(10, length(dup))))))
  }
  data.frame(
    basename = basename,
    path = file.path(dirname(grn), basename),
    dir = dirname(grn),
    stringsAsFactors = FALSE)
}

# Sample names from GEO-style headers: columns named <sample><pattern>. The
# unmethylated and methylated column sets must name the same samples.
parse_sigmatrix_samples <- function(input, unmeth_pattern, meth_pattern, detp_pattern) {
  header <- read_header(input)
  u <- header_samples(header, unmeth_pattern)
  m <- header_samples(header, meth_pattern)
  if (length(u) == 0) {
    rlang::abort(c(
      sprintf("No columns matching `unmeth_pattern` (\"%s\") in %s.", unmeth_pattern, input),
      i = sprintf("First columns: %s", paste(utils::head(header, 5), collapse = ", "))))
  }
  if (!setequal(u, m)) {
    rlang::abort(c(
      "Unmethylated and methylated columns do not name the same samples.",
      i = sprintf("Unmethylated only: %s", paste(utils::head(setdiff(u, m), 5), collapse = ", ")),
      i = sprintf("Methylated only: %s", paste(utils::head(setdiff(m, u), 5), collapse = ", "))))
  }
  if (!is.null(detp_pattern)) {
    p <- header_samples(header, detp_pattern)
    if (length(p) > 0 && !setequal(u, p)) {
      rlang::abort("Detection p-value columns do not name the same samples as the signal columns.")
    }
  }
  data.frame(basename = u, path = input, stringsAsFactors = FALSE)
}

# Column names carrying `pattern`, with the pattern and any separator
# characters around it stripped.
header_samples <- function(header, pattern) {
  hits <- header[grepl(pattern, header, fixed = TRUE)]
  s <- sub(pattern, "", hits, fixed = TRUE)
  s <- gsub("^[[:space:][:punct:]]+|[[:space:][:punct:]]+$", "", s)
  s
}

# Beta matrix: column names other than the first (probe IDs) for text files;
# colnames() of the object for .rds/.Rdata.
parse_beta_samples <- function(input) {
  ext <- tolower(tools::file_ext(sub("\\.gz$", "", input)))
  samples <- if (ext == "rds") {
    colnames(readRDS(input))
  } else if (ext %in% c("rdata", "rda")) {
    colnames(load_single_object(input))
  } else {
    read_header(input)[-1]
  }
  if (is.null(samples) || length(samples) == 0) {
    rlang::abort(sprintf("No sample columns found in %s.", input))
  }
  data.frame(basename = samples, path = input, stringsAsFactors = FALSE)
}

# First line of a delimited text file (optionally gzipped), split on tab or
# comma, quotes removed.
read_header <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  line <- readLines(con, n = 1, warn = FALSE)
  if (length(line) == 0) rlang::abort(sprintf("%s is empty.", path))
  sep <- if (grepl("\t", line)) "\t" else ","
  header <- strsplit(line, sep, fixed = TRUE)[[1]]
  gsub("^\"|\"$", "", trimws(header))
}

# ---- pheno / exclude ---------------------------------------------------------

# Read a pheno table from csv / tsv / txt / xlsx / rds / Rdata and check it.
read_pheno <- function(path) {
  ext <- tolower(tools::file_ext(path))
  df <- switch(ext,
    csv = utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE),
    tsv = ,
    txt = utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE),
    xlsx = ,
    xls = {
      rlang::check_installed("readxl", reason = "to read .xlsx pheno files")
      as.data.frame(readxl::read_excel(path), stringsAsFactors = FALSE)
    },
    rds = readRDS(path),
    rdata = ,
    rda = load_single_object(path),
    rlang::abort(sprintf(
      "Unsupported `pheno` file type \".%s\"; use .csv, .tsv, .txt, .xlsx, .rds or .Rdata.", ext)))
  check_pheno_df(df, sprintf("pheno (%s)", basename(path)))
}

# A pheno data frame must have a unique, complete `basename` column
# (`Basename` etc. are accepted and renamed). An `exclude` column is coerced
# to logical; `path`/`dir` columns are renamed `pheno_path`/`pheno_dir`.
check_pheno_df <- function(df, what) {
  if (!is.data.frame(df)) {
    rlang::abort(sprintf("`%s` must be a data frame.", what))
  }
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  if (!"basename" %in% names(df)) {
    alt <- which(tolower(names(df)) == "basename")
    if (length(alt) == 1) {
      names(df)[alt] <- "basename"
    } else {
      rlang::abort(c(
        sprintf("`%s` has no `basename` column.", what),
        i = sprintf("Columns: %s", paste(names(df), collapse = ", "))))
    }
  }
  df$basename <- as.character(df$basename)
  if (anyNA(df$basename) || any(!nzchar(df$basename))) {
    rlang::abort(sprintf("`%s` has missing values in `basename`.", what))
  }
  dup <- unique(df$basename[duplicated(df$basename)])
  if (length(dup) > 0) {
    rlang::abort(c(
      sprintf("`%s` has duplicated basenames:", what),
      stats::setNames(utils::head(dup, 10), rep("x", min(10, length(dup))))))
  }
  # `exclude` is honoured (logical); `batch` is kept when it is the selected
  # batch mode (handled by the caller); `path`/`dir` clash with our columns.
  if ("exclude" %in% names(df)) {
    ex <- df$exclude
    if (is.character(ex)) ex <- tolower(trimws(ex)) %in% c("true", "t", "yes", "y", "1")
    if (!is.logical(ex)) ex <- as.logical(ex)
    if (anyNA(ex)) ex[is.na(ex)] <- FALSE
    df$exclude <- ex
  }
  for (col in intersect(names(df), c("path", "dir"))) {
    names(df)[names(df) == col] <- paste0("pheno_", col)
  }
  df
}

# Basenames to exclude: NULL, a character vector, or a file (one per line, or
# a .csv with a basename column).
read_exclude <- function(exclude) {
  if (is.null(exclude)) return(character())
  if (length(exclude) == 1 && file.exists(exclude)) {
    ext <- tolower(tools::file_ext(exclude))
    if (ext == "csv") {
      df <- utils::read.csv(exclude, stringsAsFactors = FALSE)
      return(as.character(check_pheno_df(df, "exclude")$basename))
    }
    lines <- trimws(readLines(exclude, warn = FALSE))
    return(lines[nzchar(lines) & !startsWith(lines, "#")])
  }
  exclude
}

# The single object in an .Rdata file (or its first data frame / matrix).
load_single_object <- function(path) {
  e <- new.env()
  objs <- load(path, envir = e)
  if (length(objs) == 1) return(get(objs, e))
  ok <- objs[vapply(objs, function(o) is.data.frame(get(o, e)) || is.matrix(get(o, e)), logical(1))]
  if (length(ok) == 0) {
    rlang::abort(sprintf("%s contains no data frame or matrix.", path))
  }
  if (length(ok) > 1) {
    rlang::inform(sprintf("%s contains several objects; using `%s`.", basename(path), ok[1]))
  }
  get(ok[1], e)
}

# ---- batches -----------------------------------------------------------------

# Batch label per sample, safe to use as a directory name.
resolve_batch <- function(samples, config) {
  mode <- config$batch
  batch <- if (mode == "none") {
    rep("all", nrow(samples))
  } else if (mode == "dir") {
    rel <- substring(samples$dir, nchar(config$input) + 2)
    top <- vapply(strsplit(rel, "/", fixed = TRUE), function(p) p[1], character(1))
    if (anyNA(top) || any(!nzchar(top))) {
      rlang::abort(c(
        "`batch = \"dir\"` needs every IDAT in a subdirectory of `input`, but some are directly in it.",
        i = "Move them into a subdirectory or use `batch = \"none\"` / `\"chip\"`."))
    }
    top
  } else if (mode == "chip") {
    ok <- grepl("^[0-9]+_R[0-9]{2}C[0-9]{2}$", samples$basename)
    if (!all(ok)) {
      rlang::abort(c(
        "`batch = \"chip\"` needs basenames of the form <SentrixID>_R01C01; these do not match:",
        stats::setNames(utils::head(samples$basename[!ok], 10), rep("x", min(10, sum(!ok))))))
    }
    sub("_R[0-9]{2}C[0-9]{2}$", "", samples$basename)
  } else {
    if (!mode %in% names(samples)) {
      rlang::abort(sprintf("`batch` column \"%s\" not found in `pheno`.", mode))
    }
    if (mode == "batch" && is.logical(samples$batch)) {
      rlang::abort("`pheno` column \"batch\" is logical; it must hold batch labels.")
    }
    b <- as.character(samples[[mode]])
    if (anyNA(b) || any(!nzchar(b))) {
      rlang::abort(sprintf("`batch` column \"%s\" has missing values.", mode))
    }
    b
  }
  safe <- gsub("[^A-Za-z0-9._-]+", "_", batch)
  if (any(safe != batch)) {
    rlang::inform("Batch labels contain characters not safe for directory names; replaced with \"_\".")
  }
  safe
}
