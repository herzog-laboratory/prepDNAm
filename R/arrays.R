#' Supported Illumina methylation arrays
#'
#' Every array- or species-specific fact used by the package (annotation
#' package names, probes per array, samples per BeadChip, which probe-removal
#' lists apply, whether SNP probes / sex estimation / epigenetic clocks are
#' available) comes from one internal lookup table. `array_info()` returns the
#' row for one array; `supported_arrays()` lists the canonical names.
#'
#' Array names are matched case-insensitively and ignoring punctuation, so
#' `"EPIC v2"`, `"epicv2"` and `"IlluminaHumanMethylationEPICv2"` all resolve
#' to `"EPICv2"`. A `minfi` object (`RGChannelSet`, `MethylSet`, ...) is
#' resolved from its `annotation` slot.
#'
#' @param array Array name (see `supported_arrays()` and the aliases below), or
#'   a `minfi` object whose annotation identifies the array.
#' @return `array_info()`: a named list with elements
#'   \describe{
#'     \item{`array`}{canonical name: `"450k"`, `"EPIC"`, `"EPICv2"`, `"mouse"`}
#'     \item{`label`}{Illumina product name}
#'     \item{`species`, `genome`}{`"human"`/`"mouse"`;
#'       `"hg19"`/`"hg38"`/`"mm10"`}
#'     \item{`minfi_array`, `minfi_annotation`}{the two elements of a
#'       `minfi` object's `annotation` slot}
#'     \item{`manifest_pkg`, `anno_pkg`}{Bioconductor/GitHub package names}
#'     \item{`n_probes`}{CpG probes on the array (type I + II); `NA` if the
#'       manifest was not available when the table was built}
#'     \item{`samples_per_chip`}{8 or 12}
#'     \item{`probe_id_suffix`}{probe IDs carry a design suffix
#'       (`cg00000029_TC21`) instead of bare `cg` IDs}
#'     \item{`has_snp_probes`}{rs probes readable with [minfi::getSnpBeta()]}
#'     \item{`has_sex_estimation`}{[minfi::getSex()] applicable}
#'     \item{`has_clocks`}{immune-cell fractions, WID clocks and smoking index
#'       available (human only)}
#'     \item{`needs_compatible_subset`}{indices trained on EPIC v1 need the
#'       v1/v2 compatible subset ([subset_versionshared_CpGs()])}
#'     \item{`probe_lists`}{named character vector of `data/` objects with
#'       probes to remove, by reason (`chrY`, `non_cpg`, `snp`, `zhou`);
#'       empty if none are available for the array}
#'   }
#'   `supported_arrays()`: a character vector of canonical names.
#' @examples
#' supported_arrays()
#' array_info("EPIC v2")$anno_pkg
#' array_info("450k")$samples_per_chip
#' @export
array_info <- function(array) {
  i <- resolve_array(array)
  row <- lapply(arrays[i, ], function(col) if (is.list(col)) col[[1]] else col)
  row$aliases <- NULL
  row
}

#' @rdname array_info
#' @export
supported_arrays <- function() {
  arrays$array
}

# Row index into `arrays` for a user-supplied name or a minfi object.
resolve_array <- function(array) {
  minfi_classes <- c("RGChannelSet", "MethylSet", "GenomicMethylSet",
                     "RatioSet", "GenomicRatioSet")
  is_minfi <- vapply(minfi_classes, function(cl) methods::is(array, cl),
                     logical(1))
  if (any(is_minfi)) {
    minfi_array <- unname(minfi::annotation(array)["array"])
    i <- match(minfi_array, arrays$minfi_array)
    if (is.na(i)) {
      supported <- paste(arrays$minfi_array, collapse = ", ")
      rlang::abort(c(
        sprintf("Unsupported array annotation \"%s\".", minfi_array),
        i = sprintf("Supported: %s.", supported)))
    }
    return(i)
  }
  if (!is.character(array) || length(array) != 1 || is.na(array)) {
    rlang::abort("`array` must be a single string or a minfi object.")
  }
  key <- normalise_array_name(array)
  hits <- vapply(arrays$aliases,
                 function(a) key %in% normalise_array_name(a), logical(1))
  if (!any(hits)) {
    rlang::abort(c(
      sprintf("Unknown array \"%s\".", array),
      i = sprintf("Supported: %s.", paste(arrays$array, collapse = ", "))))
  }
  which(hits)
}

normalise_array_name <- function(x) {
  tolower(gsub("[^A-Za-z0-9]", "", x))
}
