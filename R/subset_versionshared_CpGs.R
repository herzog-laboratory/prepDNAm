#' Subset a beta matrix to EPIC v1/v2 compatible probes
#'
#' Restricts a beta matrix to the CpG probes whose sequence is identical on
#' EPIC v1 and EPIC v2 (see [mapping_compatible_v1v2]) so that results from
#' both array versions can be combined.
#'
#' For an EPIC v2 matrix, probes are renamed to their v1 ID and v2 replicate
#' probes that map to the same v1 locus are averaged (mean of the non-missing
#' replicates; `NA` if all replicates are missing). For any other array the
#' matrix is subset to the loci present in the mapping; rownames are unchanged.
#'
#' @param beta Numeric matrix of beta values, probes in rows (rownames are
#'   probe IDs) and samples in columns.
#' @param array Array type; anything matching `"v2"` (case-insensitive) is
#'   treated as EPIC v2, everything else as v1-style IDs.
#' @return A numeric matrix with EPIC v1 probe IDs as rownames, restricted to
#'   the v1/v2 compatible probes. Rows are in order of first appearance in
#'   `beta`.
#' @examples
#' beta <- matrix(runif(6), nrow = 3, dimnames = list(
#'   c("cg20029347_TC11", "cg20029347_BC11", "cg00000029_TC21"), c("s1", "s2")))
#' subset_versionshared_CpGs(beta, array = "EPICv2")
#' @export
subset_versionshared_CpGs <- function(beta, array) {
  if (!is.matrix(beta) || !is.numeric(beta)) {
    rlang::abort("`beta` must be a numeric matrix.")
  }
  if (is.null(rownames(beta))) {
    rlang::abort("`beta` must have probe IDs as rownames.")
  }
  map <- mapping_compatible_v1v2

  if (grepl("v2", array, ignore.case = TRUE)) {
    v1_id <- map$EPICv1_Loci[match(rownames(beta), map$IlmnID)]
    keep <- !is.na(v1_id)
    beta <- beta[keep, , drop = FALSE]
    v1_id <- v1_id[keep]

    # mean over replicate v2 probes mapping to the same v1 locus
    n <- rowsum(1 - is.na(beta), group = v1_id, reorder = FALSE)
    beta <- rowsum(beta, group = v1_id, na.rm = TRUE, reorder = FALSE) / n
    beta[n == 0] <- NA_real_
  } else {
    beta <- beta[rownames(beta) %in% map$EPICv1_Loci, , drop = FALSE]
  }

  beta
}
