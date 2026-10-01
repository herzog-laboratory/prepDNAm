# Datasets in data/ and the internal lookup table in R/sysdata.rda are loaded
# into the namespace; declare them so R CMD check / lintr do not flag their use
# in package code.
utils::globalVariables(c("mapping_compatible_v1v2", "arrays"))

#' @importFrom rlang .data
NULL
