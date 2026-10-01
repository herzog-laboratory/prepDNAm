# Build the internal array lookup table (R/sysdata.rda).
#
# One row per supported array. Everything that is array- or species-specific
# in package code and report pages must come from this table; no literal
# probe counts, chip sizes, annotation package names or `array == "mouse"`
# branches anywhere else.
#
# Re-run after installing an annotation package or adding an array:
#   source("data-raw/arrays.R")

count_probes <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) return(NA_integer_)
  e <- new.env()
  utils::data(list = pkg, package = pkg, envir = e)
  m <- get(pkg, e)
  nrow(minfi::getProbeInfo(m, "I")) + nrow(minfi::getProbeInfo(m, "II"))
}

human_lists <- c(chrY = "chrY_names", non_cpg = "non_CpG_names",
                 snp = "snp_names", zhou = "zhou_list")
v2_lists <- c(chrY = "chrY_0_M_names_v2", non_cpg = "non_CpG_names_v2",
              snp = "snp_names_v2", zhou = "zhou_list_v2")

arrays <- tibble::tibble(
  array   = c("450k", "EPIC", "EPICv2", "mouse"),
  label   = c("Infinium HumanMethylation450", "Infinium MethylationEPIC v1.0",
              "Infinium MethylationEPIC v2.0", "Infinium Mouse Methylation"),
  aliases = list(
    c("450k", "hm450", "HumanMethylation450", "IlluminaHumanMethylation450k"),
    c("EPIC", "EPICv1", "EPIC v1", "850k", "hm850", "MethylationEPIC",
      "IlluminaHumanMethylationEPIC"),
    c("EPICv2", "EPIC v2", "v2", "935k", "EPIC2", "MethylationEPICv2",
      "IlluminaHumanMethylationEPICv2"),
    c("mouse", "MM285", "MouseMethylation", "IlluminaMouseMethylation")),
  species = c("human", "human", "human", "mouse"),
  genome  = c("hg19", "hg19", "hg38", "mm10"),
  # minfi annotation slot: c(array = minfi_array, annotation = minfi_annotation)
  minfi_array      = c("IlluminaHumanMethylation450k", "IlluminaHumanMethylationEPIC",
                       "IlluminaHumanMethylationEPICv2", "IlluminaMouseMethylation"),
  minfi_annotation = c("ilmn12.hg19", "ilm10b4.hg19", "20a1.hg38", "12.v1.mm10"),
  manifest_pkg = paste0(minfi_array, "manifest"),
  anno_pkg     = paste0(minfi_array, "anno.", minfi_annotation),
  # chip layout; Sentrix positions are R01C01..R0<rows>C0<cols>
  samples_per_chip = c(12L, 8L, 8L, 12L),
  # probe IDs carry a design suffix (cg00000029_TC21) rather than bare cg IDs
  probe_id_suffix = c(FALSE, FALSE, TRUE, TRUE),
  # rs probes readable with minfi::getSnpBeta(); used for genotype checks
  has_snp_probes = c(TRUE, TRUE, TRUE, FALSE),
  # minfi::getSex() (X/Y median intensities)
  has_sex_estimation = c(TRUE, TRUE, TRUE, FALSE),
  # EpiDISH immune-cell fractions, WID clocks, smoking index (human only)
  has_clocks = c(TRUE, TRUE, TRUE, FALSE),
  # clocks and downstream indices trained on v1 need the v1/v2 compatible subset
  needs_compatible_subset = c(FALSE, FALSE, TRUE, FALSE),
  # names of the data/ objects with probes to remove, by reason
  probe_lists = list(human_lists, human_lists, v2_lists, character(0))
)
arrays$n_probes <- unname(vapply(arrays$manifest_pkg, count_probes, integer(1)))

missing <- arrays$array[is.na(arrays$n_probes)]
if (length(missing) > 0) {
  message("n_probes left NA (manifest package not installed): ",
          paste(missing, collapse = ", "))
}

usethis::use_data(arrays, internal = TRUE, overwrite = TRUE)
