test_that("the lookup table has one complete row per supported array", {
  expect_identical(supported_arrays(), c("450k", "EPIC", "EPICv2", "mouse"))
  expect_false(anyDuplicated(arrays$array) > 0)
  # every field is filled, except n_probes which may be NA if the manifest was
  # not installed when the table was built
  for (col in setdiff(names(arrays), c("n_probes", "aliases", "probe_lists"))) {
    expect_false(anyNA(arrays[[col]]), info = col)
  }
  expect_true(all(lengths(arrays$aliases) > 0))
  # no alias resolves to more than one array
  alias_map <- unique(data.frame(
    alias = unlist(lapply(arrays$aliases, normalise_array_name)),
    array = rep(arrays$array, lengths(arrays$aliases))))
  expect_false(anyDuplicated(alias_map$alias) > 0)
})

test_that("probe lists point at existing data objects with the right reasons", {
  for (i in seq_len(nrow(arrays))) {
    lists <- arrays$probe_lists[[i]]
    if (length(lists) == 0) next
    expect_named(lists, c("chrY", "non_cpg", "snp", "zhou"))
    for (obj in lists) {
      expect_true(exists(obj, envir = asNamespace("prepDNAm")), info = obj)
    }
  }
  # v2 lists are the v2 versions, human ones are not
  expect_true(all(grepl("_v2$", array_info("EPICv2")$probe_lists)))
  expect_false(any(grepl("_v2$", array_info("EPIC")$probe_lists)))
  expect_identical(array_info("mouse")$probe_lists, character(0))
})

test_that("known constants are right", {
  expect_identical(array_info("450k")$n_probes, 485512L)
  expect_identical(array_info("EPIC")$n_probes, 866836L)
  expect_identical(array_info("EPICv2")$n_probes, 936990L)
  expect_identical(array_info("mouse")$n_probes, 285698L)
  expect_identical(array_info("450k")$samples_per_chip, 12L)
  expect_identical(array_info("EPIC")$samples_per_chip, 8L)
  expect_identical(array_info("EPICv2")$anno_pkg, "IlluminaHumanMethylationEPICv2anno.20a1.hg38")
  expect_identical(array_info("EPIC")$anno_pkg, "IlluminaHumanMethylationEPICanno.ilm10b4.hg19")
  expect_identical(array_info("450k")$anno_pkg, "IlluminaHumanMethylation450kanno.ilmn12.hg19")
  expect_true(array_info("EPICv2")$needs_compatible_subset)
  expect_false(array_info("mouse")$has_clocks)
  expect_false(array_info("mouse")$has_snp_probes)
  expect_identical(array_info("mouse")$species, "mouse")
})

test_that("array_info() returns a flat named list", {
  info <- array_info("EPIC")
  expect_type(info, "list")
  expect_null(info$aliases)
  expect_type(info$probe_lists, "character")
  expect_identical(info$array, "EPIC")
})

test_that("names resolve case-insensitively and via aliases", {
  expect_identical(array_info("epic")$array, "EPIC")
  expect_identical(array_info("EPIC v1")$array, "EPIC")
  expect_identical(array_info("850k")$array, "EPIC")
  expect_identical(array_info("EPIC v2")$array, "EPICv2")
  expect_identical(array_info("epicv2")$array, "EPICv2")
  expect_identical(array_info("v2")$array, "EPICv2")
  expect_identical(array_info("IlluminaHumanMethylationEPICv2")$array, "EPICv2")
  expect_identical(array_info("450K")$array, "450k")
  expect_identical(array_info("HumanMethylation450")$array, "450k")
  expect_identical(array_info("Mouse")$array, "mouse")
  expect_identical(array_info("MM285")$array, "mouse")
})

test_that("unknown or malformed names give an informative error", {
  expect_error(array_info("27k"), "Unknown array \"27k\"")
  expect_error(array_info("27k"), "Supported: 450k, EPIC, EPICv2, mouse")
  expect_error(array_info(c("EPIC", "450k")), "single string")
  expect_error(array_info(NA_character_), "single string")
  expect_error(array_info(450), "single string")
})

test_that("minfi objects resolve from their annotation slot", {
  skip_if_not_installed("minfiData")
  RGsetEx <- MsetEx <- NULL
  utils::data("RGsetEx", package = "minfiData", envir = environment())
  utils::data("MsetEx", package = "minfiData", envir = environment())
  expect_identical(array_info(RGsetEx)$array, "450k")
  expect_identical(array_info(MsetEx)$array, "450k")
  # the table's minfi fields match what minfi itself records
  expect_identical(unname(minfi::annotation(RGsetEx)),
                   c(array_info("450k")$minfi_array, array_info("450k")$minfi_annotation))
})

test_that("EPIC minfi objects resolve too", {
  skip_if_not_installed("minfiDataEPIC")
  RGsetEPIC <- NULL
  utils::data("RGsetEPIC", package = "minfiDataEPIC", envir = environment())
  expect_identical(array_info(RGsetEPIC)$array, "EPIC")
})
