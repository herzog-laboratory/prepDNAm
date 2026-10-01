test_that("bundled probe lists are non-empty character vectors", {
  lists <- c("zhou_list", "zhou_list_v2", "snp_names", "snp_names_v2",
             "non_CpG_names", "non_CpG_names_v2", "chrY_names",
             "chrY_0_M_names_v2")
  for (nm in lists) {
    x <- get(nm, envir = asNamespace("prepDNAm"))
    expect_type(x, "character")
    expect_gt(length(x), 0)
    expect_false(anyNA(x))
  }
})

test_that("v1/v2 mapping has the expected columns", {
  m <- get("mapping_compatible_v1v2", envir = asNamespace("prepDNAm"))
  expect_s3_class(m, "data.frame")
  expect_true(all(c("IlmnID", "EPICv1_Loci") %in% names(m)))
  expect_gt(nrow(m), 0)
})
