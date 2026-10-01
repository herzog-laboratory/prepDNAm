# Legacy (eutopsQC 1.2.2) implementation, kept as reference for the v2 branch.
legacy_subset_v2 <- function(beta) {
  map <- prepDNAm::mapping_compatible_v1v2
  b <- beta |>
    as.data.frame() |>
    tibble::rownames_to_column(var = "IlmnID") |>
    dplyr::filter(IlmnID %in% map$IlmnID) |>
    dplyr::left_join(map, by = "IlmnID") |>
    dplyr::select(!IlmnID)
  dupes <- unique(janitor::get_dupes(b, EPICv1_Loci)$EPICv1_Loci)
  nd <- dplyr::filter(b, !EPICv1_Loci %in% dupes)
  d <- b |>
    dplyr::filter(EPICv1_Loci %in% dupes) |>
    dplyr::group_by(EPICv1_Loci) |>
    dplyr::reframe(dplyr::across(tidyselect::everything(), ~ mean(.x))) |>
    dplyr::ungroup()
  dplyr::bind_rows(d, nd) |>
    tibble::column_to_rownames("EPICv1_Loci") |>
    as.matrix()
}

# Small v2 fixture built from the real mapping: 2 loci with replicates,
# 3 single probes, 2 probes absent from the mapping.
make_v2_fixture <- function() {
  map <- prepDNAm::mapping_compatible_v1v2
  reps <- names(which(table(map$EPICv1_Loci) > 1))[1:2]
  singles <- setdiff(unique(map$EPICv1_Loci), reps)[1:3]
  ids <- c(map$IlmnID[map$EPICv1_Loci %in% reps],
           map$IlmnID[match(singles, map$EPICv1_Loci)],
           "cg99999999_TC11", "rs0000001")
  set.seed(1)
  beta <- matrix(runif(length(ids) * 3), nrow = length(ids),
                 dimnames = list(ids, paste0("s", 1:3)))
  list(beta = beta, reps = reps, singles = singles)
}

test_that("v2: keeps only mapped probes, renames to v1 IDs, averages replicates", {
  fx <- make_v2_fixture()
  out <- subset_versionshared_CpGs(fx$beta, array = "EPICv2")

  expect_true(is.matrix(out))
  expect_equal(colnames(out), colnames(fx$beta))
  expect_setequal(rownames(out), c(fx$reps, fx$singles))
  expect_false(any(duplicated(rownames(out))))

  map <- prepDNAm::mapping_compatible_v1v2
  for (r in fx$reps) {
    ids <- map$IlmnID[map$EPICv1_Loci == r]
    expect_equal(out[r, ], colMeans(fx$beta[ids, , drop = FALSE]))
  }
  for (s in fx$singles) {
    id <- map$IlmnID[match(s, map$EPICv1_Loci)]
    expect_equal(out[s, ], fx$beta[id, ])
  }
})

test_that("v2: matches the legacy implementation", {
  fx <- make_v2_fixture()
  new <- subset_versionshared_CpGs(fx$beta, array = "EPICv2")
  old <- legacy_subset_v2(fx$beta)
  expect_equal(new[sort(rownames(new)), ], old[sort(rownames(old)), ])
})

test_that("v2: NA handling averages available replicates", {
  fx <- make_v2_fixture()
  map <- prepDNAm::mapping_compatible_v1v2
  ids <- map$IlmnID[map$EPICv1_Loci == fx$reps[1]]
  beta <- fx$beta
  beta[ids[1], "s1"] <- NA          # one replicate missing
  beta[ids, "s2"] <- NA             # all replicates missing
  out <- subset_versionshared_CpGs(beta, array = "EPICv2")
  expect_equal(out[fx$reps[1], "s1"], mean(beta[ids[-1], "s1"]))
  expect_true(is.na(out[fx$reps[1], "s2"]))
  expect_equal(out[fx$reps[1], "s3"], mean(beta[ids, "s3"]))
})

test_that("v2: works with a single sample", {
  fx <- make_v2_fixture()
  out <- subset_versionshared_CpGs(fx$beta[, 1, drop = FALSE], array = "EPICv2")
  expect_equal(dim(out), c(length(c(fx$reps, fx$singles)), 1L))
})

test_that("v1: subsets to compatible loci without renaming", {
  map <- prepDNAm::mapping_compatible_v1v2
  ids <- c(unique(map$EPICv1_Loci)[1:4], "cg99999999", "rs0000001")
  beta <- matrix(runif(length(ids) * 2), nrow = length(ids),
                 dimnames = list(ids, c("a", "b")))
  out <- subset_versionshared_CpGs(beta, array = "EPIC")
  expect_equal(rownames(out), ids[1:4])
  expect_equal(out, beta[1:4, ])
})

test_that("input validation", {
  expect_error(subset_versionshared_CpGs(data.frame(a = 1), "EPICv2"), "numeric matrix")
  expect_error(subset_versionshared_CpGs(matrix(letters[1:4], 2), "EPICv2"), "numeric matrix")
  expect_error(subset_versionshared_CpGs(matrix(runif(4), 2), "EPICv2"), "rownames")
})
