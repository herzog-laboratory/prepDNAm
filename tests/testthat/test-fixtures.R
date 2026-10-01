# Smoke tests for the fixture helpers themselves.

test_that("450k IDAT fixtures exist in every layout", {
  skip_if_no_fixture_data()
  nested <- list.files(fixture_idat_450k("nested"), pattern = "_Grn\\.idat$", recursive = TRUE)
  flat <- list.files(fixture_idat_450k("flat"), pattern = "_Grn\\.idat$")
  gz <- list.files(fixture_idat_450k("gz"), pattern = "_Grn\\.idat\\.gz$")
  expect_length(nested, 6)
  expect_length(flat, 6)
  expect_length(gz, 6)
  expect_setequal(basename(dirname(list.files(fixture_idat_450k("nested"), pattern = "idat$", recursive = TRUE, full.names = TRUE))),
                  c("5723646052", "5723646053"))
  # gzipped copies read back identically
  f <- list.files(fixture_idat_450k("flat"), pattern = "_Grn\\.idat$", full.names = TRUE)[1]
  g <- file.path(fixture_idat_450k("gz"), paste0(basename(f), ".gz"))
  con <- gzfile(g, "rb")
  on.exit(close(con))
  expect_identical(readBin(con, "raw", file.size(f) + 1), readBin(f, "raw", file.size(f)))
})

test_that("EPIC IDAT fixture exists", {
  skip_if_no_fixture_data(epic = TRUE)
  expect_length(list.files(fixture_idat_epic(), pattern = "_Grn\\.idat$", recursive = TRUE), 3)
})

test_that("pheno fixtures round-trip in every format", {
  skip_if_no_fixture_data()
  df <- fixture_pheno_df()
  expect_named(df, c("basename", "plate", "sex", "age", "group", "person"))
  expect_identical(df$basename, fixture_basenames_450k())
  expect_identical(sort(unique(df$plate)), c("5723646052", "5723646053"))

  expect_identical(utils::read.csv(fixture_pheno("csv"), stringsAsFactors = FALSE)$basename, df$basename)
  expect_identical(utils::read.delim(fixture_pheno("tsv"), stringsAsFactors = FALSE)$basename, df$basename)
  expect_identical(utils::read.delim(fixture_pheno("txt"), stringsAsFactors = FALSE)$basename, df$basename)
  expect_identical(readRDS(fixture_pheno("rds")), df)
  e <- new.env(); load(fixture_pheno("Rdata"), envir = e)
  expect_identical(e$pheno, df)
  skip_if_not_installed("readxl")
  skip_if_not_installed("writexl")
  expect_identical(readxl::read_xlsx(fixture_pheno("xlsx"))$basename, df$basename)
})

test_that("pheno variants differ as intended", {
  skip_if_no_fixture_data()
  ok <- fixture_pheno_df("ok")
  expect_identical(nrow(fixture_pheno_df("missing")), nrow(ok) + 1L)
  expect_false(fixture_pheno_df("missing")$basename[7] %in% ok$basename)
  expect_identical(nrow(fixture_pheno_df("extra")), nrow(ok) - 1L)
  expect_false("basename" %in% names(fixture_pheno_df("nobasename")))
  expect_true(readLines(fixture_exclude()) %in% ok$basename)
})

test_that("signal matrix fixture has GEO-style columns for every sample", {
  skip_if_no_fixture_data()
  hdr <- strsplit(readLines(fixture_sigmatrix(), n = 1), "\t")[[1]]
  expect_identical(hdr[1], "ID_REF")
  expect_length(hdr, 1 + 3 * 6)
  for (b in fixture_basenames_450k()) {
    expect_true(all(paste(b, c("Unmethylated signal", "Methylated signal", "Detection Pval")) %in% hdr))
  }
  first <- strsplit(readLines(fixture_sigmatrix(), n = 2)[2], "\t")[[1]]
  expect_match(first[1], "^cg")
  expect_false(anyNA(as.numeric(first[-1])))
})

test_that("beta fixture has NAs and both file formats", {
  skip_if_no_fixture_data()
  fx <- fixture_beta()
  expect_identical(dim(fx$beta), c(485512L, 6L))
  expect_identical(colnames(fx$beta), fixture_basenames_450k())
  expect_true(anyNA(fx$beta))
  expect_true(all(fx$beta >= 0 & fx$beta <= 1, na.rm = TRUE))
  expect_identical(readRDS(fx$rds), fx$beta)
  csv <- utils::read.csv(fx$csv, row.names = 1, check.names = FALSE)
  expect_identical(dim(csv), dim(fx$beta))
  expect_identical(colnames(csv), colnames(fx$beta))
  expect_equal(as.matrix(csv)[1:100, ], fx$beta[1:100, ], tolerance = 1e-12)
})
