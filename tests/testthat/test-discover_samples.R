cfg_450k <- function(...) {
  suppressMessages(qc_config(
    input = fixture_idat_450k(), run_dir = file.path(fixture_root(), "run"),
    array = "450k", ...))
}

test_that("IDATs are discovered with minfi basenames, ordered by basename", {
  skip_if_no_fixture_data()
  s <- discover_samples(cfg_450k(batch = "none"))
  expect_named(s, c("basename", "batch", "exclude", "path"))
  expect_identical(s$basename, sort(fixture_basenames_450k()))
  expect_true(all(file.exists(paste0(s$path, "_Grn.idat"))))
  expect_true(all(file.exists(paste0(s$path, "_Red.idat"))))
  expect_identical(basename(s$path), s$basename)
  expect_identical(unique(s$batch), "all")
  expect_false(any(s$exclude))
})

test_that("flat and gzipped layouts are discovered too", {
  skip_if_no_fixture_data()
  for (layout in c("flat", "gz")) {
    cfg <- suppressMessages(qc_config(input = fixture_idat_450k(layout), run_dir = fixture_root(),
                                      array = "450k", batch = "none"))
    s <- discover_samples(cfg)
    expect_identical(s$basename, sort(fixture_basenames_450k()), info = layout)
    expect_identical(dirname(s$path), rep(normalizePath(fixture_idat_450k(layout)), 6), info = layout)
  }
})

test_that("a _Grn without _Red and duplicated basenames are errors", {
  skip_if_no_fixture_data()
  dir <- withr::local_tempdir()
  src <- list.files(fixture_idat_450k("flat"), full.names = TRUE)
  file.copy(src, dir)
  file.remove(file.path(dir, "5723646052_R02C02_Red.idat"))
  cfg <- suppressMessages(qc_config(input = dir, run_dir = fixture_root(), array = "450k", batch = "none"))
  expect_error(discover_samples(cfg), "no _Red IDAT")

  dir2 <- withr::local_tempdir()
  dir.create(file.path(dir2, "a")); dir.create(file.path(dir2, "b"))
  file.copy(src, file.path(dir2, "a")); file.copy(src, file.path(dir2, "b"))
  cfg <- suppressMessages(qc_config(input = dir2, run_dir = fixture_root(), array = "450k", batch = "dir"))
  expect_error(discover_samples(cfg), "occur more than once")
})

test_that("batch modes: dir, chip, pheno column", {
  skip_if_no_fixture_data()
  by_dir <- discover_samples(cfg_450k(batch = "dir"))
  expect_identical(by_dir$batch, sub("_R..C..$", "", by_dir$basename))
  by_chip <- discover_samples(cfg_450k(batch = "chip"))
  expect_identical(by_chip$batch, by_dir$batch)
  by_col <- discover_samples(cfg_450k(batch = "plate", pheno = fixture_pheno_df()))
  expect_identical(by_col$batch, by_col$plate)
  # a pheno column named `batch` is used as such
  pheno <- fixture_pheno_df(); pheno$batch <- rep(c("x", "y"), 3)
  s <- discover_samples(cfg_450k(batch = "batch", pheno = pheno))
  expect_identical(s$batch, rep(c("x", "y"), 3))
  expect_false("pheno_batch" %in% names(s))
  # ... and kept as pheno_batch when another mode is chosen
  expect_message(s <- discover_samples(cfg_450k(batch = "chip", pheno = pheno)), "pheno_batch")
  expect_identical(s$pheno_batch, rep(c("x", "y"), 3))
  expect_identical(s$batch, sub("_R..C..$", "", s$basename))
})

test_that("batch = 'dir' fails when IDATs sit directly in input", {
  skip_if_no_fixture_data()
  dir <- withr::local_tempdir()
  file.copy(list.files(fixture_idat_450k("flat"), full.names = TRUE), dir)
  dir.create(file.path(dir, "sub"))
  cfg <- suppressMessages(qc_config(input = dir, run_dir = fixture_root(), array = "450k", batch = "dir"))
  expect_error(discover_samples(cfg), "some are directly in it")
})

test_that("batch = 'chip' needs Sentrix-style basenames", {
  skip_if_no_fixture_data()
  cfg <- suppressMessages(qc_config(input = fixture_beta()$rds, run_dir = fixture_root(),
                                    array = "450k", batch = "chip"))
  expect_identical(unique(discover_samples(cfg)$batch), c("5723646052", "5723646053"))
  beta <- fixture_beta()$beta[1:10, ]
  colnames(beta) <- paste0("sample", 1:6)
  f <- withr::local_tempfile(fileext = ".rds"); saveRDS(beta, f)
  cfg <- suppressMessages(qc_config(input = f, run_dir = fixture_root(), array = "450k", batch = "chip"))
  expect_error(discover_samples(cfg), "SentrixID")
})

test_that("batch labels are made safe for directory names", {
  skip_if_no_fixture_data()
  pheno <- fixture_pheno_df(); pheno$site <- rep(c("Site A/1", "Site B (2)"), 3)
  expect_message(s <- discover_samples(cfg_450k(batch = "site", pheno = pheno)), "not safe")
  expect_identical(unique(s$batch), c("Site_A_1", "Site_B_2_"))
  pheno$site[1] <- NA
  expect_error(discover_samples(cfg_450k(batch = "site", pheno = pheno)), "missing values")
})

test_that("pheno is joined by exact basename in pheno order, every format", {
  skip_if_no_fixture_data()
  for (fmt in c("csv", "tsv", "txt", "rds", "Rdata")) {
    s <- discover_samples(cfg_450k(batch = "none", pheno = fixture_pheno(fmt)))
    expect_identical(s$basename, fixture_basenames_450k(), info = fmt)
    expect_named(s, c("basename", "batch", "exclude", "path", "plate", "sex", "age", "group", "person"), info = fmt)
    expect_identical(s$age, fixture_pheno_df()$age, info = fmt)
  }
  skip_if_not_installed("writexl")
  s <- discover_samples(cfg_450k(batch = "none", pheno = fixture_pheno("xlsx")))
  expect_identical(s$basename, fixture_basenames_450k())
})

test_that("pheno sample not on disk is an error; extra IDATs are dropped with a message", {
  skip_if_no_fixture_data()
  expect_error(discover_samples(cfg_450k(batch = "none", pheno = fixture_pheno_df("missing"))),
               "not found in `input`")
  expect_error(discover_samples(cfg_450k(batch = "none", pheno = fixture_pheno_df("missing"))),
               "5723646053_R01C01")
  expect_message(s <- discover_samples(cfg_450k(batch = "none", pheno = fixture_pheno_df("extra"))),
                 "will be dropped: 5723646053_R06C02")
  expect_identical(nrow(s), 5L)
})

test_that("pheno validation: Basename alias, duplicates, path/dir renamed", {
  skip_if_no_fixture_data()
  pheno <- fixture_pheno_df(); names(pheno)[1] <- "Basename"
  expect_identical(discover_samples(cfg_450k(batch = "none", pheno = pheno))$basename,
                   fixture_basenames_450k())
  dup <- rbind(pheno, pheno[1, ])
  expect_error(cfg_450k(batch = "none", pheno = dup), "duplicated basenames")
  pheno <- fixture_pheno_df(); pheno$path <- "x"; pheno$dir <- "y"
  s <- discover_samples(cfg_450k(batch = "none", pheno = pheno))
  expect_true(all(c("pheno_path", "pheno_dir") %in% names(s)))
  expect_identical(basename(s$path), s$basename)
})

test_that("exclusions come from a vector, a file, or a pheno column", {
  skip_if_no_fixture_data()
  b <- fixture_basenames_450k()
  s <- discover_samples(cfg_450k(batch = "none", exclude = b[1:2]))
  expect_identical(s$basename[s$exclude], sort(b[1:2]))
  s <- discover_samples(cfg_450k(batch = "none", exclude = fixture_exclude()))
  expect_identical(s$basename[s$exclude], b[2])
  pheno <- fixture_pheno_df(); pheno$exclude <- c(TRUE, FALSE, FALSE, FALSE, FALSE, NA)
  s <- discover_samples(cfg_450k(batch = "none", pheno = pheno))
  expect_identical(s$basename[s$exclude], b[1])
  pheno$exclude <- c("yes", "no", "no", "no", "no", "TRUE")
  s <- discover_samples(cfg_450k(batch = "none", pheno = pheno, exclude = b[3]))
  expect_identical(s$basename[s$exclude], b[c(1, 3, 6)])
  expect_warning(discover_samples(cfg_450k(batch = "none", exclude = "nope")),
                 "not among the samples")
})

test_that("signal-matrix samples are parsed from GEO-style headers", {
  skip_if_no_fixture_data()
  cfg <- suppressMessages(qc_config(input = fixture_sigmatrix(), run_dir = fixture_root(),
                                    array = "450k"))
  s <- discover_samples(cfg)
  expect_identical(s$basename, sort(fixture_basenames_450k()))
  expect_identical(unique(s$path), normalizePath(fixture_sigmatrix()))
  expect_identical(unique(s$batch), "all")
  # wrong pattern
  cfg <- suppressMessages(qc_config(input = fixture_sigmatrix(), run_dir = fixture_root(),
                                    array = "450k", input_type = "sigmatrix",
                                    unmeth_pattern = "Unmeth Signal"))
  expect_error(discover_samples(cfg), "No columns matching `unmeth_pattern`")
})

test_that("header_samples strips the pattern and separators", {
  expect_identical(header_samples(c("ID_REF", "GSM1 Unmethylated signal", "GSM1 Methylated signal",
                                    "GSM2_Unmethylated signal", "GSM3.Unmethylated signal"),
                                  "Unmethylated signal"),
                   c("GSM1", "GSM2", "GSM3"))
  expect_identical(header_samples(c("Unmethylated signal GSM1"), "Unmethylated signal"), "GSM1")
})

test_that("mismatched signal columns are an error", {
  f <- withr::local_tempfile(fileext = ".tsv")
  writeLines(c("ID_REF\ts1 Unmethylated signal\ts1 Methylated signal\ts2 Unmethylated signal",
               "cg1\t1\t2\t3"), f)
  cfg <- suppressMessages(qc_config(input = f, run_dir = fixture_root(), array = "450k"))
  expect_error(discover_samples(cfg), "do not name the same samples")
})

test_that("beta-matrix samples come from the column names (csv and rds)", {
  skip_if_no_fixture_data()
  for (path in c(fixture_beta()$csv, fixture_beta()$rds)) {
    cfg <- suppressMessages(qc_config(input = path, run_dir = fixture_root(), array = "450k"))
    s <- discover_samples(cfg)
    expect_identical(s$basename, sort(fixture_basenames_450k()), info = path)
    expect_identical(unique(s$path), normalizePath(path), info = path)
  }
  s <- discover_samples(suppressMessages(qc_config(input = fixture_beta()$rds, run_dir = fixture_root(),
                                                   array = "450k", pheno = fixture_pheno_df(), batch = "plate")))
  expect_identical(s$basename, fixture_basenames_450k())
  expect_identical(s$batch, s$plate)
})

test_that("discover_samples() needs a config", {
  expect_error(discover_samples(list()), "must be a `prepDNAm_config`")
})
