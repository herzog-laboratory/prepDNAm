# A valid config for the 450k IDAT fixture with all messages silenced.
cfg_450k <- function(...) {
  suppressMessages(qc_config(
    input = fixture_idat_450k(), run_dir = file.path(fixture_root(), "run"),
    array = "450k", ...))
}

test_that("a minimal IDAT config is valid, normalised and printable", {
  skip_if_no_fixture_data()
  expect_message(
    cfg <- qc_config(input = fixture_idat_450k(),
                     run_dir = file.path(fixture_root(), "run"), array = "450K"),
    "`batch` not set; using \"dir\"")
  expect_s3_class(cfg, "prepDNAm_config")
  expect_identical(cfg$input_type, "idat")
  expect_identical(cfg$array, "450k")
  expect_identical(cfg$batch, "dir")
  expect_identical(cfg$run_name, "extdata")
  expect_identical(cfg$run_dir, file.path(normalizePath(fixture_root()), "run"))
  expect_false(cfg$compatible_subset)
  expect_true(cfg$save_sex)
  expect_true(cfg$save_rs)
  expect_false(cfg$export_csv)
  expect_true(cfg$estimates)
  expect_identical(cfg$export_digits, 3)
  expect_identical(cfg$cores, 1L)
  expect_identical(cfg$seed, 1L)
  expect_identical(cfg$prepDNAm_version, as.character(utils::packageVersion("prepDNAm")))
  expect_output(print(cfg), "<prepDNAm_config>")
  expect_output(print(cfg), "input_type +idat")
})

test_that("batch default is 'none' for flat IDAT dirs and matrix input", {
  skip_if_no_fixture_data()
  expect_message(
    cfg <- qc_config(input = fixture_idat_450k("flat"),
                     run_dir = file.path(fixture_root(), "run"), array = "450k"),
    "using \"none\"")
  expect_identical(cfg$batch, "none")
  cfg <- suppressMessages(qc_config(input = fixture_beta()$rds,
                                    run_dir = file.path(fixture_root(), "run"), array = "450k"))
  expect_identical(cfg$batch, "none")
  # save_rs defaults to off, silently, where there are no SNP probes
  expect_false(cfg$save_rs)
  expect_no_message(qc_config(input = fixture_beta()$rds, run_dir = file.path(fixture_root(), "run"),
                              array = "450k", batch = "none", save_sex = FALSE), message = "save_rs")
})

test_that("batch defaults to a pheno 'batch' column if present", {
  skip_if_no_fixture_data()
  pheno <- fixture_pheno_df()
  names(pheno)[names(pheno) == "plate"] <- "batch"
  expect_message(
    cfg <- qc_config(input = fixture_idat_450k(), run_dir = file.path(fixture_root(), "run"),
                     array = "450k", pheno = pheno),
    "using \"batch\"")
  expect_identical(cfg$batch, "batch")
})

test_that("input type is inferred from the input", {
  skip_if_no_fixture_data()
  expect_identical(cfg_450k(batch = "none")$input_type, "idat")
  expect_identical(suppressMessages(qc_config(
    input = fixture_sigmatrix(), run_dir = file.path(fixture_root(), "run"),
    array = "450k"))$input_type, "sigmatrix")
  expect_identical(suppressMessages(qc_config(
    input = fixture_beta()$csv, run_dir = file.path(fixture_root(), "run"),
    array = "450k"))$input_type, "beta")
  expect_identical(suppressMessages(qc_config(
    input = fixture_beta()$rds, run_dir = file.path(fixture_root(), "run"),
    array = "450k"))$input_type, "beta")
  # run_name for a file input is the file name without extension
  expect_identical(suppressMessages(qc_config(
    input = fixture_beta()$rds, run_dir = file.path(fixture_root(), "run"),
    array = "450k"))$run_name, "beta")
})

test_that("explicit input_type must match the input", {
  skip_if_no_fixture_data()
  expect_error(cfg_450k(input_type = "beta"), "needs `input` to be a file")
  expect_error(suppressMessages(qc_config(
    input = fixture_beta()$rds, run_dir = file.path(fixture_root(), "run"),
    array = "450k", input_type = "idat")), "needs `input` to be a directory")
  expect_error(cfg_450k(input_type = "foo"), "must be one of")
})

test_that("paths are validated", {
  skip_if_no_fixture_data()
  expect_error(qc_config(input = "/nonexistent/idats", run_dir = fixture_root(), array = "450k"),
               "`input` does not exist")
  expect_error(qc_config(input = fixture_idat_450k(), run_dir = "/nonexistent/parent/run",
                         array = "450k", batch = "none"),
               "Parent directory of `run_dir` does not exist")
  expect_error(qc_config(input = fixture_idat_450k(), run_dir = fixture_root(), array = "27k"),
               "Unknown array")
})

test_that("pheno is validated at config time", {
  skip_if_no_fixture_data()
  expect_error(cfg_450k(pheno = "/nonexistent/pheno.csv"), "`pheno` file does not exist")
  expect_error(cfg_450k(pheno = fixture_pheno_df("nobasename")), "no `basename` column")
  expect_error(cfg_450k(pheno = 1:3), "must be a single non-empty string")
  cfg <- cfg_450k(pheno = fixture_pheno("csv"), batch = "plate")
  expect_identical(cfg$batch, "plate")
  expect_identical(cfg$pheno, normalizePath(fixture_pheno("csv")))
  cfg <- cfg_450k(pheno = fixture_pheno_df(), batch = "plate")
  expect_true(is.data.frame(cfg$pheno))
})

test_that("batch mode is validated against input type and pheno", {
  skip_if_no_fixture_data()
  expect_error(suppressMessages(qc_config(
    input = fixture_beta()$rds, run_dir = file.path(fixture_root(), "run"),
    array = "450k", batch = "dir")), "only available for IDAT input")
  expect_error(cfg_450k(batch = "plate"), "no `pheno` was given")
  expect_error(cfg_450k(batch = "nosuchcolumn", pheno = fixture_pheno_df()),
               "neither a batching mode nor a `pheno` column")
  expect_identical(cfg_450k(batch = "chip")$batch, "chip")
})

test_that("array-dependent features are switched off with a message", {
  skip_if_no_fixture_data()
  expect_true(cfg_450k(batch = "none", save_rs = TRUE)$save_rs)
  expect_false(cfg_450k(batch = "none", save_rs = FALSE)$save_rs)
  # mouse: no rs probes, no sex, no probe lists
  msgs <- character()
  withCallingHandlers(
    qc_config(input = fixture_idat_450k(), run_dir = fixture_root(), array = "mouse",
              batch = "none", save_rs = TRUE),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
  expect_match(msgs, "`save_rs` ignored", all = FALSE)
  expect_match(msgs, "`save_sex` ignored", all = FALSE)
  expect_match(msgs, "No probe-removal lists", all = FALSE)
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = fixture_root(),
                                    array = "mouse", batch = "none", save_rs = TRUE))
  expect_false(cfg$save_rs)
  expect_false(cfg$save_sex)
  expect_false(cfg$remove_chrY)
  expect_false(cfg$remove_zhou)
  # beta input: no rs / sex either
  cfg <- suppressMessages(qc_config(input = fixture_beta()$rds, run_dir = fixture_root(),
                                    array = "EPIC", save_rs = TRUE, save_sex = TRUE))
  expect_false(cfg$save_rs)
  expect_false(cfg$save_sex)
})

test_that("compatible_subset defaults per array and is EPICv2-only", {
  skip_if_no_fixture_data()
  expect_false(cfg_450k(batch = "none")$compatible_subset)
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = fixture_root(),
                                    array = "EPICv2", batch = "none"))
  expect_true(cfg$compatible_subset)
  expect_error(cfg_450k(batch = "none", compatible_subset = TRUE), "only applies to EPIC v2")
})

test_that("thresholds, flags and compute settings are validated", {
  skip_if_no_fixture_data()
  expect_error(cfg_450k(batch = "none", detp_threshold = 2), "in \\[0, 1\\]")
  expect_error(cfg_450k(batch = "none", intensity_threshold = -1), ">= 0")
  expect_error(cfg_450k(batch = "none", impute = "yes"), "must be TRUE or FALSE")
  expect_error(cfg_450k(batch = "none", normalise = "maybe"), "\"auto\", TRUE or FALSE")
  expect_identical(cfg_450k(batch = "none", normalise = FALSE)$normalise, FALSE)
  expect_error(cfg_450k(batch = "none", cores = 0), ">= 1")
  expect_message(cfg <- qc_config(input = fixture_idat_450k(), run_dir = fixture_root(),
                                  array = "450k", batch = "none", cores = 10000),
                 "capped at")
  expect_lte(cfg$cores, parallelly::availableCores())
  expect_identical(cfg_450k(batch = "none", seed = 42)$seed, 42L)
})

test_that("exclude accepts a vector or a file", {
  skip_if_no_fixture_data()
  expect_identical(cfg_450k(batch = "none", exclude = c("a", "b"))$exclude, c("a", "b"))
  expect_identical(cfg_450k(batch = "none", exclude = fixture_exclude())$exclude,
                   normalizePath(fixture_exclude()))
  expect_error(cfg_450k(batch = "none", exclude = 1), "character vector of basenames")
})
