sheet <- function(...) data.frame(basename = paste0("s", 1:4), ..., stringsAsFactors = FALSE)

test_that("sex values: F/M and spellings accepted, unknowns allowed", {
  df <- sheet(sex = c("F", "male", "Female", "unknown"))
  expect_identical(check_sample_metadata(df), df)
  expect_identical(normalise_sex(df$sex), c("F", "M", "F", NA))
  expect_no_error(check_sample_metadata(sheet(sex = c("", NA, "NA", "m"))))
})

test_that("unrecognised and numeric sex codes are an error with guidance", {
  expect_error(check_sample_metadata(sheet(sex = c("F", "M", "x", "F"))), "s3: \"x\"")
  expect_error(check_sample_metadata(sheet(sex = c(0, 1, 1, 0))), "Numeric codes")
  expect_error(check_sample_metadata(sheet(sex = c("1", "2", NA, "1"))), "Numeric codes")
})

test_that("age must be numeric years; text numbers are converted", {
  out <- check_sample_metadata(sheet(age = c("34.5", " 51", "", "unknown")))
  expect_identical(out$age, c(34.5, 51, NA, NA))
  expect_error(check_sample_metadata(sheet(age = c("45y", "30", "40-50", NA))), "s1: \"45y\"")
  expect_error(check_sample_metadata(sheet(age = c(-1, 30, 40, 50))), "negative")
  expect_warning(check_sample_metadata(sheet(age = c(30, 400, 40, 50))), "above 120")
})

test_that("alternative column names are used with a message", {
  expect_message(out <- check_sample_metadata(sheet(Gender = c("F", "M", "F", "M"))), "name it `sex`")
  expect_message(check_sample_metadata(sheet(gender = c("F", "M", "q", "M"))), "gender") |>
    expect_error("unrecognised")
  expect_message(check_sample_metadata(sheet(chron_age = c(1, 2, 3, 4))), "name it `age`")
  expect_identical(find_sex_col(sheet(sex = "F", gender = "M")), "sex")
  expect_identical(find_age_col(sheet(age = 1), "age"), "age")
  expect_error(find_age_col(sheet(age = 1), "nope"), "No column")
  expect_null(find_sex_col(sheet(x = 1)))
})

test_that("discover_samples validates the pheno sex/age columns", {
  skip_if_no_fixture_data()
  pheno <- fixture_pheno_df(); pheno$sex[2] <- "0"
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = file.path(fixture_root(), "run"),
                                    array = "450k", batch = "none", pheno = pheno))
  expect_error(discover_samples(cfg), "unrecognised value in the sex column")
  pheno <- fixture_pheno_df(); pheno$age <- as.character(pheno$age)
  cfg <- suppressMessages(qc_config(input = fixture_idat_450k(), run_dir = file.path(fixture_root(), "run"),
                                    array = "450k", batch = "none", pheno = pheno))
  expect_true(is.numeric(suppressMessages(discover_samples(cfg))$age))
})
