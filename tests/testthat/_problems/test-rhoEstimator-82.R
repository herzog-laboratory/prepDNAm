# Extracted from test-rhoEstimator.R:82

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "prepDNAm", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
simulate_beta <- function(rho, n = 20000, beta_e = 0.3309135,
                          beta_s1 = 0.019570, beta_s2 = 0.970650, sd = 0.02) {
  signal <- sample(c(beta_s1, beta_s2), n, replace = TRUE)
  obs <- (rho * signal + beta_e) / (1 + rho) + stats::rnorm(n, sd = sd)
  pmin(pmax(obs, 0), 1)
}

# test -------------------------------------------------------------------------
skip_if_not_installed("IlluminaHumanMethylationEPICv2manifest")
skip_if_not_installed("IlluminaHumanMethylationEPICv2anno.20a1.hg38")
suppressWarnings(suppressMessages(requireNamespace("minfi", quietly = TRUE)))
types <- manifest_types_for("IlluminaHumanMethylationEPICv2manifest")
ann <- minfi::getAnnotation(IlluminaHumanMethylationEPICv2anno.20a1.hg38::IlluminaHumanMethylationEPICv2anno.20a1.hg38,
                              what = "Manifest")
