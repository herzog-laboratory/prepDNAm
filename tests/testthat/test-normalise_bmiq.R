# Simulate a sample from BMIQ's own model: a three-state (U/H/M) beta mixture
# for each probe type, with the type II U and M components pulled towards the
# middle as on a real array. BMIQ should move them onto the type I components.
simulate_biased_beta <- function(n1 = 20000, n2 = 60000) {
  w <- c(0.45, 0.1, 0.45)
  draw <- function(n, shapes) {
    state <- sample(1:3, n, replace = TRUE, prob = w)
    stats::rbeta(n, shapes[state, 1], shapes[state, 2])
  }
  shapes_I <- rbind(c(2, 40), c(10, 10), c(40, 2))   # peaks ~0.025, 0.5, 0.975
  shapes_II <- rbind(c(4, 20), c(10, 10), c(20, 4))  # peaks ~0.14, 0.5, 0.86
  list(beta = c(draw(n1, shapes_I), draw(n2, shapes_II)),
       probe_type = rep(c("I", "II"), c(n1, n2)))
}

mode_of <- function(x) {
  d <- stats::density(x, bw = 0.01, from = 0, to = 1, na.rm = TRUE)
  d$x[which.max(d$y)]
}

test_that("type II peaks are moved onto the type I peaks", {
  skip_if_not_installed("ChAMP")
  set.seed(7)
  sim <- simulate_biased_beta()
  beta <- cbind(s1 = sim$beta, s2 = sim$beta)
  is_II <- sim$probe_type == "II"
  lo <- function(b) b < 0.3
  hi <- function(b) b > 0.7

  # sanity: before normalisation the peaks differ by ~0.1
  expect_gt(mode_of(beta[is_II & lo(beta[, 1]), 1]) -
              mode_of(beta[!is_II & lo(beta[, 1]), 1]), 0.08)

  res <- normalise_bmiq(beta, probe_type = sim$probe_type)

  expect_identical(dim(res), dim(beta))
  expect_identical(dimnames(res), dimnames(beta))
  expect_true(all(res >= 0 & res <= 1))
  # type I probes are untouched
  expect_equal(res[!is_II, ], beta[!is_II, ])
  # type II low and high peaks now coincide with the type I peaks
  expect_lt(abs(mode_of(res[is_II & lo(res[, 1]), 1]) -
                  mode_of(res[!is_II & lo(res[, 1]), 1])), 0.01)
  expect_lt(abs(mode_of(res[is_II & hi(res[, 1]), 1]) -
                  mode_of(res[!is_II & hi(res[, 1]), 1])), 0.01)
  # identical samples give identical output (internal seed)
  expect_identical(res[, 1], res[, 2])
})

test_that("NAs are excluded from the fit and returned as NA; 0/1 are clamped", {
  skip_if_not_installed("ChAMP")
  set.seed(8)
  sim <- simulate_biased_beta(n1 = 5000, n2 = 15000)
  beta <- matrix(sim$beta, ncol = 1, dimnames = list(NULL, "s1"))
  na_idx <- c(1:50, 19000:19050)
  beta[na_idx, 1] <- NA
  beta[100, 1] <- 0
  beta[101, 1] <- 1

  expect_message(res <- normalise_bmiq(beta, probe_type = sim$probe_type, nfit = 2000),
                 "replaced 1 zeros and 1 ones")
  expect_true(all(is.na(res[na_idx, 1])))
  expect_false(anyNA(res[-na_idx, 1]))
  expect_true(all(res[-na_idx, 1] > 0 & res[-na_idx, 1] < 1))
})

test_that("nfit larger than the number of probes of a type is capped", {
  skip_if_not_installed("ChAMP")
  set.seed(10)
  sim <- simulate_biased_beta(n1 = 3000, n2 = 9000)
  beta <- matrix(sim$beta, ncol = 1, dimnames = list(NULL, "s1"))
  expect_no_error(normalise_bmiq(beta, probe_type = sim$probe_type, nfit = 10000))
})

test_that("parallel and serial give identical results", {
  skip_if_not_installed("ChAMP")
  skip_on_os("windows")
  set.seed(9)
  sim_a <- simulate_biased_beta()
  sim_b <- simulate_biased_beta()
  beta <- cbind(a = sim_a$beta, b = sim_b$beta)
  pt <- sim_a$probe_type
  expect_identical(normalise_bmiq(beta, pt, cores = 1),
                   normalise_bmiq(beta, pt, cores = 2))
})

test_that("input validation", {
  beta <- matrix(stats::runif(20), ncol = 2)
  expect_error(normalise_bmiq(beta), "probe_type")
  expect_error(normalise_bmiq(beta, probe_type = rep("II", 3)), "one entry per row")
  expect_error(normalise_bmiq(beta, probe_type = rep("III", 10)), "\"I\" and \"II\"")
  expect_error(normalise_bmiq(as.data.frame(beta), probe_type = rep("II", 10)), "MethylSet")
})

test_that("a sample without a beta distribution fails with its name", {
  skip_if_not_installed("ChAMP")
  beta <- matrix(rep(0.5, 2000), ncol = 1, dimnames = list(NULL, "flat"))
  expect_error(normalise_bmiq(beta, probe_type = rep(c("I", "II"), 1000)),
               "BMIQ failed for sample flat")
})

test_that("a degenerate fit warns once, naming the sample", {
  skip_if_not_installed("ChAMP")
  set.seed(12)
  sim <- simulate_biased_beta()
  # probe types misaligned with the values -> empty mixture components
  beta <- matrix(rev(sim$beta), ncol = 1, dimnames = list(NULL, "scrambled"))
  expect_warning(normalise_bmiq(beta, probe_type = sim$probe_type),
                 "degenerate mixture fit for sample scrambled")
})

test_that("the caller's RNG state is left untouched", {
  skip_if_not_installed("ChAMP")
  set.seed(11)
  sim <- simulate_biased_beta(n1 = 5000, n2 = 15000)
  beta <- matrix(sim$beta, ncol = 1, dimnames = list(NULL, "s1"))
  set.seed(123)
  expected <- stats::runif(3)
  set.seed(123)
  normalise_bmiq(beta, probe_type = sim$probe_type, nfit = 2000)
  expect_identical(stats::runif(3), expected)
})

test_that("MethylSet input reproduces ChAMP::champ.norm() on Noob-processed minfiData 450k", {
  skip_if_not_installed("ChAMP")
  skip_if_not_installed("minfiData")
  skip_if_not_installed("IlluminaHumanMethylation450kmanifest")
  RGsetEx <- NULL
  utils::data("RGsetEx", package = "minfiData", envir = environment())
  # as in eutopsQC: BMIQ on the betas of the single-sample Noob MethylSet
  mset <- minfi::preprocessNoob(RGsetEx[, 1:2], dyeCorr = TRUE, dyeMethod = "single")

  res <- normalise_bmiq(mset)

  # champ.norm() needs ChAMP attached for its bare data(probeInfoALL.lv) call
  # and writes into resultsDir
  withr::local_package("ChAMP")
  beta <- minfi::getBeta(mset)
  withr::with_tempdir(
    suppressMessages(utils::capture.output(
      out <- ChAMP::champ.norm(beta = beta, method = "BMIQ", arraytype = "450K",
                               cores = 1, resultsDir = "champ")))
  )
  expect_identical(dimnames(res), dimnames(beta))
  expect_equal(res, out, tolerance = 1e-12)
})
