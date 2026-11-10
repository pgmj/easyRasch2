# Regression tests for the shared estimation helpers in R/utils-theta.R.
# Findings are numbered as in dev/code-review/utils-theta.md.

pcm_data <- function() {
  skip_if_not_installed("eRm")
  e <- new.env()
  utils::data(pcmdat2, package = "eRm", envir = e)
  as.matrix(e$pcmdat2)
}

test_that("WLE clamps to the boundary on the root's side (finding 1)", {
  # A single very hard item answered 0: the WLE root is above +6.
  expect_equal(unname(.theta_wle(0, list(7.5), c(-6, 6))["theta"]), 6)
  expect_equal(unname(.theta_wle(c(1, 0), list(9, 9), c(-6, 6))["theta"]), 6)
  # Mirror image: very easy items, root below -6.
  expect_equal(unname(.theta_wle(1, list(-7.5), c(-6, 6))["theta"]), -6)
  # A root inside the range is unaffected.
  expect_true(abs(.theta_wle(0, list(7.5), c(-10, 10))["theta"]) < 10)
})

test_that("item shapes psychotools would reshape stop with the item named (findings 2, 2b, 11)", {
  x <- pcm_data()
  interior <- x
  interior[interior[, 1] == 1, 1] <- 2
  expect_error(.fit_cml_thresholds(interior), "I1 \\(empty category 1\\)")

  lowest <- x
  lowest[lowest[, 1] == 0, 1] <- 1
  expect_error(.fit_cml_thresholds(lowest), "I1 \\(empty category 0\\)")

  constant <- x
  constant[, 2] <- 0
  expect_error(.fit_cml_thresholds(constant), "I2 \\(only category 0\\)")
})

test_that("an empty top category is allowed and keeps alignment", {
  x <- pcm_data()
  x[x[, 1] == 2, 1] <- 1
  thr <- suppressWarnings(.fit_cml_thresholds(x))
  expect_length(thr, ncol(x))
  expect_equal(unname(vapply(thr, length, integer(1))), c(1L, 2L, 2L, 2L))
})

test_that("all-NA respondents are dropped before the CML fit (finding 12)", {
  x <- pcm_data()
  ref <- .fit_cml_thresholds(x)
  x_na <- rbind(x, NA)
  expect_equal(.fit_cml_thresholds(x_na), ref)

  dich <- (x > 0) * 1
  ref_d <- .fit_cml_thresholds(dich)
  # raschmodel() segfaults on an all-NA row, so this would crash without the fix.
  expect_equal(.fit_cml_thresholds(rbind(NA, dich)), ref_d)
})

test_that("EAP returns NA for respondents with no answered items (finding 3)", {
  x <- pcm_data()
  thr <- .fit_cml_thresholds(x)
  x_na <- rbind(NA, x)
  est <- .estimate_thetas(x_na, thr, method = "EAP")
  expect_true(is.na(est$theta[1]))
  expect_true(is.na(est$sem[1]))
  # The empty row does not enter the prior-SD estimate.
  ref <- .estimate_thetas(x, thr, method = "EAP")
  expect_equal(attr(est, "prior_sd"), attr(ref, "prior_sd"))
  expect_equal(est$theta[-1], ref$theta)
})

test_that("zero-row input returns an empty data frame (finding 10)", {
  out <- .estimate_thetas(matrix(numeric(0), 0, 3), list(-1, 0, 1))
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 0L)
  expect_named(out, c("theta", "sem"))
})

test_that("latent moments recover an off-target normal distribution", {
  skip_if_not_installed("eRm")
  set.seed(20261110)
  d <- as.matrix(eRm::sim.rasch(stats::rnorm(2000, 3, 0.6),
                                seq(-2, 2, length.out = 20)))
  lat <- .latent_moments(d, .fit_cml_thresholds(d))
  expect_equal(lat$mean, 3, tolerance = 0.1)
  expect_equal(lat$sd, 0.6, tolerance = 0.1)
})
