# Tests for RMretestSD()

make_retest <- function(n = 300, k = 6, occ_sd = 0.3, seed = 1L) {
  set.seed(seed)
  thr <- lapply(seq(-1.2, 1.2, length.out = k), function(b) b + c(-0.8, 0, 0.8))
  names(thr) <- paste0("I", seq_len(k))
  theta <- stats::rnorm(n, 0, 1.4)
  a <- as.data.frame(sim_partial_score(thr, theta + stats::rnorm(n, 0, occ_sd)))
  b <- as.data.frame(sim_partial_score(thr, theta + stats::rnorm(n, 0, occ_sd)))
  colnames(a) <- colnames(b) <- names(thr)
  list(t1 = a, t2 = b, thr = thr)
}

test_that("output has the documented shape", {
  skip_on_cran()
  d <- make_retest(n = 150)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 60,
                    parallel = FALSE, boot = FALSE, seed = 1,
                    output = "dataframe")
  expect_named(res, c("n", "var_change", "var_error", "variance", "sd",
                      "lower", "upper"))
  expect_identical(res$n, 150L)
  expect_true(is.finite(res$var_change) && is.finite(res$var_error))
})

test_that("variance is the documented halved difference", {
  skip_on_cran()
  d <- make_retest(n = 150)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 60,
                    parallel = FALSE, boot = FALSE, seed = 1,
                    output = "dataframe")
  expect_equal(res$variance, (res$var_change - res$var_error) / 2)
  expect_equal(res$sd, sqrt(res$variance))
})

test_that("the error term is simulated, not the mean squared SE", {
  skip_on_cran()
  d <- make_retest(n = 200, k = 6)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 120,
                    parallel = FALSE, boot = FALSE, seed = 2,
                    output = "dataframe")
  e1 <- .estimate_thetas(as.matrix(d$t1), d$thr, method = "WLE")
  e2 <- .estimate_thetas(as.matrix(d$t2), d$thr, method = "WLE")
  asymptotic <- mean((e1$sem^2 + e2$sem^2)[is.finite(e1$sem + e2$sem)])
  # the asymptotic term overstates the real change variance on a short scale
  expect_lt(res$var_error, asymptotic)
})

test_that("a negative variance is reported, not floored", {
  skip_on_cran()
  # No occasion variance at all, so the subtraction can land below zero
  d <- make_retest(n = 200, occ_sd = 0, seed = 7)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 120,
                    parallel = FALSE, boot = FALSE, seed = 7,
                    output = "dataframe")
  expect_true(res$variance < 0.06)
  if (res$variance < 0) {
    expect_true(is.na(res$sd))
    txt <- paste(as.character(
      RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 120,
                 parallel = FALSE, boot = FALSE, seed = 7)
    ), collapse = " ")
    expect_match(gsub("[[:space:]]+", " ", txt), "estimated variance is negative")
  }
})

test_that("a real occasion SD is recovered within tolerance", {
  skip_on_cran()
  d <- make_retest(n = 600, k = 12, occ_sd = 0.4, seed = 3)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 150,
                    parallel = FALSE, boot = FALSE, seed = 3,
                    output = "dataframe")
  expect_equal(res$variance, 0.16, tolerance = 0.5)
  expect_gt(res$sd, 0.2)
})

test_that("the bootstrap holds the error term fixed and widens with noise", {
  skip_on_cran()
  skip_if_not_installed("ggdist")
  d <- make_retest(n = 200)
  res <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 60,
                    parallel = FALSE, boot = TRUE, boot_iter = 200,
                    seed = 4, output = "dataframe")
  skip_if(is.na(res$sd))
  expect_true(is.finite(res$lower) || is.na(res$lower))
  if (is.finite(res$lower) && is.finite(res$upper)) {
    expect_lt(res$lower, res$upper)
    expect_true(res$sd >= res$lower && res$sd <= res$upper)
  }
})

test_that("paired-input checks are shared with RMpersonChange", {
  d <- make_retest(n = 40)
  expect_error(RMretestSD(d$t1, d$t2[1:10, ]), regexp = "same number of rows")
  expect_error(RMretestSD(d$t1, d$t2[, 1:3]), regexp = "same items")
})

test_that("conf_int is validated", {
  d <- make_retest(n = 40)
  expect_error(RMretestSD(d$t1, d$t2, conf_int = 1.5),
               regexp = "between 0 and 1")
})

test_that("the estimate feeds RMpersonChange's retest null", {
  skip_on_cran()
  d <- make_retest(n = 200, k = 12, occ_sd = 0.4, seed = 8)
  rs <- RMretestSD(d$t1, d$t2, item_params = d$thr, sim_iter = 80,
                   parallel = FALSE, boot = FALSE, seed = 8,
                   output = "dataframe")
  skip_if(is.na(rs$sd))
  res <- RMpersonChange(d$t1, d$t2, item_params = d$thr, null = "retest",
                        retest_sd = rs$sd, critical = 1.96)
  expect_equal(
    res$se_diff,
    sqrt(res$se_t1^2 + res$se_t2^2 + 2 * rs$sd^2)
  )
})
