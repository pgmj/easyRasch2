# Tests for RMreliability() and RMUreliability()

make_dichotomous <- function(n = 200, k = 10, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# RMreliability() -- input validation
# ---------------------------------------------------------------------
test_that("RMreliability errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L
  expect_error(RMreliability(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# RMreliability() -- output structure (no bootstrap path)
# ---------------------------------------------------------------------
test_that("RMreliability output = 'dataframe' returns four reliability rows", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  res <- RMreliability(df, draws = 50, rmu_iter = 20,
                       output = "dataframe", seed = 1)
  expect_s3_class(res, "data.frame")
  # Four metrics: alpha, PSI, Empirical, RMU
  expect_true(nrow(res) >= 4L)
  expect_true(all(c("metric", "estimate") %in% names(res)))
  # All point estimates should be finite numerics
  expect_true(all(is.finite(res$estimate)))
  # All estimates should be in plausible reliability range
  expect_true(all(res$estimate >= -1 & res$estimate <= 1))
})

test_that("RMreliability output = 'kable' returns knitr_kable", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  out <- RMreliability(df, draws = 30, rmu_iter = 10,
                       output = "kable", seed = 1)
  expect_s3_class(out, "knitr_kable")
})

# ---------------------------------------------------------------------
# RMUreliability() -- input validation
# ---------------------------------------------------------------------
test_that("RMUreliability errors when fewer than 2 columns", {
  skip_if_not_installed("ggdist")
  m <- matrix(rnorm(50), nrow = 50, ncol = 1)
  expect_error(RMUreliability(m), regexp = "at least 2 columns")
})

test_that("RMUreliability errors when level is out of range", {
  skip_if_not_installed("ggdist")
  m <- matrix(rnorm(200), nrow = 50, ncol = 4)
  expect_error(RMUreliability(m, level = 0), regexp = "between 0 and 1")
  expect_error(RMUreliability(m, level = 1), regexp = "between 0 and 1")
})

# ---------------------------------------------------------------------
# RMUreliability() -- output structure
# ---------------------------------------------------------------------
test_that("RMUreliability on synthetic high-reliability draws returns plausible structure", {
  skip_if_not_installed("ggdist")
  # Build draws with a strong common signal -> high reliability expected
  set.seed(42)
  n_subj  <- 100
  n_draws <- 20
  theta   <- stats::rnorm(n_subj, mean = 0, sd = 1)
  draws   <- matrix(NA_real_, nrow = n_subj, ncol = n_draws)
  for (j in seq_len(n_draws)) {
    draws[, j] <- theta + stats::rnorm(n_subj, sd = 0.3)
  }
  res <- RMUreliability(draws)
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 1L)
  expect_true(all(c("rmu_estimate", "hdci_lowerbound",
                    "hdci_upperbound") %in% names(res)))
  expect_true(is.finite(res$rmu_estimate))
  # Strong signal should give reliability well above 0.5
  expect_gt(res$rmu_estimate, 0.5)
})
