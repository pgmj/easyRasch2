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

# ---------------------------------------------------------------------
# RMreliability() -- bootstrap path (sequential)
# ---------------------------------------------------------------------
test_that("RMreliability with boot = TRUE (sequential) returns finite CIs", {
  skip_on_cran()
  skip_if_not_installed("eRm")
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous(n = 150, k = 8)
  res <- RMreliability(df, draws = 30, rmu_iter = 5,
                       boot = TRUE, boot_iter = 10, parallel = FALSE,
                       seed = 42, output = "dataframe")
  expect_s3_class(res, "data.frame")
  psi_row  <- res[res$metric == "PSI", ]
  marg_row <- res[res$metric == "Marginal", ]
  expect_true(is.finite(psi_row$lower)  && is.finite(psi_row$upper))
  expect_true(is.finite(marg_row$lower) && is.finite(marg_row$upper))
  expect_match(psi_row$notes, "bootstrap resamples")
})

# ---------------------------------------------------------------------
# RMreliability() -- bootstrap path (parallel via mirai)
# ---------------------------------------------------------------------
test_that("RMreliability with boot = TRUE (parallel) runs via mirai", {
  skip_on_cran()
  skip_if_not_installed("eRm")
  skip_if_not_installed("mirt")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("mirai")
  df <- make_dichotomous(n = 120, k = 6)
  res <- RMreliability(df, draws = 20, rmu_iter = 3,
                       boot = TRUE, boot_iter = 4, parallel = TRUE,
                       n_cores = 2, seed = 7, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 4L)
})

# ---------------------------------------------------------------------
# cronbach_alpha() -- internal edge cases
# ---------------------------------------------------------------------
test_that("cronbach_alpha returns NA for <2 items or zero total variance", {
  one_col <- data.frame(a = sample(0:1, 20, replace = TRUE))
  expect_true(is.na(easyRasch2:::cronbach_alpha(one_col)))
  const <- as.data.frame(matrix(1L, nrow = 30, ncol = 5))
  expect_true(is.na(easyRasch2:::cronbach_alpha(const)))
})

# ---------------------------------------------------------------------
# RMUreliability() -- warnings / verbose / degenerate columns
# ---------------------------------------------------------------------
test_that("RMUreliability warns on zero-SD columns", {
  skip_if_not_installed("ggdist")
  set.seed(1)
  m <- cbind(rnorm(60), rnorm(60), rnorm(60), 5)  # last column constant
  expect_warning(RMUreliability(m), regexp = "zero standard deviation")
})

test_that("RMUreliability warns on NA values", {
  skip_if_not_installed("ggdist")
  set.seed(2)
  m <- matrix(rnorm(60 * 4), nrow = 60, ncol = 4)
  m[1, 1] <- NA
  expect_warning(RMUreliability(m), regexp = "NA value")
})

test_that("RMUreliability prints input summary when verbose = TRUE", {
  skip_if_not_installed("ggdist")
  set.seed(3)
  m <- matrix(rnorm(60 * 4), nrow = 60, ncol = 4)
  expect_message(RMUreliability(m, verbose = TRUE), regexp = "Subjects")
})
