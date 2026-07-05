# Tests for RMitemInfitCutoff() branches not exercised in
# test-conditional_infit.R (which covers the dichotomous
# quantile/hdci/conditional/resample paths): the polytomous DGP paths, the
# parallel runner, verbose progress, and error/fallback branches.

make_poly <- function(n = 200, k = 5, seed = 3L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dich <- function(n = 200, k = 5, seed = 3L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

test_that("polytomous data (resample DGP) returns per-item cutoffs", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  res <- RMitemInfitCutoff(make_poly(), iterations = 15, parallel = FALSE,
                           seed = 1, cutoff_method = "quantile")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_equal(res$actual_iterations, 15L)
  expect_true(all(res$item_cutoffs$infit_low < res$item_cutoffs$infit_high))
  expect_true(all(res$item_cutoffs$outfit_low < res$item_cutoffs$outfit_high))
})

test_that("polytomous data with conditional DGP runs", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  res <- RMitemInfitCutoff(make_poly(), iterations = 15, parallel = FALSE,
                           seed = 1, dgp = "conditional",
                           cutoff_method = "quantile")
  expect_identical(res$dgp, "conditional")
  expect_equal(nrow(res$item_cutoffs), 5L)
})

test_that("verbose = TRUE runs the sequential progress path", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  suppressMessages(invisible(utils::capture.output(
    res <- RMitemInfitCutoff(make_dich(), iterations = 5, parallel = FALSE,
                             seed = 1, cutoff_method = "quantile",
                             verbose = TRUE)
  )))
  expect_equal(res$actual_iterations, 5L)
})

test_that("errors when no complete cases remain after na.omit", {
  skip_if_not_installed("iarm")
  d <- data.frame(a = c(0, NA, 1), b = c(NA, 1, NA))
  expect_snapshot(
    RMitemInfitCutoff(d, iterations = 3, parallel = FALSE),
    error = TRUE
  )
})

test_that("parallel runner produces the same-shaped result", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("mirai")
  res <- suppressMessages(
    RMitemInfitCutoff(make_dich(), iterations = 6, parallel = TRUE,
                      n_cores = 2, seed = 1, cutoff_method = "quantile")
  )
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_equal(res$actual_iterations, 6L)
})

test_that("parallel = TRUE without n_cores or mc.cores warns and falls back", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("mirai")
  withr::local_options(mc.cores = NULL)
  expect_snapshot(
    res <- RMitemInfitCutoff(make_dich(), iterations = 4, parallel = TRUE,
                             seed = 1, cutoff_method = "quantile")
  )
  expect_equal(res$actual_iterations, 4L)
})
