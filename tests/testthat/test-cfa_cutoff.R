# Tests for RMcfaCutoff() / RMcfaPlot()
#
# Iterations are intentionally very small (10-20) so the file runs in
# a few seconds. Each test uses parallel = FALSE because parallel
# requires mirai daemons and complicates CRAN test runs.

# ---------------------------------------------------------------------
# Validation errors (no lavaan needed for these)
# ---------------------------------------------------------------------
test_that("RMcfaCutoff errors when lavaan is not installed", {
  skip_if(requireNamespace("lavaan", quietly = TRUE),
          "lavaan is installed; skipping missing-package error path test")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  expect_error(RMcfaCutoff(df), regexp = "lavaan")
})

test_that("RMcfaCutoff errors when data has non-zero minimum", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(1:3, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  expect_error(RMcfaCutoff(df), regexp = "scored starting at 0")
})

test_that("RMcfaCutoff errors when fewer than 3 items", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 100, replace = TRUE),
                             nrow = 50, ncol = 2))
  expect_error(RMcfaCutoff(df), regexp = "at least 3 items")
})

test_that("RMcfaCutoff rejects ML-family estimators", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  for (est in c("ML", "MLR", "MLM", "MLF")) {
    expect_error(
      RMcfaCutoff(df, estimator = est, iterations = 5,
                  parallel = FALSE),
      regexp = "not appropriate for ordinal-CFA"
    )
  }
})

test_that("RMcfaCutoff rejects out-of-range percentile", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  for (bad in c(0, 50, 100, 101, NA, -1)) {
    expect_error(
      RMcfaCutoff(df, percentile = bad,
                  iterations = 5, parallel = FALSE),
      regexp = "percentile"
    )
  }
})

# ---------------------------------------------------------------------
# Default output = "kable" with attr "result"
# ---------------------------------------------------------------------
test_that("RMcfaCutoff default output is a knitr_kable with the result list attached", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res <- RMcfaCutoff(
    raschdat1[, 1:8],
    iterations = 10,
    parallel   = FALSE,
    seed       = 1
  )

  expect_s3_class(res, "knitr_kable")
  result_list <- attr(res, "result")
  expect_type(result_list, "list")
  expect_equal(result_list$percentile, 99)
})

# ---------------------------------------------------------------------
# Output structure (output = "list")
# ---------------------------------------------------------------------
test_that("RMcfaCutoff output = 'list' returns the expected structure on dichotomous data", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res <- RMcfaCutoff(
    raschdat1[, 1:8],
    iterations = 10,
    parallel   = FALSE,
    seed       = 1,
    output     = "list"
  )

  expect_type(res, "list")
  expected_names <- c("observed", "simulated", "percentile", "cutoffs",
                      "flagged", "actual_iterations", "sample_n",
                      "n_items", "item_names", "is_polytomous",
                      "estimator")
  expect_true(all(expected_names %in% names(res)))

  expect_named(res$observed, c("cfi", "rmsea", "srmr"))
  expect_s3_class(res$simulated, "data.frame")
  expect_named(res$simulated, c("iteration", "cfi", "rmsea", "srmr"))
  expect_true(nrow(res$simulated) > 0L)
  expect_true(nrow(res$simulated) <= 10L)

  expect_named(res$cutoffs, c("cfi", "rmsea", "srmr"))
  expect_named(res$flagged, c("cfi", "rmsea", "srmr"))
  expect_type(res$flagged, "logical")
  expect_equal(res$percentile, 99)

  expect_false(res$is_polytomous)
  expect_equal(res$estimator, "WLSMV")
  expect_equal(res$n_items, 8L)
})

test_that("Cutoffs change with different percentile values", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res99 <- RMcfaCutoff(raschdat1[, 1:8], iterations = 30,
                       parallel = FALSE, seed = 1,
                       percentile = 99, output = "list")
  res95 <- RMcfaCutoff(raschdat1[, 1:8], iterations = 30,
                       parallel = FALSE, seed = 1,
                       percentile = 95, output = "list")

  # Same simulated distribution (same seed), different cutoffs
  expect_equal(res99$simulated, res95$simulated)
  # 99th percentile is more extreme than 95th
  expect_gte(res99$cutoffs[["rmsea"]], res95$cutoffs[["rmsea"]])
  expect_gte(res99$cutoffs[["srmr"]],  res95$cutoffs[["srmr"]])
  expect_lte(res99$cutoffs[["cfi"]],   res95$cutoffs[["cfi"]])
})

# ---------------------------------------------------------------------
# Output structure (PCM/polytomous on pcmdat2)
# ---------------------------------------------------------------------
test_that("RMcfaCutoff auto-picks PCM on polytomous data", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("pcmdat2", package = "eRm")

  res <- RMcfaCutoff(
    pcmdat2[1:100, ],
    iterations = 10,
    parallel   = FALSE,
    seed       = 2,
    output     = "list"
  )
  expect_true(res$is_polytomous)
  expect_true(res$actual_iterations >= 1L)
})

# ---------------------------------------------------------------------
# Reproducibility via seed
# ---------------------------------------------------------------------
test_that("Same seed produces identical simulated distributions", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res1 <- RMcfaCutoff(raschdat1[, 1:8], iterations = 5,
                      parallel = FALSE, seed = 42, output = "list")
  res2 <- RMcfaCutoff(raschdat1[, 1:8], iterations = 5,
                      parallel = FALSE, seed = 42, output = "list")
  expect_equal(res1$simulated, res2$simulated)
  expect_equal(res1$observed,  res2$observed)
})

# ---------------------------------------------------------------------
# Plot
# ---------------------------------------------------------------------
test_that("RMcfaPlot returns a ggplot object from a list result", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  data("raschdat1", package = "eRm")

  res <- RMcfaCutoff(raschdat1[, 1:8], iterations = 5,
                     parallel = FALSE, seed = 1, output = "list")
  p <- RMcfaPlot(res)
  expect_s3_class(p, "ggplot")
})

test_that("RMcfaPlot accepts the kable output (reads result via attr)", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  data("raschdat1", package = "eRm")

  res_kbl <- RMcfaCutoff(raschdat1[, 1:8], iterations = 5,
                         parallel = FALSE, seed = 1)  # default kable
  p <- RMcfaPlot(res_kbl)
  expect_s3_class(p, "ggplot")
})

test_that("RMcfaPlot percentile override recomputes cutoffs without re-simulating", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  data("raschdat1", package = "eRm")

  res <- RMcfaCutoff(raschdat1[, 1:8], iterations = 30,
                     parallel = FALSE, seed = 1, output = "list")
  expect_s3_class(RMcfaPlot(res, percentile = 95), "ggplot")
  expect_s3_class(RMcfaPlot(res, percentile = 99.5), "ggplot")
  expect_error(RMcfaPlot(res, percentile = 50),
               regexp = "percentile")
})

test_that("RMcfaPlot rejects non-RMcfaCutoff input", {
  skip_if_not_installed("ggplot2")
  expect_error(RMcfaPlot(list(foo = 1)),
               regexp = "must be the result returned by RMcfaCutoff")
})
