test_that("RMlocdepQ3 works without cutoff (returns knitr_kable)", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df)
  expect_s3_class(result, "knitr_kable")
})

test_that("RMlocdepQ3 with cutoff = NULL returns knitr_kable without cutoff annotation", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL)
  expect_s3_class(result, "knitr_kable")
})

test_that("RMlocdepQ3 errors when data has non-zero minimum", {
  df <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3 errors when data is not a data.frame or matrix", {
  expect_error(RMlocdepQ3(list(a = 1:5)), regexp = "data.frame or matrix")
})

test_that("RMlocdepQ3 returns a data.frame when output = 'dataframe' and no cutoff", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, output = "dataframe")
  expect_s3_class(result, "data.frame")
  # Lower triangle: off-diagonal lower should be numeric, upper + diag should be NA
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 returns a data.frame when output = 'dataframe' with cutoff", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "dataframe")
  expect_s3_class(result, "data.frame")
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 returns a knitr_kable object when output = 'kable' with cutoff", {
  skip_if_not_installed("mirt")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "kable")
  expect_s3_class(result, "knitr_kable")
})

test_that("RMlocdepQ3 errors on invalid cutoff", {
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMlocdepQ3(df, cutoff = "abc"), regexp = "single numeric value")
  expect_error(RMlocdepQ3(df, cutoff = c(0.1, 0.2)), regexp = "single numeric value")
})

test_that("RMlocdepQ3cutoff errors on bad data", {
  expect_error(RMlocdepQ3cutoff(list(a = 1:5)), regexp = "data.frame or matrix")
  df_bad <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3cutoff(df_bad), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3cutoff returns expected structure for dichotomous data", {
  skip_on_cran()
  skip_if_not_installed("eRm")
  skip_if_not_installed("psychotools")
  skip_if_not_installed("foreach")
  skip_if_not_installed("doParallel")
  skip_if_not_installed("doRNG")
  skip_if_not_installed("mirt")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 400, replace = TRUE), nrow = 80, ncol = 5))
  result <- suppressMessages(RMlocdepQ3cutoff(df, iterations = 10, cpu = 1, seed = 1))
  expect_type(result, "list")
  expect_true(all(c("results", "actual_iterations", "sample_n", "sample_summary",
                    "max_diff", "sd_diff", "p95", "p99", "p995", "p999",
                    "suggested_cutoff") %in% names(result)))
  expect_s3_class(result$results, "data.frame")
  expect_true(all(c("mean", "max", "diff") %in% names(result$results)))
  expect_equal(result$sample_n, 80L)
  expect_true(is.numeric(result$suggested_cutoff))
  expect_length(result$suggested_cutoff, 1L)
})

