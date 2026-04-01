test_that("RMlocdepQ3 errors when data has non-zero minimum", {
  df <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3 errors when data is not a data.frame or matrix", {
  expect_error(RMlocdepQ3(list(a = 1:5)), regexp = "data.frame or matrix")
})

test_that("RMlocdepQ3 errors when cutoff is not a single numeric", {
  df <- as.data.frame(matrix(sample(0:1, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df, cutoff = "bad"), regexp = "single numeric")
  expect_error(RMlocdepQ3(df, cutoff = c(0.1, 0.2)), regexp = "single numeric")
})

test_that("RMlocdepQ3 with NULL cutoff returns data.frame with Q3 values", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL, output = "dataframe")
  expect_s3_class(result, "data.frame")
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 default (no cutoff) returns data.frame correctly", {
  skip_if_not_installed("mirt")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, output = "dataframe")
  expect_s3_class(result, "data.frame")
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

test_that("RMlocdepQ3 returns a knitr_kable object without cutoff", {
  skip_if_not_installed("mirt")
  set.seed(3)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, output = "kable")
  expect_s3_class(result, "knitr_kable")
})

test_that("RMlocdepQ3 returns a knitr_kable object with cutoff", {
  skip_if_not_installed("mirt")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "kable")
  expect_s3_class(result, "knitr_kable")
  # Caption should include dynamic cut-off info
  expect_true(grepl("Dynamic cut-off", attr(result, "caption")))
})

test_that("RMlocdepQ3 kable without cutoff has informative caption", {
  skip_if_not_installed("mirt")
  set.seed(4)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, output = "kable")
  expect_true(grepl("No cut-off", attr(result, "caption")))
})

test_that("RMlocdepQ3cutoff errors on bad data", {
  expect_error(RMlocdepQ3cutoff(list(a = 1:5)), regexp = "data.frame or matrix")
  df_bad <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3cutoff(df_bad), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3cutoff returns expected structure", {
  skip_on_cran()
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("psychotools")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 400, replace = TRUE), nrow = 80, ncol = 5))
  result <- RMlocdepQ3cutoff(df, iterations = 20, seed = 1, progress = FALSE)
  expect_type(result, "list")
  expect_true(all(c("results", "actual_iterations", "sample_n", "sample_summary",
                     "max_diff", "sd_diff", "p95", "p99", "p995", "p999",
                     "suggested_cutoff") %in% names(result)))
  expect_s3_class(result$results, "data.frame")
  expect_true(all(c("mean", "max", "diff") %in% colnames(result$results)))
  expect_equal(result$sample_n, 80L)
  expect_identical(result$suggested_cutoff, result$p99)
})
