test_that("RMlocdepQ3 errors when cutoff is missing", {
  df <- as.data.frame(matrix(sample(0:1, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df), regexp = "cutoff")
})

test_that("RMlocdepQ3 errors when data has non-zero minimum", {
  df <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df, cutoff = 0.2), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3 errors when data is not a data.frame or matrix", {
  expect_error(RMlocdepQ3(list(a = 1:5), cutoff = 0.2), regexp = "data.frame or matrix")
})

test_that("RMlocdepQ3 returns a data.frame when output = 'dataframe'", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "dataframe")
  expect_s3_class(result, "data.frame")
  # Lower triangle: off-diagonal lower should be numeric, upper + diag should be NA
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 returns a knitr_kable object when output = 'kable'", {
  skip_if_not_installed("mirt")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "kable")
  expect_s3_class(result, "knitr_kable")
})
