test_that("RMlocdepQ3 no longer errors when cutoff is missing", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  # NULL cutoff should return a knitr_kable, not an error
  expect_s3_class(RMlocdepQ3(df), "knitr_kable")
})

test_that("RMlocdepQ3 with cutoff = NULL returns raw Q3 kable", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL)
  expect_s3_class(result, "knitr_kable")
  # Caption should mention "Raw Q3"
  cap <- attr(result, "caption")
  if (is.null(cap)) cap <- paste(as.character(result), collapse = "\n")
  expect_true(grepl("Raw Q3", cap))
})

test_that("RMlocdepQ3 with cutoff = NULL returns raw Q3 dataframe", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL, output = "dataframe")
  expect_s3_class(result, "data.frame")
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
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
  # Caption should mention dynamic cut-off
  cap <- attr(result, "caption")
  if (is.null(cap)) cap <- paste(as.character(result), collapse = "\n")
  expect_true(grepl("Dynamic cut-off", cap))
})
