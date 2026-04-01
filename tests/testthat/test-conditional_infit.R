test_that("RMiteminfit errors when iarm is not installed", {
  skip_if(
    requireNamespace("iarm", quietly = TRUE),
    "iarm is installed; skipping missing-package error path test"
  )
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMiteminfit(df), regexp = "iarm")
})

test_that("RMiteminfit errors when data does not start at 0", {
  df <- as.data.frame(matrix(sample(1:3, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMiteminfit(df), regexp = "scored starting at 0")
})

test_that("RMiteminfit errors when data is all NA", {
  df <- as.data.frame(matrix(NA_integer_, nrow = 10, ncol = 4))
  expect_error(RMiteminfit(df), regexp = "No complete|no non-missing")
})

test_that("RMiteminfit errors when data is not a data.frame or matrix", {
  expect_error(RMiteminfit(list(a = 1:5)), regexp = "data.frame or matrix")
})

test_that("RMiteminfit errors when no complete cases", {
  skip_if_not_installed("iarm")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  # Make all rows have at least one NA by cycling through columns
  df[cbind(seq_len(nrow(df)), ((seq_len(nrow(df)) - 1L) %% ncol(df)) + 1L)] <- NA
  expect_error(RMiteminfit(df), regexp = "No complete cases")
})

test_that("RMiteminfit output = 'dataframe' returns correct structure (dichotomous)", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMiteminfit(df, output = "dataframe")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5L)
  expect_equal(names(result), c("Item", "Infit_MSQ", "Relative_location"))
  expect_equal(result$Item, colnames(df))
  expect_true(all(result$Infit_MSQ > 0))
})

test_that("RMiteminfit output = 'kable' returns knitr_kable with complete-cases caption", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMiteminfit(df, output = "kable")
  expect_s3_class(result, "knitr_kable")
  cap <- attr(result, "caption")
  expect_true(grepl("complete cases", cap))
})

test_that("RMiteminfit sort = 'infit' sorts by Infit_MSQ descending", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result_sorted   <- RMiteminfit(df, output = "dataframe", sort = "infit")
  result_unsorted <- RMiteminfit(df, output = "dataframe")

  expect_true(
    all(diff(result_sorted$Infit_MSQ) <= 0),
    info = "Rows should be in descending order of Infit_MSQ"
  )
  expect_setequal(result_sorted$Item, result_unsorted$Item)
})
