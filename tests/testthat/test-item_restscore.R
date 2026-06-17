test_that("RMitemRestscore errors when iarm is not installed", {
  # This test is only meaningful when iarm is NOT installed.
  skip_if(
    requireNamespace("iarm", quietly = TRUE),
    "iarm is installed; skipping missing-package error path test"
  )
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMitemRestscore(df), regexp = "iarm")
})

test_that("RMitemRestscore errors when data does not start at 0", {
  df <- as.data.frame(matrix(sample(1:3, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMitemRestscore(df), regexp = "scored starting at 0")
})

test_that("RMitemRestscore errors when data is all NA", {
  df <- as.data.frame(matrix(NA_integer_, nrow = 10, ncol = 4))
  expect_error(RMitemRestscore(df), regexp = "no non-missing|No complete|non-missing")
})

test_that("RMitemRestscore errors when no complete cases", {
  skip_if_not_installed("iarm")
  # Create data where every row has at least one NA
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  # Make all rows have at least one NA by cycling through columns
  df[cbind(seq_len(nrow(df)), ((seq_len(nrow(df)) - 1L) %% ncol(df)) + 1L)] <- NA
  expect_error(RMitemRestscore(df), regexp = "No complete cases")
})

test_that("RMitemRestscore output = 'dataframe' returns correct structure (dichotomous)", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMitemRestscore(df, output = "dataframe")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5L)
  expected_cols <- c(
    "Item", "Observed", "Expected", "Difference",
    "p_adjusted", "Flagged", "Relative_location"
  )
  expect_equal(names(result), expected_cols)
  expect_true(all(result$Flagged %in% c("overfit", "underfit", "")))
  # overfit <=> observed above expected, underfit <=> below (when flagged)
  flg <- result$Flagged != ""
  expect_true(all((result$Difference[flg] > 0) == (result$Flagged[flg] == "overfit")))
  expect_equal(result$Item, colnames(df))
  # Signed difference equals Observed - Expected (allowing rounding noise)
  expect_equal(result$Difference,
               round(result$Observed - result$Expected, 3),
               tolerance = 1e-3)
})

test_that("RMitemRestscore output = 'kable' returns knitr_kable object", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMitemRestscore(df, output = "kable")
  expect_s3_class(result, "knitr_kable")
})

test_that("RMitemRestscore sort = 'diff' sorts by absolute Difference descending", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result_sorted   <- RMitemRestscore(df, output = "dataframe", sort = "diff")
  result_unsorted <- RMitemRestscore(df, output = "dataframe")

  expect_true(
    all(diff(abs(result_sorted$Difference)) <= 0),
    info = "Rows should be in descending order of |Difference|"
  )
  # Unsorted should contain the same rows, just in a different order
  expect_setequal(result_sorted$Item, result_unsorted$Item)
})
