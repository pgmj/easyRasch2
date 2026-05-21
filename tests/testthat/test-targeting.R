# Tests for RMtargeting()

make_polytomous <- function(n = 200, k = 6, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 200, k = 8, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMtargeting errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L
  expect_error(RMtargeting(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMtargeting default returns a patchwork/ggplot figure on polytomous data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  p  <- RMtargeting(df)
  # patchwork-composed plots also inherit ggplot
  expect_s3_class(p, "ggplot")
})

test_that("RMtargeting default returns a figure on dichotomous data", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  p  <- RMtargeting(df)
  expect_s3_class(p, "ggplot")
})

test_that("RMtargeting output = 'list' returns sub-plots", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  res <- RMtargeting(df, output = "list")
  expect_type(res, "list")
  expect_true(length(res) >= 1L)
})

test_that("RMtargeting accepts robust = TRUE and sort_items = 'location'", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  p <- RMtargeting(df, robust = TRUE, sort_items = "location",
                   ci_level = 0.84)
  expect_s3_class(p, "ggplot")
})
