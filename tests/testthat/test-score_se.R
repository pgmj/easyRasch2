# Tests for RMscoreSE()

make_polytomous <- function(n = 150, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 150, k = 8, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMscoreSE errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L
  expect_error(RMscoreSE(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMscoreSE output = 'dataframe' returns one row per raw score", {
  skip_if_not_installed("eRm")
  df  <- make_polytomous()
  res <- RMscoreSE(df, output = "dataframe")
  expect_s3_class(res, "data.frame")
  # k = 5, max category = 3 -> max raw score = 15 -> 16 rows
  expect_equal(nrow(res), ncol(df) * 3L + 1L)
  expect_named(res, c("raw_score", "logit_score", "logit_se"))
  expect_equal(res$raw_score, 0:15)
})

test_that("RMscoreSE WLE on dichotomous data produces valid table", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  res <- RMscoreSE(df, method = "WLE", output = "dataframe")
  expect_equal(nrow(res), ncol(df) + 1L)
  expect_true(all(is.finite(res$logit_score[2:(ncol(df))])))
})

test_that("RMscoreSE output = 'kable' returns knitr_kable", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  df <- make_polytomous()
  out <- RMscoreSE(df, output = "kable")
  expect_s3_class(out, "knitr_kable")
})

test_that("RMscoreSE output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMscoreSE(df, output = "ggplot")
  expect_s3_class(p, "ggplot")
})

# EAP method goes through mirt; skip if not available
test_that("RMscoreSE method = 'EAP' produces a valid table when mirt is installed", {
  skip_if_not_installed("mirt")
  df  <- make_polytomous()
  res <- RMscoreSE(df, method = "EAP", output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("raw_score", "logit_score", "logit_se"))
  expect_true(nrow(res) > 0L)
})
