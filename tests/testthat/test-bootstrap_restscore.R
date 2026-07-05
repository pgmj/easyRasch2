# Tests for RMitemRestscoreBoot()

make_dichotomous <- function(n = 200, k = 8, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("Item", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMitemRestscoreBoot errors when data has non-zero minimum", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous() + 1L
  expect_error(RMitemRestscoreBoot(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMitemRestscoreBoot output = 'dataframe' returns one row per item × classification", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  res <- RMitemRestscoreBoot(df, iterations = 10L, samplesize = 100L,
                         parallel = FALSE, seed = 1L,
                         output = "dataframe")
  expect_s3_class(res, "data.frame")
  # 8 items × 3 classifications (overfit / underfit / no misfit) = 24 rows
  expect_equal(nrow(res), ncol(df) * 3L)
  expect_true(all(c("Item", "item_restscore", "n", "percent",
                    "Infit_MSQ", "Relative_location") %in% names(res)))
})

test_that("RMitemRestscoreBoot output = 'raw' returns per-iteration data", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  res <- RMitemRestscoreBoot(df, iterations = 10L, samplesize = 100L,
                         parallel = FALSE, seed = 1L,
                         output = "raw")
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) > 0L)
})

test_that("RMitemRestscoreBoot output = 'kable' returns a knitr_kable", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  out <- RMitemRestscoreBoot(df, iterations = 10L, samplesize = 100L,
                         parallel = FALSE, seed = 1L,
                         output = "kable")
  expect_s3_class(out, "knitr_kable")
})

test_that("RMitemRestscoreBoot is reproducible with the same seed", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  r1 <- RMitemRestscoreBoot(df, iterations = 5L, samplesize = 100L,
                        parallel = FALSE, seed = 42L,
                        output = "dataframe")
  r2 <- RMitemRestscoreBoot(df, iterations = 5L, samplesize = 100L,
                        parallel = FALSE, seed = 42L,
                        output = "dataframe")
  expect_equal(r1, r2)
})
