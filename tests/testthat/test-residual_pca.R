# Tests for RMdimResidualPCA() and RMdimResidualPCACutoff()

make_dichotomous <- function(n = 200, k = 10, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_polytomous <- function(n = 200, k = 6, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# RMdimResidualPCA -- input validation
# ---------------------------------------------------------------------
test_that("RMdimResidualPCA errors when data has non-zero minimum", {
  df <- make_dichotomous() + 1L
  expect_error(RMdimResidualPCA(df), regexp = "scored starting at 0")
})

test_that("RMdimResidualPCA errors when cutoff is the wrong type", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  expect_error(RMdimResidualPCA(df, cutoff = "abc"),
               regexp = "cutoff")
})

# ---------------------------------------------------------------------
# RMdimResidualPCA -- output structure
# ---------------------------------------------------------------------
test_that("RMdimResidualPCA output = 'dataframe' returns eigenvalue table + variance_partition attr", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  res <- RMdimResidualPCA(df, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("Component", "Eigenvalue", "Proportion_of_variance"))
  expect_equal(nrow(res), 5L)  # default n_components
  expect_true(all(res$Eigenvalue > 0))
  expect_true(all(res$Proportion_of_variance > 0 & res$Proportion_of_variance <= 1))

  vp <- attr(res, "variance_partition")
  expect_type(vp, "list")
  expect_true(all(c("total", "explained", "unexplained",
                    "pct_explained", "pct_unexplained",
                    "n_persons") %in% names(vp)))
})

test_that("RMdimResidualPCA output = 'kable' returns knitr_kable", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  out <- RMdimResidualPCA(df, output = "kable")
  expect_s3_class(out, "knitr_kable")
})

test_that("RMdimResidualPCA output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  p  <- RMdimResidualPCA(df, output = "ggplot")
  expect_s3_class(p, "ggplot")
})

test_that("RMdimResidualPCA output = 'loadings' is a backward-compatible alias for 'ggplot'", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  p  <- RMdimResidualPCA(df, output = "loadings")
  expect_s3_class(p, "ggplot")
})

test_that("RMdimResidualPCA auto-picks PCM on polytomous data", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  res <- RMdimResidualPCA(df, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_true(nrow(res) >= 1L)
})

test_that("RMdimResidualPCA accepts numeric cutoff and adds Flagged column", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  res <- RMdimResidualPCA(df, cutoff = 2.0, output = "dataframe")
  expect_true("Flagged" %in% names(res))
  expect_type(res$Flagged, "logical")
})

test_that("RMdimResidualPCA n_components controls number of rows", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  res <- RMdimResidualPCA(df, n_components = 3L, output = "dataframe")
  expect_equal(nrow(res), 3L)
})

# ---------------------------------------------------------------------
# RMdimResidualPCACutoff -- small iterations to keep tests fast
# ---------------------------------------------------------------------
test_that("RMdimResidualPCACutoff returns a list with simulated eigenvalues + percentile cutoffs", {
  skip_if_not_installed("eRm")
  df  <- make_dichotomous()
  res <- RMdimResidualPCACutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(c("results", "p95", "p99", "p995", "p999",
                    "suggested_cutoff", "actual_iterations",
                    "sample_n", "item_names") %in% names(res)))
  expect_s3_class(res$results, "data.frame")
  expect_true(res$actual_iterations >= 1L)
  expect_true(res$actual_iterations <= 10L)
  expect_true(is.finite(res$suggested_cutoff))
})

test_that("RMdimResidualPCA accepts an RMdimResidualPCACutoff result for cutoff arg", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  bound <- RMdimResidualPCACutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  res   <- RMdimResidualPCA(df, cutoff = bound, output = "dataframe")
  expect_true("Flagged" %in% names(res))
})

test_that("RMdimResidualPCACutoff is reproducible with the same seed", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  r1 <- RMdimResidualPCACutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMdimResidualPCACutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$results, r2$results)
})

test_that("RMdimResidualPCA warns on sparse/zero-variance items (CML)", {
  skip_if_not_installed("eRm")
  df <- make_dichotomous()
  df[[1]] <- 0L   # constant item -> destabilises CML; should warn, not fail silently
  expect_warning(try(RMdimResidualPCA(df, output = "dataframe"), silent = TRUE),
                 regexp = "[Ss]parse|zero-variance")
})
