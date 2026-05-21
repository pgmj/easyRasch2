# Tests for RMmartinLof() and RMmartinLofResiduals()

make_polytomous <- function(n = 200, k = 8, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 300, k = 8, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMmartinLof errors when partition has fewer than 2 subscales", {
  df <- make_polytomous()
  expect_error(
    RMmartinLof(df, partition = rep(1L, ncol(df)),
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = ">= 2 subscales"
  )
})

test_that("RMmartinLof errors when a subscale has fewer than 2 items", {
  df <- make_polytomous()
  partition <- c(1L, rep(2L, ncol(df) - 1L))
  expect_error(
    RMmartinLof(df, partition = partition,
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = "at least 2 items"
  )
})

# ---------------------------------------------------------------------
# Output structure (polytomous)
# ---------------------------------------------------------------------
test_that("RMmartinLof returns expected fields on polytomous data", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMmartinLof(
    df, partition = list(c("I1","I2","I3","I4"), c("I5","I6","I7","I8")),
    iterations = 20L, parallel = FALSE, seed = 1L,
    stopping = "sequential", h = 10L
  )
  expect_type(res, "list")
  expect_true(all(c("T_obs", "p_value", "actual_iterations", "rejected",
                    "partition", "n_subscales", "is_polytomous",
                    "sample_n", "n_items", "stopping", "h",
                    "T_rep", "wle_scores", "wle_correlation") %in% names(res)))
  expect_true(res$is_polytomous)
  expect_equal(res$n_subscales, 2L)
  expect_true(is.finite(res$T_obs) && res$T_obs >= 0)
  expect_true(res$p_value >= 0 && res$p_value <= 1)
})

# ---------------------------------------------------------------------
# Output structure (dichotomous)
# ---------------------------------------------------------------------
test_that("RMmartinLof works on dichotomous data (Rasch model path)", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  res <- RMmartinLof(
    df, partition = c(1,1,1,1,2,2,2,2),
    iterations = 20L, parallel = FALSE, seed = 1L,
    stopping = "sequential", h = 10L
  )
  expect_type(res, "list")
  expect_false(res$is_polytomous)
  expect_true(is.finite(res$T_obs))
})

test_that("RMmartinLof p_value is reproducible with the same seed", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous(n = 150)
  r1 <- RMmartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  r2 <- RMmartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  expect_equal(r1$T_obs,  r2$T_obs)
  expect_equal(r1$T_rep,  r2$T_rep)
})

# ---------------------------------------------------------------------
# RMmartinLofResiduals
# ---------------------------------------------------------------------
test_that("RMmartinLofResiduals output = 'dataframe' returns expected columns", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMmartinLofResiduals(
    df, partition = list(c("I1","I2","I3","I4"), c("I5","I6","I7","I8")),
    output = "dataframe"
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("t1", "t2", "total", "observed", "expected",
                    "residual", "flagged") %in% names(res)))
  expect_type(res$flagged, "logical")
})

test_that("RMmartinLofResiduals output = 'kable' returns a kable-like object", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("knitr")
  df  <- make_polytomous()
  out <- RMmartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "kable"
  )
  expect_true(inherits(out, "knitr_kable") || inherits(out, "knit_asis"))
})

test_that("RMmartinLofResiduals output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMmartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "ggplot"
  )
  expect_s3_class(p, "ggplot")
})
