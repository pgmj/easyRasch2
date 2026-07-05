# Tests for RMdimMartinLof() and RMdimMartinLofResiduals()

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
test_that("RMdimMartinLof errors when partition has fewer than 2 subscales", {
  df <- make_polytomous()
  expect_error(
    RMdimMartinLof(df, partition = rep(1L, ncol(df)),
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = ">= 2 subscales"
  )
})

test_that("RMdimMartinLof errors when a subscale has fewer than 2 items", {
  df <- make_polytomous()
  partition <- c(1L, rep(2L, ncol(df) - 1L))
  expect_error(
    RMdimMartinLof(df, partition = partition,
                iterations = 5L, parallel = FALSE, seed = 1L),
    regexp = "at least 2 items"
  )
})

# ---------------------------------------------------------------------
# Output structure (polytomous)
# ---------------------------------------------------------------------
test_that("RMdimMartinLof returns expected fields on polytomous data", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMdimMartinLof(
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
test_that("RMdimMartinLof works on dichotomous data (Rasch model path)", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  res <- RMdimMartinLof(
    df, partition = c(1,1,1,1,2,2,2,2),
    iterations = 20L, parallel = FALSE, seed = 1L,
    stopping = "sequential", h = 10L
  )
  expect_type(res, "list")
  expect_false(res$is_polytomous)
  expect_true(is.finite(res$T_obs))
})

test_that("RMdimMartinLof p_value is reproducible with the same seed", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_polytomous(n = 150)
  r1 <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  r2 <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                    iterations = 10L, parallel = FALSE, seed = 42L,
                    stopping = "sequential", h = 5L)
  expect_equal(r1$T_obs,  r2$T_obs)
  expect_equal(r1$T_rep,  r2$T_rep)
})

# ---------------------------------------------------------------------
# RMdimMartinLofResiduals
# ---------------------------------------------------------------------
test_that("RMdimMartinLofResiduals output = 'dataframe' returns expected columns", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  res <- RMdimMartinLofResiduals(
    df, partition = list(c("I1","I2","I3","I4"), c("I5","I6","I7","I8")),
    output = "dataframe"
  )
  expect_s3_class(res, "data.frame")
  expect_true(all(c("t1", "t2", "total", "observed", "expected",
                    "residual", "flagged") %in% names(res)))
  expect_type(res$flagged, "logical")
})

test_that("RMdimMartinLofResiduals output = 'kable' returns a kable-like object", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("knitr")
  df  <- make_polytomous()
  out <- RMdimMartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "kable"
  )
  expect_true(inherits(out, "knitr_kable") || inherits(out, "knit_asis"))
})

test_that("RMdimMartinLofResiduals output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMdimMartinLofResiduals(
    df, partition = c(1,1,1,1,2,2,2,2),
    output = "ggplot"
  )
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# stopping = "none" Monte Carlo runner (the default; existing tests above
# all use stopping = "sequential")
# ---------------------------------------------------------------------
test_that("RMdimMartinLof stopping = 'none' returns a valid p-value", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                        iterations = 20L, stopping = "none",
                        parallel = FALSE, seed = 1L)
  expect_equal(res$stopping, "none")
  expect_equal(res$actual_iterations, 20L)
  expect_true(res$p_value > 0 && res$p_value <= 1)
})

test_that("RMdimMartinLof verbose = TRUE runs the sequential progress path", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  suppressMessages(invisible(utils::capture.output(
    res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                          iterations = 10L, stopping = "none",
                          parallel = FALSE, seed = 1L, verbose = TRUE)
  )))
  expect_true(res$p_value <= 1)
})

test_that("RMdimMartinLof parallel path returns a valid result", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  skip_if_not_installed("mirai")
  df <- make_dichotomous()
  res <- suppressMessages(
    RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                   iterations = 8L, stopping = "none",
                   parallel = TRUE, n_cores = 2L, seed = 1L)
  )
  expect_true(res$p_value > 0 && res$p_value <= 1)
})

test_that("RMdimMartinLof parallel without n_cores/mc.cores warns and falls back", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  skip_if_not_installed("mirai")
  withr::local_options(mc.cores = NULL)
  df <- make_dichotomous()
  expect_warning(
    res <- RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                          iterations = 6L, stopping = "none",
                          parallel = TRUE, seed = 1L),
    regexp = "Falling back to sequential"
  )
  expect_true(res$p_value <= 1)
})

test_that("RMdimMartinLof errors with fewer than 30 complete cases", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()[1:20, ]
  expect_error(
    RMdimMartinLof(df, partition = c(1,1,1,1,2,2,2,2),
                   iterations = 5L, parallel = FALSE),
    regexp = "at least 30 complete cases"
  )
})

# ---------------------------------------------------------------------
# Partition normalisation branches
# ---------------------------------------------------------------------
test_that("RMdimMartinLof warns and drops items not assigned to a subscale", {
  skip_on_cran()
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()  # 8 items
  expect_warning(
    res <- RMdimMartinLof(df, partition = list(c(1, 2, 3), c(4, 5, 6)),
                          iterations = 10L, parallel = FALSE, seed = 1L),
    regexp = "not assigned"
  )
  expect_equal(res$n_items, 6L)
})

test_that("RMdimMartinLof errors on overlapping partition", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  expect_error(
    RMdimMartinLof(df, partition = list(c(1, 2, 3), c(3, 4, 5)),
                   iterations = 5L, parallel = FALSE),
    regexp = "overlapping"
  )
})

test_that("RMdimMartinLof errors on an invalid partition format", {
  skip_if_not_installed("psychotools")
  df <- make_dichotomous()
  expect_error(
    RMdimMartinLof(df, partition = TRUE, iterations = 5L, parallel = FALSE),
    regexp = "Invalid partition format"
  )
})

# ---------------------------------------------------------------------
# RMdimMartinLofResiduals: option branches
# ---------------------------------------------------------------------
test_that("RMdimMartinLofResiduals min_expected blanks sparse cells", {
  skip_if_not_installed("psychotools")
  df  <- make_polytomous()
  res <- RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                                 output = "dataframe", min_expected = 1)
  expect_true(all(is.na(res$residual[res$expected < 1])))
})

test_that("RMdimMartinLofResiduals ggplot color_by = 'n' returns a ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p  <- RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                                output = "ggplot", color_by = "n",
                                color_limits = c(0, 50))
  expect_s3_class(p, "ggplot")
})

test_that("RMdimMartinLofResiduals supports D = 3 for kable and ggplot", {
  skip_if_not_installed("psychotools")
  skip_if_not_installed("knitr")
  df <- make_polytomous(k = 6)
  k  <- RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3), output = "kable")
  expect_true(inherits(k, "knitr_kable") || inherits(k, "knit_asis"))
  skip_if_not_installed("ggplot2")
  p <- RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3), output = "ggplot")
  expect_s3_class(p, "ggplot")
})

test_that("RMdimMartinLofResiduals rejects bad min_expected / color_limits / D>3 ggplot", {
  skip_if_not_installed("psychotools")
  df <- make_polytomous()
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                            min_expected = -1),
    regexp = "min_expected"
  )
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,1,1,2,2,2,2),
                            color_limits = c(5, 1)),
    regexp = "color_limits"
  )
  expect_error(
    RMdimMartinLofResiduals(df, partition = c(1,1,2,2,3,3,4,4),
                            output = "ggplot"),
    regexp = "D = 2 or 3"
  )
})
