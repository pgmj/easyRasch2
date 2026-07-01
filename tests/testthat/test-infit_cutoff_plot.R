# Tests for RMitemInfitPlot() (and its deprecated alias RMitemInfitCutoffPlot())

make_dichotomous <- function(n = 200, k = 8, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Output structures: simulated-only distribution
# ---------------------------------------------------------------------
test_that("RMitemInfitPlot returns a ggplot from a simfit result (no data overlay)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df     <- make_dichotomous()
  simfit <- RMitemInfitCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  p <- RMitemInfitPlot(simfit)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# Output structures: with observed-data overlay
# ---------------------------------------------------------------------
test_that("RMitemInfitPlot with observed data overlay returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df     <- make_dichotomous()
  simfit <- RMitemInfitCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  p <- RMitemInfitPlot(simfit, df)  # default statistic = "infit"
  expect_s3_class(p, "ggplot")
})

test_that("RMitemInfitPlot statistic = 'outfit' returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df     <- make_dichotomous()
  simfit <- RMitemInfitCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  p <- RMitemInfitPlot(simfit, df, statistic = "outfit")
  expect_s3_class(p, "ggplot")
})

test_that("RMitemInfitPlot statistic = 'both' returns a (patchwork) ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("patchwork")
  df     <- make_dichotomous()
  simfit <- RMitemInfitCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  p <- RMitemInfitPlot(simfit, df, statistic = "both")
  # patchwork compositions also inherit ggplot
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# Deprecated alias
# ---------------------------------------------------------------------
test_that("RMitemInfitCutoffPlot is a deprecated alias forwarding to RMitemInfitPlot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df     <- make_dichotomous()
  simfit <- RMitemInfitCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  expect_snapshot(p <- RMitemInfitCutoffPlot(simfit))
  expect_s3_class(p, "ggplot")
})
