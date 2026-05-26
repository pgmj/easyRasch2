# Tests for RMdifGamma(), RMdifGammaCutoff(), and RMdifGammaPlot()

make_dichotomous <- function(n = 200, k = 10, seed = 42L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("Item", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMdifGamma errors when iarm is not installed", {
  skip_if(requireNamespace("iarm", quietly = TRUE),
          "iarm is installed; skipping missing-package error path")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  expect_error(RMdifGamma(df, grp), regexp = "iarm")
})

test_that("RMdifGamma errors when dif_var length is wrong", {
  skip_if_not_installed("iarm")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df) - 1L, replace = TRUE))
  expect_error(RMdifGamma(df, grp),
               regexp = "same length")
})

test_that("RMdifGamma errors when dif_var has fewer than 2 levels", {
  skip_if_not_installed("iarm")
  df  <- make_dichotomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMdifGamma(df, grp), regexp = "at least 2")
})

# ---------------------------------------------------------------------
# RMdifGamma output structure
# ---------------------------------------------------------------------
test_that("RMdifGamma output = 'dataframe' returns one row per item", {
  skip_if_not_installed("iarm")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  res <- RMdifGamma(df, grp, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), ncol(df))
  expect_true(all(c("Item", "gamma", "se", "lower", "upper",
                    "padj_bh", "Significance") %in% names(res)))
  expect_true(all(res$Significance %in% c("", ".", "*", "**", "***")))
})

test_that("RMdifGamma output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  out <- RMdifGamma(df, grp, output = "kable")
  expect_s3_class(out, "knitr_kable")
})

# ---------------------------------------------------------------------
# RMdifGammaCutoff -- small iterations
# ---------------------------------------------------------------------
test_that("RMdifGammaCutoff returns a list with item_cutoffs + actual_iterations", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  res <- RMdifGammaCutoff(df, grp, iterations = 10L,
                       parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true("item_cutoffs" %in% names(res))
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), ncol(df))
  expect_true(all(c("Item", "gamma_low", "gamma_high") %in%
                  names(res$item_cutoffs)))
  expect_true(res$actual_iterations >= 1L)
})

test_that("RMdifGamma accepts an RMdifGammaCutoff result and adds Flagged column", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(df, grp, iterations = 5L,
                        parallel = FALSE, seed = 1L)
  res  <- RMdifGamma(df, grp, cutoff = cuts, output = "dataframe")
  expect_true("flagged" %in% names(res))
  expect_type(res$flagged, "logical")
})

test_that("RMdifGammaCutoff is reproducible with the same seed", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  r1 <- RMdifGammaCutoff(df, grp, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMdifGammaCutoff(df, grp, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$item_cutoffs, r2$item_cutoffs)
})

# ---------------------------------------------------------------------
# RMdifGammaPlot
# ---------------------------------------------------------------------
test_that("RMdifGammaPlot returns a ggplot from an RMdifGammaCutoff result", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(df, grp, iterations = 5L,
                        parallel = FALSE, seed = 1L)
  p <- RMdifGammaPlot(cuts)
  expect_s3_class(p, "ggplot")
})

test_that("RMdifGammaPlot with observed data overlay returns a ggplot", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df  <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(df, grp, iterations = 5L,
                        parallel = FALSE, seed = 1L)
  p <- RMdifGammaPlot(cuts, data = df, dif_var = grp)
  expect_s3_class(p, "ggplot")
})
