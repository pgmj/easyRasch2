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
  skip_if(
    requireNamespace("iarm", quietly = TRUE),
    "iarm is installed; skipping missing-package error path"
  )
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  expect_error(RMdifGamma(df, grp), regexp = "iarm")
})

test_that("RMdifGamma errors when dif_var length is wrong", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df) - 1L, replace = TRUE))
  expect_error(RMdifGamma(df, grp), regexp = "same length")
})

test_that("RMdifGamma errors when dif_var has fewer than 2 levels", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMdifGamma(df, grp), regexp = "at least 2")
})

# ---------------------------------------------------------------------
# RMdifGamma output structure
# ---------------------------------------------------------------------
test_that("RMdifGamma output = 'dataframe' returns one row per item", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  res <- RMdifGamma(df, grp, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), ncol(df))
  expect_true(all(
    c("Item", "gamma", "se", "lower", "upper", "padj_bh", "Significance") %in%
      names(res)
  ))
  expect_true(all(res$Significance %in% c("", ".", "*", "**", "***")))
})

test_that("RMdifGamma output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
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
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  res <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 10L,
    parallel = FALSE,
    seed = 1L
  )
  expect_type(res, "list")
  expect_true("item_cutoffs" %in% names(res))
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), ncol(df))
  expect_true(all(
    c("Item", "gamma_low", "gamma_high") %in%
      names(res$item_cutoffs)
  ))
  expect_true(res$actual_iterations >= 1L)
})

test_that("RMdifGamma accepts an RMdifGammaCutoff result and adds Flagged column", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L
  )
  res <- RMdifGamma(df, grp, cutoff = cuts, output = "dataframe")
  expect_true("flagged" %in% names(res))
  expect_type(res$flagged, "logical")
})

test_that("RMdifGammaCutoff is reproducible with the same seed", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
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
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L
  )
  p <- RMdifGammaPlot(cuts)
  expect_s3_class(p, "ggplot")
})

test_that("RMdifGammaPlot with observed data overlay returns a ggplot", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L
  )
  p <- RMdifGammaPlot(cuts, data = df, dif_var = grp)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# RMdifGamma -- cutoff validation and band flagging
# ---------------------------------------------------------------------
test_that("RMdifGamma rejects malformed cutoff arguments", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  expect_error(RMdifGamma(df, grp, cutoff = "abc"), regexp = "cutoff")
  expect_error(
    RMdifGamma(df, grp, cutoff = data.frame(Item = "Item1", gamma_low = -1)),
    regexp = "missing required columns"
  )
})

test_that("RMdifGamma errors when cutoff item names do not match data", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L
  )
  bad <- cuts$item_cutoffs
  bad$Item[1L] <- "NotAnItem"
  expect_error(RMdifGamma(df, grp, cutoff = bad), regexp = "do not match")
})

test_that("RMdifGamma band flag matches gamma vs bounds", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 10L,
    parallel = FALSE,
    seed = 1L
  )
  res <- RMdifGamma(df, grp, cutoff = cuts, output = "dataframe")
  expect_equal(
    res$flagged,
    res$gamma < res$gamma_low | res$gamma > res$gamma_high
  )
})

test_that("RMdifGamma handles missing values in data and dif_var", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  df[1:3, 1] <- NA
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  grp[4:5] <- NA
  res <- RMdifGamma(df, grp, output = "dataframe")
  expect_equal(nrow(res), ncol(df))
  kbl <- RMdifGamma(df, grp)
  expect_match(paste(kbl, collapse = "\n"), "complete cases")
})

# ---------------------------------------------------------------------
# RMdifGamma -- bootstrap p-values
# ---------------------------------------------------------------------
test_that("RMdifGamma p_value requires the full cutoff object", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  expect_error(
    RMdifGamma(df, grp, p_value = TRUE),
    regexp = "full RMdifGammaCutoff"
  )
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L
  )
  expect_error(
    RMdifGamma(df, grp, cutoff = cuts$item_cutoffs, p_value = TRUE),
    regexp = "full RMdifGammaCutoff"
  )
})

test_that("RMdifGamma p_value replaces asymptotic p columns and flags at alpha", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 10L,
    parallel = FALSE,
    seed = 1L
  )
  expect_warning(
    res <- RMdifGamma(
      df,
      grp,
      cutoff = cuts,
      p_value = TRUE,
      output = "dataframe"
    ),
    regexp = "only 10 simulation iterations"
  )
  expect_named(
    res,
    c(
      "Item",
      "gamma",
      "se",
      "lower",
      "upper",
      "gamma_low",
      "gamma_high",
      "p_gamma",
      "padj_gamma",
      "flagged"
    )
  )
  expect_false(any(c("padj_bh", "Significance") %in% names(res)))
  # p is rounded to 4 dp, so compare against the rounded floor
  B <- cuts$actual_iterations
  p_floor <- round(1 / (B + 1), 4)
  expect_true(all(res$p_gamma >= p_floor & res$p_gamma <= 1, na.rm = TRUE))
  expect_true(all(res$padj_gamma >= res$p_gamma, na.rm = TRUE))
  expect_equal(res$flagged, !is.na(res$padj_gamma) & res$padj_gamma < 0.05)

  # kable renders with the correction label
  suppressWarnings(
    kbl <- RMdifGamma(df, grp, cutoff = cuts, p_value = TRUE)
  )
  expect_s3_class(kbl, "knitr_kable")
  expect_match(paste(kbl, collapse = "\n"), "Westfall-Young")
})

test_that("RMdifGamma p_value correction variants run and 'none' is marginal", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  cuts <- RMdifGammaCutoff(
    df,
    grp,
    iterations = 10L,
    parallel = FALSE,
    seed = 1L
  )
  suppressWarnings({
    none <- RMdifGamma(
      df,
      grp,
      cutoff = cuts,
      p_value = TRUE,
      correction = "none",
      output = "dataframe"
    )
    bh <- RMdifGamma(
      df,
      grp,
      cutoff = cuts,
      p_value = TRUE,
      correction = "fdr_bh",
      output = "dataframe"
    )
  })
  expect_equal(none$p_gamma, none$padj_gamma)
  expect_true(all(bh$padj_gamma >= bh$p_gamma, na.rm = TRUE))
})

test_that("RMdifGamma validates alpha", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  expect_error(RMdifGamma(df, grp, alpha = 0), regexp = "alpha")
  expect_error(RMdifGamma(df, grp, alpha = 1.2), regexp = "alpha")
})
