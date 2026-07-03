# Tests for RMdimCFACutoff() / RMdimCFA() / RMdimCFAPlot()
#
# Iterations are intentionally very small (5-30) so the file runs in
# a few seconds. Each test uses parallel = FALSE because parallel
# requires mirai daemons and complicates CRAN test runs.

# ---------------------------------------------------------------------
# Validation errors
# ---------------------------------------------------------------------
test_that("RMdimCFACutoff errors when lavaan is not installed", {
  skip_if(requireNamespace("lavaan", quietly = TRUE),
          "lavaan is installed; skipping missing-package error path test")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  expect_error(RMdimCFACutoff(df), regexp = "lavaan")
})

test_that("RMdimCFACutoff errors when data has non-zero minimum", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(1:3, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  expect_error(RMdimCFACutoff(df), regexp = "scored starting at 0")
})

test_that("RMdimCFACutoff errors when fewer than 3 items", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 100, replace = TRUE),
                             nrow = 50, ncol = 2))
  expect_error(RMdimCFACutoff(df), regexp = "at least 3 items")
})

test_that("RMdimCFACutoff rejects ML-family estimators", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  for (est in c("ML", "MLR", "MLM", "MLF")) {
    expect_error(
      RMdimCFACutoff(df, estimator = est, iterations = 5, parallel = FALSE),
      regexp = "not appropriate for ordinal-CFA"
    )
  }
})

test_that("RMdimCFACutoff rejects out-of-range percentile", {
  skip_if_not_installed("lavaan")
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE),
                             nrow = 40, ncol = 5))
  for (bad in c(0, 50, 100, 101, NA, -1)) {
    expect_error(
      RMdimCFACutoff(df, percentile = bad, iterations = 5, parallel = FALSE),
      regexp = "percentile"
    )
  }
})

test_that("RMdimCFACutoff(output = 'kable') errors with a migration message", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")
  expect_error(
    RMdimCFACutoff(raschdat1[, 1:8], output = "kable",
                   iterations = 5, parallel = FALSE, seed = 1),
    regexp = "RMdimCFA"
  )
})

# ---------------------------------------------------------------------
# RMdimCFACutoff output structure (simulate-only list)
# ---------------------------------------------------------------------
test_that("RMdimCFACutoff returns the simulate-only list on dichotomous data", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 10,
                        parallel = FALSE, seed = 1)

  expect_type(res, "list")
  expected_names <- c("simulated", "simulated_loadings", "percentile",
                      "cutoffs", "loading_cutoffs", "actual_iterations",
                      "sample_n", "n_items", "item_names", "is_polytomous",
                      "estimator")
  expect_true(all(expected_names %in% names(res)))
  # No observed/flagged any more -- those moved to RMdimCFA()
  expect_false(any(c("observed", "flagged") %in% names(res)))

  expect_s3_class(res$simulated, "data.frame")
  expect_named(res$simulated, c("iteration", "cfi", "rmsea", "srmr"))
  expect_s3_class(res$simulated_loadings, "data.frame")
  expect_equal(names(res$simulated_loadings),
               c("iteration", res$item_names))
  expect_equal(nrow(res$simulated), nrow(res$simulated_loadings))

  expect_named(res$cutoffs, c("cfi", "rmsea", "srmr"))
  expect_s3_class(res$loading_cutoffs, "data.frame")
  expect_named(res$loading_cutoffs, c("Item", "low", "high"))
  expect_equal(nrow(res$loading_cutoffs), 8L)
  expect_true(all(res$loading_cutoffs$low <= res$loading_cutoffs$high))

  expect_equal(res$percentile, 99)
  expect_false(res$is_polytomous)
  expect_equal(res$estimator, "WLSMV")
  expect_equal(res$n_items, 8L)
})

test_that("Cutoffs change with different percentile values, same simulation", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res99 <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 30, parallel = FALSE,
                          seed = 1, percentile = 99)
  res95 <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 30, parallel = FALSE,
                          seed = 1, percentile = 95)

  expect_equal(res99$simulated, res95$simulated)
  expect_equal(res99$simulated_loadings, res95$simulated_loadings)
  expect_gte(res99$cutoffs[["rmsea"]], res95$cutoffs[["rmsea"]])
  expect_lte(res99$cutoffs[["cfi"]],   res95$cutoffs[["cfi"]])
  # Wider percentile => wider expected loading interval
  expect_true(all(res99$loading_cutoffs$low  <= res95$loading_cutoffs$low))
  expect_true(all(res99$loading_cutoffs$high >= res95$loading_cutoffs$high))
})

test_that("RMdimCFACutoff auto-picks PCM on polytomous data", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("pcmdat2", package = "eRm")

  res <- RMdimCFACutoff(pcmdat2[1:100, ], iterations = 10,
                        parallel = FALSE, seed = 2)
  expect_true(res$is_polytomous)
  expect_true(res$actual_iterations >= 1L)
})

test_that("Same seed produces identical simulated distributions", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  res1 <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 5,
                         parallel = FALSE, seed = 42)
  res2 <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 5,
                         parallel = FALSE, seed = 42)
  expect_equal(res1$simulated, res2$simulated)
  expect_equal(res1$simulated_loadings, res2$simulated_loadings)
})

# ---------------------------------------------------------------------
# RMdimCFA (tables)
# ---------------------------------------------------------------------
test_that("RMdimCFA requires a cutoff object", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")
  expect_error(RMdimCFA(raschdat1[, 1:8]), regexp = "cutoff")
  expect_error(RMdimCFA(raschdat1[, 1:8], cutoff = list(foo = 1)),
               regexp = "RMdimCFACutoff")
})

test_that("RMdimCFA returns a list of two kables by default", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("knitr")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 20,
                        parallel = FALSE, seed = 1)
  tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim)

  expect_type(tabs, "list")
  expect_named(tabs, c("fit", "loadings"))
  expect_s3_class(tabs$fit, "knitr_kable")
  expect_s3_class(tabs$loadings, "knitr_kable")
})

test_that("RMdimCFA output = 'dataframe' returns the two tables as data.frames", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  sim  <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 20,
                         parallel = FALSE, seed = 1)
  tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim, output = "dataframe")

  expect_named(tabs, c("fit", "loadings"))
  expect_s3_class(tabs$fit, "data.frame")
  expect_named(tabs$fit, c("Index", "Observed", "Cutoff", "Direction",
                           "Flagged"))
  expect_equal(tabs$fit$Index, c("CFI", "RMSEA", "SRMR"))

  expect_s3_class(tabs$loadings, "data.frame")
  expect_named(tabs$loadings, c("Item", "Observed", "Expected_low",
                                "Expected_high", "Flagged"))
  expect_equal(nrow(tabs$loadings), 8L)
  expect_true(all(tabs$loadings$Flagged %in% c("below", "above", "")))
})

test_that("RMdimCFA p_value adds p/padj and redefines Flagged at alpha", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 20,
                        parallel = FALSE, seed = 1)
  expect_warning(
    tabs <- RMdimCFA(raschdat1[, 1:8], cutoff = sim, p_value = TRUE,
                     output = "dataframe"),
    regexp = "only 20 simulation iterations"
  )

  # p is rounded to 4 dp in the table, so compare against the rounded floor
  B <- sim$actual_iterations
  p_floor <- round(1 / (B + 1), 4)
  expect_named(tabs$fit, c("Index", "Observed", "Cutoff", "Direction",
                           "p", "padj", "Flagged"))
  expect_true(all(tabs$fit$p >= p_floor & tabs$fit$p <= 1, na.rm = TRUE))
  expect_true(all(tabs$fit$padj >= tabs$fit$p, na.rm = TRUE))

  expect_named(tabs$loadings, c("Item", "Observed", "Expected_low",
                                "Expected_high", "p_loading",
                                "padj_loading", "Flagged"))
  expect_true(all(tabs$loadings$p_loading >= p_floor &
                    tabs$loadings$p_loading <= 1, na.rm = TRUE))
  expect_true(all(tabs$loadings$Flagged %in% c("below", "above", "")))

  # kable output renders with the extra columns and correction label
  suppressWarnings(
    kbls <- RMdimCFA(raschdat1[, 1:8], cutoff = sim, p_value = TRUE)
  )
  expect_s3_class(kbls$fit, "knitr_kable")
  expect_s3_class(kbls$loadings, "knitr_kable")
  expect_match(paste(kbls$fit, collapse = "\n"), "Westfall-Young")
})

test_that("RMdimCFA p_value correction variants run and 'none' is marginal", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 20,
                        parallel = FALSE, seed = 1)
  suppressWarnings({
    none <- RMdimCFA(raschdat1[, 1:8], cutoff = sim, p_value = TRUE,
                     correction = "none", output = "dataframe")
    bh <- RMdimCFA(raschdat1[, 1:8], cutoff = sim, p_value = TRUE,
                   correction = "fdr_bh", output = "dataframe")
  })
  expect_equal(none$fit$p, none$fit$padj)
  expect_true(all(bh$loadings$padj_loading >= bh$loadings$p_loading,
                  na.rm = TRUE))
})

# ---------------------------------------------------------------------
# RMdimCFAPlot (two-plot list)
# ---------------------------------------------------------------------
test_that("RMdimCFAPlot returns a list of two ggplots", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 20,
                        parallel = FALSE, seed = 1)
  p <- RMdimCFAPlot(sim, data = raschdat1[, 1:8])

  expect_type(p, "list")
  expect_named(p, c("loadings", "fit"))
  expect_s3_class(p$loadings, "ggplot")
  expect_s3_class(p$fit, "ggplot")
})

test_that("RMdimCFAPlot requires data", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 5,
                        parallel = FALSE, seed = 1)
  expect_error(RMdimCFAPlot(sim), regexp = "data")
})

test_that("RMdimCFAPlot percentile override recomputes without re-simulating", {
  skip_if_not_installed("lavaan")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  data("raschdat1", package = "eRm")

  sim <- RMdimCFACutoff(raschdat1[, 1:8], iterations = 30,
                        parallel = FALSE, seed = 1)
  expect_s3_class(RMdimCFAPlot(sim, data = raschdat1[, 1:8],
                               percentile = 95)$fit, "ggplot")
  expect_error(RMdimCFAPlot(sim, data = raschdat1[, 1:8], percentile = 50),
               regexp = "percentile")
})

test_that("RMdimCFAPlot rejects non-RMdimCFACutoff input", {
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("ggdist")
  expect_error(RMdimCFAPlot(list(foo = 1), data = data.frame(a = 1)),
               regexp = "must be the result returned by RMdimCFACutoff")
})
