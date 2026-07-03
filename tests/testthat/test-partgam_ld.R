# Tests for RMlocdepGamma(), RMlocdepGammaCutoff(), and RMlocdepGammaPlot()

make_dichotomous <- function(n = 200, k = 10, seed = 42L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("Item", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMlocdepGamma errors when iarm is not installed", {
  skip_if(
    requireNamespace("iarm", quietly = TRUE),
    "iarm is installed; skipping missing-package error path"
  )
  df <- make_dichotomous()
  expect_error(RMlocdepGamma(df), regexp = "iarm")
})

test_that("RMlocdepGamma errors when data has non-zero minimum", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous() + 1L
  expect_error(RMlocdepGamma(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# RMlocdepGamma output structure
# ---------------------------------------------------------------------
test_that("RMlocdepGamma output = 'dataframe' returns a list of two long-format tables", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  res <- RMlocdepGamma(df, output = "dataframe")
  # Two perspectives: each item paired with every other item (both directions)
  expect_type(res, "list")
  expect_true(length(res) >= 1L)
  expect_s3_class(res[[1L]], "data.frame")
  expect_true(all(
    c("Item1", "Item2", "gamma", "padj_bh", "Significance") %in%
      names(res[[1L]])
  ))
  # Significance is a star-string from iarm; should be one of the valid markers
  expect_true(all(
    res[[1L]]$Significance %in%
      c("", ".", "*", "**", "***")
  ))
  # 10 items -> 45 unique pairs
  expect_equal(nrow(res[[1L]]), choose(ncol(df), 2L))
})

test_that("RMlocdepGamma output = 'kable' returns a custom-class list", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  out <- RMlocdepGamma(df, output = "kable")
  expect_s3_class(out, "RMlocdepGamma")
  expect_true(all(c("direction1", "direction2", ".combined") %in% names(out)))
  expect_s3_class(out$direction1, "knitr_kable")
  expect_s3_class(out$direction2, "knitr_kable")
})

test_that("RMlocdepGamma knit_print emits two clean pipe tables", {
  # Regression test for the old asis-blob bug: an earlier implementation
  # used `paste(kable1, "\n\n", kable2)` which silently interleaved the
  # two multi-line character vectors row-by-row. The combined asis string
  # is now stored in $.combined and surfaced via knit_print.RMlocdepGamma().
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  out <- RMlocdepGamma(df, output = "kable")
  asis <- knitr::knit_print(out)
  expect_s3_class(asis, "knit_asis")
  txt <- as.character(asis)
  # Exactly two "Table:" captions (one per rest-score direction)
  expect_equal(length(gregexpr("\nTable: ", paste0("\n", txt))[[1]]), 2L)
  # Exactly two header rows "|Item 1 |Item 2 ..."
  expect_equal(length(gregexpr("|Item 1 |Item 2", txt, fixed = TRUE)[[1]]), 2L)
  # No row should contain the column header twice (interleaving symptom)
  expect_false(any(grepl("Item 1.*Item 1", strsplit(txt, "\n")[[1]])))
})

test_that("RMlocdepGamma print method emits both tables", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df <- make_dichotomous()
  out <- RMlocdepGamma(df, output = "kable")
  printed <- capture.output(print(out))
  joined <- paste(printed, collapse = "\n")
  expect_equal(length(gregexpr("\nTable: ", paste0("\n", joined))[[1]]), 2L)
  expect_equal(
    length(gregexpr("|Item 1 |Item 2", joined, fixed = TRUE)[[1]]),
    2L
  )
})

# ---------------------------------------------------------------------
# RMlocdepGammaCutoff -- small iterations
# ---------------------------------------------------------------------
test_that("RMlocdepGammaCutoff returns a list with pair_cutoffs + actual_iterations", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  res <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true("pair_cutoffs" %in% names(res))
  expect_s3_class(res$pair_cutoffs, "data.frame")
  expect_equal(nrow(res$pair_cutoffs), choose(ncol(df), 2L))
  expect_true(all(
    c("Item1", "Item2", "gamma_low", "gamma_high") %in% names(res$pair_cutoffs)
  ))
  expect_true(res$actual_iterations >= 1L)
})

test_that("RMlocdepGamma n_pairs trims to top-N by |gamma| per direction", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous() # 10 items -> 45 unique pairs per direction
  res_all <- RMlocdepGamma(df, output = "dataframe")
  expect_equal(nrow(res_all[[1L]]), choose(ncol(df), 2L))

  res_top <- RMlocdepGamma(df, output = "dataframe", n_pairs = 5L)
  expect_equal(nrow(res_top[[1L]]), 5L)
  expect_equal(nrow(res_top[[2L]]), 5L)

  # Each direction should contain exactly the 5 pairs with the largest |gamma|
  expected_top <- res_all[[1L]][
    order(abs(res_all[[1L]]$gamma), decreasing = TRUE)[seq_len(5L)],
    c("Item1", "Item2"),
    drop = FALSE
  ]
  expect_setequal(
    paste(res_top[[1L]]$Item1, res_top[[1L]]$Item2),
    paste(expected_top$Item1, expected_top$Item2)
  )

  # |gamma| within the returned top should be non-increasing
  expect_true(all(diff(abs(res_top[[1L]]$gamma)) <= 0))
})

test_that("RMlocdepGamma n_pairs > total pairs silently returns all pairs", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  res <- RMlocdepGamma(df, output = "dataframe", n_pairs = 10000L)
  expect_equal(nrow(res[[1L]]), choose(ncol(df), 2L))
})

test_that("RMlocdepGamma validates n_pairs", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  expect_error(RMlocdepGamma(df, n_pairs = 0L), regexp = "positive integer")
  expect_error(RMlocdepGamma(df, n_pairs = -3L), regexp = "positive integer")
  expect_error(RMlocdepGamma(df, n_pairs = 1.5), regexp = "positive integer")
  expect_error(
    RMlocdepGamma(df, n_pairs = c(1, 2)),
    regexp = "positive integer"
  )
})

test_that("RMlocdepGamma accepts an RMlocdepGammaCutoff result and adds flagged column", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  res <- RMlocdepGamma(df, cutoff = cuts, output = "dataframe")
  expect_type(res, "list")
  expect_true("flagged" %in% names(res[[1L]]))
  expect_type(res[[1L]]$flagged, "logical")
})

test_that("RMlocdepGammaCutoff is reproducible with the same seed", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  r1 <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$pair_cutoffs, r2$pair_cutoffs)
})

# ---------------------------------------------------------------------
# RMlocdepGammaPlot
# ---------------------------------------------------------------------
test_that("RMlocdepGammaPlot returns a ggplot from an RMlocdepGammaCutoff result", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  p <- RMlocdepGammaPlot(cuts)
  expect_s3_class(p, "ggplot")
})

test_that("RMlocdepGammaPlot with observed data overlay returns a ggplot", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  p <- RMlocdepGammaPlot(cuts, data = df)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# RMlocdepGamma -- cutoff validation
# ---------------------------------------------------------------------
test_that("RMlocdepGamma rejects malformed cutoff arguments", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  expect_error(RMlocdepGamma(df, cutoff = "abc"), regexp = "cutoff")
  expect_error(
    RMlocdepGamma(df, cutoff = data.frame(Item1 = "a", Item2 = "b")),
    regexp = "missing required columns"
  )
})

test_that("RMlocdepGamma band flag matches gamma vs bounds in both directions", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  res <- RMlocdepGamma(df, cutoff = cuts, output = "dataframe")
  for (dir in res) {
    expect_equal(
      dir$flagged,
      !is.na(dir$gamma_low) &
        (dir$gamma < dir$gamma_low | dir$gamma > dir$gamma_high)
    )
  }
})

test_that("RMlocdepGamma handles missing values in data", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  df[1:3, 1] <- NA
  res <- RMlocdepGamma(df, output = "dataframe")
  expect_equal(nrow(res$direction1), choose(ncol(df), 2))
  kbl <- RMlocdepGamma(df)
  expect_match(kbl$.combined, "complete cases")
})

# ---------------------------------------------------------------------
# RMlocdepGamma -- bootstrap p-values
# ---------------------------------------------------------------------
test_that("RMlocdepGamma p_value requires the full cutoff object", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  expect_error(
    RMlocdepGamma(df, p_value = TRUE),
    regexp = "full RMlocdepGammaCutoff"
  )
  cuts <- RMlocdepGammaCutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_error(
    RMlocdepGamma(df, cutoff = cuts$pair_cutoffs, p_value = TRUE),
    regexp = "full RMlocdepGammaCutoff"
  )
})

test_that("RMlocdepGamma p_value adds one-sided p per pair, mirrored across directions", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  expect_warning(
    res <- RMlocdepGamma(
      df,
      cutoff = cuts,
      p_value = TRUE,
      output = "dataframe"
    ),
    regexp = "only 10 simulation iterations"
  )
  expected_cols <- c(
    "Item1",
    "Item2",
    "gamma",
    "gamma_low",
    "gamma_high",
    "p_gamma",
    "padj_gamma",
    "flagged"
  )
  expect_named(res$direction1, expected_cols)
  expect_named(res$direction2, expected_cols)
  expect_false(any(c("padj_bh", "Significance") %in% names(res$direction1)))

  # p is rounded to 4 dp, so compare against the rounded floor
  B <- cuts$actual_iterations
  p_floor <- round(1 / (B + 1), 4)
  expect_true(all(
    res$direction1$p_gamma >= p_floor &
      res$direction1$p_gamma <= 1,
    na.rm = TRUE
  ))
  expect_true(all(
    res$direction1$padj_gamma >= res$direction1$p_gamma,
    na.rm = TRUE
  ))

  # One test per pair: p/padj identical across the two direction tables
  key <- function(d) paste(pmin(d$Item1, d$Item2), pmax(d$Item1, d$Item2))
  m <- match(key(res$direction2), key(res$direction1))
  expect_equal(res$direction2$p_gamma, res$direction1$p_gamma[m])
  expect_equal(res$direction2$padj_gamma, res$direction1$padj_gamma[m])
  expect_equal(res$direction2$flagged, res$direction1$flagged[m])

  # kable renders with the correction label and canonical-direction note
  suppressWarnings(kbl <- RMlocdepGamma(df, cutoff = cuts, p_value = TRUE))
  expect_s3_class(kbl, "RMlocdepGamma")
  expect_match(kbl$.combined, "Westfall-Young")
  expect_match(kbl$.combined, "canonical direction")
})

test_that("RMlocdepGamma p_value correction runs on the full family before n_pairs", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  suppressWarnings({
    full <- RMlocdepGamma(
      df,
      cutoff = cuts,
      p_value = TRUE,
      output = "dataframe"
    )
    top3 <- RMlocdepGamma(
      df,
      cutoff = cuts,
      p_value = TRUE,
      n_pairs = 3,
      output = "dataframe"
    )
  })
  expect_equal(nrow(top3$direction1), 3L)
  key <- function(d) paste(pmin(d$Item1, d$Item2), pmax(d$Item1, d$Item2))
  m <- match(key(top3$direction1), key(full$direction1))
  expect_equal(top3$direction1$padj_gamma, full$direction1$padj_gamma[m])
})

test_that("RMlocdepGamma p_value correction variants run and 'none' is marginal", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  cuts <- RMlocdepGammaCutoff(df, iterations = 10L, parallel = FALSE, seed = 1L)
  suppressWarnings({
    none <- RMlocdepGamma(
      df,
      cutoff = cuts,
      p_value = TRUE,
      correction = "none",
      output = "dataframe"
    )
    by <- RMlocdepGamma(
      df,
      cutoff = cuts,
      p_value = TRUE,
      correction = "fdr_by",
      output = "dataframe"
    )
  })
  expect_equal(none$direction1$p_gamma, none$direction1$padj_gamma)
  expect_true(all(
    by$direction1$padj_gamma >= by$direction1$p_gamma,
    na.rm = TRUE
  ))
})

test_that("RMlocdepGamma validates alpha", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous()
  expect_error(RMlocdepGamma(df, alpha = -1), regexp = "alpha")
  expect_error(RMlocdepGamma(df, alpha = "x"), regexp = "alpha")
})
