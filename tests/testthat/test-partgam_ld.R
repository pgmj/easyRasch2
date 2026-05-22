# Tests for RMpartgamLD(), RMpgLDcutoff(), and RMpgLDplot()

make_dichotomous <- function(n = 200, k = 10, seed = 42L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("Item", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMpartgamLD errors when iarm is not installed", {
  skip_if(requireNamespace("iarm", quietly = TRUE),
          "iarm is installed; skipping missing-package error path")
  df <- make_dichotomous()
  expect_error(RMpartgamLD(df), regexp = "iarm")
})

test_that("RMpartgamLD errors when data has non-zero minimum", {
  skip_if_not_installed("iarm")
  df <- make_dichotomous() + 1L
  expect_error(RMpartgamLD(df), regexp = "scored starting at 0")
})

# ---------------------------------------------------------------------
# RMpartgamLD output structure
# ---------------------------------------------------------------------
test_that("RMpartgamLD output = 'dataframe' returns a list of two long-format tables", {
  skip_if_not_installed("iarm")
  df  <- make_dichotomous()
  res <- RMpartgamLD(df, output = "dataframe")
  # Two perspectives: each item paired with every other item (both directions)
  expect_type(res, "list")
  expect_true(length(res) >= 1L)
  expect_s3_class(res[[1L]], "data.frame")
  expect_true(all(c("Item1", "Item2", "gamma", "padj_bh") %in%
                  names(res[[1L]])))
  # 10 items -> 45 unique pairs
  expect_equal(nrow(res[[1L]]), choose(ncol(df), 2L))
})

test_that("RMpartgamLD output = 'kable' returns a kable-like object", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  out <- RMpartgamLD(df, output = "kable")
  # RMpartgamLD glues two kable tables into one knit_asis string.
  expect_true(inherits(out, "knitr_kable") || inherits(out, "knit_asis"))
})

test_that("RMpartgamLD kable output contains two clean pipe tables", {
  # Regression test: an earlier implementation used
  #   paste(kable1, "\n\n", kable2)
  # which silently interleaved the two multi-line character vectors
  # row-by-row, producing an unreadable corrupt table in vignettes.
  # We now collapse each table first, then join with a blank line.
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df  <- make_dichotomous()
  out <- RMpartgamLD(df, output = "kable")
  txt <- as.character(out)
  # Exactly two "Table:" captions (one per rest-score direction)
  expect_equal(length(gregexpr("\nTable: ", paste0("\n", txt))[[1]]), 2L)
  # Exactly two header rows "|Item 1 |Item 2 ..."
  expect_equal(length(gregexpr("|Item 1 |Item 2", txt,
                               fixed = TRUE)[[1]]), 2L)
  # No row should contain the column header twice (interleaving symptom)
  expect_false(any(grepl("Item 1.*Item 1", strsplit(txt, "\n")[[1]])))
})

# ---------------------------------------------------------------------
# RMpgLDcutoff -- small iterations
# ---------------------------------------------------------------------
test_that("RMpgLDcutoff returns a list with pair_cutoffs + actual_iterations", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df  <- make_dichotomous()
  res <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true("pair_cutoffs" %in% names(res))
  expect_s3_class(res$pair_cutoffs, "data.frame")
  expect_equal(nrow(res$pair_cutoffs), choose(ncol(df), 2L))
  expect_true(all(c("Item1", "Item2", "gamma_low",
                    "gamma_high") %in% names(res$pair_cutoffs)))
  expect_true(res$actual_iterations >= 1L)
})

test_that("RMpartgamLD accepts an RMpgLDcutoff result and adds flagged column", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df  <- make_dichotomous()
  cuts <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  res  <- RMpartgamLD(df, cutoff = cuts, output = "dataframe")
  expect_type(res, "list")
  expect_true("flagged" %in% names(res[[1L]]))
  expect_type(res[[1L]]$flagged, "logical")
})

test_that("RMpgLDcutoff is reproducible with the same seed", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  df <- make_dichotomous()
  r1 <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$pair_cutoffs, r2$pair_cutoffs)
})

# ---------------------------------------------------------------------
# RMpgLDplot
# ---------------------------------------------------------------------
test_that("RMpgLDplot returns a ggplot from an RMpgLDcutoff result", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df   <- make_dichotomous()
  cuts <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  p    <- RMpgLDplot(cuts)
  expect_s3_class(p, "ggplot")
})

test_that("RMpgLDplot with observed data overlay returns a ggplot", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  df   <- make_dichotomous()
  cuts <- RMpgLDcutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  p    <- RMpgLDplot(cuts, data = df)
  expect_s3_class(p, "ggplot")
})
