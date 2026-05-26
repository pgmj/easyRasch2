test_that("RMlocdepQ3 no longer errors when cutoff is missing", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  # NULL cutoff should return a knitr_kable, not an error
  expect_s3_class(RMlocdepQ3(df), "knitr_kable")
})

test_that("RMlocdepQ3 with cutoff = NULL returns raw Q3 kable", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL)
  expect_s3_class(result, "knitr_kable")
  # Caption should mention "Raw Q3"
  cap <- attr(result, "caption")
  if (is.null(cap)) cap <- paste(as.character(result), collapse = "\n")
  expect_true(grepl("Raw Q3", cap))
})

test_that("RMlocdepQ3 with cutoff = NULL returns raw Q3 dataframe", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = NULL, output = "dataframe")
  expect_s3_class(result, "data.frame")
  expect_true(all(is.na(diag(as.matrix(result)))))
  expect_true(all(is.na(result[upper.tri(result)])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 errors when data has non-zero minimum", {
  df <- as.data.frame(matrix(sample(1:3, 100, replace = TRUE), nrow = 20, ncol = 5))
  expect_error(RMlocdepQ3(df, cutoff = 0.2), regexp = "scored starting at 0")
})

test_that("RMlocdepQ3 errors when data is not a data.frame or matrix", {
  expect_error(RMlocdepQ3(list(a = 1:5), cutoff = 0.2), regexp = "data.frame or matrix")
})

test_that("RMlocdepQ3 returns a data.frame when output = 'dataframe'", {
  skip_if_not_installed("mirt")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "dataframe")
  expect_s3_class(result, "data.frame")
  # Lower triangle: off-diagonal lower should be numeric, upper + diag should be NA
  item_cols <- setdiff(names(result), "above_cutoff")
  expect_true(all(is.na(diag(as.matrix(result[item_cols])))))
  expect_true(all(is.na(result[item_cols][upper.tri(result[item_cols])])))
  expect_false(all(is.na(result[lower.tri(result)])))
})

test_that("RMlocdepQ3 returns a knitr_kable object when output = 'kable'", {
  skip_if_not_installed("mirt")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  result <- RMlocdepQ3(df, cutoff = 0.2, output = "kable")
  expect_s3_class(result, "knitr_kable")
  # Caption should mention dynamic cut-off
  cap <- attr(result, "caption")
  if (is.null(cap)) cap <- paste(as.character(result), collapse = "\n")
  expect_true(grepl("Dynamic cut-off", cap))
})

# ---------------------------------------------------------------------
# RMlocdepQ3cutoff -- small iterations
# ---------------------------------------------------------------------
test_that("RMlocdepQ3cutoff returns the expected list structure", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  res <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(c("results", "p95", "p99", "p995", "p999",
                    "suggested_cutoff", "actual_iterations",
                    "sample_n") %in% names(res)))
  expect_true(res$actual_iterations >= 1L)
  expect_true(is.finite(res$suggested_cutoff))
})

test_that("RMlocdepQ3cutoff result is consumable by RMlocdepQ3", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  cu  <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  out <- RMlocdepQ3(df, cutoff = cu$suggested_cutoff)
  expect_s3_class(out, "knitr_kable")
})

test_that("RMlocdepQ3cutoff is reproducible with the same seed", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  r1 <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$results, r2$results)
})

# ---------------------------------------------------------------------
# Per-pair additions (pair_results / pair_cutoffs / item_names)
# ---------------------------------------------------------------------
test_that("RMlocdepQ3cutoff returns pair_results / pair_cutoffs / item_names", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  res <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)

  # New top-level slots
  expect_true(all(c("pair_results", "pair_cutoffs", "item_names",
                    "cutoff_method", "hdci_width") %in% names(res)))

  # pair_results: long format, one row per pair per successful iteration
  expect_s3_class(res$pair_results, "data.frame")
  expect_true(all(c("Item1", "Item2", "Q3", "iteration") %in%
                  names(res$pair_results)))
  expect_equal(nrow(res$pair_results),
               choose(ncol(df), 2L) * res$actual_iterations)
  expect_setequal(unique(c(res$pair_results$Item1, res$pair_results$Item2)),
                  colnames(df))

  # pair_cutoffs: one row per unordered pair
  expect_s3_class(res$pair_cutoffs, "data.frame")
  expect_true(all(c("Item1", "Item2", "Q3_low", "Q3_high") %in%
                  names(res$pair_cutoffs)))
  expect_equal(nrow(res$pair_cutoffs), choose(ncol(df), 2L))
  expect_true(all(res$pair_cutoffs$Q3_low <= res$pair_cutoffs$Q3_high))

  # item_names preserves the data colnames
  expect_equal(res$item_names, colnames(df))
})

test_that("RMlocdepQ3cutoff cutoff_method = 'quantile' works without ggdist", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  res <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L,
                          cutoff_method = "quantile")
  expect_equal(res$cutoff_method, "quantile")
  expect_true(all(c("Q3_low", "Q3_high") %in% names(res$pair_cutoffs)))
})

# ---------------------------------------------------------------------
# RMlocdepQ3 accepts the full simfit object
# ---------------------------------------------------------------------
test_that("RMlocdepQ3 extracts $suggested_cutoff when given full simfit", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu  <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  # Pass the FULL list — should match passing just the scalar
  a <- RMlocdepQ3(df, cutoff = cu,                  output = "dataframe")
  b <- RMlocdepQ3(df, cutoff = cu$suggested_cutoff, output = "dataframe")
  expect_equal(a, b)
})

# ---------------------------------------------------------------------
# RMlocdepQ3plot
# ---------------------------------------------------------------------
test_that("RMlocdepQ3plot returns a ggplot in both cases", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  # Case 1: no data
  p1 <- RMlocdepQ3plot(cu)
  expect_s3_class(p1, "ggplot")

  # Case 2: observed data overlay
  p2 <- RMlocdepQ3plot(cu, data = df)
  expect_s3_class(p2, "ggplot")
})

test_that("RMlocdepQ3plot n_pairs trims and orders by deviation", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  p <- RMlocdepQ3plot(cu, data = df, n_pairs = 3L)
  expect_s3_class(p, "ggplot")
  # Pair factor should have exactly 3 levels
  expect_equal(nlevels(p$data$Pair), 3L)
})

test_that("RMlocdepQ3plot validates inputs", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  # Missing required simfit slots
  expect_error(RMlocdepQ3plot(list(results = data.frame())),
               regexp = "missing required components")
  # n_pairs validation
  expect_error(RMlocdepQ3plot(cu, n_pairs = 0L),  regexp = "positive integer")
  expect_error(RMlocdepQ3plot(cu, n_pairs = 1.5), regexp = "positive integer")
  # items validation
  expect_error(RMlocdepQ3plot(cu, items = "no_such_item"),
               regexp = "Unknown item")
  expect_error(RMlocdepQ3plot(cu, items = c("I1")),
               regexp = "at least 2 item names")
})
