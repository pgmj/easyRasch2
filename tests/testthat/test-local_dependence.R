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
# RMlocdepQ3Cutoff -- small iterations
# ---------------------------------------------------------------------
test_that("RMlocdepQ3Cutoff returns the expected list structure", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  res <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(c("results", "p95", "p99", "p995", "p999",
                    "suggested_cutoff", "actual_iterations",
                    "sample_n") %in% names(res)))
  expect_true(res$actual_iterations >= 1L)
  expect_true(is.finite(res$suggested_cutoff))
})

test_that("RMlocdepQ3Cutoff result is consumable by RMlocdepQ3", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  cu  <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  out <- RMlocdepQ3(df, cutoff = cu$suggested_cutoff)
  expect_s3_class(out, "knitr_kable")
})

test_that("RMlocdepQ3Cutoff is reproducible with the same seed", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  r1 <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$results, r2$results)
})

# ---------------------------------------------------------------------
# Per-pair additions (pair_results / pair_cutoffs / item_names)
# ---------------------------------------------------------------------
test_that("RMlocdepQ3Cutoff returns pair_results / pair_cutoffs / item_names", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 8, replace = TRUE), 200, 8))
  colnames(df) <- paste0("I", 1:8)
  res <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 42L)

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

test_that("RMlocdepQ3Cutoff cutoff_method = 'quantile' works without ggdist", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  res <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L,
                          cutoff_method = "quantile")
  expect_equal(res$cutoff_method, "quantile")
  expect_true(all(c("Q3_low", "Q3_high") %in% names(res$pair_cutoffs)))
})

# ---------------------------------------------------------------------
# RMlocdepQ3 accepts the full simfit object
# ---------------------------------------------------------------------
test_that("full simfit returns list whose $matrix matches the scalar cutoff", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu  <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)
  # Full object -> list($matrix, $pairs); $matrix matches the scalar-cutoff matrix
  a <- RMlocdepQ3(df, cutoff = cu,                  output = "dataframe")
  b <- RMlocdepQ3(df, cutoff = cu$suggested_cutoff, output = "dataframe")
  expect_named(a, c("matrix", "pairs"))
  expect_equal(a$matrix, b)
})

# ---------------------------------------------------------------------
# RMlocdepQ3Plot
# ---------------------------------------------------------------------
test_that("RMlocdepQ3Plot returns a ggplot in both cases", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  # Case 1: no data
  p1 <- RMlocdepQ3Plot(cu)
  expect_s3_class(p1, "ggplot")

  # Case 2: observed data overlay
  p2 <- RMlocdepQ3Plot(cu, data = df)
  expect_s3_class(p2, "ggplot")
})

test_that("RMlocdepQ3Plot n_pairs trims and orders by deviation", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("scales")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  p <- RMlocdepQ3Plot(cu, data = df, n_pairs = 3L)
  expect_s3_class(p, "ggplot")
  # Pair factor should have exactly 3 levels
  expect_equal(nlevels(p$data$Pair), 3L)
})

test_that("RMlocdepQ3Plot validates inputs", {
  skip_if_not_installed("mirt")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 6, replace = TRUE), 200, 6))
  colnames(df) <- paste0("I", 1:6)
  cu <- RMlocdepQ3Cutoff(df, iterations = 5L, parallel = FALSE, seed = 1L)

  # Missing required simfit slots
  expect_error(RMlocdepQ3Plot(list(results = data.frame())),
               regexp = "missing required components")
  # n_pairs validation
  expect_error(RMlocdepQ3Plot(cu, n_pairs = 0L),  regexp = "positive integer")
  expect_error(RMlocdepQ3Plot(cu, n_pairs = 1.5), regexp = "positive integer")
  # items validation
  expect_error(RMlocdepQ3Plot(cu, items = "no_such_item"),
               regexp = "Unknown item")
  expect_error(RMlocdepQ3Plot(cu, items = c("I1")),
               regexp = "at least 2 item names")
})

# ---------------------------------------------------------------------
# RMlocdepQ3() bootstrap p-values (per item pair)
# ---------------------------------------------------------------------

q3_null_data <- function(n = 300, J = 7, seed = 11L) {
  set.seed(seed)
  theta <- rnorm(n, 0, 1.3); beta <- seq(-1.5, 1.5, length.out = J)
  df <- as.data.frame(sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j]))))
  colnames(df) <- paste0("I", seq_len(J)); df
}

test_that("full cutoff object returns list($matrix, $pairs)", {
  skip_if_not_installed("mirt"); skip_if_not_installed("ggdist")
  df  <- q3_null_data()
  sim <- RMlocdepQ3Cutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  res <- RMlocdepQ3(df, cutoff = sim, output = "dataframe")
  expect_named(res, c("matrix", "pairs"))
  expect_named(res$pairs, c("Item1", "Item2", "Observed", "Low", "High", "Flagged"))
  expect_equal(nrow(res$pairs), choose(ncol(df), 2L))          # one row per pair
  expect_true(all(res$pairs$Flagged %in% c("above", "below", "")))
})

test_that("n_pairs truncates the $pairs table", {
  skip_if_not_installed("mirt"); skip_if_not_installed("ggdist")
  df  <- q3_null_data()
  sim <- RMlocdepQ3Cutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  res <- RMlocdepQ3(df, cutoff = sim, n_pairs = 3, output = "dataframe")
  expect_equal(nrow(res$pairs), 3L)
})

test_that("p_value = TRUE adds p columns to $pairs and flags on padj", {
  skip_if_not_installed("mirt"); skip_if_not_installed("ggdist")
  df  <- q3_null_data()
  sim <- RMlocdepQ3Cutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  res <- suppressWarnings(
    RMlocdepQ3(df, cutoff = sim, p_value = TRUE, output = "dataframe")
  )
  expect_true(all(c("p_q3", "padj_q3") %in% names(res$pairs)))
  expect_true(all(res$pairs$p_q3 >= 0 & res$pairs$p_q3 <= 1, na.rm = TRUE))
  expect_true(all(res$pairs$padj_q3 >= res$pairs$p_q3 - 1e-9, na.rm = TRUE))
})

test_that("RMlocdepQ3 p_value = TRUE errors without the full cutoff object", {
  skip_if_not_installed("mirt"); skip_if_not_installed("ggdist")
  df  <- q3_null_data()
  sim <- RMlocdepQ3Cutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  expect_error(RMlocdepQ3(df, p_value = TRUE), regexp = "full RMlocdepQ3Cutoff")
  expect_error(RMlocdepQ3(df, cutoff = sim$suggested_cutoff, p_value = TRUE),
               regexp = "full RMlocdepQ3Cutoff")
})

test_that("NULL cutoff is unchanged (single square matrix)", {
  skip_if_not_installed("mirt")
  df  <- q3_null_data()
  res <- RMlocdepQ3(df, output = "dataframe")
  expect_false(is.list(res) && all(c("matrix", "pairs") %in% names(res)))
  expect_equal(nrow(res), ncol(df))                            # square matrix
})

test_that("RMlocdepQ3 detects an injected locally dependent pair", {
  skip_if_not_installed("mirt"); skip_if_not_installed("ggdist")
  set.seed(3); n <- 400; J <- 7
  theta <- rnorm(n, 0, 1.3); beta <- seq(-1.5, 1.5, length.out = J)
  m <- sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j])))
  cp <- sample(n, 0.7 * n); m[cp, 2] <- m[cp, 1]              # I1, I2 dependent
  df <- as.data.frame(m); colnames(df) <- paste0("I", seq_len(J))
  sim <- RMlocdepQ3Cutoff(df, iterations = 500, parallel = FALSE, seed = 4)
  res <- RMlocdepQ3(df, cutoff = sim, output = "dataframe")$pairs
  dep <- res[(res$Item1 == "I1" & res$Item2 == "I2") |
             (res$Item1 == "I2" & res$Item2 == "I1"), ]
  expect_identical(dep$Flagged, "above")
  # the injected pair has the largest departure -> sorted to the top row
  expect_true((res$Item1[1] == "I1" && res$Item2[1] == "I2") ||
              (res$Item1[1] == "I2" && res$Item2[1] == "I1"))
})
