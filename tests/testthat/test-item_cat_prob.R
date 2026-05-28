# Tests for RMitemCatProb()

make_polytomous <- function(n = 200, k = 4, seed = 1L, max_score = 3L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:max_score, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_dichotomous <- function(n = 200, k = 5, seed = 2L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMitemCatProb errors on min != 0", {
  skip_if_not_installed("eRm")
  df <- make_polytomous() + 1L
  expect_error(RMitemCatProb(df), regexp = "scored starting at 0")
})

test_that("RMitemCatProb errors when fewer than 2 items", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()[, 1L, drop = FALSE]
  expect_error(RMitemCatProb(df), regexp = "at least 2 columns")
})

test_that("RMitemCatProb errors on invalid theta_range", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  expect_error(RMitemCatProb(df, theta_range = c(2, -2)),
               regexp = "theta_range\\[1\\] < theta_range\\[2\\]")
  expect_error(RMitemCatProb(df, theta_range = 1),
               regexp = "length 2")
})

test_that("RMitemCatProb errors on bad n_points", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  expect_error(RMitemCatProb(df, n_points = 1L), regexp = "integer >= 2")
  expect_error(RMitemCatProb(df, n_points = 50.5), regexp = "integer >= 2")
})

test_that("RMitemCatProb errors on bad viridis_end", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  expect_error(RMitemCatProb(df, viridis_end = 0),  regexp = "in \\(0, 1\\]")
  expect_error(RMitemCatProb(df, viridis_end = 1.5), regexp = "in \\(0, 1\\]")
})

test_that("RMitemCatProb errors when item_labels length mismatches", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  expect_error(RMitemCatProb(df, item_labels = c("a", "b")),
               regexp = "same length")
})

test_that("RMitemCatProb errors when category_labels length mismatches", {
  skip_if_not_installed("eRm")
  df <- make_polytomous(max_score = 3L)  # max=3 -> 4 categories
  expect_error(RMitemCatProb(df, category_labels = c("low", "high")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# Polytomous output
# ---------------------------------------------------------------------
test_that("RMitemCatProb returns a ggplot for polytomous data", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemCatProb(df)
  expect_s3_class(p, "ggplot")
})

test_that("RMitemCatProb dataframe output has correct structure + sums to 1", {
  skip_if_not_installed("eRm")
  df <- make_polytomous(k = 4L, max_score = 3L)
  res <- RMitemCatProb(df, n_points = 50L, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_named(res, c("Item", "Category", "Theta", "Probability"))
  # 4 items * 4 categories (0..3) * 50 theta points = 800 rows
  expect_equal(nrow(res), 4L * 4L * 50L)
  expect_s3_class(res$Item, "factor")
  # Per-(Item, Theta) cell, probabilities sum to exactly 1
  sums <- tapply(res$Probability, list(res$Item, res$Theta), sum)
  expect_true(all(abs(sums - 1) < 1e-10))
  # Probabilities in [0, 1]
  expect_true(all(res$Probability >= 0 & res$Probability <= 1))
})

test_that("RMitemCatProb item_labels are applied to facet titles", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()
  labs <- paste("Custom", seq_len(ncol(df)))
  res <- RMitemCatProb(df, item_labels = labs, output = "dataframe")
  # Each item facet level should contain the custom label
  uniq_levels <- levels(res$Item)
  expect_true(all(grepl("Custom", uniq_levels)))
})

# ---------------------------------------------------------------------
# Dichotomous fallback (RM)
# ---------------------------------------------------------------------
test_that("RMitemCatProb handles dichotomous data via RM", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  p <- RMitemCatProb(df)
  expect_s3_class(p, "ggplot")
  res <- RMitemCatProb(df, n_points = 30L, output = "dataframe")
  # 5 items * 2 categories * 30 theta points = 300 rows
  expect_equal(nrow(res), 5L * 2L * 30L)
  # Sums to 1 per (Item, Theta)
  sums <- tapply(res$Probability, list(res$Item, res$Theta), sum)
  expect_true(all(abs(sums - 1) < 1e-10))
})

# ---------------------------------------------------------------------
# Monotonic boundaries: at the extreme low theta, P(Category 0) -> 1
# At the extreme high theta, P(Category max) -> 1
# ---------------------------------------------------------------------
# ---------------------------------------------------------------------
# Path mode (geomtextpath single-item branch)
# ---------------------------------------------------------------------
test_that("RMitemCatProb path mode requires the `item` argument", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous()
  expect_error(
    RMitemCatProb(df, label_curves = "path"),
    regexp = "requires an `item` argument"
  )
})

test_that("RMitemCatProb path mode rejects unknown item name", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous()
  expect_error(
    RMitemCatProb(df, label_curves = "path", item = "Z9"),
    regexp = "not found in data"
  )
})

test_that("RMitemCatProb path mode rejects out-of-range item index", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous()
  expect_error(
    RMitemCatProb(df, label_curves = "path", item = 999L),
    regexp = "index must be an integer"
  )
  expect_error(
    RMitemCatProb(df, label_curves = "path", item = 1.5),
    regexp = "index must be an integer"
  )
})

test_that("RMitemCatProb path mode returns a ggplot for polytomous data", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous()
  p <- RMitemCatProb(df, label_curves = "path", item = "I1")
  expect_s3_class(p, "ggplot")
})

test_that("RMitemCatProb path mode also works for dichotomous (no middle tier)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("geomtextpath")
  df <- make_dichotomous()
  p <- RMitemCatProb(df, label_curves = "path", item = "I1")
  expect_s3_class(p, "ggplot")
})

test_that("RMitemCatProb path mode dataframe output is subset to the chosen item", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df  <- make_polytomous(k = 4L, max_score = 3L)
  res <- RMitemCatProb(df, label_curves = "path", item = "I2",
                       n_points = 30L, output = "dataframe")
  # 1 item * 4 categories * 30 theta points
  expect_equal(nrow(res), 1L * 4L * 30L)
  expect_true(all(grepl("^I2", as.character(res$Item))))
})

test_that("RMitemCatProb path mode default theta_range widens to (-5, 5)", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df  <- make_polytomous()
  res <- RMitemCatProb(df, label_curves = "path", item = "I1",
                       output = "dataframe")
  expect_equal(range(res$Theta), c(-5, 5))

  # User-supplied theta_range should still be honoured
  res2 <- RMitemCatProb(df, label_curves = "path", item = "I1",
                        theta_range = c(-3, 3), output = "dataframe")
  expect_equal(range(res2$Theta), c(-3, 3))
})

test_that("RMitemCatProb path mode validates text_size", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous()
  expect_error(RMitemCatProb(df, label_curves = "path", item = "I1",
                             text_size = -1),
               regexp = "positive number")
  expect_error(RMitemCatProb(df, label_curves = "path", item = "I1",
                             text_size = c(3, 4)),
               regexp = "positive number")
})

test_that("RMitemCatProb path mode produces one geom_textpath layer per category", {
  # Sanity check the new per-category hjust strategy: there should be
  # exactly (max_score + 1) geom_textpath layers, one per response
  # category. This catches regressions if the layer loop is ever
  # changed back to a single layer.
  skip_if_not_installed("eRm")
  skip_if_not_installed("geomtextpath")
  df <- make_polytomous(max_score = 3L)  # categories 0..3 -> 4 layers
  p  <- RMitemCatProb(df, label_curves = "path", item = "I1")
  geom_classes <- vapply(
    p$layers,
    function(l) class(l$geom)[1L],
    character(1L)
  )
  expect_equal(sum(geom_classes == "GeomTextpath"), 4L)
})

test_that("RMitemCatProb boundary probabilities behave correctly", {
  skip_if_not_installed("eRm")
  df  <- make_polytomous(max_score = 3L)
  res <- RMitemCatProb(df, theta_range = c(-10, 10), n_points = 100L,
                       output = "dataframe")
  # At theta = -10, category 0 should dominate; at theta = +10, category 3 should
  low  <- res[res$Theta == min(res$Theta), ]
  high <- res[res$Theta == max(res$Theta), ]
  # Per item: argmax of probability at low end should be category 0
  for (it in unique(low$Item)) {
    sub <- low[low$Item == it, ]
    expect_equal(sub$Category[which.max(sub$Probability)], 0L,
                 info = paste("Item:", it))
  }
  for (it in unique(high$Item)) {
    sub <- high[high$Item == it, ]
    expect_equal(sub$Category[which.max(sub$Probability)], 3L,
                 info = paste("Item:", it))
  }
})
