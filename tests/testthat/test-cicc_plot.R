# Tests for RMitemICCPlot()

# Polytomous data gives enough populated total-score levels for class
# intervals and DIF groups at modest n.
make_polytomous <- function(n = 200, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation (reachable without the plotting packages)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot errors when data has non-zero minimum", {
  df <- make_polytomous() + 1L
  expect_error(RMitemICCPlot(df), regexp = "scored starting at 0")
})

test_that("RMitemICCPlot errors when fewer than 2 items", {
  df <- make_polytomous()[, 1L, drop = FALSE]
  expect_error(RMitemICCPlot(df), regexp = "at least 2 items")
})

test_that("RMitemICCPlot errors when class_intervals < 2", {
  df <- make_polytomous()
  expect_error(RMitemICCPlot(df, class_intervals = 1),
               regexp = "class_intervals")
})

test_that("RMitemICCPlot errors when conf_level is out of range", {
  df <- make_polytomous()
  expect_error(RMitemICCPlot(df, conf_level = 1.2), regexp = "conf_level")
})

test_that("RMitemICCPlot errors when dif_var length mismatches nrow(data)", {
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMitemICCPlot(df, dif_var = grp), regexp = "same length")
})

test_that("RMitemICCPlot errors when dif_var has fewer than 2 levels", {
  df  <- make_polytomous()
  grp <- factor(rep("A", nrow(df)))
  expect_error(RMitemICCPlot(df, dif_var = grp), regexp = "at least 2 distinct")
})

# ---------------------------------------------------------------------
# Output structures (no iarm needed for the non-DIF CML plot)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot default output is a patchwork/ggplot composite", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemICCPlot(df, class_intervals = 4L)
  expect_s3_class(p, "patchwork")
  expect_s3_class(p, "ggplot")
})

test_that("RMitemICCPlot output = 'list' returns one ggplot per item, named", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  plist <- RMitemICCPlot(df, class_intervals = 4L, output = "list")
  expect_type(plist, "list")
  expect_equal(length(plist), ncol(df))
  expect_equal(names(plist), names(df))
  for (p in plist) expect_s3_class(p, "ggplot")
})

test_that("RMitemICCPlot method = 'score' and error_band return a composite", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  expect_s3_class(RMitemICCPlot(df, method = "score"), "patchwork")
  expect_s3_class(RMitemICCPlot(df, error_band = TRUE, ci = FALSE), "patchwork")
})

test_that("RMitemICCPlot `items` subsets the rendered panels", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  plist <- RMitemICCPlot(df, items = c("I2", "I4"), output = "list")
  expect_equal(names(plist), c("I2", "I4"))
})

# ---------------------------------------------------------------------
# DIF mode (needs iarm for the partial-gamma magnitude)
# ---------------------------------------------------------------------
test_that("RMitemICCPlot DIF mode returns a composite with a gamma table attr", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous(n = 300)
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  p <- RMitemICCPlot(df, class_intervals = 4L, dif_var = grp)
  expect_s3_class(p, "patchwork")
  g <- attr(p, "dif_gamma")
  expect_s3_class(g, "data.frame")
  expect_true(all(c("Item", "gamma", "lower", "upper", "padj_bh") %in% names(g)))
})

# ---------------------------------------------------------------------
# NA handling
# ---------------------------------------------------------------------
test_that("RMitemICCPlot drops rows with NA in data or dif_var", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("patchwork")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous(n = 300)
  df[1:5, 1] <- NA_integer_
  grp <- factor(sample(c("A", "B"), nrow(df), replace = TRUE))
  grp[6:10] <- NA
  expect_s3_class(
    RMitemICCPlot(df, class_intervals = 4L, dif_var = grp), "patchwork"
  )
})
