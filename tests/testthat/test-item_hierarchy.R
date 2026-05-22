# Tests for RMitemHierarchy()

make_polytomous <- function(n = 200, k = 5, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMitemHierarchy errors when data has non-zero minimum", {
  df <- make_polytomous() + 1L
  expect_error(RMitemHierarchy(df), regexp = "scored starting at 0")
})

test_that("RMitemHierarchy errors when fewer than 2 items", {
  df <- make_polytomous()[, 1L, drop = FALSE]
  expect_error(RMitemHierarchy(df), regexp = "at least 2 items")
})

test_that("RMitemHierarchy errors on dichotomous data with helpful pointer", {
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:1, 200 * 5, replace = TRUE), 200, 5))
  colnames(df) <- paste0("I", 1:5)
  expect_error(RMitemHierarchy(df),
               regexp = "polytomous data")
})

test_that("RMitemHierarchy errors on invalid sem_multiplier", {
  df <- make_polytomous()
  expect_error(RMitemHierarchy(df, sem_multiplier = -1),
               regexp = "positive numeric")
  expect_error(RMitemHierarchy(df, sem_multiplier = c(1, 2)),
               regexp = "positive numeric")
})

test_that("RMitemHierarchy errors when item_labels length mismatches ncol(data)", {
  df <- make_polytomous()
  expect_error(RMitemHierarchy(df, item_labels = c("a", "b")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# Output structures
# ---------------------------------------------------------------------
test_that("RMitemHierarchy default returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemHierarchy(df)
  expect_s3_class(p, "ggplot")
})

test_that("RMitemHierarchy show_numbers = FALSE returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemHierarchy(df, show_numbers = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("RMitemHierarchy with custom sem_multiplier returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemHierarchy(df, sem_multiplier = 1.96)
  expect_s3_class(p, "ggplot")
})

test_that("RMitemHierarchy with custom item_labels returns a ggplot", {
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMitemHierarchy(df, item_labels = paste("Item", seq_len(ncol(df))))
  expect_s3_class(p, "ggplot")
})

test_that("RMitemHierarchy actually applies item_labels to the x-axis", {
  # Regression test: scale_x_discrete(labels = ...) does a named lookup
  # against the factor levels, so axis_labels must be named by the
  # column names (item_order). Otherwise the labels silently revert to
  # the raw column names. (#item_labels-bug, 2026)
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  custom <- paste("Label", seq_len(ncol(df)))
  p <- RMitemHierarchy(df, item_labels = custom)

  sc <- ggplot2::ggplot_build(p)$plot$scales$get_scales("x")
  expect_true(all(names(sc$labels) %in% names(df)))
  # Every displayed label must contain its custom suffix
  for (i in seq_along(custom)) {
    expect_true(any(grepl(custom[i], sc$labels, fixed = TRUE)),
                info = paste("custom label not found:", custom[i]))
  }
})

# ---------------------------------------------------------------------
# Dataframe output
# ---------------------------------------------------------------------
test_that("RMitemHierarchy output = 'dataframe' returns long-format data", {
  skip_if_not_installed("eRm")
  df <- make_polytomous()  # 5 items × 3 thresholds (categories 0..3)
  res <- RMitemHierarchy(df, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), ncol(df) * 3L)  # one row per item × threshold
  expect_named(res, c("Item", "ItemLabel", "Threshold",
                      "ThresholdLocation", "ThresholdSE", "ItemLocation"))
  # Item is a factor ordered by ItemLocation
  expect_s3_class(res$Item, "factor")
  # ThresholdLocations centred at the grand mean -> mean should be ~ 0
  expect_lt(abs(mean(res$ThresholdLocation, na.rm = TRUE)), 1e-8)
  # SE values should all be positive and finite
  expect_true(all(res$ThresholdSE > 0 & is.finite(res$ThresholdSE)))
  # Per-item Location = mean of that item's centred thresholds
  per_item <- tapply(res$ThresholdLocation, res$Item, mean, na.rm = TRUE)
  expect_equal(
    as.numeric(per_item),
    as.numeric(tapply(res$ItemLocation, res$Item, unique)),
    tolerance = 1e-8
  )
})

test_that("RMitemHierarchy item factor levels are ordered by ItemLocation", {
  skip_if_not_installed("eRm")
  df  <- make_polytomous()
  res <- RMitemHierarchy(df, output = "dataframe")
  # Pull the unique (Item, ItemLocation) pairs in factor-level order
  ord_locs <- tapply(res$ItemLocation, res$Item, unique)
  expect_true(all(diff(ord_locs) >= 0),
              info = "Item factor levels should be in ascending ItemLocation order")
})
