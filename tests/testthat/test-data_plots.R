# Tests for RMplotBar() and RMplotStackedbar()

make_polytomous <- function(n = 200, k = 4, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# RMplotBar --- input validation
# ---------------------------------------------------------------------
test_that("RMplotBar errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotBar(matrix(0:3, 4, 4)), regexp = "data.frame")
})

test_that("RMplotBar errors when fewer than 2 items", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotBar(make_polytomous(k = 1L)),
               regexp = "at least 2 columns")
})

test_that("RMplotBar errors on non-numeric columns", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df$group <- letters[seq_len(nrow(df))]
  expect_error(RMplotBar(df), regexp = "must be numeric")
})

test_that("RMplotBar errors on non-integer values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1, 1] <- 0.5
  expect_error(RMplotBar(df), regexp = "integers")
})

test_that("RMplotBar errors when item_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotBar(make_polytomous(), item_labels = c("a", "b")),
               regexp = "same length")
})

test_that("RMplotBar errors when category_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotBar(make_polytomous(),
                         category_labels = c("Low", "High")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# RMplotBar --- output
# ---------------------------------------------------------------------
test_that("RMplotBar returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(RMplotBar(make_polytomous()), "ggplot")
})

test_that("RMplotBar accepts custom item_labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMplotBar(df, item_labels = paste("Item", seq_len(ncol(df))))
  expect_s3_class(p, "ggplot")
})

test_that("RMplotBar accepts custom category_labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMplotBar(df, category_labels = c("None", "Some", "More", "Most"))
  expect_s3_class(p, "ggplot")
})

test_that("RMplotBar handles 0/1 dichotomous data", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200 * 5, replace = TRUE), 200, 5))
  colnames(df) <- paste0("I", 1:5)
  expect_s3_class(RMplotBar(df), "ggplot")
})

test_that("RMplotBar tolerates NA values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1:5, 1] <- NA
  expect_s3_class(RMplotBar(df), "ggplot")
})

# ---------------------------------------------------------------------
# RMplotStackedbar --- input validation
# ---------------------------------------------------------------------
test_that("RMplotStackedbar errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotStackedbar(matrix(0:3, 4, 4)), regexp = "data.frame")
})

test_that("RMplotStackedbar errors when fewer than 2 items", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotStackedbar(make_polytomous(k = 1L)),
               regexp = "at least 2 columns")
})

test_that("RMplotStackedbar errors on non-numeric columns", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df$group <- letters[seq_len(nrow(df))]
  expect_error(RMplotStackedbar(df), regexp = "must be numeric")
})

test_that("RMplotStackedbar errors on non-integer values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1, 1] <- 0.5
  expect_error(RMplotStackedbar(df), regexp = "integers")
})

test_that("RMplotStackedbar errors when item_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotStackedbar(make_polytomous(),
                                item_labels = c("a", "b")),
               regexp = "same length")
})

test_that("RMplotStackedbar errors when category_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMplotStackedbar(make_polytomous(),
                                category_labels = c("Low", "High")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# RMplotStackedbar --- output
# ---------------------------------------------------------------------
test_that("RMplotStackedbar returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(RMplotStackedbar(make_polytomous()), "ggplot")
})

test_that("RMplotStackedbar accepts custom item & category labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMplotStackedbar(
    df,
    item_labels     = paste("Item", seq_len(ncol(df))),
    category_labels = c("None", "Some", "More", "Most")
  )
  expect_s3_class(p, "ggplot")
})

test_that("RMplotStackedbar supports show_percent and min_label_n", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMplotStackedbar(df, show_percent = TRUE, show_n = FALSE,
                        min_label_n = 5L)
  expect_s3_class(p, "ggplot")
})

test_that("RMplotStackedbar can suppress all labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMplotStackedbar(df, show_n = FALSE, show_percent = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("RMplotStackedbar handles 0/1 dichotomous data", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  df <- as.data.frame(matrix(sample(0:1, 200 * 5, replace = TRUE), 200, 5))
  colnames(df) <- paste0("I", 1:5)
  expect_s3_class(RMplotStackedbar(df), "ggplot")
})
