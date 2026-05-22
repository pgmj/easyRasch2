# Tests for RMbarplot() and RMstackedbarplot()

make_polytomous <- function(n = 200, k = 4, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# RMbarplot --- input validation
# ---------------------------------------------------------------------
test_that("RMbarplot errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  expect_error(RMbarplot(matrix(0:3, 4, 4)), regexp = "data.frame")
})

test_that("RMbarplot errors when fewer than 2 items", {
  skip_if_not_installed("ggplot2")
  expect_error(RMbarplot(make_polytomous(k = 1L)),
               regexp = "at least 2 columns")
})

test_that("RMbarplot errors on non-numeric columns", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df$group <- letters[seq_len(nrow(df))]
  expect_error(RMbarplot(df), regexp = "must be numeric")
})

test_that("RMbarplot errors on non-integer values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1, 1] <- 0.5
  expect_error(RMbarplot(df), regexp = "integers")
})

test_that("RMbarplot errors when item_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMbarplot(make_polytomous(), item_labels = c("a", "b")),
               regexp = "same length")
})

test_that("RMbarplot errors when category_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMbarplot(make_polytomous(),
                         category_labels = c("Low", "High")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# RMbarplot --- output
# ---------------------------------------------------------------------
test_that("RMbarplot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(RMbarplot(make_polytomous()), "ggplot")
})

test_that("RMbarplot accepts custom item_labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMbarplot(df, item_labels = paste("Item", seq_len(ncol(df))))
  expect_s3_class(p, "ggplot")
})

test_that("RMbarplot accepts custom category_labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMbarplot(df, category_labels = c("None", "Some", "More", "Most"))
  expect_s3_class(p, "ggplot")
})

test_that("RMbarplot handles 0/1 dichotomous data", {
  skip_if_not_installed("ggplot2")
  set.seed(2)
  df <- as.data.frame(matrix(sample(0:1, 200 * 5, replace = TRUE), 200, 5))
  colnames(df) <- paste0("I", 1:5)
  expect_s3_class(RMbarplot(df), "ggplot")
})

test_that("RMbarplot tolerates NA values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1:5, 1] <- NA
  expect_s3_class(RMbarplot(df), "ggplot")
})

# ---------------------------------------------------------------------
# RMstackedbarplot --- input validation
# ---------------------------------------------------------------------
test_that("RMstackedbarplot errors on non-data.frame input", {
  skip_if_not_installed("ggplot2")
  expect_error(RMstackedbarplot(matrix(0:3, 4, 4)), regexp = "data.frame")
})

test_that("RMstackedbarplot errors when fewer than 2 items", {
  skip_if_not_installed("ggplot2")
  expect_error(RMstackedbarplot(make_polytomous(k = 1L)),
               regexp = "at least 2 columns")
})

test_that("RMstackedbarplot errors on non-numeric columns", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df$group <- letters[seq_len(nrow(df))]
  expect_error(RMstackedbarplot(df), regexp = "must be numeric")
})

test_that("RMstackedbarplot errors on non-integer values", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1, 1] <- 0.5
  expect_error(RMstackedbarplot(df), regexp = "integers")
})

test_that("RMstackedbarplot errors when item_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMstackedbarplot(make_polytomous(),
                                item_labels = c("a", "b")),
               regexp = "same length")
})

test_that("RMstackedbarplot errors when category_labels length mismatches", {
  skip_if_not_installed("ggplot2")
  expect_error(RMstackedbarplot(make_polytomous(),
                                category_labels = c("Low", "High")),
               regexp = "same length")
})

# ---------------------------------------------------------------------
# RMstackedbarplot --- output
# ---------------------------------------------------------------------
test_that("RMstackedbarplot returns a ggplot", {
  skip_if_not_installed("ggplot2")
  expect_s3_class(RMstackedbarplot(make_polytomous()), "ggplot")
})

test_that("RMstackedbarplot accepts custom item & category labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMstackedbarplot(
    df,
    item_labels     = paste("Item", seq_len(ncol(df))),
    category_labels = c("None", "Some", "More", "Most")
  )
  expect_s3_class(p, "ggplot")
})

test_that("RMstackedbarplot supports show_percent and min_label_n", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMstackedbarplot(df, show_percent = TRUE, show_n = FALSE,
                        min_label_n = 5L)
  expect_s3_class(p, "ggplot")
})

test_that("RMstackedbarplot can suppress all labels", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  p <- RMstackedbarplot(df, show_n = FALSE, show_percent = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("RMstackedbarplot handles 0/1 dichotomous data", {
  skip_if_not_installed("ggplot2")
  set.seed(3)
  df <- as.data.frame(matrix(sample(0:1, 200 * 5, replace = TRUE), 200, 5))
  colnames(df) <- paste0("I", 1:5)
  expect_s3_class(RMstackedbarplot(df), "ggplot")
})
