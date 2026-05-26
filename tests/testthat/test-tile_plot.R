# Tests for RMplotTile()

make_dichotomous <- function(n = 80, k = 4, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:1, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

make_polytomous <- function(n = 80, k = 4, seed = 1L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:3, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df
}

# ---------------------------------------------------------------------
# Input validation
# ---------------------------------------------------------------------
test_that("RMplotTile errors when data is not a data.frame", {
  m <- matrix(0:1, nrow = 4, ncol = 4)
  expect_error(RMplotTile(m), regexp = "data\\.frame")
})

test_that("RMplotTile errors when fewer than 2 columns", {
  df <- data.frame(I1 = sample(0:1, 20, replace = TRUE))
  expect_error(RMplotTile(df), regexp = "at least 2 columns")
})

test_that("RMplotTile errors when all values are NA", {
  df <- data.frame(I1 = NA_integer_, I2 = NA_integer_)
  expect_error(RMplotTile(df), regexp = "only NA")
})

test_that("RMplotTile errors when non-numeric column present", {
  df <- data.frame(I1 = sample(0:1, 20, replace = TRUE),
                   I2 = letters[1:20])
  expect_error(RMplotTile(df), regexp = "numeric")
})

test_that("RMplotTile errors when group length mismatches nrow(data)", {
  df  <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMplotTile(df, group = grp),
               regexp = "same length as nrow")
})

test_that("RMplotTile errors when item_labels length mismatches ncol(data)", {
  df <- make_dichotomous()
  expect_error(RMplotTile(df, item_labels = c("a", "b")),
               regexp = "same length as the number of")
})

# ---------------------------------------------------------------------
# Output: dataframe
# ---------------------------------------------------------------------
test_that("RMplotTile output = 'dataframe' returns one row per item x category", {
  df <- make_polytomous()  # values 0..3 -> 4 categories
  out <- RMplotTile(df, output = "dataframe")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(df) * 4L)
  expect_true(all(c("item", "category", "n", "percentage") %in% names(out)))
  expect_true(all(out$n >= 0))
  expect_equal(sum(out$n), nrow(df) * ncol(df))
})

test_that("RMplotTile dataframe has group column when group is supplied", {
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  out <- RMplotTile(df, group = grp, output = "dataframe")
  expect_true("group" %in% names(out))
  # 4 items × 4 categories × 2 groups = 32 rows
  expect_equal(nrow(out), 4L * 4L * 2L)
})

# ---------------------------------------------------------------------
# Output: ggplot
# ---------------------------------------------------------------------
test_that("RMplotTile output = 'ggplot' returns a ggplot", {
  skip_if_not_installed("ggplot2")
  df <- make_dichotomous()
  p  <- RMplotTile(df)  # default output
  expect_s3_class(p, "ggplot")
})

test_that("RMplotTile with group + percent + custom cutoff returns ggplot", {
  skip_if_not_installed("ggplot2")
  df  <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  p   <- RMplotTile(df, group = grp, percent = TRUE, cutoff = 5)
  expect_s3_class(p, "ggplot")
})
