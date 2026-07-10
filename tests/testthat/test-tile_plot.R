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
  df <- data.frame(I1 = sample(0:1, 20, replace = TRUE), I2 = letters[1:20])
  expect_error(RMplotTile(df), regexp = "numeric")
})

test_that("RMplotTile errors when group length mismatches nrow(data)", {
  df <- make_dichotomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df) - 1L))
  expect_error(RMplotTile(df, group = grp), regexp = "same length as nrow")
})

test_that("RMplotTile errors when item_labels length mismatches ncol(data)", {
  df <- make_dichotomous()
  expect_error(
    RMplotTile(df, item_labels = c("a", "b")),
    regexp = "same length as the number of"
  )
})

# ---------------------------------------------------------------------
# Output: dataframe
# ---------------------------------------------------------------------
test_that("RMplotTile output = 'dataframe' returns one row per item x category", {
  df <- make_polytomous() # values 0..3 -> 4 categories
  out <- RMplotTile(df, output = "dataframe")
  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), ncol(df) * 4L)
  expect_true(all(c("item", "category", "n", "percentage") %in% names(out)))
  expect_true(all(out$n >= 0))
  expect_equal(sum(out$n), nrow(df) * ncol(df))
})

test_that("RMplotTile dataframe has group column when group is supplied", {
  df <- make_polytomous()
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
  p <- RMplotTile(df) # default output
  expect_s3_class(p, "ggplot")
})

test_that("RMplotTile with group + percent + custom cutoff returns ggplot", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  grp <- factor(rep(c("A", "B"), length.out = nrow(df)))
  p <- RMplotTile(df, group = grp, percent = TRUE, cutoff = 5)
  expect_s3_class(p, "ggplot")
})

# ---------------------------------------------------------------------
# Missingness reporting (message + caption)
# ---------------------------------------------------------------------
test_that("RMplotTile messages when group has NA and reports n in caption", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  grp[1:5] <- NA

  expect_message(
    p <- RMplotTile(df, group = grp),
    regexp = "5 row\\(s\\) with NA in `group` dropped"
  )
  expect_match(p$labels$caption, "n = 75 of 80 respondents")
  expect_match(p$labels$caption, "complete group data")
})

test_that("RMplotTile is silent with complete group and plain caption", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))

  expect_no_message(p <- RMplotTile(df, group = grp))
  expect_match(p$labels$caption, "n = 80 respondents")
  expect_no_match(p$labels$caption, "of 80")
})

test_that("RMplotTile caption notes retained incomplete responses", {
  skip_if_not_installed("ggplot2")
  df <- make_polytomous()
  df[1:3, 1] <- NA
  p <- RMplotTile(df)
  expect_match(p$labels$caption, "n = 80 respondents")
  expect_match(p$labels$caption, "incomplete responses retained")
})

test_that("RMplotTile group-NA message also fires for dataframe output", {
  df <- make_polytomous()
  grp <- rep(c("A", "B"), length.out = nrow(df))
  grp[1L] <- NA
  expect_message(
    RMplotTile(df, group = grp, output = "dataframe"),
    regexp = "1 row\\(s\\) with NA in `group` dropped"
  )
})

test_that("zero cells are taken off the fill scale and text_size applies", {
  skip_if_not_installed("ggplot2")
  set.seed(1)
  df <- as.data.frame(matrix(sample(0:2, 60 * 3, replace = TRUE), nrow = 60))
  colnames(df) <- paste0("i", 1:3)
  df$i1[df$i1 == 2] <- 1   # force an empty category on item 1

  p <- RMplotTile(df, text_size = 5)
  expect_s3_class(p, "ggplot")
  # zero cells have NA fill (mapped to zero_fill via na.value)
  expect_true(anyNA(p$data$fill_n))
  expect_true(all(is.na(p$data$fill_n[p$data$n == 0])))
  # caption mentions the convention
  expect_match(p$labels$caption, "no responses")

  # opting out restores the old fill mapping
  p0 <- RMplotTile(df, zero_fill = NULL)
  expect_false(anyNA(p0$data$fill_n))
})
