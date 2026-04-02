test_that("RMiteminfit errors when iarm is not installed", {
  skip_if(
    requireNamespace("iarm", quietly = TRUE),
    "iarm is installed; skipping missing-package error path test"
  )
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMiteminfit(df), regexp = "iarm")
})

test_that("RMiteminfit errors when data does not start at 0", {
  df <- as.data.frame(matrix(sample(1:3, 200, replace = TRUE), nrow = 40, ncol = 5))
  expect_error(RMiteminfit(df), regexp = "scored starting at 0")
})

test_that("RMiteminfit errors when data is all NA", {
  df <- as.data.frame(matrix(NA_integer_, nrow = 10, ncol = 4))
  expect_error(RMiteminfit(df), regexp = "No complete|no non-missing")
})

test_that("RMiteminfit errors when data is not a data.frame or matrix", {
  expect_error(RMiteminfit(list(a = 1:5)), regexp = "data.frame or matrix")
})

test_that("RMiteminfit errors when no complete cases", {
  skip_if_not_installed("iarm")
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  # Make all rows have at least one NA by cycling through columns
  df[cbind(seq_len(nrow(df)), ((seq_len(nrow(df)) - 1L) %% ncol(df)) + 1L)] <- NA
  expect_error(RMiteminfit(df), regexp = "No complete cases")
})

test_that("RMiteminfit output = 'dataframe' returns correct structure (dichotomous)", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMiteminfit(df, output = "dataframe")

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5L)
  expect_equal(names(result), c("Item", "Infit_MSQ", "Relative_location"))
  expect_equal(result$Item, colnames(df))
  expect_true(all(result$Infit_MSQ > 0))
})

test_that("RMiteminfit output = 'kable' returns knitr_kable with complete-cases caption", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMiteminfit(df, output = "kable")
  expect_s3_class(result, "knitr_kable")
  cap <- attr(result, "caption")
  expect_true(grepl("complete cases", cap))
})

test_that("RMiteminfit sort = 'infit' sorts by Infit_MSQ descending", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result_sorted   <- RMiteminfit(df, output = "dataframe", sort = "infit")
  result_unsorted <- RMiteminfit(df, output = "dataframe")

  expect_true(
    all(diff(result_sorted$Infit_MSQ) <= 0),
    info = "Rows should be in descending order of Infit_MSQ"
  )
  expect_setequal(result_sorted$Item, result_unsorted$Item)
})

# ---------------------------------------------------------------------------
# cutoff parameter tests
# ---------------------------------------------------------------------------

test_that("RMiteminfit cutoff = NULL (default) behaves as before", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  result <- RMiteminfit(df, cutoff = NULL, output = "dataframe")
  expect_equal(names(result), c("Item", "Infit_MSQ", "Relative_location"))
})

test_that("RMiteminfit with mock cutoff data.frame adds expected columns", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  mock_cutoff <- data.frame(
    Item        = paste0("Item", 1:5),
    infit_low   = rep(0.7, 5),
    infit_high  = rep(1.3, 5),
    outfit_low  = rep(0.7, 5),
    outfit_high = rep(1.3, 5),
    stringsAsFactors = FALSE
  )

  result <- RMiteminfit(df, cutoff = mock_cutoff, output = "dataframe")
  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 5L)
  expect_equal(
    names(result),
    c("Item", "Infit_MSQ", "Infit_low", "Infit_high", "Flagged", "Relative_location")
  )
  expect_type(result$Flagged, "logical")
  expect_true(all(result$Infit_low  == 0.7))
  expect_true(all(result$Infit_high == 1.3))
})

test_that("RMiteminfit with cutoff list (from RMinfitcutoff shape) extracts item_cutoffs", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  mock_cutoff_list <- list(
    item_cutoffs = data.frame(
      Item        = paste0("Item", 1:5),
      infit_low   = rep(0.7, 5),
      infit_high  = rep(1.3, 5),
      outfit_low  = rep(0.7, 5),
      outfit_high = rep(1.3, 5),
      stringsAsFactors = FALSE
    ),
    actual_iterations = 100L
  )

  result_df <- RMiteminfit(df, cutoff = mock_cutoff_list, output = "dataframe")
  expect_equal(
    names(result_df),
    c("Item", "Infit_MSQ", "Infit_low", "Infit_high", "Flagged", "Relative_location")
  )

  result_kable <- RMiteminfit(df, cutoff = mock_cutoff_list, output = "kable")
  expect_s3_class(result_kable, "knitr_kable")
  cap <- attr(result_kable, "caption")
  expect_true(grepl("100 simulation iterations", cap))
})

test_that("RMiteminfit caption includes HDI info when cutoff_method = 'hdi'", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  mock_cutoff_hdi <- list(
    item_cutoffs = data.frame(
      Item        = paste0("Item", 1:5),
      infit_low   = rep(0.7, 5),
      infit_high  = rep(1.3, 5),
      outfit_low  = rep(0.7, 5),
      outfit_high = rep(1.3, 5),
      stringsAsFactors = FALSE
    ),
    actual_iterations = 200L,
    cutoff_method     = "hdi",
    hdi_width         = 0.999
  )

  result_kable <- RMiteminfit(df, cutoff = mock_cutoff_hdi, output = "kable")
  cap <- attr(result_kable, "caption")
  expect_true(grepl("200 simulation iterations", cap))
  expect_true(grepl("HDI", cap))
  expect_true(grepl("99.9", cap))
})

test_that("RMiteminfit caption includes quantile info when cutoff_method = 'quantile'", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  mock_cutoff_quantile <- list(
    item_cutoffs = data.frame(
      Item        = paste0("Item", 1:5),
      infit_low   = rep(0.7, 5),
      infit_high  = rep(1.3, 5),
      outfit_low  = rep(0.7, 5),
      outfit_high = rep(1.3, 5),
      stringsAsFactors = FALSE
    ),
    actual_iterations = 150L,
    cutoff_method     = "quantile",
    hdi_width         = 0.999
  )

  result_kable <- RMiteminfit(df, cutoff = mock_cutoff_quantile, output = "kable")
  cap <- attr(result_kable, "caption")
  expect_true(grepl("150 simulation iterations", cap))
  expect_true(grepl("percentile", cap))
})

test_that("RMinfitcutoff errors when cutoff_method = 'hdi' but ggdist is not installed", {
  skip_if(
    requireNamespace("ggdist", quietly = TRUE),
    "ggdist is installed; skipping missing-package error path test"
  )
  set.seed(42)
  df <- as.data.frame(matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5))
  colnames(df) <- paste0("Item", 1:5)
  expect_error(
    RMinfitcutoff(df, iterations = 5, parallel = FALSE, seed = 42,
                  cutoff_method = "hdi"),
    regexp = "ggdist"
  )
})

test_that("RMinfitcutoff cutoff_method = 'quantile' returns valid cutoffs", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  res <- RMinfitcutoff(df, iterations = 20, parallel = FALSE, seed = 42,
                       cutoff_method = "quantile")

  expect_equal(res$cutoff_method, "quantile")
  expect_true(is.numeric(res$hdi_width))
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(res$item_cutoffs$infit_low < res$item_cutoffs$infit_high))
  expect_true(all(res$item_cutoffs$outfit_low < res$item_cutoffs$outfit_high))
})

test_that("RMinfitcutoff cutoff_method = 'hdi' returns valid cutoffs", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  skip_if_not_installed("ggdist")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  res <- RMinfitcutoff(df, iterations = 20, parallel = FALSE, seed = 42,
                       cutoff_method = "hdi", hdi_width = 0.999)

  expect_equal(res$cutoff_method, "hdi")
  expect_equal(res$hdi_width, 0.999)
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(res$item_cutoffs$infit_low < res$item_cutoffs$infit_high))
  expect_true(all(res$item_cutoffs$outfit_low < res$item_cutoffs$outfit_high))
})

test_that("RMiteminfit with cutoff data.frame (no actual_iterations) uses generic caption", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  mock_cutoff <- data.frame(
    Item        = paste0("Item", 1:5),
    infit_low   = rep(0.7, 5),
    infit_high  = rep(1.3, 5),
    stringsAsFactors = FALSE
  )

  result_kable <- RMiteminfit(df, cutoff = mock_cutoff, output = "kable")
  expect_s3_class(result_kable, "knitr_kable")
  cap <- attr(result_kable, "caption")
  expect_true(grepl("Simulation-based cutoff values applied", cap))
})

test_that("RMiteminfit errors when cutoff item names do not match data", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  bad_cutoff <- data.frame(
    Item        = paste0("X", 1:5),
    infit_low   = rep(0.7, 5),
    infit_high  = rep(1.3, 5),
    stringsAsFactors = FALSE
  )

  expect_error(
    RMiteminfit(df, cutoff = bad_cutoff),
    regexp = "Item names in `cutoff` do not match"
  )
})

test_that("RMiteminfit errors when cutoff data.frame has missing required columns", {
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  bad_cutoff <- data.frame(
    Item      = paste0("Item", 1:5),
    infit_low = rep(0.7, 5),
    stringsAsFactors = FALSE
  )

  expect_error(
    RMiteminfit(df, cutoff = bad_cutoff),
    regexp = "missing required columns"
  )
})

test_that("RMiteminfit Flagged column correctly identifies values outside cutoff range", {
  skip_if_not_installed("iarm")
  skip_if_not_installed("eRm")
  set.seed(42)
  df <- as.data.frame(
    matrix(sample(0:1, 200, replace = TRUE), nrow = 40, ncol = 5)
  )
  colnames(df) <- paste0("Item", 1:5)

  # Compute actual infit values to build tight vs. loose cutoffs
  base_result <- RMiteminfit(df, output = "dataframe")
  infit_vals  <- base_result$Infit_MSQ

  # Tight cutoff: all items should be flagged (narrow window around each)
  tight_cutoff <- data.frame(
    Item        = paste0("Item", 1:5),
    infit_low   = infit_vals + 0.1,
    infit_high  = infit_vals + 0.2,
    stringsAsFactors = FALSE
  )
  tight_result <- RMiteminfit(df, cutoff = tight_cutoff, output = "dataframe")
  expect_true(all(tight_result$Flagged))

  # Loose cutoff: no items should be flagged (window encompasses all values)
  loose_cutoff <- data.frame(
    Item        = paste0("Item", 1:5),
    infit_low   = rep(0.0, 5),
    infit_high  = rep(10.0, 5),
    stringsAsFactors = FALSE
  )
  loose_result <- RMiteminfit(df, cutoff = loose_cutoff, output = "dataframe")
  expect_false(any(loose_result$Flagged))
})
