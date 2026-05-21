# Tests for RMinfitcutoff_mi() and RMiteminfit_mi()

# Build a small `mids` (multiply-imputed dataset) object once for the file.
make_mids <- function(n = 200, k = 5, prop_na = 0.05, seed = 1L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, n * k, replace = TRUE), n, k)
  mat[sample(length(mat), round(prop_na * length(mat)))] <- NA
  df <- as.data.frame(mat)
  colnames(df) <- paste0("I", seq_len(k))
  suppressWarnings(mice::mice(df, m = 3L, method = "pmm",
                               seed = seed, printFlag = FALSE))
}

# ---------------------------------------------------------------------
# RMinfitcutoff_mi
# ---------------------------------------------------------------------
test_that("RMinfitcutoff_mi returns the expected list structure", {
  skip_if_not_installed("mice")
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  res <- RMinfitcutoff_mi(imp, iterations = 5L,
                          parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(c("results", "item_cutoffs", "actual_iterations",
                    "sample_n", "item_names", "cutoff_method",
                    "hdci_width", "n_imputations") %in% names(res)))
  expect_equal(res$n_imputations, 3L)
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(c("Item", "infit_low", "infit_high") %in%
                  names(res$item_cutoffs)))
})

test_that("RMinfitcutoff_mi is reproducible with the same seed", {
  skip_if_not_installed("mice")
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  r1 <- RMinfitcutoff_mi(imp, iterations = 5L,
                         parallel = FALSE, seed = 42L)
  r2 <- RMinfitcutoff_mi(imp, iterations = 5L,
                         parallel = FALSE, seed = 42L)
  expect_equal(r1$item_cutoffs, r2$item_cutoffs)
})

# ---------------------------------------------------------------------
# RMiteminfit_mi
# ---------------------------------------------------------------------
test_that("RMiteminfit_mi output = 'dataframe' returns pooled infit per item", {
  skip_if_not_installed("mice")
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  res <- RMiteminfit_mi(imp, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_true(all(c("Item", "Infit_MSQ") %in% names(res)))
})

test_that("RMiteminfit_mi accepts an RMinfitcutoff_mi result and adds Flagged column", {
  skip_if_not_installed("mice")
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp  <- make_mids()
  cuts <- RMinfitcutoff_mi(imp, iterations = 5L,
                           parallel = FALSE, seed = 1L)
  res  <- RMiteminfit_mi(imp, cutoff = cuts, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_true("Flagged" %in% names(res))
  expect_type(res$Flagged, "logical")
})

test_that("RMiteminfit_mi output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("mice")
  skip_if_not_installed("eRm")
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")

  imp <- make_mids()
  out <- RMiteminfit_mi(imp, output = "kable")
  expect_s3_class(out, "knitr_kable")
})
