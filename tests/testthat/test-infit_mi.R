# Tests for RMitemInfitCutoffMI() and RMitemInfitMI()

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
# RMitemInfitCutoffMI
# ---------------------------------------------------------------------
test_that("RMitemInfitCutoffMI returns the expected list structure", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  res <- RMitemInfitCutoffMI(imp, iterations = 5L,
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

test_that("RMitemInfitCutoffMI is reproducible with the same seed", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  r1 <- RMitemInfitCutoffMI(imp, iterations = 5L,
                         parallel = FALSE, seed = 42L)
  r2 <- RMitemInfitCutoffMI(imp, iterations = 5L,
                         parallel = FALSE, seed = 42L)
  expect_equal(r1$item_cutoffs, r2$item_cutoffs)
})

test_that("RMitemInfitCutoffMI cutoff_method = 'quantile' aggregates without ggdist", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  res <- RMitemInfitCutoffMI(imp, iterations = 6L, parallel = FALSE, seed = 1L,
                             cutoff_method = "quantile")
  expect_equal(res$cutoff_method, "quantile")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(res$item_cutoffs$infit_low < res$item_cutoffs$infit_high))
})

test_that("RMitemInfitCutoffMI verbose = TRUE runs the per-imputation progress path", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  suppressMessages(invisible(utils::capture.output(
    res <- RMitemInfitCutoffMI(imp, iterations = 6L, parallel = FALSE, seed = 1L,
                               cutoff_method = "quantile", verbose = TRUE)
  )))
  expect_equal(res$n_imputations, 3L)
})

test_that("RMitemInfitCutoffMI errors on a non-mids object", {
  skip_if_not_installed("mice")
  expect_snapshot(
    RMitemInfitCutoffMI(data.frame(a = 0:1, b = 1:0)),
    error = TRUE
  )
})

# ---------------------------------------------------------------------
# RMitemInfitMI
# ---------------------------------------------------------------------
test_that("RMitemInfitMI output = 'dataframe' returns pooled infit per item", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  res <- RMitemInfitMI(imp, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_equal(nrow(res), 5L)
  expect_true(all(c("Item", "Infit_MSQ") %in% names(res)))
})

test_that("RMitemInfitMI accepts an RMitemInfitCutoffMI result and adds Flagged column", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp  <- make_mids()
  cuts <- RMitemInfitCutoffMI(imp, iterations = 5L,
                           parallel = FALSE, seed = 1L)
  res  <- RMitemInfitMI(imp, cutoff = cuts, output = "dataframe")
  expect_s3_class(res, "data.frame")
  expect_true("Flagged" %in% names(res))
  expect_type(res$Flagged, "character")
  expect_true(all(res$Flagged %in% c("overfit", "underfit", "")))
})

test_that("RMitemInfitMI output = 'kable' returns a knitr_kable", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")

  imp <- make_mids()
  out <- RMitemInfitMI(imp, output = "kable")
  expect_s3_class(out, "knitr_kable")
})
