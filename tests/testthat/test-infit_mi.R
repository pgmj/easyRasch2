# Tests for RMitemInfitCutoffMI() and RMitemInfitMI()

# Build a small `mids` (multiply-imputed dataset) object once for the file.
make_mids <- function(n = 200, k = 5, prop_na = 0.05, seed = 1L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, n * k, replace = TRUE), n, k)
  mat[sample(length(mat), round(prop_na * length(mat)))] <- NA
  df <- as.data.frame(mat)
  colnames(df) <- paste0("I", seq_len(k))
  suppressWarnings(mice::mice(
    df,
    m = 3L,
    method = "pmm",
    seed = seed,
    printFlag = FALSE
  ))
}

# A `mids` built from ordered factors + mice's `polr` method, as recommended
# for ordinal items. mice::complete() then returns ordered factors, which the
# MI functions must coerce back to numeric internally.
make_mids_polr <- function(n = 200, k = 5, prop_na = 0.05, seed = 1L) {
  set.seed(seed)
  mat <- matrix(sample(0:1, n * k, replace = TRUE), n, k)
  mat[sample(length(mat), round(prop_na * length(mat)))] <- NA
  df <- as.data.frame(mat)
  colnames(df) <- paste0("I", seq_len(k))
  df[] <- lapply(df, function(x) factor(x, ordered = TRUE))
  suppressWarnings(mice::mice(
    df,
    m = 3L,
    method = "polr",
    seed = seed,
    printFlag = FALSE
  ))
}

# ---------------------------------------------------------------------
# RMitemInfitCutoffMI
# ---------------------------------------------------------------------
test_that("RMitemInfitCutoffMI accepts a polr/ordered-factor mids (coerced to numeric)", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("MASS")
  imp <- make_mids_polr()
  res <- RMitemInfitCutoffMI(
    imp,
    iterations = 5L,
    parallel = FALSE,
    seed = 1L,
    cutoff_method = "quantile"
  )
  expect_equal(res$n_imputations, 3L)
  expect_equal(nrow(res$item_cutoffs), 5L)
})

test_that("RMitemInfitCutoffMI returns the expected list structure", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  res <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 1L)
  expect_type(res, "list")
  expect_true(all(
    c(
      "results",
      "item_cutoffs",
      "actual_iterations",
      "sample_n",
      "item_names",
      "cutoff_method",
      "hdci_width",
      "n_imputations"
    ) %in%
      names(res)
  ))
  expect_equal(res$n_imputations, 3L)
  expect_s3_class(res$item_cutoffs, "data.frame")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(
    c("Item", "infit_low", "infit_high") %in%
      names(res$item_cutoffs)
  ))
})

test_that("RMitemInfitCutoffMI is reproducible with the same seed", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")

  imp <- make_mids()
  r1 <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 42L)
  r2 <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 42L)
  expect_equal(r1$item_cutoffs, r2$item_cutoffs)
})

test_that("RMitemInfitCutoffMI cutoff_method = 'quantile' aggregates without ggdist", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  res <- RMitemInfitCutoffMI(
    imp,
    iterations = 6L,
    parallel = FALSE,
    seed = 1L,
    cutoff_method = "quantile"
  )
  expect_equal(res$cutoff_method, "quantile")
  expect_equal(nrow(res$item_cutoffs), 5L)
  expect_true(all(res$item_cutoffs$infit_low < res$item_cutoffs$infit_high))
})

test_that("RMitemInfitCutoffMI verbose = TRUE runs the per-imputation progress path", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")

  imp <- make_mids()
  suppressMessages(invisible(utils::capture.output(
    res <- RMitemInfitCutoffMI(
      imp,
      iterations = 6L,
      parallel = FALSE,
      seed = 1L,
      cutoff_method = "quantile",
      verbose = TRUE
    )
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

  imp <- make_mids()
  cuts <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 1L)
  res <- RMitemInfitMI(imp, cutoff = cuts, output = "dataframe")
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

test_that("RMitemInfitMI errors on a non-mids object", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  expect_snapshot(
    RMitemInfitMI(data.frame(a = 0:1, b = 1:0)),
    error = TRUE
  )
})

test_that("RMitemInfitMI rejects malformed cutoff arguments", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  imp <- make_mids()
  expect_snapshot(
    RMitemInfitMI(imp, cutoff = "abc"),
    error = TRUE
  )
  expect_snapshot(
    RMitemInfitMI(imp, cutoff = data.frame(Item = "I1", infit_low = 0.8)),
    error = TRUE
  )
})

test_that("RMitemInfitMI errors when cutoff item names do not match", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  imp <- make_mids()
  cuts <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 1L)
  bad <- cuts$item_cutoffs
  bad$Item[1L] <- "NotAnItem"
  expect_error(RMitemInfitMI(imp, cutoff = bad), regexp = "do not match")
})

test_that("RMitemInfitMI accepts the bare $item_cutoffs data.frame", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  imp <- make_mids()
  cuts <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 1L)
  res <- RMitemInfitMI(imp, cutoff = cuts$item_cutoffs, output = "dataframe")
  expect_true("Flagged" %in% names(res))
  # Without the full object, the caption falls back to the generic wording
  kbl <- RMitemInfitMI(imp, cutoff = cuts$item_cutoffs)
  expect_match(
    paste(kbl, collapse = "\n"),
    "Simulation-based cutoff values applied"
  )
})

test_that("RMitemInfitMI sort = 'infit' sorts descending", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  imp <- make_mids()
  res <- RMitemInfitMI(imp, output = "dataframe", sort = "infit")
  expect_false(is.unsorted(rev(res$Infit_MSQ)))
})

test_that("RMitemInfitMI kable captions report pooling and cutoff metadata", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  skip_if_not_installed("ggdist")
  imp <- make_mids()

  plain <- paste(RMitemInfitMI(imp), collapse = "\n")
  expect_match(plain, "Pooled MSQ values from 3 imputations")
  expect_match(plain, "missing values imputed, m = 3")
  expect_match(plain, "per imputed dataset")

  cuts <- RMitemInfitCutoffMI(imp, iterations = 5L, parallel = FALSE, seed = 1L)
  with_cut <- paste(RMitemInfitMI(imp, cutoff = cuts), collapse = "\n")
  expect_match(with_cut, "across 3 imputations")
  expect_match(with_cut, "overfit = infit below range")
})

# ---------------------------------------------------------------------
# RMitemInfitMI -- imputation-failure paths
# ---------------------------------------------------------------------

# Tamper with a mids object so that the completed dataset for the given
# imputation(s) contains a negative response, which fails
# validate_response_data() inside the per-imputation loop.
break_imputations <- function(imp, which_m) {
  col_with_na <- names(Filter(function(x) nrow(x) > 0L, imp$imp))[1L]
  for (i in which_m) {
    imp$imp[[col_with_na]][1L, i] <- -1
  }
  imp
}

test_that("RMitemInfitMI warns about failed imputations and pools the rest", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  imp <- break_imputations(make_mids(), which_m = 1L)
  expect_warning(
    expect_warning(
      res <- RMitemInfitMI(imp, output = "dataframe"),
      regexp = "failed for imputation 1"
    ),
    regexp = "1 of 3 imputed datasets failed"
  )
  expect_equal(nrow(res), 5L)
  expect_true(all(is.finite(res$Infit_MSQ)))
})

test_that("RMitemInfitMI errors when fewer than 2 imputations succeed", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  imp <- break_imputations(make_mids(), which_m = 1:2)
  expect_error(
    suppressWarnings(RMitemInfitMI(imp, output = "dataframe")),
    regexp = "At least 2 are required"
  )
})

test_that("RMitemInfitMI errors when all imputations fail", {
  skip_if_not_installed("mice")
  skip_if_not_installed("iarm")
  imp <- break_imputations(make_mids(), which_m = 1:3)
  expect_error(
    suppressWarnings(RMitemInfitMI(imp, output = "dataframe")),
    regexp = "failed for all 3 imputed datasets"
  )
})
