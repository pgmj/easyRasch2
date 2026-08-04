# Tests for the bootstrap p-value layer of RMitemInfit() and the shared
# multiplicity helper .bootstrap_pvalues().

sim_rasch_null <- function(n = 300, J = 8, seed = 1L) {
  set.seed(seed)
  theta <- rnorm(n, 0, 1.2); beta <- seq(-1.5, 1.5, length.out = J)
  m <- sapply(seq_len(J), function(j) rbinom(n, 1, plogis(theta - beta[j])))
  df <- as.data.frame(m); colnames(df) <- paste0("I", seq_len(J))
  df
}

# ---------------------------------------------------------------------
# .bootstrap_pvalues() unit tests (no model fitting needed)
# ---------------------------------------------------------------------

test_that(".bootstrap_pvalues returns valid p-values and FWER >= marginal", {
  set.seed(1)
  k <- 6; B <- 2000
  sim <- matrix(rnorm(B * k, mean = 1, sd = 0.1), nrow = B)
  colnames(sim) <- paste0("I", seq_len(k))
  observed <- setNames(c(1.0, 1.0, 1.0, 1.0, 1.0, 1.6), colnames(sim))  # I6 extreme
  res <- easyRasch2:::.bootstrap_pvalues(observed, sim, correction = "fwer")
  expect_named(res, c("name", "p", "padj"))
  expect_true(all(res$p >= 0 & res$p <= 1))
  expect_true(all(res$padj >= res$p - 1e-9))          # FWER no smaller than marginal
  expect_lt(res$padj[res$name == "I6"], 0.05)         # the extreme item is detected
  expect_true(all(res$padj[res$name != "I6"] > 0.05)) # the null items are not
})

test_that(".bootstrap_pvalues FDR detects strong signals and is valid", {
  set.seed(2)
  k <- 8; B <- 2000
  sim <- matrix(rnorm(B * k, 1, 0.1), nrow = B); colnames(sim) <- paste0("I", 1:k)
  observed <- setNames(c(1.5, 1.45, rep(1.0, k - 2)), colnames(sim))
  bh <- easyRasch2:::.bootstrap_pvalues(observed, sim, "fdr_bh")
  expect_true(all(bh$padj >= 0 & bh$padj <= 1))
  expect_true(all(bh$padj[1:2] < 0.05))               # the two real effects detected
  expect_true(all(bh$padj[3:k] > 0.05))               # null comparisons not flagged
})

test_that(".bootstrap_pvalues handles zero-variance comparisons as NA", {
  sim <- cbind(I1 = rnorm(500, 1, 0.1), I2 = rep(1, 500))
  observed <- setNames(c(1.3, 1.0), c("I1", "I2"))
  res <- easyRasch2:::.bootstrap_pvalues(observed, sim, "fwer")
  expect_true(is.na(res$p[res$name == "I2"]))
  expect_false(is.na(res$p[res$name == "I1"]))
})

# ---------------------------------------------------------------------
# RMitemInfit() integration
# ---------------------------------------------------------------------

test_that("p_value = TRUE adds p columns and bases Flagged on padj", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  res <- suppressWarnings(
    RMitemInfit(df, cutoff = sim, p_value = TRUE, output = "dataframe")
  )
  expect_true(all(c("p_infit", "padj_infit", "Flagged") %in% names(res)))
  expect_true(all(res$p_infit >= 0 & res$p_infit <= 1, na.rm = TRUE))
  expect_true(all(res$padj_infit >= res$p_infit - 1e-9, na.rm = TRUE))
})

test_that("p_value = TRUE errors without the full cutoff object", {
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  expect_error(RMitemInfit(df, p_value = TRUE), regexp = "full RMitemInfitCutoff")
  expect_error(RMitemInfit(df, cutoff = sim$item_cutoffs, p_value = TRUE),
               regexp = "full RMitemInfitCutoff")
})

test_that("p_value = FALSE leaves the output unchanged (backward compatible)", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  res <- RMitemInfit(df, cutoff = sim, output = "dataframe")
  expect_named(res, c("Item", "Infit_MSQ", "Infit_low", "Infit_high",
                      "Flagged", "Relative_location"))
})

test_that("all correction methods run and return the expected columns", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  for (corr in c("fwer", "fdr_bh", "fdr_by", "none")) {
    res <- suppressWarnings(
      RMitemInfit(df, cutoff = sim, p_value = TRUE, correction = corr,
                  output = "dataframe")
    )
    expect_true(all(res$padj_infit >= 0 & res$padj_infit <= 1, na.rm = TRUE))
  }
})

test_that("low iteration count triggers a warning", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  expect_warning(RMitemInfit(df, cutoff = sim, p_value = TRUE,
                             output = "dataframe"),
                 regexp = "iterations")
})

test_that("invalid alpha is rejected", {
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df <- sim_rasch_null()
  expect_error(RMitemInfit(df, alpha = 1.5), regexp = "alpha")
})

test_that("kable output works with p-values", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  expect_s3_class(
    suppressWarnings(RMitemInfit(df, cutoff = sim, p_value = TRUE)),
    "knitr_kable"
  )
})

# ---------------------------------------------------------------------
# statistic = "outfit"
# ---------------------------------------------------------------------

test_that("statistic = 'outfit' names columns for outfit and uses OutfitMSQ", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 300, parallel = FALSE, seed = 1)
  res <- suppressWarnings(
    RMitemInfit(df, cutoff = sim, p_value = TRUE, statistic = "outfit",
                output = "dataframe")
  )
  expect_named(res, c("Item", "Outfit_MSQ", "Outfit_low", "Outfit_high",
                      "p_outfit", "padj_outfit", "Flagged",
                      "Relative_location"))
  expect_true(all(res$p_outfit >= 0 & res$p_outfit <= 1, na.rm = TRUE))
  expect_true(all(res$padj_outfit >= res$p_outfit - 1e-9, na.rm = TRUE))
  # the reported statistic is iarm's outfit, not infit
  cfit <- iarm::out_infit(psychotools::pcmodel(df))
  expect_equal(res$Outfit_MSQ, as.numeric(cfit$Outfit))
  # and the p-values come from the outfit column of the simulated null
  sim_mat <- tapply(sim$results$OutfitMSQ,
                    list(sim$results$iteration, sim$results$Item),
                    function(x) x[1L])
  expected <- easyRasch2:::.bootstrap_pvalues(
    setNames(as.numeric(cfit$Outfit), names(df)), sim_mat, "fwer"
  )
  expect_equal(res$p_outfit, expected$p[match(res$Item, expected$name)])
})

test_that("statistic = 'outfit' works without p-values and against a bare cutoff df", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  res <- RMitemInfit(df, cutoff = sim$item_cutoffs, statistic = "outfit",
                     output = "dataframe")
  expect_named(res, c("Item", "Outfit_MSQ", "Outfit_low", "Outfit_high",
                      "Flagged", "Relative_location"))
  expect_equal(res$Outfit_low, sim$item_cutoffs$outfit_low)
  # a cutoff df carrying only infit bounds is rejected
  expect_error(
    RMitemInfit(df, cutoff = sim$item_cutoffs[, c("Item", "infit_low", "infit_high")],
                statistic = "outfit"),
    regexp = "outfit_low, outfit_high"
  )
})

test_that("statistic = 'infit' is the default and leaves output unchanged", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  expect_identical(
    RMitemInfit(df, cutoff = sim, output = "dataframe"),
    RMitemInfit(df, cutoff = sim, statistic = "infit", output = "dataframe")
  )
  expect_error(RMitemInfit(df, statistic = "both"), regexp = "should be one of")
})
