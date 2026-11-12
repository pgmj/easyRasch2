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

test_that(".bootstrap_pvalues joint moments make the step-down exact", {
  # Rotate each of the B + 1 rows into the observed slot. With joint moments
  # every rotation scores all rows on the same scale, so the smallest
  # adjusted p-values are the ranks 1..(B + 1) over B + 1. Moments from the
  # simulated values alone break this.
  set.seed(3)
  B <- 59; k <- 6
  x <- matrix(rnorm((B + 1) * k), B + 1, k,
              dimnames = list(NULL, paste0("I", seq_len(k))))
  min_padj <- function(tail, moments) {
    vapply(seq_len(B + 1), function(r) {
      min(easyRasch2:::.bootstrap_pvalues(
        x[r, ], x[-r, , drop = FALSE], "fwer", tail = tail, moments = moments
      )$padj)
    }, numeric(1))
  }
  ranks <- seq_len(B + 1) / (B + 1)
  for (tail in c("two.sided", "upper")) {
    expect_equal(sort(min_padj(tail, "joint")), ranks)
    expect_false(isTRUE(all.equal(sort(min_padj(tail, "simulated")), ranks)))
  }
})

test_that(".bootstrap_pvalues upper-tail marginal p is the raw rank under both moments", {
  set.seed(4)
  sim <- matrix(rnorm(300 * 4, mean = 1:4), ncol = 4, byrow = TRUE,
                dimnames = list(NULL, paste0("I", 1:4)))
  observed <- setNames(c(1.5, 2, 4.2, 3), colnames(sim))
  raw <- vapply(1:4, function(i) (1 + sum(sim[, i] >= observed[i])) / 301,
                numeric(1))
  j <- easyRasch2:::.bootstrap_pvalues(observed, sim, "none", tail = "upper")
  s <- easyRasch2:::.bootstrap_pvalues(observed, sim, "none", tail = "upper",
                                       moments = "simulated")
  expect_equal(j$p, raw)
  expect_equal(s$p, raw)
})

test_that(".bootstrap_pvalues leaves missing simulated values out (UM1)", {
  set.seed(2)
  sim <- matrix(rnorm(400), 200, 2, dimnames = list(NULL, c("a", "b")))
  obs <- c(a = 1.6, b = 0)
  na_sim <- sim
  na_sim[1:100, "a"] <- NA
  for (tail in c("two.sided", "upper")) {
    p_na <- easyRasch2:::.bootstrap_pvalues(obs, na_sim, "none", tail = tail,
                                            moments = "simulated")
    p_drop <- easyRasch2:::.bootstrap_pvalues(obs["a"], sim[101:200, "a", drop = FALSE],
                                              "none", tail = tail, moments = "simulated")
    # same as dropping the missing rows, not counted as non-exceedances
    expect_equal(p_na$p[1], p_drop$p)
    # the complete column keeps its own denominator
    expect_equal(p_na$p[2], easyRasch2:::.bootstrap_pvalues(
      obs["b"], sim[, "b", drop = FALSE], "none", tail = tail, moments = "simulated")$p)
  }
})

test_that("the step-down uses complete iterations only, without warnings (UM2)", {
  set.seed(5)
  sim <- matrix(rnorm(300 * 3), 300, 3, dimnames = list(NULL, c("a", "b", "c")))
  obs <- c(a = 2.5, b = 0.3, c = -1)
  part <- sim
  part[sample(300, 40), "b"] <- NA
  part[7, ] <- NA # an all-missing iteration
  cc <- stats::complete.cases(part)
  # on fixed statistics, an incomplete iteration is simply left out
  expect_no_warning(wy <- easyRasch2:::.wy_stepdown(abs(obs), abs(part)))
  expect_equal(wy, easyRasch2:::.wy_stepdown(abs(obs), abs(part[cc, ])))
  expect_no_warning(easyRasch2:::.bootstrap_pvalues(obs, part, "fwer"))
  # no complete iteration at all: no step-down, no warning
  alt <- sim[, 1:2]
  alt[c(TRUE, FALSE), "a"] <- NA
  alt[c(FALSE, TRUE), "b"] <- NA
  expect_no_warning(r2 <- easyRasch2:::.bootstrap_pvalues(obs[1:2], alt, "fwer"))
  expect_true(all(is.na(r2$padj)))
  expect_true(all(is.finite(r2$p)))
})

test_that("a constant simulated column gives NA under joint moments too (UM3)", {
  set.seed(6)
  sim <- cbind(a = rep(1, 200), b = rnorm(200))
  obs <- c(a = 1.01, b = 0)
  for (mo in c("joint", "simulated")) {
    res <- easyRasch2:::.bootstrap_pvalues(obs, sim, "fwer", moments = mo)
    expect_true(is.na(res$p[1]))
    expect_true(is.na(res$padj[1]))
    expect_false(is.na(res$p[2]))
  }
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

test_that("explicit p_value = FALSE gives the pre-1.2.0 columns", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  # Since 1.2.0 the default (p_value = NULL) would add the p-value columns
  # here, because `sim` carries the simulations. Opting out restores the
  # interval-based output.
  res <- suppressMessages(
    RMitemInfit(df, cutoff = sim, p_value = FALSE, output = "dataframe")
  )
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

test_that("low iteration count triggers the calibration message", {
  skip_on_cran()
  skip_if_not_installed("eRm"); skip_if_not_installed("iarm")
  df  <- sim_rasch_null()
  sim <- RMitemInfitCutoff(df, iterations = 200, parallel = FALSE, seed = 1)
  # Since 1.2.0 this is a once-per-session message rather than a warning, and
  # the threshold is the recommended 400 rather than 1000.
  for (cr in c("fwer", "fdr_bh", "fdr_by", "none")) rlang::reset_message_verbosity(paste0("easyRasch2_low_iterations_", cr))
  expect_message(RMitemInfit(df, cutoff = sim, p_value = TRUE,
                             output = "dataframe"),
                 regexp = "below the recommended 400")
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

test_that("small numbers in captions are plain decimals, never 7e-04", {
  expect_identical(.fmt_small(1 / 1501), "0.00067")
  expect_identical(.fmt_small(1 / 401), "0.0025")
  expect_identical(.fmt_small(1 / 41), "0.024")
  expect_identical(.fmt_small(1 / 100001), "0.00001")
  expect_identical(.fmt_small(1 / 5001, 3L), "0.0002")
  w <- tryCatch(.warn_fdr_floor(5000, 2000, "fdr_by"), warning = conditionMessage)
  expect_no_match(w, "[0-9]e-[0-9]")
})

test_that("the RMitemInfit caption states the p-value floor as a plain decimal", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  skip_if_not_installed("knitr")
  set.seed(8)
  thr <- lapply(seq(-1, 1, length.out = 5), function(b) b + c(-0.6, 0.6))
  d <- as.data.frame(sim_partial_score(thr, stats::rnorm(150)))
  names(d) <- paste0("I", 1:5)
  cu <- suppressMessages(RMitemInfitCutoff(d, iterations = 30, parallel = FALSE, seed = 1))
  # only the caption reads the iteration count, so pretend it was 1500
  cu$actual_iterations <- 1500L
  cu$requested_iterations <- 1500L
  cap <- paste(as.character(suppressMessages(RMitemInfit(d, cutoff = cu))), collapse = " ")
  expect_match(cap, "1/(1500+1) = 0.00067", fixed = TRUE)
  expect_no_match(cap, "[0-9]e-[0-9]")
})

test_that(".bootstrap_pvalues stops on unnamed or mismatched observed names (UM10)", {
  sim <- matrix(stats::rnorm(40), 20, 2, dimnames = list(NULL, c("a", "b")))
  msg <- "needs unique names that are all column names of `sim`"
  expect_error(.bootstrap_pvalues(c(1, 2), sim), msg, fixed = TRUE)
  expect_error(.bootstrap_pvalues(c(a = 1, a = 2), sim), msg, fixed = TRUE)
  expect_error(.bootstrap_pvalues(c(a = 1, z = 2), sim), msg, fixed = TRUE)
  expect_error(.bootstrap_pvalues(c(a = 1), unname(sim)), msg, fixed = TRUE)
  # a subset of the columns, in any order, is fine
  expect_identical(.bootstrap_pvalues(c(b = 1, a = 0), sim)$name, c("b", "a"))
})

test_that(".wy_stepdown equals the direct step-down definition (UM6)", {
  # Reference: at each step a fresh row-wise max over the remaining columns.
  ref_wy <- function(obs, sim) {
    B <- nrow(sim)
    ord <- order(obs, decreasing = TRUE)
    padj <- numeric(length(obs))
    prev <- 0
    for (j in seq_along(ord)) {
      mx <- apply(sim[, ord[j:length(ord)], drop = FALSE], 1L, max)
      prev <- max(prev, (1 + sum(mx >= obs[ord[j]])) / (B + 1))
      padj[ord[j]] <- min(prev, 1)
    }
    padj
  }
  set.seed(9)
  for (r in 1:30) {
    k <- sample(c(1, 2, 6, 25), 1)
    B <- sample(c(19, 99, 250), 1)
    sim <- matrix(round(abs(stats::rnorm(B * k)), sample(c(1, 6), 1)), B, k)
    obs <- round(abs(stats::rnorm(k, 0, 2)), 1)
    if (k > 2) obs[2] <- obs[1] # tied observed values
    expect_identical(as.vector(.wy_stepdown(obs, sim)), ref_wy(obs, sim))
  }
})

test_that("Westfall-Young padj is never below the marginal p (UM12)", {
  set.seed(1)
  for (r in 1:200) {
    K <- sample(2:8, 1)
    sim <- matrix(stats::rnorm(200 * K), 200, K,
                  dimnames = list(NULL, letters[1:K]))
    sim[matrix(stats::runif(200 * K) < stats::runif(1, 0, 0.2), 200, K)] <- NA
    obs <- stats::setNames(stats::rnorm(K, 0, 2), letters[1:K])
    res <- suppressMessages(.bootstrap_pvalues(obs, sim, correction = "fwer"))
    expect_true(all(res$padj >= res$p, na.rm = TRUE))
    expect_null(attributes(res$padj))
  }
})

test_that("A step-down that loses iterations says so (UM11)", {
  set.seed(2)
  sim <- matrix(stats::rnorm(300), 100, 3, dimnames = list(NULL, c("a", "b", "c")))
  obs <- c(a = 3, b = 0, c = 0)
  expect_no_message(.bootstrap_pvalues(obs, sim, correction = "fwer"))
  sim[1:4, "b"] <- NA # 96 of 100 complete: below the 5% threshold
  expect_no_message(.bootstrap_pvalues(obs, sim, correction = "fwer"))
  sim[1:20, "b"] <- NA
  expect_message(
    res <- .bootstrap_pvalues(obs, sim, correction = "fwer"),
    "used the 80 of 100 simulated datasets in which every comparison could be computed. Corrected p-values cannot be smaller than 1/(80+1) = 0.012.",
    fixed = TRUE
  )
  expect_gte(res$padj[1], 1 / 81)
  # The marginal corrections do not use the step-down
  expect_no_message(.bootstrap_pvalues(obs, sim, correction = "fdr_bh"))
})

test_that("The low-iterations notice is shown once per correction (UM16)", {
  id <- "easyRasch2_test_um16"
  for (cr in c("fwer", "none")) rlang::reset_message_verbosity(paste0(id, "_", cr))
  expect_message(.notify_low_iterations(150, NULL, id = id, correction = "fwer"))
  expect_no_message(.notify_low_iterations(150, NULL, id = id, correction = "fwer"))
  expect_message(.notify_low_iterations(150, NULL, id = id, correction = "none"))
})
