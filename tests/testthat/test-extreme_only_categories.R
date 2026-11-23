# A category used only by respondents with a zero or perfect total score is
# empty in the data CML is estimated from: its threshold diverges and a
# parametric bootstrap never generates it again.

make_extreme_top <- function(n = 200, k = 6, seed = 5L) {
  set.seed(seed)
  df <- as.data.frame(matrix(sample(0:2, n * k, replace = TRUE), n, k))
  colnames(df) <- paste0("I", seq_len(k))
  df$I6[df$I6 == 2] <- 1L
  df[1:2, ] <- 2L # the only respondents in I6's top category, both perfect
  df
}

test_that(".extreme_only_categories finds a category used only by perfect scorers", {
  ec <- easyRasch2:::.extreme_only_categories(make_extreme_top())
  expect_equal(ec, data.frame(Item = "I6", Category = 2L))
})

test_that(".extreme_only_categories returns no rows for ordinary data", {
  set.seed(6)
  df <- as.data.frame(matrix(sample(0:2, 300 * 5, replace = TRUE), 300, 5))
  expect_equal(nrow(easyRasch2:::.extreme_only_categories(df)), 0L)
})

test_that(".extreme_only_categories handles dichotomous items, zero scorers and NA", {
  set.seed(7)
  df <- as.data.frame(matrix(stats::rbinom(200 * 5, 1, 0.5), 200, 5))
  df[, 1] <- 1L
  df[1:2, ] <- 0L # item 1 is 0 only for zero scorers
  expect_equal(easyRasch2:::.extreme_only_categories(df),
               data.frame(Item = "V1", Category = 0L))
  # A respondent who answered two items, both at 0, is also a zero scorer.
  df[3, ] <- c(0L, 0L, NA, NA, NA)
  expect_equal(easyRasch2:::.extreme_only_categories(df),
               data.frame(Item = "V1", Category = 0L))
})

test_that("the message names each item and category", {
  msg <- easyRasch2:::.extreme_only_message(
    data.frame(Item = c("A", "A", "B"), Category = c(0L, 3L, 4L))
  )
  expect_match(msg, "A (categories 0, 3), B (category 4)", fixed = TRUE)
  expect_match(msg, "merging", fixed = TRUE)
})

test_that("RMitemParameters warns about a category used only by extreme scorers", {
  expect_warning(
    RMitemParameters(make_extreme_top(), output = "dataframe"),
    "zero or perfect total score: I6 \\(category 2\\)"
  )
})

test_that("simulation cutoffs stop before simulating on an extreme-only category", {
  # Before, RMitemInfitCutoff() simulated first; when a few refits survived it
  # returned cutoffs from those few without a word.
  d <- make_extreme_top()
  pat <- "zero or perfect total score: I6 \\(category 2\\).*Nothing was simulated"
  expect_error(RMitemInfitCutoff(d, iterations = 5, parallel = FALSE, seed = 1),
               pat, class = "easyRasch2_category_error")
  expect_error(RMitemRestscoreCutoff(d, iterations = 5, parallel = FALSE, seed = 1),
               pat)
  expect_error(RMdimResidualPCACutoff(d, iterations = 5, parallel = FALSE, seed = 1),
               pat)
  expect_error(RMlocdepQ3Cutoff(d, iterations = 5, parallel = FALSE, seed = 1),
               pat)
  expect_error(RMlocdepGammaCutoff(d, iterations = 5, parallel = FALSE, seed = 1),
               pat)
  grp <- rep(c("a", "b"), length.out = nrow(d))
  expect_error(RMdifGammaCutoff(d, dif_var = grp, iterations = 5,
                                parallel = FALSE, seed = 1, dgp = "conditional"),
               pat)
})

test_that("the DIF gamma permutation DGP needs no model and still runs", {
  d <- make_extreme_top()
  grp <- rep(c("a", "b"), length.out = nrow(d))
  res <- RMdifGammaCutoff(d, dif_var = grp, iterations = 5, parallel = FALSE,
                          seed = 1, dgp = "permutation")
  expect_equal(res$actual_iterations, 5L)
})

test_that("RMitemRestscoreBoot stops before resampling on an extreme-only category", {
  expect_error(
    RMitemRestscoreBoot(make_extreme_top(), iterations = 5, samplesize = 100,
                        parallel = FALSE, seed = 1),
    "zero or perfect total score: I6 \\(category 2\\).*Nothing was resampled",
    class = "easyRasch2_category_error"
  )
})

test_that(".all_sims_failed_message names what failed", {
  d <- make_extreme_top()
  expect_match(easyRasch2:::.all_sims_failed_message(d, what = "bootstrap"),
               "^All bootstrap iterations failed\\. Response categories")
  data("pcmdat2", package = "eRm")
  expect_equal(easyRasch2:::.all_sims_failed_message(pcmdat2),
               "All simulation iterations failed. Check your data.")
})
