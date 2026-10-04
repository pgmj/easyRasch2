# Tests for .check_cutoff_sample(), the warning that a cutoff object was
# simulated for a different sample than the data it is applied to.

test_that("silent when the samples agree or cannot be compared", {
  expect_silent(.check_cutoff_sample(200, 200, "RMitemInfitCutoff()"))
  expect_silent(.check_cutoff_sample(NULL, 200, "RMitemInfitCutoff()"))
  expect_silent(.check_cutoff_sample(NA_real_, 200, "RMitemInfitCutoff()"))
  expect_silent(.check_cutoff_sample(
    200, 200, "RMdifGammaCutoff()",
    groups_cutoff = c(100L, 100L), groups_used = c(100L, 100L)
  ))
  # Group sizes absent (older DIF objects) skip the group comparison
  expect_silent(.check_cutoff_sample(
    200, 200, "RMdifGammaCutoff()", groups_used = c(150L, 50L)
  ))
})

test_that("warns on a different sample size, naming the policy", {
  expect_warning(
    .check_cutoff_sample(200, 400, "RMitemInfitCutoff()"),
    "simulated for n = 200 but `data` has n = 400 \\(complete cases\\)"
  )
  expect_warning(
    .check_cutoff_sample(200, 190, "RMlocdepQ3Cutoff()",
                         policy = "respondents with at least one response"),
    "at least one response"
  )
})

test_that("warns on different DIF group sizes at the same n", {
  expect_warning(
    .check_cutoff_sample(
      200, 200, "RMdifGammaCutoff()",
      groups_cutoff = c(100L, 100L), groups_used = c(150L, 50L)
    ),
    "group sizes 100/100 but `dif_var` has 150/50"
  )
})

test_that("consumers warn when handed a cutoff from another sample", {
  skip_on_cran()
  skip_if_not_installed("iarm")
  set.seed(1)
  mk <- function(n) {
    d <- as.data.frame(matrix(sample(0:1, n * 5, replace = TRUE), n, 5))
    colnames(d) <- paste0("I", 1:5)
    d
  }
  small <- mk(150)
  big <- mk(300)
  ci <- RMitemInfitCutoff(small, iterations = 5, parallel = FALSE, seed = 1,
                          cutoff_method = "quantile")
  expect_warning(suppressMessages(RMitemInfit(big, cutoff = ci)),
                 "simulated for n = 150")
  expect_no_warning(suppressMessages(RMitemInfit(small, cutoff = ci)))

  cr <- RMitemRestscoreCutoff(small, iterations = 5, parallel = FALSE,
                              seed = 1, cutoff_method = "quantile")
  expect_warning(suppressMessages(RMitemRestscore(big, cutoff = cr)),
                 "simulated for n = 150")

  # A complete-case cutoff matches data with missing values when the
  # complete-case counts agree
  small_na <- small
  small_na[1:5, 2] <- NA
  cr_na <- RMitemRestscoreCutoff(small_na, iterations = 5, parallel = FALSE,
                                 seed = 1, cutoff_method = "quantile")
  expect_no_warning(suppressMessages(RMitemRestscore(small_na, cutoff = cr_na)))
})
