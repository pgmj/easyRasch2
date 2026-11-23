# `iterations` is checked the same way in every simulation function (F1).

test_that(".check_iterations accepts whole numbers at or above the minimum", {
  chk <- easyRasch2:::.check_iterations
  expect_invisible(chk(1))
  expect_invisible(chk(400L))
  expect_invisible(chk(0, min = 0L))
  expect_invisible(chk(20, min = 20L))
})

test_that(".check_iterations rejects everything else", {
  chk <- easyRasch2:::.check_iterations
  msg <- "`iterations` must be a single integer >= 1."
  for (bad in list(0, -5, 2.5, NA_real_, Inf, "10", c(10, 20), NULL)) {
    expect_error(chk(bad), msg, fixed = TRUE)
  }
  expect_error(chk(19, min = 20L), ">= 20", fixed = TRUE)
})

test_that("simulation functions stop on invalid iterations before doing work", {
  d <- as.data.frame(matrix(sample(0:2, 400, replace = TRUE), 100, 4))
  colnames(d) <- paste0("I", 1:4)
  pat <- "`iterations` must be a single integer >= 1"
  expect_error(RMitemInfitCutoff(d, iterations = 0), pat, fixed = TRUE)
  expect_error(RMitemRestscoreCutoff(d, iterations = -1), pat, fixed = TRUE)
  expect_error(RMitemRestscoreBoot(d, iterations = 2.5), pat, fixed = TRUE)
  expect_error(RMdimResidualPCACutoff(d, iterations = "a"), pat, fixed = TRUE)
  expect_error(RMdimCFACutoff(d, iterations = 0), pat, fixed = TRUE)
  expect_error(RMlocdepQ3Cutoff(d, iterations = 0), pat, fixed = TRUE)
  expect_error(RMlocdepGammaCutoff(d, iterations = 0), pat, fixed = TRUE)
  expect_error(RMdifGammaCutoff(d, dif_var = rep(1:2, 50), iterations = 0),
               pat, fixed = TRUE)
  expect_error(RMdimMartinLof(d, partition = list(1:2, 3:4), iterations = 0),
               pat, fixed = TRUE)
  # RMpersonFit: 0 means "statistics only", so only negatives stop
  expect_error(RMpersonFit(d, iterations = -1), ">= 0", fixed = TRUE)
  expect_error(RMitemCatFit(d, iterations = 10), ">= 20", fixed = TRUE)
})
