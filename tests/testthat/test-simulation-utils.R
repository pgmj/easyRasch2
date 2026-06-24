test_that("sim_poly_item returns integer vector of correct length", {
  deltas <- c(-1, 0, 1)
  thetas <- rnorm(50)
  result <- easyRasch2:::sim_poly_item(deltas, thetas)
  expect_length(result, 50)
  expect_true(all(result >= 0L & result <= length(deltas)))
})

test_that("sim_partial_score returns matrix with correct dimensions", {
  set.seed(1)
  deltaslist <- list(c(-1, 1), c(-0.5, 0.5), c(-2, 0, 2))
  thetavec <- rnorm(30)
  result <- easyRasch2:::sim_partial_score(deltaslist, thetavec)
  expect_equal(dim(result), c(30L, 3L))
})

test_that("sim_poly_item responses are within valid range", {
  set.seed(42)
  deltas <- c(-1, 1)   # 3-category item
  thetas <- rnorm(100)
  result <- easyRasch2:::sim_poly_item(deltas, thetas)
  expect_true(all(result %in% 0:2))
})

test_that(".fit_cml_thresholds returns a per-item threshold list (poly)", {
  set.seed(7)
  poly_data <- as.data.frame(
    matrix(sample(0:2, 120, replace = TRUE), nrow = 30, ncol = 4)
  )
  thr <- easyRasch2:::.fit_cml_thresholds(poly_data)
  expect_type(thr, "list")
  # One element per item, each with (categories - 1) = 2 thresholds
  expect_length(thr, 4L)
  expect_true(all(vapply(thr, length, integer(1L)) == 2L))
  # Grand-mean centred
  expect_equal(mean(unlist(thr)), 0, tolerance = 1e-8)
})
