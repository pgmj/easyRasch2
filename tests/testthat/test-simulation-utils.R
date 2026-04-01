test_that("sim_poly_item returns vector of correct length and valid range", {
  deltas <- c(-1.5, 0, 1.5)
  thetas <- seq(-2, 2, length.out = 50)
  result <- easyRasch2:::sim_poly_item(deltas, thetas)
  expect_length(result, 50)
  expect_true(all(result >= 0L))
  expect_true(all(result <= length(deltas)))
})

test_that("sim_partial_score returns matrix with correct dimensions", {
  deltas1 <- c(-1, 1)
  deltas2 <- c(-0.5, 0.5)
  deltaslist <- list(list(deltas1), list(deltas2))
  thetas <- seq(-2, 2, length.out = 30)
  result <- easyRasch2:::sim_partial_score(deltaslist, thetas)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 30)
  expect_equal(ncol(result), 2)
  expect_true(all(result >= 0))
})

test_that("extract_item_thresholds returns a matrix", {
  skip_if_not_installed("eRm")
  set.seed(10)
  # Create simple polytomous data (3 categories, 0-2)
  df <- as.data.frame(
    matrix(sample(0:2, 200, replace = TRUE), nrow = 50, ncol = 4)
  )
  result <- easyRasch2:::extract_item_thresholds(df)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4)
})
