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

test_that("extract_item_thresholds works with polytomous data", {
  skip_if_not_installed("eRm")
  set.seed(7)
  # Simulate polytomous data (0-2)
  poly_data <- as.data.frame(
    matrix(sample(0:2, 120, replace = TRUE), nrow = 30, ncol = 4)
  )
  thresh <- easyRasch2:::extract_item_thresholds(poly_data)
  expect_true(is.matrix(thresh))
  # Should have 4 rows (one per item)
  expect_equal(nrow(thresh), 4L)
  # Location column should not be present
  expect_false("Location" %in% colnames(thresh))
})
