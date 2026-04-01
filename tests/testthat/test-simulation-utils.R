test_that("sim_poly_item returns vector of correct length with values in valid range", {
  thetas <- rnorm(50)
  deltas <- c(-1, 0, 1)  # 4-category item
  result <- easyRasch2:::sim_poly_item(deltas, thetas)
  expect_length(result, 50L)
  expect_true(all(result %in% 0L:3L))
})

test_that("sim_poly_item handles single-threshold (dichotomous-like) item", {
  thetas <- rnorm(30)
  deltas <- 0.5  # 2-category item
  result <- easyRasch2:::sim_poly_item(deltas, thetas)
  expect_length(result, 30L)
  expect_true(all(result %in% 0L:1L))
})

test_that("sim_partial_score returns matrix of correct dimensions", {
  thetas <- rnorm(40)
  itemlist <- list(
    list(c(-1, 0)),
    list(c(-0.5, 0.5)),
    list(c(0, 1))
  )
  result <- easyRasch2:::sim_partial_score(deltaslist = itemlist, thetavec = thetas)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 40L)
  expect_equal(ncol(result), 3L)
})

test_that("extract_item_thresholds returns a matrix", {
  skip_on_cran()
  skip_if_not_installed("eRm")
  set.seed(1)
  # Create minimal polytomous data (0-2 scoring)
  df <- as.data.frame(matrix(
    sample(0:2, 200, replace = TRUE, prob = c(0.3, 0.4, 0.3)),
    nrow = 50, ncol = 4
  ))
  result <- easyRasch2:::extract_item_thresholds(df)
  expect_true(is.matrix(result))
  expect_equal(nrow(result), 4L)  # one row per item
})
