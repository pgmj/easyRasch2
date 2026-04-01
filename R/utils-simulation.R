# Internal helpers for simulation used by RMlocdepQ3cutoff()

#' Simulate responses for one polytomous item
#'
#' Given a vector of threshold parameters and a vector of person abilities,
#' returns a vector of simulated item responses under the Partial Credit Model.
#'
#' @param deltas Numeric vector of item threshold parameters (length = k - 1
#'   where k is the number of response categories).
#' @param thetas Numeric vector of person ability parameters.
#'
#' @return Integer vector of length `length(thetas)` with values in
#'   `0:(length(deltas))`.
#'
#' @noRd
sim_poly_item <- function(deltas, thetas) {
  k <- length(deltas) + 1L
  n <- length(thetas)
  Y <- outer(thetas, deltas, "-")
  cumsums <- t(rbind(rep(0, times = n), apply(X = Y, MARGIN = 1, FUN = cumsum)))
  expcumsums <- exp(cumsums)
  norms <- apply(X = expcumsums, MARGIN = 1, FUN = sum)
  z <- expcumsums / norms
  vapply(X = seq_len(n), FUN = function(j) {
    sample(x = 0L:(k - 1L), size = 1L, replace = TRUE, prob = z[j, ])
  }, FUN.VALUE = 1L)
}

#' Simulate a full polytomous dataset
#'
#' Applies `sim_poly_item()` to each item in `deltaslist` to produce a
#' simulated response matrix.
#'
#' @param deltaslist A list of length equal to the number of items. Each
#'   element is itself a list (or vector when unlisted) of threshold parameters
#'   for that item.
#' @param thetavec Numeric vector of person ability parameters.
#'
#' @return A matrix with `length(thetavec)` rows and `length(deltaslist)`
#'   columns.
#'
#' @noRd
sim_partial_score <- function(deltaslist, thetavec) {
  deltaslist <- lapply(X = deltaslist, FUN = unlist)
  sapply(X = deltaslist, FUN = sim_poly_item, thetas = thetavec)
}

#' Extract item threshold matrix from polytomous data
#'
#' Fits a Partial Credit Model via `eRm::PCM()` and extracts the item
#' threshold parameters using `eRm::thresholds()`.
#'
#' @param data A data.frame or matrix of polytomous item responses scored
#'   starting at 0.
#'
#' @return A numeric matrix with one row per item and one column per threshold.
#'   Items with fewer categories than the maximum will have `NA` in unused
#'   threshold columns.
#'
#' @noRd
extract_item_thresholds <- function(data) {
  pcm_fit <- eRm::PCM(data)
  thresh_obj <- eRm::thresholds(pcm_fit)
  # thresh_obj$threshpar is a matrix: items x thresholds
  thresh_obj$threshpar
}
