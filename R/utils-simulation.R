# Internal simulation helpers for parametric bootstrap functions

#' Simulate responses for a single polytomous item
#'
#' @param deltas Numeric vector of threshold parameters.
#' @param thetas Numeric vector of person parameters.
#' @return Integer vector of simulated responses (0-indexed).
#' @noRd
sim_poly_item <- function(deltas, thetas) {
  k <- length(deltas) + 1L
  n <- length(thetas)
  Y <- outer(thetas, deltas, "-")
  cumsums <- t(rbind(rep(0, times = n), apply(X = Y, MARGIN = 1, FUN = cumsum)))
  expcumsums <- exp(cumsums)
  norms <- apply(X = expcumsums, MARGIN = 1, FUN = sum)
  z <- expcumsums / norms
  vapply(X = seq_len(n), FUN = function(x) {
    sample(x = 0L:(k - 1L), size = 1L, replace = TRUE, prob = z[x, ])
  }, FUN.VALUE = 1L)
}

#' Simulate partial score (polytomous) response matrix
#'
#' @param deltaslist List of threshold parameter vectors, one per item.
#' @param thetavec Numeric vector of person parameters.
#' @return Integer matrix of simulated responses (rows = persons, cols = items).
#' @noRd
sim_partial_score <- function(deltaslist, thetavec) {
  deltaslist <- lapply(X = deltaslist, FUN = unlist)
  sapply(X = deltaslist, FUN = sim_poly_item, thetas = thetavec)
}

#' Extract item threshold parameters from polytomous data via CML
#'
#' Fits a Partial Credit Model using `eRm::PCM()` and returns the threshold
#' matrix (rows = items, cols = thresholds), without the location column.
#'
#' @param data A data.frame or matrix of polytomous item responses (min = 0).
#' @return Numeric matrix of threshold parameters.
#' @noRd
extract_item_thresholds <- function(data) {
  pcm_fit <- eRm::PCM(data)
  thresh <- eRm::thresholds(pcm_fit)
  thresh_mat <- thresh$threshtable$`1`
  if ("Location" %in% colnames(thresh_mat)) {
    thresh_mat <- thresh_mat[, colnames(thresh_mat) != "Location", drop = FALSE]
  }
  thresh_mat
}
