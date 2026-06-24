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

# CML item thresholds are now obtained via `.fit_cml_thresholds()` (psychotools,
# grand-mean centred; see utils-theta.R), which returns the per-item threshold
# list the bootstrap generators consume directly. The former eRm-based
# `extract_item_thresholds()` matrix helper has been removed.
