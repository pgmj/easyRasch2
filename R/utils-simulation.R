# Internal simulation helpers for parametric bootstrap functions

#' Simulate responses for one polytomous item
#'
#' @param deltas Numeric vector of threshold (delta) parameters for the item.
#' @param thetas Numeric vector of person ability parameters.
#'
#' @return Integer vector of length `length(thetas)` with simulated responses
#'   in the range `0:(length(deltas))`.
#'
#' @noRd
sim_poly_item <- function(deltas, thetas) {
  k <- length(deltas) + 1L
  n <- length(thetas)
  Y <- outer(thetas, deltas, "-")
  cumsums <- t(rbind(rep(0, times = n), apply(X = Y, MARGIN = 1L, FUN = cumsum)))
  expcumsums <- exp(cumsums)
  norms <- apply(X = expcumsums, MARGIN = 1L, FUN = sum)
  z <- expcumsums / norms
  vapply(X = seq_len(n), FUN = function(i) {
    sample(x = 0L:(k - 1L), size = 1L, replace = TRUE, prob = z[i, ])
  }, FUN.VALUE = integer(1L))
}

#' Simulate a full polytomous dataset
#'
#' @param deltaslist A list of length equal to the number of items. Each
#'   element should be a list containing a single numeric vector of threshold
#'   parameters for that item (as returned by [extract_item_thresholds()]).
#' @param thetavec Numeric vector of person ability parameters.
#'
#' @return An integer matrix with `length(thetavec)` rows and
#'   `length(deltaslist)` columns.
#'
#' @noRd
sim_partial_score <- function(deltaslist, thetavec) {
  deltaslist <- lapply(X = deltaslist, FUN = function(x) unlist(x))
  sapply(X = deltaslist, FUN = sim_poly_item, thetas = thetavec)
}

#' Extract item threshold matrix from polytomous data
#'
#' Uses an existing `eRm::PCM()` fit (or fits one if not provided) and returns
#' a matrix of item thresholds (one row per item, one column per threshold;
#' shorter items are padded with `NA`).
#'
#' @param data A data.frame or matrix of polytomous item responses scored from
#'   0. Only used when `pcm_fit` is `NULL`.
#' @param pcm_fit An optional `eRm` object from `eRm::PCM()`. When supplied,
#'   `data` is ignored.
#'
#' @return A numeric matrix with `ncol(data)` rows. Column \eqn{j} contains
#'   the \eqn{j}th threshold for each item, padded with `NA` when an item has
#'   fewer thresholds.
#'
#' @noRd
extract_item_thresholds <- function(data, pcm_fit = NULL) {
  if (is.null(pcm_fit)) pcm_fit <- eRm::PCM(data)
  thresh_out <- eRm::thresholds(pcm_fit)
  # thresh_out$threshtable is a list of matrices, one per item grouping.
  # For a standard PCM the relevant thresholds are in thresh_out$threshtable[[1]].
  thresh_mat_raw <- thresh_out$threshtable[[1L]]
  # thresh_mat_raw is a matrix: rows = items, cols = thresholds (may include Location)
  # Remove "Location" column if present
  loc_col <- which(colnames(thresh_mat_raw) == "Location")
  if (length(loc_col) > 0L) {
    thresh_mat_raw <- thresh_mat_raw[, -loc_col, drop = FALSE]
  }
  thresh_mat_raw
}
