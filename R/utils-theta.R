# Internal helpers for item- and person-parameter estimation.
#
# These low-level routines are shared by RMitemParameters() and
# RMpersonParameters(). They operate on a normalised item-parameter
# representation: a list `thr_list` with one element per item, each a
# numeric vector of Andrich thresholds (the category-boundary locations
# on the logit difficulty scale). A dichotomous item is simply an item
# with a single threshold equal to its difficulty, so one set of
# polytomous routines handles both the Rasch and Partial Credit cases.

# ---------------------------------------------------------------------
# Category-response probabilities (Partial Credit / Rasch)
# ---------------------------------------------------------------------

#' Category probabilities for one item at a given theta
#'
#' @param theta Numeric person location.
#' @param thr Numeric vector of Andrich thresholds (length = n_cat - 1).
#' @return Numeric vector of length `length(thr) + 1` giving the
#'   probability of scoring in categories 0, 1, ..., m. Computed with a
#'   log-sum-exp shift for numerical stability.
#' @keywords internal
#' @noRd
.pcm_cat_probs <- function(theta, thr) {
  cum_eta <- c(0, cumsum(theta - thr))
  shift   <- max(cum_eta)
  exp(cum_eta - (shift + log(sum(exp(cum_eta - shift)))))
}

#' Centre item thresholds to grand mean zero
#'
#' The standard Rasch identification: the mean of all thresholds is set
#' to zero so that person locations are on a common, comparable scale
#' across functions (matching the convention in [RMitemParameters()] and
#' the score-to-logit table).
#'
#' @param thr_list List of threshold vectors, one per item.
#' @return The list with every threshold shifted by the grand mean.
#' @keywords internal
#' @noRd
.center_thresholds <- function(thr_list) {
  shift <- mean(unlist(thr_list))
  lapply(thr_list, function(x) x - shift)
}

# ---------------------------------------------------------------------
# Per-pattern WLE (Warm, 1989) with analytic SEM
# ---------------------------------------------------------------------

#' Warm's weighted likelihood estimate for a single response pattern
#'
#' Newton-Raphson with adaptive step damping on the Warm-corrected score
#' equation. Only the items the person actually answered (non-`NA`)
#' contribute, so partial missingness is handled directly on the response
#' pattern rather than via a sum-score lookup.
#'
#' @param resp Numeric response vector for one person (may contain `NA`).
#' @param thr_list List of threshold vectors, one per item, aligned with
#'   `resp`.
#' @param theta_range Length-2 numeric search/boundary range.
#' @return Named numeric vector `c(theta, sem)`. Minimum and maximum
#'   possible scores yield finite WLE estimates (Warm's bias correction
#'   keeps the weighted likelihood equation solvable at the boundaries,
#'   unlike MLE); only if the root falls outside `theta_range` is the
#'   boundary returned with `NA` SEM.
#' @keywords internal
#' @noRd
.theta_wle <- function(resp, thr_list, theta_range, tol = 1e-6) {
  ok <- !is.na(resp)
  if (!any(ok)) return(c(theta = NA_real_, sem = NA_real_))

  resp     <- resp[ok]
  thr_list <- thr_list[ok]
  r         <- sum(resp)

  # Information at theta (also used for the SEM).
  info_at <- function(theta) {
    s <- 0
    for (i in seq_along(thr_list)) {
      cats <- 0:length(thr_list[[i]])
      p    <- .pcm_cat_probs(theta, thr_list[[i]])
      E    <- sum(cats * p)
      s    <- s + sum((cats - E)^2 * p)
    }
    s
  }

  # Warm's weighted-likelihood score function: l'(theta) + J(theta) / (2 I(theta)),
  # with J the summed third central moment of the item scores. Solving this for
  # zero gives the WLE, which is finite even at extreme scores.
  score <- function(theta) {
    dll <- r; info <- 0; third <- 0
    for (i in seq_along(thr_list)) {
      cats <- 0:length(thr_list[[i]])
      p    <- .pcm_cat_probs(theta, thr_list[[i]])
      E    <- sum(cats * p)
      dll   <- dll  - E
      info  <- info + sum((cats - E)^2 * p)
      third <- third + sum((cats - E)^3 * p)
    }
    dll + third / (2 * info)
  }

  lo <- theta_range[1L]; hi <- theta_range[2L]
  g_lo <- score(lo); g_hi <- score(hi)

  if (is.finite(g_lo) && is.finite(g_hi) && g_lo * g_hi < 0) {
    theta <- stats::uniroot(score, c(lo, hi), tol = tol)$root
    info  <- info_at(theta)
    sem   <- if (info < 1e-12) NA_real_ else 1 / sqrt(info)
  } else {
    # Root outside the search range: clamp to the nearer boundary.
    theta <- if (r <= sum(vapply(thr_list, length, integer(1))) / 2) lo else hi
    sem   <- NA_real_
  }
  c(theta = theta, sem = sem)
}

# ---------------------------------------------------------------------
# Quadrature machinery for EAP and prior-SD estimation
# ---------------------------------------------------------------------

#' Per-item log-probability tables over a quadrature grid
#'
#' Precomputes, for each item, a `length(grid) x n_cat` matrix of
#' log category probabilities. This lets per-person EAP and the marginal
#' likelihood reuse the same tables instead of recomputing probabilities.
#'
#' @param thr_list List of threshold vectors, one per item.
#' @param grid Numeric vector of quadrature nodes.
#' @return List of log-probability matrices, one per item.
#' @keywords internal
#' @noRd
.logp_tables <- function(thr_list, grid) {
  lapply(thr_list, function(thr) {
    tab <- vapply(grid, function(th) log(.pcm_cat_probs(th, thr)),
                  numeric(length(thr) + 1L))
    t(tab)   # rows = grid nodes, cols = categories 0..m
  })
}

#' Person-by-node log-likelihood matrix
#'
#' Sums the per-item log probabilities of each person's observed
#' responses across the quadrature grid, skipping missing responses.
#'
#' @param data_mat Numeric response matrix (persons x items).
#' @param logp_tabs List of per-item log-probability tables from
#'   [.logp_tables()].
#' @param grid Numeric vector of quadrature nodes.
#' @return Numeric matrix (persons x nodes) of log-likelihoods.
#' @keywords internal
#' @noRd
.grid_loglik <- function(data_mat, logp_tabs, grid) {
  n_persons <- nrow(data_mat)
  loglik    <- matrix(0, n_persons, length(grid))
  for (i in seq_along(logp_tabs)) {
    x   <- data_mat[, i]
    obs <- !is.na(x)
    if (!any(obs)) next
    # category index (1-based) for each observed response on item i
    cat_idx <- x[obs] + 1L
    loglik[obs, ] <- loglik[obs, ] + t(logp_tabs[[i]][, cat_idx, drop = FALSE])
  }
  loglik
}

#' Marginal-maximum-likelihood estimate of the prior SD
#'
#' Holds the item parameters fixed and finds the normal-prior SD that
#' maximises the marginal likelihood of the data over the quadrature
#' grid. Used as the default EAP prior when the user does not supply one.
#'
#' @param loglik Person-by-node log-likelihood matrix from [.grid_loglik()].
#' @param grid Numeric quadrature nodes.
#' @param prior_mean Numeric prior mean.
#' @return Numeric scalar: the estimated prior SD.
#' @keywords internal
#' @noRd
.estimate_prior_sd <- function(loglik, grid, prior_mean) {
  dgrid <- if (length(grid) > 1L) mean(diff(grid)) else 1
  neg_marg_ll <- function(sigma) {
    lprior <- stats::dnorm(grid, mean = prior_mean, sd = sigma, log = TRUE)
    # log integral over nodes for each person, via log-sum-exp
    M  <- sweep(loglik, 2L, lprior, `+`)
    mx <- apply(M, 1L, max)
    person_ll <- mx + log(rowSums(exp(M - mx)) * dgrid)
    -sum(person_ll[is.finite(person_ll)])
  }
  stats::optimize(neg_marg_ll, interval = c(0.05, 5))$minimum
}

#' EAP point estimates and posterior SDs from a grid log-likelihood
#'
#' @param loglik Person-by-node log-likelihood matrix.
#' @param grid Numeric quadrature nodes.
#' @param prior_mean,prior_sd Normal-prior parameters.
#' @return data.frame with columns `theta` (EAP) and `sem` (posterior SD).
#'   Persons with no answered items return `NA`.
#' @keywords internal
#' @noRd
.theta_eap <- function(loglik, grid, prior_mean, prior_sd) {
  lprior <- stats::dnorm(grid, mean = prior_mean, sd = prior_sd, log = TRUE)
  M      <- sweep(loglik, 2L, lprior, `+`)
  mx     <- apply(M, 1L, max)
  w      <- exp(M - mx)
  wsum   <- rowSums(w)

  eap  <- as.numeric((w %*% grid) / wsum)
  eap2 <- as.numeric((w %*% grid^2) / wsum)
  sd   <- sqrt(pmax(eap2 - eap^2, 0))
  # rows with no information (all-missing) -> NA
  empty <- !is.finite(eap)
  eap[empty] <- NA_real_
  sd[empty]  <- NA_real_
  data.frame(theta = eap, sem = sd)
}
