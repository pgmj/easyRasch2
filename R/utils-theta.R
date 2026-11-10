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
#' Solves the Warm-corrected score equation with `stats::uniroot()` on
#' `theta_range`. Only the items the person actually answered (non-`NA`)
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
    info  <- .test_information(thr_list, theta)
    sem   <- if (info < 1e-12) NA_real_ else 1 / sqrt(info)
  } else {
    # Root outside the search range: clamp to the boundary on the root's side.
    # The score function decreases in theta, so a positive value at both ends
    # puts the root above `hi`. The raw-score rule is only a fallback when the
    # score function is not finite at a boundary: it picks the wrong side when
    # the answered items all sit far from the person.
    theta <- if (is.finite(g_lo) && is.finite(g_hi)) {
      if (g_lo > 0) hi else lo
    } else if (r <= sum(vapply(thr_list, length, integer(1))) / 2) {
      lo
    } else {
      hi
    }
    sem <- NA_real_
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
#'   `.logp_tables()`.
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
#' @param loglik Person-by-node log-likelihood matrix from `.grid_loglik()`.
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
#'   A person with no answered items has a log-likelihood of zero at every
#'   node and so gets the prior mean and SD here. `.estimate_thetas()` sets
#'   those rows to `NA`, which this function cannot do from `loglik` alone.
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
  # Rows whose posterior cannot be formed (non-finite log-likelihood) -> NA.
  # All-missing rows are finite here and are handled in .estimate_thetas().
  empty <- !is.finite(eap)
  eap[empty] <- NA_real_
  sd[empty]  <- NA_real_
  data.frame(theta = eap, sem = sd)
}

#' Test information at one or more theta values
#'
#' Fisher test information for a Rasch / partial credit model,
#' \eqn{I(\theta) = \sum_i \mathrm{Var}_i(\mathrm{score} \mid \theta)}, summed
#' over items from the centred Andrich thresholds. The conditional standard
#' error of measurement follows as \eqn{1/\sqrt{I(\theta)}}.
#'
#' Extracted from `.marginal_rxx()` so that the same information curve backs
#' the marginal coefficient, `RMreliabilityCurve()`, and (later) the
#' conditional critical values in the change functions.
#'
#' @param thr_list List of centred Andrich threshold vectors, one per item.
#' @param theta Numeric vector of latent locations.
#' @return Numeric vector of test information, same length as `theta`.
#' @keywords internal
#' @noRd
.test_information <- function(thr_list, theta) {
  vapply(
    theta,
    function(th) {
      sum(vapply(
        thr_list,
        function(thr) {
          cats <- 0:length(thr)
          P <- .pcm_cat_probs(th, thr)
          E <- sum(cats * P)
          sum((cats - E)^2 * P)
        },
        numeric(1L)
      ))
    },
    numeric(1L)
  )
}

#' Marginal-maximum-likelihood estimate of the latent distribution
#'
#' Estimates **both** the mean and the SD of the latent distribution by
#' marginal maximum likelihood with the item parameters held fixed.
#'
#' Before 1.3.1 only the SD was estimated, with the mean held at 0. That is an
#' assumption rather than a consequence of centring the thresholds: grand-mean
#' centring fixes the *item* mean at 0 and says nothing about where the
#' respondents sit, and the distance between the two is targeting. With the
#' mean held at 0 a normal density has to reach an off-target sample by
#' widening, so \eqn{\sigma} absorbed the mistargeting and came out too large,
#' which inflated every quantity resting on it. Verified on simulated data
#' with a true SD of 0.90: shifting the sample 2 logits off target returned
#' 2.02 with the mean fixed and 0.96 with it free.
#'
#' The grid is built in two passes. The first is wide and fixed, to find the
#' mean without assuming where it is; the second re-centres on that estimate so
#' the quadrature covers the distribution it is integrating. One pass on a
#' fixed grid cannot do both.
#'
#' @param data Numeric response matrix or data.frame (items from 0).
#' @param thr_list List of centred Andrich threshold vectors.
#' @param n_nodes Number of quadrature nodes for the grid.
#' @return A list with `mean` and `sd`, or both `NA_real_` when the
#'   optimisation does not return usable values.
#' @keywords internal
#' @noRd
.latent_moments <- function(data, thr_list, n_nodes = 81L) {
  data_mat <- as.matrix(data)

  fit_on <- function(grid, start) {
    loglik <- .grid_loglik(data_mat, .logp_tables(thr_list, grid), grid)
    .estimate_prior_moments(loglik, grid, start = start)
  }

  # The first pass only has to find the mean, so it runs on a coarser grid:
  # its cost is linear in the node count and its answer is refined anyway.
  first <- fit_on(seq(-10, 10, length.out = max(41L, n_nodes %/% 2L)),
                  start = c(0, 0))
  if (!is.finite(first$mean) || !is.finite(first$sd) || first$sd <= 0) {
    return(list(mean = NA_real_, sd = NA_real_))
  }

  # The second pass starts from the first-pass answer, which lies inside its
  # grid by construction and saves objective evaluations in bootstrap loops.
  span <- max(6 * first$sd, 3)
  second <- fit_on(seq(first$mean - span, first$mean + span,
                       length.out = n_nodes),
                   start = c(first$mean, log(first$sd)))
  if (!is.finite(second$mean) || !is.finite(second$sd) || second$sd <= 0) {
    return(first)
  }
  second
}

#' Marginal-maximum-likelihood estimate of the latent mean and SD
#'
#' The two-parameter counterpart of `.estimate_prior_sd()`, which optimises the
#' SD alone against a fixed mean and is still what the EAP prior uses.
#'
#' The SD is optimised on the log scale and bounded to the same `[0.05, 5]`
#' interval `.estimate_prior_sd()` searches. The floor matters: data with no
#' latent variance at all, such as responses drawn at random, drives an
#' unbounded search to a SD of zero and every quantity resting on it to `NaN`.
#' The old estimator floored it implicitly through `optimize()`'s interval, and
#' keeping the same bound keeps degenerate data behaving as it did. The mean is
#' bounded to the quadrature grid, outside which the integral is not evaluated.
#'
#' Two parameters mean many more objective evaluations than the SD-only search
#' needed, so the person-by-node likelihood is exponentiated **once** per call
#' rather than once per evaluation. Writing `E = exp(loglik - offset)`, each
#' evaluation is then the matrix-vector product `E %*% w` with
#' `w = dnorm(grid, mu, sigma)`, which is one BLAS call instead of a `sweep()`,
#' a full `exp()` and a row-wise `apply()`. The offset is a fixed upper bound on
#' each row's maximum, so the log-sum-exp stays stable without recomputing the
#' row maxima for every candidate. On long scales a narrow candidate prior far
#' from a respondent can still underflow every product to 0; those rows are
#' recomputed in log space, so the objective equals the log-space one.
#'
#' The prior weights are normalised to sum to one over the grid. Without that,
#' a SD near the 0.05 floor on a coarse grid lets the optimiser put the mean on
#' a node and collect a weight far above one, so the first pass of
#' `.latent_moments()` (node spacing 0.5) could snap to a spurious narrow fit.
#' On the fine second-pass grid the normalising constant is one to within
#' rounding, so reported moments are essentially unchanged.
#'
#' @param loglik Person-by-node log-likelihood matrix from `.grid_loglik()`.
#' @param grid Numeric quadrature nodes.
#' @param start Starting values `c(mean, log(sd))`, moved inside the bounds.
#' @return A list with `mean` and `sd`.
#' @keywords internal
#' @noRd
.estimate_prior_moments <- function(loglik, grid, start = c(0, 0)) {

  # Row maxima of the log-likelihood, plus the largest log-density any
  # candidate prior can contribute, bound the row maxima of loglik + lprior
  # from above. Subtracting that keeps every exponential at or below 1.
  row_max_raw <- do.call(pmax, as.data.frame(loglik))
  row_max <- row_max_raw
  row_max[!is.finite(row_max)] <- 0
  max_lprior <- stats::dnorm(0, 0, 0.05, log = TRUE)
  offset <- row_max + max_lprior
  expl <- exp(loglik - offset)
  expl[!is.finite(expl)] <- 0

  neg_marg_ll <- function(par) {
    lprior <- stats::dnorm(grid, mean = par[1L], sd = exp(par[2L]), log = TRUE)
    # log of the grid sum of the prior weights, so the weights sum to one
    lmax <- max(lprior)
    log_norm <- lmax + log(sum(exp(lprior - lmax)))
    w <- exp(lprior)
    person_ll <- offset + log(as.numeric(expl %*% w)) - log_norm
    # A narrow candidate prior far from a respondent can make every product
    # underflow to 0, giving -Inf where the log-space sum is merely very
    # negative. Dropping such rows would make the candidate look better, so
    # they are recomputed in log space, exactly as before the speed-up.
    under <- !is.finite(person_ll) & is.finite(row_max_raw)
    if (any(under)) {
      M  <- sweep(loglik[under, , drop = FALSE], 2L, lprior, `+`)
      mx <- apply(M, 1L, max)
      person_ll[under] <- mx + log(rowSums(exp(M - mx))) - log_norm
    }
    out <- -sum(person_ll[is.finite(person_ll)])
    if (!is.finite(out)) .Machine$double.xmax else out
  }

  lo <- c(min(grid), log(0.05))
  hi <- c(max(grid), log(5))
  opt <- try(
    stats::optim(pmin(pmax(start, lo), hi), neg_marg_ll, method = "L-BFGS-B",
                 lower = lo, upper = hi),
    silent = TRUE
  )
  if (inherits(opt, "try-error") || !all(is.finite(opt$par))) {
    return(list(mean = NA_real_, sd = NA_real_))
  }
  list(mean = opt$par[1L], sd = exp(opt$par[2L]))
}

#' Latent-density-weighted marginal summaries of the precision curve
#'
#' The single place the marginal reliability scalars are formed, so that
#' `.marginal_rxx()` (behind `RMreliability()`) and `RMreliabilityCurve()`
#' cannot drift apart. All three are integrated over
#' \eqn{N(\mu, \sigma^2)} on a \eqn{\mu \pm 6\sigma} grid. Before 1.3.1
#' \eqn{\mu} was held at 0 rather than estimated; see `.latent_moments()`.
#'
#' @param thr_list List of Andrich threshold vectors.
#' @param sigma Latent SD, from `.latent_moments()`.
#' @param mu Latent mean, from `.latent_moments()`.
#' @param n_nodes Number of quadrature nodes.
#' @return A list with `sem_average` (root mean error variance), `ratio` (the
#'   density-weighted mean of the bounded reliability curve, which is what
#'   both functions report) and `green` (the superseded subtractive
#'   coefficient, retained for comparison).
#' @keywords internal
#' @noRd
.marginal_summaries <- function(thr_list, sigma, mu = 0, n_nodes = 161L) {
  g <- seq(mu - 6 * sigma, mu + 6 * sigma, length.out = n_nodes)
  w <- stats::dnorm(g, mu, sigma)
  w <- w / sum(w)
  sem2 <- 1 / .test_information(thr_list, g)
  sem2_bar <- sum(w * sem2)
  list(
    sem_average = sqrt(sem2_bar),
    ratio = sum(w * (sigma^2 / (sigma^2 + sem2))),
    green = 1 - sem2_bar / sigma^2
  )
}

# ---------------------------------------------------------------------
# Shared CML/WLE estimation engine (used by the Q3 local-dependence
# functions; intended as the common entry point for migrating other
# functions off mirt MML/EAP and eRm MLE).
# ---------------------------------------------------------------------

#' Conditional-maximum-likelihood item thresholds via psychotools
#'
#' Fits a Rasch model (dichotomous) or Partial Credit Model (polytomous)
#' by CML using `psychotools::raschmodel()` / `psychotools::pcmodel()`
#' and returns a list of grand-mean-centred Andrich threshold vectors,
#' one per item --- the normalised representation consumed by
#' `.pcm_cat_probs()`, `.theta_wle()`, and the residual routines.
#' `psychotools` is used (rather than `eRm`) because its CML fitter is
#' an order of magnitude faster, handles missing data without the
#' per-NA-pattern slowdown of `eRm`, and returns numerically identical
#' parameters after centring.
#'
#' Rows with no answered items are dropped before fitting: they carry no
#' information for CML, `pcmodel()` errors on them and `raschmodel()`
#' (psychotools 0.7.7) crashes the R session. `.check_cml_categories()` then
#' stops on items psychotools would silently reshape, and the fitted list is
#' checked against the data so a threshold vector can never be paired with
#' the wrong item.
#'
#' @param data Numeric response matrix or data.frame (items scored from 0).
#' @return List of centred Andrich threshold vectors, one per item.
#' @keywords internal
#' @noRd
.fit_cml_thresholds <- function(data) {
  data <- as.matrix(data)
  data <- data[rowSums(!is.na(data)) > 0L, , drop = FALSE]
  .check_cml_categories(data)
  is_poly <- max(data, na.rm = TRUE) > 1L
  # A dichotomous item is a 2-category PCM, so pcmodel() handles both types and
  # gives the identical CML estimate (same log-likelihood, difficulties match to
  # 0). raschmodel() is used for the dichotomous case only because it is ~2x
  # faster, which matters in the per-iteration bootstrap loops that call this.
  # The one-time `.rasch_fit_cml()` (item_parameters.R) uses pcmodel() for both,
  # where that speed difference is immaterial.
  fit <- if (is_poly) psychotools::pcmodel(data) else psychotools::raschmodel(data)
  thr_list <- lapply(psychotools::threshpar(fit), as.numeric)
  n_thr <- vapply(thr_list, length, integer(1L))
  max_obs <- apply(data, 2L, max, na.rm = TRUE)
  if (length(thr_list) != ncol(data) || any(n_thr != max_obs)) {
    stop(
      "The CML fit returned thresholds that do not match the observed ",
      "response categories, so item parameters cannot be paired with items.",
      call. = FALSE
    )
  }
  .center_thresholds(thr_list)
}

#' Stop on items that psychotools would silently reshape
#'
#' `psychotools::pcmodel()` and `raschmodel()` handle three item shapes in
#' ways that break the threshold list's alignment with the data, which every
#' consumer of `.fit_cml_thresholds()` assumes:
#' - an item with one observed category is dropped from the list,
#' - an empty lowest category makes psychotools rescale the item to start at
#'   0, while the data keep the original codes,
#' - an empty interior category loses its threshold, while the data keep the
#'   categories above it.
#' An empty top category is harmless (the item gets one threshold fewer and
#' the data agree), so it is not flagged. A category used only by extreme
#' scorers keeps its threshold in psychotools and is reported elsewhere by
#' `.extreme_only_categories()`.
#'
#' @param data Numeric response matrix, all-`NA` rows already removed.
#' @return Invisibly `TRUE`; stops naming each problem item otherwise.
#' @keywords internal
#' @noRd
.check_cml_categories <- function(data) {
  items <- colnames(data)
  if (is.null(items)) items <- paste0("item ", seq_len(ncol(data)))
  problems <- character(0)
  for (j in seq_len(ncol(data))) {
    x <- data[!is.na(data[, j]), j]
    if (length(x) == 0L) {
      problems <- c(problems, paste0(items[j], " (no responses)"))
      next
    }
    if (length(unique(x)) == 1L) {
      problems <- c(problems, paste0(items[j], " (only category ", x[1L], ")"))
      next
    }
    empty <- which(tabulate(x + 1L, nbins = max(x) + 1L) == 0L) - 1L
    if (length(empty) > 0L) {
      problems <- c(problems, paste0(
        items[j],
        if (length(empty) == 1L) " (empty category " else " (empty categories ",
        paste(empty, collapse = ", "), ")"
      ))
    }
  }
  if (length(problems) > 0L) {
    stop(
      "Item thresholds cannot be estimated for: ",
      paste(problems, collapse = ", "),
      ". Remove items with a single response category, and merge an empty ",
      "category with an adjacent one (an empty lowest category means ",
      "recoding the item to start at 0).",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Estimate person locations under a fixed item model
#'
#' Shared person-estimation entry point. Given centred item thresholds,
#' returns a per-respondent location estimate (and its SEM) by Warm's
#' weighted likelihood (`"WLE"`, the default) or expected a posteriori
#' (`"EAP"`). WLE depends only on the sufficient statistic --- the
#' answered-item set and the partial raw score over it --- so estimates
#' are solved once per distinct (answered-set, score) key and mapped
#' back to respondents. On complete data this collapses to one solve
#' per distinct raw score.
#'
#' @param data Numeric response matrix or data.frame (items from 0; `NA`
#'   allowed).
#' @param thr_list List of centred Andrich threshold vectors, one per
#'   item, aligned with the columns of `data` (e.g. from
#'   `.fit_cml_thresholds()`).
#' @param method `"WLE"` (default) or `"EAP"`.
#' @param theta_range Length-2 search/boundary range for WLE and the EAP
#'   quadrature grid. Person-reporting callers (`RMpersonParameters()`,
#'   `RMpersonFit()`) pass the wider `c(-10, 10)` to locate extreme scorers;
#'   the default `c(-6, 6)` is used by the parametric-bootstrap theta pools
#'   (via `.wle_theta_pool()`), where it is the validated/calibrated range and
#'   only the resampled extremes are affected.
#' @param prior_mean,prior_sd Normal-prior parameters for `"EAP"`. When
#'   `prior_sd` is `NULL` it is estimated from the data by marginal
#'   maximum likelihood.
#' @param n_nodes Number of EAP quadrature nodes (default 121; the single
#'   source for every EAP caller so the same data yields the same estimate).
#' @return data.frame with columns `theta` and `sem`, one row per
#'   respondent (in input order). Respondents with no answered items get
#'   `NA` under both methods, and are left out when `prior_sd` is estimated.
#'   For `method = "EAP"` the result carries
#'   `attr(., "prior_mean")` and `attr(., "prior_sd")` (the prior actually
#'   used).
#' @keywords internal
#' @noRd
.estimate_thetas <- function(data, thr_list, method = c("WLE", "EAP"),
                             theta_range = c(-6, 6), prior_mean = 0,
                             prior_sd = NULL, n_nodes = 121L) {
  method <- match.arg(method)
  data <- as.matrix(data)

  if (method == "EAP") {
    grid   <- seq(theta_range[1L], theta_range[2L], length.out = n_nodes)
    logp   <- .logp_tables(thr_list, grid)
    loglik <- .grid_loglik(data, logp, grid)
    answered <- rowSums(!is.na(data)) > 0L
    if (is.null(prior_sd)) {
      prior_sd <- .estimate_prior_sd(loglik[answered, , drop = FALSE], grid,
                                     prior_mean)
    }
    out <- .theta_eap(loglik, grid, prior_mean, prior_sd)
    # An all-missing row would otherwise report the prior as its estimate.
    out$theta[!answered] <- NA_real_
    out$sem[!answered]   <- NA_real_
    # Expose the prior used so callers (e.g. RMpersonParameters) can report it.
    attr(out, "prior_mean") <- prior_mean
    attr(out, "prior_sd")   <- prior_sd
    return(out)
  }

  # WLE: cache by sufficient statistic (answered-set, partial score).
  pat <- apply(!is.na(data), 1L, function(z) paste0(which(z), collapse = ","))
  rs  <- rowSums(data, na.rm = TRUE)
  key <- paste(pat, rs, sep = "|")
  rep_idx <- which(!duplicated(key))
  est <- vapply(rep_idx, function(i) .theta_wle(data[i, ], thr_list, theta_range),
                c(theta = 0, sem = 0))             # rows: theta, sem; cols: reps
  m <- match(key, key[rep_idx])
  data.frame(theta = est["theta", m], sem = est["sem", m])
}

#' Finite WLE person-location pool for parametric-bootstrap DGPs
#'
#' Shared setup for the simulation-cutoff functions: fits CML item thresholds
#' and returns them together with the finite WLE person locations used as the
#' resampling pool. Centralises the block repeated across the `*Cutoff()`
#' functions (Q3, infit, CFA, gamma DIF/LD).
#'
#' @param data Numeric response matrix or data.frame (items from 0).
#' @return A list with `thr_list` (centred CML thresholds) and `thetas` (the
#'   finite WLE person locations).
#' @keywords internal
#' @noRd
.wle_theta_pool <- function(data) {
  data     <- as.matrix(data)
  thr_list <- .fit_cml_thresholds(data)
  thetas   <- .estimate_thetas(data, thr_list, method = "WLE")$theta
  list(thr_list = thr_list, thetas = thetas[is.finite(thetas)])
}

#' Standardized Rasch/PCM residual matrix (CML items + chosen person estimator)
#'
#' Computes the model standardized residuals \eqn{(x - E)/\sqrt{Var}} for
#' every person-by-item cell, using CML item thresholds and WLE (or EAP)
#' person locations. This is the residual object underlying Yen's Q3
#' (its column-wise correlation matrix) and is reusable for residual-based
#' diagnostics. Non-finite residuals (e.g. zero-variance cells) are set
#' to `NA`.
#'
#' @param data Numeric response matrix or data.frame (items from 0; `NA`
#'   allowed).
#' @param method Person estimator passed to `.estimate_thetas()`.
#' @return Numeric matrix (persons x items) of standardized residuals,
#'   with the item names of `data` as column names.
#' @keywords internal
#' @noRd
.rasch_std_residuals <- function(data, method = "WLE") {
  data <- as.matrix(data)
  thr_list <- .fit_cml_thresholds(data)
  th <- .estimate_thetas(data, thr_list, method = method)$theta
  N <- nrow(data); K <- ncol(data)
  E <- W <- matrix(NA_real_, N, K)
  for (i in seq_len(K)) {
    cats <- 0:length(thr_list[[i]])
    P  <- vapply(th, function(z) .pcm_cat_probs(z, thr_list[[i]]),
                 numeric(length(cats)))            # (n_cat x N)
    Ei <- as.numeric(crossprod(cats, P))
    W[, i] <- as.numeric(crossprod(cats^2, P)) - Ei^2
    E[, i] <- Ei
  }
  R <- (data - E) / sqrt(W)
  R[!is.finite(R)] <- NA
  colnames(R) <- colnames(data)
  R
}
