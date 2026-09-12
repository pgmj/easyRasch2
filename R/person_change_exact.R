# Exact null distribution of the Rasch change index.
#
# With complete data the sum score is a sufficient statistic, so theta-hat and
# its standard error are deterministic functions of the score and the whole
# null distribution of the RCI can be enumerated rather than simulated. See
# dev/rci_exact_null.qmd for the derivation and the validation against
# simulation.

#' Lord-Wingersky recursion: exact sum-score distribution given theta
#'
#' Convolves the per-item category probabilities one item at a time.
#'
#' @param thr_list List of Andrich threshold vectors for the answered items.
#' @param theta Numeric scalar.
#' @return Numeric vector of probabilities, index `score + 1`.
#' @keywords internal
#' @noRd
.pc_score_dist <- function(thr_list, theta) {
  d <- 1
  for (thr in thr_list) {
    p <- .pcm_cat_probs(theta, thr)
    nd <- numeric(length(d) + length(p) - 1L)
    for (a in seq_along(d)) {
      idx <- a + seq_along(p) - 1L
      nd[idx] <- nd[idx] + d[a] * p
    }
    d <- nd
  }
  d
}

#' Score distribution averaged over the occasion deviation
#'
#' Under `null = "retest"` each occasion carries an independent deviation
#' \eqn{u \sim N(0, \sigma^2)}. Because the two deviations are independent, the
#' joint over the two scores factorises into the product of two identical
#' averaged distributions, so the retest null costs no more than the
#' measurement null.
#'
#' @keywords internal
#' @noRd
.pc_score_dist_retest <- function(thr_list, theta, sigma, n_nodes = 21L) {
  if (sigma <= 0) {
    return(.pc_score_dist(thr_list, theta))
  }
  nodes <- seq(-4 * sigma, 4 * sigma, length.out = n_nodes)
  w <- stats::dnorm(nodes, 0, sigma)
  w <- w / sum(w)
  out <- 0
  for (i in seq_along(nodes)) {
    out <- out + w[i] * .pc_score_dist(thr_list, theta + nodes[i])
  }
  out
}

#' Score-to-theta lookup for one answered-item set
#'
#' Uses the shared estimation engine, so the locations and standard errors
#' match [RMpersonParameters()] and [RMscoreSE()].
#'
#' @param thr_list Full threshold list.
#' @param answered Integer indices of the answered items.
#' @param method `"WLE"` or `"EAP"`.
#' @param theta_range Search range.
#' @param prior_sd EAP prior SD, taken from the observed fit so the lookup uses
#'   the same prior the respondents were scored under.
#' @return A list with `theta` and `sem`, indexed by `score + 1`.
#' @keywords internal
#' @noRd
.pc_score_lookup <- function(thr_list, answered, method, theta_range,
                             prior_sd = NULL) {
  steps <- vapply(thr_list, length, integer(1L))
  sub_steps <- steps[answered]
  max_r <- sum(sub_steps)
  pat <- t(vapply(0:max_r, function(r) {
    resp <- rep(NA_integer_, length(steps))
    resp[answered] <- as.integer(.score_pattern(r, sub_steps))
    resp
  }, integer(length(steps))))
  est <- .estimate_thetas(
    pat, thr_list, method = method, theta_range = theta_range,
    prior_sd = prior_sd
  )
  list(theta = est$theta, sem = est$sem)
}

#' Exact critical values and p-values for the change index
#'
#' Enumerates the null over score pairs. Respondents are grouped by their pair
#' of answered-item sets, so incomplete data is handled by enumerating within
#' each pattern rather than by falling back to simulation.
#'
#' @return A list with `crit` (a list of `lower` and `upper` vectors) and
#'   `p_value`.
#' @keywords internal
#' @noRd
.pc_exact_null <- function(
  thr_list,
  theta_null,
  mat_t1,
  mat_t2,
  rci_obs,
  sigma_retest,
  method,
  theta_range,
  prior_sd,
  probs,
  direction,
  conditional_crit,
  n_nodes = 21L
) {
  n <- length(theta_null)
  key1 <- apply(!is.na(mat_t1), 1L, function(z) paste0(which(z), collapse = ","))
  key2 <- apply(!is.na(mat_t2), 1L, function(z) paste0(which(z), collapse = ","))

  # One lookup per distinct answered-item set, shared across occasions.
  lookups <- new.env(parent = emptyenv())
  get_lookup <- function(key) {
    if (!nzchar(key)) {
      return(NULL)
    }
    if (!is.null(lookups[[key]])) {
      return(lookups[[key]])
    }
    answered <- as.integer(strsplit(key, ",", fixed = TRUE)[[1L]])
    lk <- .pc_score_lookup(thr_list, answered, method, theta_range, prior_sd)
    lookups[[key]] <- lk
    lk
  }

  # The RCI grid depends only on the pair of answered sets, so cache it too.
  grids <- new.env(parent = emptyenv())
  get_grid <- function(k1, k2) {
    gk <- paste(k1, k2, sep = "|")
    if (!is.null(grids[[gk]])) {
      return(grids[[gk]])
    }
    l1 <- get_lookup(k1)
    l2 <- get_lookup(k2)
    g <- outer(
      seq_along(l1$theta), seq_along(l2$theta),
      function(i, j) {
        (l2$theta[j] - l1$theta[i]) /
          sqrt(l1$sem[i]^2 + l2$sem[j]^2 + 2 * sigma_retest^2)
      }
    )
    grids[[gk]] <- g
    g
  }

  crit_lower <- rep(NA_real_, n)
  crit_upper <- rep(NA_real_, n)
  p_value <- rep(NA_real_, n)

  pooled_val <- vector("list", n)
  pooled_wt <- vector("list", n)

  per_person <- vector("list", n)

  for (i in seq_len(n)) {
    if (!nzchar(key1[i]) || !nzchar(key2[i]) || !is.finite(theta_null[i])) {
      next
    }
    a1 <- as.integer(strsplit(key1[i], ",", fixed = TRUE)[[1L]])
    a2 <- as.integer(strsplit(key2[i], ",", fixed = TRUE)[[1L]])
    d1 <- .pc_score_dist_retest(thr_list[a1], theta_null[i], sigma_retest,
                                n_nodes)
    d2 <- .pc_score_dist_retest(thr_list[a2], theta_null[i], sigma_retest,
                                n_nodes)
    g <- get_grid(key1[i], key2[i])
    w <- outer(d1, d2)
    per_person[[i]] <- list(v = as.numeric(g), w = as.numeric(w) / sum(w))
    pooled_val[[i]] <- as.numeric(g)
    pooled_wt[[i]] <- as.numeric(w) / sum(w)
  }

  ok <- !vapply(per_person, is.null, logical(1L))
  if (!any(ok)) {
    return(list(
      crit = list(lower = crit_lower, upper = crit_upper),
      p_value = p_value
    ))
  }

  if (isTRUE(conditional_crit)) {
    for (i in which(ok)) {
      cr <- .pc_disc_quantile(per_person[[i]]$v, per_person[[i]]$w, probs)
      crit_lower[i] <- cr[1L]
      crit_upper[i] <- cr[2L]
      p_value[i] <- .pc_disc_pvalue(
        per_person[[i]]$v, per_person[[i]]$w, rci_obs[i], direction
      )
    }
  } else {
    v <- unlist(pooled_val[ok], use.names = FALSE)
    w <- unlist(pooled_wt[ok], use.names = FALSE)
    w <- w / sum(w)
    cr <- .pc_disc_quantile(v, w, probs)
    crit_lower[ok] <- cr[1L]
    crit_upper[ok] <- cr[2L]
    p_value[ok] <- vapply(
      which(ok),
      function(i) .pc_disc_pvalue(v, w, rci_obs[i], direction),
      numeric(1L)
    )
  }

  list(
    crit = list(lower = crit_lower, upper = crit_upper),
    p_value = p_value
  )
}

#' Quantiles of a weighted discrete distribution
#'
#' Returns the smallest value whose cumulative weight reaches each probability,
#' which is the usual convention for a discrete null. An infinite probability
#' marks the unused side of a one-sided test.
#'
#' @keywords internal
#' @noRd
.pc_disc_quantile <- function(v, w, probs) {
  o <- order(v)
  vs <- v[o]
  cw <- cumsum(w[o])
  pick <- function(p) {
    if (!is.finite(p)) {
      return(if (p < 0) -Inf else Inf)
    }
    idx <- which(cw >= p - 1e-12)[1L]
    if (is.na(idx)) vs[length(vs)] else vs[idx]
  }
  c(pick(probs$lo), pick(probs$hi))
}

#' Exact p-value against a weighted discrete null
#'
#' The observed point mass is included, which is the conservative convention
#' for a discrete test and the counterpart of the `(1 + count) / (B + 1)` rule
#' used on the simulated path.
#'
#' @keywords internal
#' @noRd
.pc_disc_pvalue <- function(v, w, obs, direction) {
  if (!is.finite(obs)) {
    return(NA_real_)
  }
  tol <- 1e-12
  keep <- switch(direction,
    two.sided = abs(v) >= abs(obs) - tol,
    increase = v >= obs - tol,
    decrease = v <= obs + tol
  )
  min(sum(w[keep]), 1)
}
