# Internal engine for RMdifImpact().
#
# Ported from the development study (dev/dif_shift_cml.R,
# dev/dif_impact_theory.R, dev/dif_impact_sim.R, dev/dif_impact_ext_sim.R,
# dev/dif_impact_ext_sim_3b.R; design and validation in
# dev/dif-impact-design.md, Part B, and dev/dif_impact.qmd). The intervals
# that the function reports are the ones validated there.
#
# Notation. Two groups, reference (1) and focal (2). Threshold lists are
# lists of Andrich threshold vectors, one per item. A split model gives each
# DIF item a uniform focal shift: focal thresholds = reference thresholds +
# shift. Its thresholds are centred so that the reference group's thresholds
# have mean zero, and every theta quantity is on that scale.

# ---------------------------------------------------------------------
# Probability machinery on a theta grid
# ---------------------------------------------------------------------

#' Category probabilities of one item on a theta grid
#' @return Matrix, one row per grid point, columns categories 0..m.
#' @keywords internal
#' @noRd
.dimp_cat_probs <- function(thr, grid) {
  m <- length(thr)
  eta <- outer(grid, 0:m) -
    matrix(c(0, cumsum(thr)), length(grid), m + 1L, byrow = TRUE)
  eta <- eta - do.call(pmax, as.data.frame(eta))
  p <- exp(eta)
  p / rowSums(p)
}

#' Convolve a raw-score distribution with one more item
#' @keywords internal
#' @noRd
.dimp_conv <- function(SP, pc) {
  m <- ncol(pc) - 1L
  w <- ncol(SP)
  new <- matrix(0, nrow(SP), w + m)
  for (c in 0:m) {
    idx <- (c + 1L):(c + w)
    new[, idx] <- new[, idx] + SP * pc[, c + 1L]
  }
  new
}

#' P(raw score = r | theta), r = 0..R, on a grid
#' @keywords internal
#' @noRd
.dimp_score_probs <- function(thr_list, grid) {
  out <- matrix(1, length(grid), 1L)
  for (thr in thr_list) out <- .dimp_conv(out, .dimp_cat_probs(thr, grid))
  out
}

#' Sums over items of the expected score, information and third central
#' moment at each theta
#' @keywords internal
#' @noRd
.dimp_moments <- function(thr_list, theta) {
  E <- I <- J <- 0
  for (thr in thr_list) {
    p <- .dimp_cat_probs(thr, theta)
    cats <- 0:length(thr)
    e <- drop(p %*% cats)
    dv <- outer(-e, cats, "+")
    E <- E + e
    I <- I + rowSums(p * dv^2)
    J <- J + rowSums(p * dv^3)
  }
  list(E = E, I = I, J = J)
}

#' Normal quadrature (equally spaced, +-7 SD)
#' @keywords internal
#' @noRd
.dimp_quad <- function(mu, sigma, n_nodes = 201L) {
  g <- seq(mu - 7 * sigma, mu + 7 * sigma, length.out = n_nodes)
  w <- stats::dnorm(g, mu, sigma)
  list(grid = g, w = w / sum(w))
}

#' Expected sufficient statistics per person for latent N(mu, sigma)
#' @return List with `pr` (raw-score distribution) and `ctot` (expected
#'   count per category above 0, items in order).
#' @keywords internal
#' @noRd
.dimp_pop_stats <- function(thr_list, mu, sigma) {
  q <- .dimp_quad(mu, sigma)
  pr <- colSums(q$w * .dimp_score_probs(thr_list, q$grid))
  ctot <- unlist(lapply(thr_list, function(thr) {
    colSums(q$w * .dimp_cat_probs(thr, q$grid))[-1L]
  }))
  list(pr = pr, ctot = ctot)
}

# ---------------------------------------------------------------------
# Conditional ML for the PCM with uniform focal shifts
# ---------------------------------------------------------------------
#
# Item j has cumulative category parameters psi_j1..psi_jm (pcmodel()'s
# parameterisation). For a DIF item, the focal group has psi_jc + c * shift,
# which moves every threshold of the item by the shift. With two groups the
# conditional likelihood is the sum of the groups' PCM conditional
# likelihoods, each from its category counts and raw-score distribution. The
# first category parameter of the first item is fixed at 0.

#' Model structure and sufficient statistics for the shift CML
#' @param m Number of thresholds per item.
#' @param dif_items Indices of the shifted items.
#' @param stats Two lists (reference, focal) with `ctot` (category counts
#'   above 0) and `ptot` (raw-score counts 0..R).
#' @keywords internal
#' @noRd
.dimp_prep <- function(m, dif_items, stats) {
  P <- sum(m)
  item_of <- rep.int(seq_along(m), m)
  cat_of <- unlist(lapply(m, seq_len))
  dif_items <- sort(as.integer(dif_items))
  S <- matrix(0, P, length(dif_items))
  for (q in seq_along(dif_items)) {
    pos <- which(item_of == dif_items[q])
    S[pos, q] <- cat_of[pos]
  }
  list(m = m, P = P, item_of = item_of, cat_of = cat_of,
       dif_items = dif_items, S = S, stats = stats)
}

#' Sufficient statistics of a complete response matrix per group
#' @param X Integer matrix, items from 0, categories 0..m_j all used.
#' @param focal Logical, focal group membership.
#' @keywords internal
#' @noRd
.dimp_data_stats <- function(X, focal, m) {
  P <- sum(m)
  lapply(list(!focal, focal), function(rows) {
    Xg <- X[rows, , drop = FALSE]
    ctot <- unlist(lapply(seq_along(m), function(j) {
      tabulate(Xg[, j], nbins = m[j])
    }))
    ptot <- tabulate(rowSums(Xg) + 1L, nbins = P + 1L)
    list(ctot = ctot, ptot = ptot)
  })
}

#' Stop when the shift of a split item has no finite estimate
#'
#' Given the raw scores, the shift of item j is estimated from each group's
#' total on the item. A person with raw score r can score from
#' max(0, r - (R - m_j)) to min(m_j, r) on the item. When every person in a
#' group is at the lower (or upper) end, the conditional likelihood keeps
#' rising as the shift grows, and the estimate diverges (with an item that
#' does not vary in a group as the common case).
#' @param labels Group labels (reference, focal), for the message.
#' @keywords internal
#' @noRd
.dimp_check_shifts <- function(X, focal, items, labels) {
  m <- apply(X, 2L, max)
  R <- sum(m)
  r <- rowSums(X)
  problems <- character(0)
  bad <- character(0)
  for (j in items) {
    lo <- pmax(0, r - (R - m[j]))
    hi <- pmin(m[j], r)
    for (gi in 1:2) {
      rows <- if (gi == 1L) !focal else focal
      s <- sum(X[rows, j])
      end <- if (s == sum(lo[rows])) "lowest" else
        if (s == sum(hi[rows])) "highest" else NA_character_
      if (!is.na(end)) {
        bad <- union(bad, colnames(X)[j])
        problems <- c(problems, paste0(
          colnames(X)[j], " (", end, " possible scores in group '",
          labels[gi], "')"
        ))
      }
    }
  }
  if (length(problems)) {
    stop(.category_error(paste0(
      "The DIF of ", paste(problems, collapse = ", "), " cannot be ",
      "estimated: every respondent in that group has the lowest (or ",
      "highest) score on the item that their raw score allows, so the shift ",
      "has no finite estimate. Combine groups, merge categories, or leave ",
      "the item out of `items`."
    ), bad))
  }
  invisible(TRUE)
}

#' @keywords internal
#' @noRd
.dimp_unpack <- function(par, prep) {
  P <- prep$P
  psi <- c(0, par[seq_len(P - 1L)])
  nd <- length(prep$dif_items)
  d <- if (nd) par[P - 1L + seq_len(nd)] else numeric(0)
  psi2 <- if (nd) psi + drop(prep$S %*% d) else psi
  list(psi = psi, d = d, psi_g = list(psi, psi2))
}

#' Negative conditional log-likelihood of the shift model
#' @keywords internal
#' @noRd
.dimp_negll <- function(par, prep) {
  u <- .dimp_unpack(par, prep)
  ll <- 0
  for (k in 1:2) {
    st <- prep$stats[[k]]
    esf <- psychotools::elementary_symmetric_functions(
      split(u$psi_g[[k]], prep$item_of), order = 0L
    )[[1L]]
    keep <- st$ptot > 0
    ll <- ll - sum(st$ctot * u$psi_g[[k]]) -
      sum(st$ptot[keep] * log(esf[keep]))
  }
  if (!is.finite(ll)) return(.Machine$double.xmax)
  -ll
}

#' Analytic gradient of `.dimp_negll()`
#'
#' Per group, with gamma0 the elementary symmetric functions and gamma1
#' their derivatives: `d(-ll)/d(psi_g) = ctot_g - sum_r ptot_gr *
#' gamma1[r, ] / gamma0[r]`.
#' @keywords internal
#' @noRd
.dimp_grad <- function(par, prep) {
  u <- .dimp_unpack(par, prep)
  gpsi <- vector("list", 2L)
  for (k in 1:2) {
    st <- prep$stats[[k]]
    esf <- psychotools::elementary_symmetric_functions(
      split(u$psi_g[[k]], prep$item_of), order = 1L
    )
    keep <- which(st$ptot > 0)
    g0 <- esf[[1L]][keep]
    g1 <- esf[[2L]][keep, , drop = FALSE]
    gpsi[[k]] <- st$ctot - colSums(st$ptot[keep] * g1 / g0)
  }
  g_psi <- (gpsi[[1L]] + gpsi[[2L]])[-1L]
  g_d <- if (length(prep$dif_items)) {
    drop(crossprod(prep$S, gpsi[[2L]]))
  } else {
    numeric(0)
  }
  c(g_psi, g_d)
}

#' Fit the shift model
#' @param start Starting values (length P - 1 + number of shifts), or
#'   `NULL` for zeros.
#' @return List with `par`, `vcov` (`NULL` without Hessian or when it is
#'   singular), `convergence`.
#' @keywords internal
#' @noRd
.dimp_fit <- function(prep, start = NULL, hessian = TRUE, reltol = 1e-10,
                      maxit = 2000L) {
  np <- prep$P - 1L + length(prep$dif_items)
  if (is.null(start)) start <- rep(0, np)
  opt <- stats::optim(start, .dimp_negll, .dimp_grad, prep = prep,
                      method = "BFGS",
                      control = list(reltol = reltol, maxit = maxit))
  vc <- NULL
  if (hessian) {
    H <- stats::optimHess(opt$par, .dimp_negll, .dimp_grad, prep = prep)
    vc <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(vc) && any(!is.finite(diag(vc)) | diag(vc) <= 0)) vc <- NULL
  }
  list(par = opt$par, vcov = vc, convergence = opt$convergence)
}

#' Andrich thresholds of both groups from shift-model parameters
#'
#' Centred so that the reference group's thresholds have mean zero.
#' @keywords internal
#' @noRd
.dimp_par_to_thr <- function(par, prep, item_names) {
  u <- .dimp_unpack(par, prep)
  thr <- .center_thresholds(
    lapply(split(u$psi, prep$item_of), function(p) diff(c(0, p)))
  )
  names(thr) <- item_names
  thr2 <- thr
  for (q in seq_along(prep$dif_items)) {
    j <- prep$dif_items[q]
    thr2[[j]] <- thr2[[j]] + u$d[q]
  }
  list(ref = thr, focal = thr2)
}

#' Set the focal thresholds of the split items to reference + shift
#' @keywords internal
#' @noRd
.dimp_with_shifts <- function(thr, items, d) {
  for (q in seq_along(items)) {
    thr$focal[[items[q]]] <- thr$ref[[items[q]]] + d[q]
  }
  thr
}

# ---------------------------------------------------------------------
# Person scores by raw score
# ---------------------------------------------------------------------

#' Theta range for the grids and solvers
#'
#' -10 to 10, widened when a threshold lies beyond +-5 so that the range
#' reaches at least 5 logits past the most extreme threshold. Items with very
#' low or high thresholds then do not push the score curves, the WLEs or a
#' latent distribution off the grid. The ends are multiples of 5, so the
#' grids stay the same when the delta method perturbs the thresholds
#' slightly.
#' @param ... Threshold lists.
#' @keywords internal
#' @noRd
.dimp_span <- function(...) {
  t <- unlist(list(...))
  c(min(-10, 5 * floor((min(t) - 5) / 5)),
    max(10, 5 * ceiling((max(t) + 5) / 5)))
}

#' WLE for every raw score 0..R, by vectorised bisection on Warm's score
#' function, searched over `.dimp_span()`
#' @keywords internal
#' @noRd
.dimp_wle_by_score <- function(thr_list, lo = .dimp_span(thr_list)[1L],
                               hi = .dimp_span(thr_list)[2L], iter = 40L) {
  r <- 0:sum(lengths(thr_list))
  a <- rep(lo, length(r))
  b <- rep(hi, length(r))
  for (it in seq_len(iter)) {
    mid <- (a + b) / 2
    mo <- .dimp_moments(thr_list, mid)
    s <- r - mo$E + mo$J / (2 * mo$I)
    up <- s > 0
    a[up] <- mid[up]
    b[!up] <- mid[!up]
  }
  (a + b) / 2
}

#' theta with expected raw score r (r strictly between 0 and R)
#' @keywords internal
#' @noRd
.dimp_inv_true_score <- function(thr_list, r,
                                 lo = .dimp_span(thr_list)[1L] - 5,
                                 hi = .dimp_span(thr_list)[2L] + 5,
                                 iter = 40L) {
  a <- rep(lo, length(r))
  b <- rep(hi, length(r))
  for (it in seq_len(iter)) {
    mid <- (a + b) / 2
    up <- .dimp_moments(thr_list, mid)$E < r
    a[up] <- mid[up]
    b[!up] <- mid[!up]
  }
  (a + b) / 2
}

#' True-score equating of focal raw scores onto the reference items,
#' T_ref(T_focal^-1(r)); 0 and R map to themselves
#' @keywords internal
#' @noRd
.dimp_equate <- function(thr_ref, thr_focal) {
  R <- sum(lengths(thr_ref))
  r <- 0:R
  out <- r
  inner <- r > 0 & r < R
  if (any(inner)) {
    th <- .dimp_inv_true_score(thr_focal, r[inner])
    out[inner] <- .dimp_moments(thr_ref, th)$E
  }
  out
}

# ---------------------------------------------------------------------
# Latent mean and SD by MML with the item parameters fixed
# ---------------------------------------------------------------------

#' MML grid over `.dimp_span()` with step 0.05
#' @keywords internal
#' @noRd
.dimp_mml_grid <- function(thr_list = NULL) {
  sp <- if (is.null(thr_list)) c(-10, 10) else .dimp_span(thr_list)
  seq(sp[1L], sp[2L], length.out = round((sp[2L] - sp[1L]) / 0.05) + 1L)
}

#' Log-likelihood of each raw score at each grid node, up to a constant
#' @keywords internal
#' @noRd
.dimp_mml_L <- function(thr_list, grid = .dimp_mml_grid(thr_list)) {
  logZ <- 0
  for (thr in thr_list) {
    m <- length(thr)
    eta <- outer(grid, 0:m) -
      matrix(c(0, cumsum(thr)), length(grid), m + 1L, byrow = TRUE)
    mx <- do.call(pmax, as.data.frame(eta))
    logZ <- logZ + mx + log(rowSums(exp(eta - mx)))
  }
  R <- sum(lengths(thr_list))
  outer(0:R, grid) - matrix(logZ, R + 1L, length(grid), byrow = TRUE)
}

#' Negative marginal log-likelihood per person, par = (mu, log sd)
#' @keywords internal
#' @noRd
.dimp_mml_negll <- function(par, L, pr, grid) {
  lp <- stats::dnorm(grid, par[1L], exp(par[2L]), log = TRUE)
  lmx <- max(lp)
  lp <- lp - (lmx + log(sum(exp(lp - lmx))))
  M <- sweep(L, 2L, lp, "+")
  mx <- do.call(pmax, as.data.frame(M))
  keep <- pr > 0
  -sum(pr[keep] * (mx + log(rowSums(exp(M - mx))))[keep])
}

#' Latent normal moments of one group from its raw-score distribution
#' @param pr Raw-score distribution (proportions).
#' @return List with `par` (mu, log sd), `mu`, `sd`, `convergence`
#'   (optim code), and with `hessian`, `vcov1` (covariance of par for one
#'   person; divide by n; `NULL` when the Hessian is not positive definite).
#' @keywords internal
#' @noRd
.dimp_mml_fit <- function(thr_list, pr, start = c(0, 0), hessian = TRUE) {
  grid <- .dimp_mml_grid(thr_list)
  L <- .dimp_mml_L(thr_list, grid)
  opt <- stats::optim(start, .dimp_mml_negll, L = L, pr = pr, grid = grid,
                      method = "BFGS",
                      control = list(reltol = 1e-15, maxit = 2000L))
  vcov1 <- NULL
  if (hessian) {
    H <- stats::optimHess(opt$par, .dimp_mml_negll, L = L, pr = pr,
                          grid = grid)
    vcov1 <- tryCatch(solve(H), error = function(e) NULL)
    # Not positive definite: no usable SE (callers stop or report NA).
    if (!is.null(vcov1) &&
          (!all(is.finite(vcov1)) || any(diag(vcov1) <= 0))) {
      vcov1 <- NULL
    }
  }
  list(par = opt$par, mu = opt$par[1L], sd = exp(opt$par[2L]),
       vcov1 = vcov1, convergence = opt$convergence)
}

# ---------------------------------------------------------------------
# Group comparison
# ---------------------------------------------------------------------

#' The group difference with three estimators, from raw-score
#' distributions
#'
#' @param thr List (ref, focal) of threshold lists used for scoring.
#' @param pr List (ref, focal) of raw-score distributions (proportions).
#' @param n Group sizes (ref, focal), for the sampling SEs.
#' @param equate Equate the focal sum score onto the reference items.
#' @param lat_start Starting values for the two MML fits.
#' @param sample_var Use sample variances (n - 1) for the WLE and sum score,
#'   as `stats::var()` on the persons. `FALSE` for population values.
#' @return List with `group` (matrix, rows latent / wle / sum, columns
#'   delta, se, smd, sd_pooled, var_ref, var_focal; var_* are per-person
#'   variances, so se = sqrt(var_ref / n_ref + var_focal / n_focal)) and
#'   `lat` (the two MML fits).
#' @keywords internal
#' @noRd
.dimp_compare <- function(thr, pr, n, equate, lat_start = NULL,
                          hessian = TRUE, sample_var = TRUE) {
  if (is.null(lat_start)) lat_start <- list(c(0, 0), c(0, 0))
  R <- length(pr[[1L]]) - 1L
  lat <- Map(.dimp_mml_fit, thr, pr, lat_start,
             MoreArgs = list(hessian = hessian))
  row <- function(m, v, sd_p) {
    m <- unname(m)
    v <- unname(v)
    sd_p <- unname(sd_p)
    c(delta = m[2] - m[1], se = sqrt(v[1] / n[1] + v[2] / n[2]),
      smd = (m[2] - m[1]) / sd_p, sd_pooled = sd_p,
      var_ref = v[1], var_focal = v[2])
  }
  # Latent mean: per-person variance of mu from the MML information.
  vmu <- vapply(lat, function(f) {
    if (is.null(f$vcov1)) NA_real_ else f$vcov1[1L, 1L]
  }, numeric(1L))
  latent <- row(c(lat[[1L]]$mu, lat[[2L]]$mu), vmu,
                sqrt((lat[[1L]]$sd^2 + lat[[2L]]$sd^2) / 2))
  # WLE mean and sum score: means and variances over the raw-score
  # distribution.
  mv <- function(score_list) {
    m <- mapply(function(s, p) sum(p * s), score_list, pr)
    v <- mapply(function(s, p, mm) sum(p * (s - mm)^2), score_list, pr, m)
    list(m = m, v = v)
  }
  w <- mv(lapply(thr, .dimp_wle_by_score))
  s <- mv(list(0:R, if (equate) .dimp_equate(thr$ref, thr$focal) else 0:R))
  # Sample variances (n - 1), as stats::var() on the persons.
  bessel <- if (sample_var) n / pmax(n - 1, 1) else c(1, 1)
  wle <- row(w$m, w$v * bessel, sqrt(mean(w$v * bessel)))
  sum_ <- row(s$m, s$v * bessel, sqrt(mean(s$v * bessel)))
  list(group = rbind(latent = latent, wle = wle, sum = sum_), lat = lat)
}

# ---------------------------------------------------------------------
# Local bias b(theta)
# ---------------------------------------------------------------------
#
# b(theta) = T_ref^-1(T_focal(theta)) - theta: how far a focal person at
# theta lands from a reference person at the same theta when both are
# scored with the reference parameters (the horizontal distance between
# the two test characteristic curves). Negative when the DIF makes the
# items harder for the focal group.

#' Central range of the two-group latent mixture
#' @keywords internal
#' @noRd
.dimp_range <- function(mu, sd, w_focal, q = c(0.025, 0.975)) {
  pmix <- function(x) {
    (1 - w_focal) * stats::pnorm(x, mu[1L], sd[1L]) +
      w_focal * stats::pnorm(x, mu[2L], sd[2L])
  }
  vapply(q, function(p) {
    stats::uniroot(function(x) pmix(x) - p,
                   c(min(-15, mu - 10 * sd), max(15, mu + 10 * sd)),
                   tol = 1e-10)$root
  }, numeric(1L))
}

#' b(theta), SEM(theta) of the reference test and b / SEM
#' @keywords internal
#' @noRd
.dimp_b_at <- function(thr, theta) {
  Tf <- .dimp_moments(thr$focal, theta)$E
  mr <- .dimp_moments(thr$ref, theta)
  b <- .dimp_inv_true_score(thr$ref, Tf) - theta
  sem <- 1 / sqrt(mr$I)
  list(b = b, sem = sem, b_sem = b / sem)
}

#' Where on a grid b(theta) and b / SEM are largest in absolute value
#' @param g Output of `.dimp_b_at()` on `grid`.
#' @keywords internal
#' @noRd
.dimp_local_where <- function(g, grid) {
  c(theta_b = grid[which.max(abs(g$b))],
    theta_b_sem = grid[which.max(abs(g$b_sem))])
}

#' Signed b and b / SEM at fixed points, and at the cut
#'
#' The largest absolute value is reported signed, at the point where it
#' occurs. Holding that point fixed when the shifts are perturbed gives the
#' same first-order derivative as re-maximising (envelope theorem), and the
#' signed value keeps a symmetric interval meaningful when b is near zero,
#' where an interval for abs(b) would extend below zero.
#' @keywords internal
#' @noRd
.dimp_local_at <- function(thr, where, theta_cut = NULL) {
  out <- c(b_max = .dimp_b_at(thr, where[["theta_b"]])$b,
           b_sem_max = .dimp_b_at(thr, where[["theta_b_sem"]])$b_sem)
  if (!is.null(theta_cut)) {
    ct <- .dimp_b_at(thr, theta_cut)
    out <- c(out, b_cut = ct$b, b_sem_cut = ct$b_sem)
  }
  out
}

# ---------------------------------------------------------------------
# Classification at a cut
# ---------------------------------------------------------------------
#
# The rule (which raw scores are classified at or above the cut, per group
# and strategy) is taken from the fit and held fixed. Its rates (share
# classified at or above, false negatives, false positives, as shares of
# the group) are then smooth in the parameters, so the delta method applies
# (dev/dif-impact-design.md, "Step 3b results").

.dimp_th_h <- 0.02
.dimp_th_grid <- function(span = c(-10, 10)) {
  seq(span[1L], span[2L], by = .dimp_th_h)
}

#' What the rates need from one set of split-model thresholds
#'
#' Score probabilities on the grid, the true-status point for a sum-score
#' cut (theta where the reference expected score is cut_raw - 0.5), and,
#' with `rule = TRUE`, what the classification rule is built from: WLE
#' tables per group and equated focal raw scores. These use the same
#' solvers as the keep WLE table and the impact (`.dimp_wle_by_score()`,
#' `.dimp_equate()`), so keep and split classify a raw score alike when
#' their parameters agree. The rule is fixed at the estimates, so the
#' perturbed calls of the delta method skip it, and pass the grid of the
#' estimates (`grid`, `NULL` for one over `.dimp_span()`).
#' @keywords internal
#' @noRd
.dimp_rate_parts <- function(thr, items, cut_raw = NULL, rule = TRUE,
                             grid = NULL) {
  if (is.null(grid)) grid <- .dimp_th_grid(.dimp_span(thr$ref, thr$focal))
  common <- setdiff(seq_along(thr$ref), items)
  R <- sum(lengths(thr$ref))
  SPc <- .dimp_score_probs(thr$ref[common], grid)
  SPr <- SPc
  SPf <- SPc
  for (j in items) {
    SPr <- .dimp_conv(SPr, .dimp_cat_probs(thr$ref[[j]], grid))
    SPf <- .dimp_conv(SPf, .dimp_cat_probs(thr$focal[[j]], grid))
  }
  out <- list(R = R, grid = grid, SP = list(ref = SPr, focal = SPf),
              theta_c_sum = if (is.null(cut_raw)) NA_real_ else
                stats::approx(.dimp_moments(thr$ref, grid)$E, grid,
                              xout = cut_raw - 0.5)$y)
  if (rule) {
    out$wle_ref <- .dimp_wle_by_score(thr$ref)
    out$wle_focal <- .dimp_wle_by_score(thr$focal)
    out$eq <- .dimp_equate(thr$ref, thr$focal)
  }
  out
}

#' Rates of one rule in one group with latent N(mu, sd)
#'
#' The true-status indicator is interpolated across the grid cell that holds
#' the cut, so the rates are smooth in mu and sd.
#' @keywords internal
#' @noRd
.dimp_rates_one <- function(SP, cls, mu, sd, theta_c, grid = NULL) {
  if (is.null(grid)) grid <- .dimp_th_grid()
  w <- stats::dnorm(grid, mu, sd) * .dimp_th_h
  above <- pmin(pmax((grid - theta_c) / .dimp_th_h + 0.5, 0), 1)
  pp <- drop(SP %*% as.numeric(cls))
  c(prev = sum(w * pp), fn = sum(w * above * (1 - pp)),
    fp = sum(w * (1 - above) * pp))
}

#' The classification rule from the fit
#' @param wk Keep WLE table on the split scale (theta cut only).
#' @return List keep / split, each a list ref / focal of logical vectors
#'   over raw scores 0..R.
#' @keywords internal
#' @noRd
.dimp_rule <- function(pt, wk, scale, cut) {
  r <- 0:pt$R
  if (scale == "theta") {
    list(keep = list(ref = wk >= cut, focal = wk >= cut),
         split = list(ref = pt$wle_ref >= cut, focal = pt$wle_focal >= cut))
  } else {
    # An equated focal score counts as at or above the cut when it is
    # nearest to a whole score at or above it (eq >= cut - 0.5), the same
    # point that defines the true status. Comparing eq with the cut itself
    # would move the focal cut by a whole score for any DIF, however small,
    # that makes the items easier for the focal group (eq just below r).
    list(keep = list(ref = r >= cut, focal = r >= cut),
         split = list(ref = r >= cut, focal = pt$eq >= cut - 0.5))
  }
}

#' Lowest raw score classified at or above the cut (R + 1 if none)
#' @keywords internal
#' @noRd
.dimp_rule_cuts <- function(rule) {
  out <- c()
  for (st in c("keep", "split")) for (g in c("ref", "focal")) {
    cl <- rule[[st]][[g]]
    out[paste(st, g, sep = ".")] <- if (any(cl)) which(cl)[1L] - 1L else
      length(cl)
  }
  out
}

#' Rates of a fixed rule: per strategy and group, and keep - split
#' @keywords internal
#' @noRd
.dimp_fixed_rates <- function(pt, rule, mu, sd, theta_c) {
  out <- c()
  for (st in c("keep", "split")) for (gi in 1:2) {
    g <- c("ref", "focal")[gi]
    v <- .dimp_rates_one(pt$SP[[g]], rule[[st]][[g]], mu[gi], sd[gi], theta_c,
                         pt$grid)
    out[paste(st, g, names(v), sep = ".")] <- v
  }
  for (g in c("ref", "focal")) for (m in c("prev", "fn", "fp")) {
    out[paste("diff", g, m, sep = ".")] <-
      out[[paste("keep", g, m, sep = ".")]] -
      out[[paste("split", g, m, sep = ".")]]
  }
  out
}

# ---------------------------------------------------------------------
# Projection to a target population
# ---------------------------------------------------------------------

#' Keep (pooled) thresholds that CML converges to when the data follow the
#' split model with latent N(mu, sd) per group and group sizes n
#' @param start Starting values for the CML fit (`NULL` for zeros).
#' @keywords internal
#' @noRd
.dimp_pop_keep <- function(thr, design, start = NULL) {
  st <- list(.dimp_pop_stats(thr$ref, design$mean[1L], design$sd[1L]),
             .dimp_pop_stats(thr$focal, design$mean[2L], design$sd[2L]))
  stats <- Map(function(s, nn) list(ctot = nn * s$ctot, ptot = nn * s$pr),
               st, design$n)
  prep <- .dimp_prep(lengths(thr$ref), integer(0), stats)
  fk <- .dimp_fit(prep, start = start, hessian = FALSE, reltol = 1e-15,
                  maxit = 10000L)
  list(thr = .dimp_par_to_thr(fk$par, prep, names(thr$ref))$ref,
       pr = list(st[[1L]]$pr, st[[2L]]$pr), par = fk$par,
       convergence = fk$convergence)
}

#' Keep impact on Delta in the target population, and the SE of Delta
#' (split) in a study of the target size
#'
#' @param start `NULL`, or the `start` element of an earlier call: warm
#'   starts for the population CML and the four MML fits, used by the
#'   perturbed calls of the delta method.
#' @param se `FALSE` skips the MML Hessians, so `se_target` has `NA` for
#'   the latent mean. The perturbed calls only need the impact.
#' @return List with `impact`, `se_target`, `thr_keep`, `start` and
#'   `convergence` (`TRUE` when every fit converged).
#' @keywords internal
#' @noRd
.dimp_project <- function(thr, target, start = NULL, se = TRUE) {
  pk <- .dimp_pop_keep(thr, target, start = start$keep)
  ck <- .dimp_compare(list(ref = pk$thr, focal = pk$thr), pk$pr, target$n,
                      equate = FALSE, lat_start = start$lat_keep,
                      hessian = FALSE, sample_var = FALSE)
  cs <- .dimp_compare(thr, pk$pr, target$n, equate = TRUE,
                      lat_start = start$lat_split, hessian = se,
                      sample_var = FALSE)
  lat_par <- function(cc) lapply(cc$lat, `[[`, "par")
  conv <- c(pk$convergence,
            vapply(c(ck$lat, cs$lat), `[[`, integer(1L), "convergence"))
  list(impact = ck$group[, "delta"] - cs$group[, "delta"],
       se_target = cs$group[, "se"], thr_keep = pk$thr,
       start = list(keep = pk$par, lat_keep = lat_par(ck),
                    lat_split = lat_par(cs)),
       convergence = all(conv == 0L))
}
