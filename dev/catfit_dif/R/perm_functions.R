# Permutation bands for group-split PAV calibration curves (prototype)
#
# Uses easyRasch2 internals (.fit_cml_thresholds, .cond_cat_probs), so source
# after devtools::load_all() from the package root. See ../PLAN.md, section
# "Validating the permutation bands".
#
# Idea. The predicted probability of event X_i >= k depends only on the total
# score r and the pooled CML item parameters. Group labels are permuted within
# total scores, which keeps each group's count at every score, so each group's
# curve has the same points (scores) in every permutation. All curves are
# computed from per-score sums, aggregated by group.

# -- Data generation ----------------------------------------------------------

#' PCM responses. theta: abilities. delta_list: threshold vectors per item.
sim_pcm <- function(theta, delta_list) {
  sapply(delta_list, function(d) {
    m <- length(d)
    lp <- sapply(0:m, function(j) j * theta - sum(d[seq_len(j)]))
    lp <- matrix(lp, nrow = length(theta))
    p <- exp(lp - apply(lp, 1, max))
    p <- p / rowSums(p)
    rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
  })
}

#' Thresholds for k items with four categories: locations from -1.5 to 1.5,
#' thresholds at location -1, 0, +1.
design_thresholds <- function(k) {
  loc <- seq(-1.5, 1.5, length.out = k)
  lapply(loc, function(l) l + c(-1, 0, 1))
}

#' Data for one replication, no DIF. Group means are spread symmetrically
#' around 0 with range `impact`. The RNG state is set by the caller.
#' props: group proportions. Returns list(X, grp).
gen_null <- function(n, k, impact, props) {
  G <- length(props)
  ng <- round(n * props)
  ng[G] <- n - sum(ng[-G])
  grp <- rep(seq_len(G), ng)
  mu <- if (G == 1L) 0 else impact * ((seq_len(G) - 1) / (G - 1) - 0.5)
  theta <- stats::rnorm(n, mu[grp], 1)
  X <- sim_pcm(theta, design_thresholds(k))
  colnames(X) <- paste0("I", seq_len(k))
  list(X = X, grp = grp)
}

# -- Preparation --------------------------------------------------------------

#' Pooled CML fit, extreme scores dropped, event indicators and predictions.
#' Returns NULL when the CML fit fails (an empty category).
prep_dif <- function(X, grp) {
  X <- as.matrix(X)
  storage.mode(X) <- "integer"
  thr <- tryCatch(.fit_cml_thresholds(X), error = function(e) NULL)
  if (is.null(thr)) return(NULL)
  P <- .cond_cat_probs(thr)
  m <- lengths(thr)
  K <- sum(m)
  r <- rowSums(X)
  keep <- r > 0 & r < K
  X <- X[keep, , drop = FALSE]
  r <- r[keep]
  grp <- as.integer(factor(grp[keep]))
  ev_item <- rep(seq_along(m), m)
  ev_k <- sequence(m)
  # Event indicators (n x E) and predictions per score (K + 1 x E).
  Y <- vapply(seq_along(ev_item), function(e) {
    as.integer(X[, ev_item[e]] >= ev_k[e])
  }, integer(nrow(X)))
  Pev <- vapply(seq_along(ev_item), function(e) {
    Pi <- P[[ev_item[e]]]
    rowSums(Pi[, (ev_k[e] + 1L):ncol(Pi), drop = FALSE])
  }, numeric(K + 1L))
  sc <- sort(unique(r))
  list(X = X, Y = Y, r = r, grp = grp, G = max(grp), sc = sc,
       s_idx = match(r, sc), Pev = Pev[sc + 1L, , drop = FALSE],
       ev_item = ev_item, ev_k = ev_k, thr = thr)
}

# -- Aggregation and PAV ------------------------------------------------------

#' Counts per group x score (G x S) and event sums per group x score
#' ((G * S) x E, row (g - 1) * S + s).
agg_counts <- function(Y, s_idx, grp, G, S) {
  key <- (grp - 1L) * S + s_idx
  n <- tabulate(key, G * S)
  sums <- matrix(0, G * S, ncol(Y))
  rs <- rowsum(Y, key, reorder = TRUE)
  sums[as.integer(rownames(rs)), ] <- rs
  list(n = matrix(n, G, S, byrow = TRUE), sums = sums)
}

#' Weighted PAV of values v (ordered by the predictor) with weights w > 0.
#' Returns the fitted value at each point.
pav_w <- function(v, w) {
  J <- length(v)
  bv <- numeric(J); bw <- numeric(J); bl <- integer(J)
  t <- 0L
  for (j in seq_len(J)) {
    t <- t + 1L
    bv[t] <- v[j]; bw[t] <- w[j]; bl[t] <- 1L
    while (t > 1L && bv[t - 1L] > bv[t]) {
      bv[t - 1L] <- (bv[t - 1L] * bw[t - 1L] + bv[t] * bw[t]) /
        (bw[t - 1L] + bw[t])
      bw[t - 1L] <- bw[t - 1L] + bw[t]
      bl[t - 1L] <- bl[t - 1L] + bl[t]
      t <- t - 1L
    }
  }
  rep(bv[seq_len(t)], bl[seq_len(t)])
}

#' Curve layout: one curve per group x event, points at the scores where the
#' group has respondents. Fixed under within-score permutation.
curve_layout <- function(pr) {
  S <- length(pr$sc)
  a <- agg_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)
  E <- ncol(pr$Y)
  do.call(rbind, lapply(seq_len(E), function(e) {
    do.call(rbind, lapply(seq_len(pr$G), function(g) {
      s <- which(a$n[g, ] > 0)
      data.frame(curve = (e - 1L) * pr$G + g, event = e,
                 item = pr$ev_item[e], k = pr$ev_k[e], group = g,
                 s = s, score = pr$sc[s], n = a$n[g, s],
                 predicted = pr$Pev[s, e])
    }))
  }))
}

#' PAV values at every point of the layout, for one group labelling.
curve_values <- function(pr, grp, lay) {
  S <- length(pr$sc)
  a <- agg_counts(pr$Y, pr$s_idx, grp, pr$G, S)
  out <- numeric(nrow(lay))
  for (cv in unique(lay$curve)) {
    i <- which(lay$curve == cv)
    g <- lay$group[i[1L]]
    e <- lay$event[i[1L]]
    s <- lay$s[i]
    w <- a$n[g, s]
    v <- a$sums[(g - 1L) * S + s, e] / w
    out[i] <- pav_w(v, w)
  }
  out
}

#' Raw per-score shares (no PAV), for the mirror check.
curve_raw <- function(pr, grp, lay) {
  S <- length(pr$sc)
  a <- agg_counts(pr$Y, pr$s_idx, grp, pr$G, S)
  a$sums[cbind((lay$group - 1L) * S + lay$s, lay$event)] /
    a$n[cbind(lay$group, lay$s)]
}

#' Permute group labels within total scores.
permute_within <- function(grp, s_idx) {
  ord <- order(s_idx)
  out <- grp
  out[ord] <- grp[ord][order(s_idx[ord], stats::runif(length(grp)))]
  out
}

# -- Bands --------------------------------------------------------------------

#' Observed curves, B permutations, pointwise and simultaneous p-values.
#' Returns list(points, curves) data.frames, plus the permutation matrix
#' when keep_perm = TRUE.
perm_bands <- function(pr, B = 999L, keep_perm = FALSE) {
  lay <- curve_layout(pr)
  obs <- curve_values(pr, pr$grp, lay)
  Pm <- vapply(seq_len(B), function(b) {
    curve_values(pr, permute_within(pr$grp, pr$s_idx), lay)
  }, numeric(nrow(lay)))
  Pm <- matrix(Pm, nrow = nrow(lay))
  # Pointwise two-sided permutation p (observed counted in).
  lo <- rowSums(Pm <= obs + 1e-12)
  hi <- rowSums(Pm >= obs - 1e-12)
  p_point <- pmin(1, 2 * pmin(lo + 1, hi + 1) / (B + 1))
  # Simultaneous per curve: max |z| over the curve's points.
  mu <- rowMeans(Pm)
  # Centred SD: the mean-of-squares form leaves rounding noise (~1e-8)
  # when every permutation gives the same value.
  sdv <- sqrt(rowMeans((Pm - mu)^2))
  # Points that never move (sd 0, e.g. certain events) get z = 0.
  z_fun <- function(x) {
    z <- abs(x - mu) / sdv
    z[sdv <= 1e-12] <- 0
    z
  }
  z_obs <- z_fun(obs)
  Z <- z_fun(Pm)
  curves <- unique(lay[, c("curve", "event", "item", "k", "group")])
  ci <- split(seq_len(nrow(lay)), lay$curve)[as.character(curves$curve)]
  T_obs <- vapply(ci, function(i) max(z_obs[i]), numeric(1))
  T_perm <- t(vapply(ci, function(i) {
    apply(Z[i, , drop = FALSE], 2L, max)
  }, numeric(B)))
  T_perm <- matrix(T_perm, nrow = length(ci))
  curves$p_sim <- (1 + rowSums(T_perm >= T_obs - 1e-12)) / (B + 1)
  curves$min_p_point <- vapply(ci, function(i) min(p_point[i]), numeric(1))
  curves$n_points <- lengths(ci)
  points <- cbind(lay, observed = obs, p_point = p_point)
  out <- list(points = points, curves = curves)
  if (keep_perm) out$perm <- Pm
  out
}

# -- Generic bands (any value function over a fixed layout) -------------------
#
# lay: data.frame with at least `curve`. value_fun(grp) returns one value per
# row of lay (NA where a point is empty under that labelling). Same p-value
# definitions as perm_bands(). Points that are NA in the observed data must be
# removed from lay by the caller.

perm_bands_generic <- function(pr, lay, value_fun, B = 999L) {
  obs <- value_fun(pr$grp)
  Pm <- vapply(seq_len(B), function(b) {
    value_fun(permute_within(pr$grp, pr$s_idx))
  }, numeric(nrow(lay)))
  Pm <- matrix(Pm, nrow = nrow(lay))
  ok <- !is.na(Pm)
  nb <- rowSums(ok)
  lo <- rowSums(ok & Pm <= obs + 1e-12)
  hi <- rowSums(ok & Pm >= obs - 1e-12)
  p_point <- pmin(1, 2 * pmin(lo + 1, hi + 1) / (nb + 1))
  P0 <- Pm
  P0[!ok] <- 0
  mu <- rowSums(P0) / pmax(nb, 1)
  D0 <- (Pm - mu)^2
  D0[!ok] <- 0
  sdv <- sqrt(rowSums(D0) / pmax(nb, 1))
  z_fun <- function(x) {
    z <- abs(x - mu) / sdv
    z[sdv <= 1e-12 | is.na(z)] <- 0
    z
  }
  z_obs <- z_fun(obs)
  Z <- z_fun(Pm)
  ci <- split(seq_len(nrow(lay)), lay$curve)
  T_obs <- vapply(ci, function(i) max(z_obs[i]), numeric(1))
  T_perm <- matrix(t(vapply(ci, function(i) {
    apply(Z[i, , drop = FALSE], 2L, max)
  }, numeric(B))), nrow = length(ci))
  first <- vapply(ci, `[`, integer(1), 1L)
  curves <- lay[first, intersect(c("curve", "item", "k", "pair", "group"),
                                 names(lay)), drop = FALSE]
  curves$p_sim <- (1 + rowSums(T_perm >= T_obs - 1e-12)) / (B + 1)
  curves$min_p_point <- vapply(ci, function(i) min(p_point[i]), numeric(1))
  curves$n_points <- lengths(ci)
  # A point "varies" when its permutation SD is positive.
  curves$n_points_var <- vapply(ci, function(i) sum(sdv[i] > 1e-12),
                                numeric(1))
  points <- cbind(lay, observed = obs, p_point = p_point,
                  varies = sdv > 1e-12)
  list(points = points, curves = curves)
}

#' Calibration through the generic code (regression check against perm_bands).
calib_setup <- function(pr) {
  lay <- curve_layout(pr)
  list(lay = lay, fun = function(grp) curve_values(pr, grp, lay))
}

#' Adjacent-category shares per group and score.
adjacent_setup <- function(pr) {
  S <- length(pr$sc)
  m <- lengths(pr$thr)
  # Category counts per item, one column per (item, category).
  cat_item <- rep(seq_along(m), m + 1L)
  cat_c <- sequence(m + 1L) - 1L
  Cm <- vapply(seq_along(cat_item), function(j) {
    as.integer(pr$X[, cat_item[j]] == cat_c[j])
  }, integer(nrow(pr$X)))
  col_of <- function(i, c) which(cat_item == i & cat_c == c)
  counts <- function(grp) {
    key <- (grp - 1L) * S + pr$s_idx
    out <- matrix(0, pr$G * S, ncol(Cm))
    rs <- rowsum(Cm, key, reorder = TRUE)
    out[as.integer(rownames(rs)), ] <- rs
    out
  }
  n0 <- agg_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n
  lay <- do.call(rbind, lapply(seq_along(m), function(i) {
    do.call(rbind, lapply(seq_len(m[i]), function(j) {
      do.call(rbind, lapply(seq_len(pr$G), function(g) {
        s <- which(n0[g, ] > 0)
        data.frame(item = i, pair = j, group = g, s = s, score = pr$sc[s],
                   lo_col = col_of(i, j - 1L), hi_col = col_of(i, j))
      }))
    }))
  }))
  lay$curve <- as.integer(factor(paste(lay$item, lay$pair, lay$group),
                                 levels = unique(paste(lay$item, lay$pair,
                                                       lay$group))))
  fun_all <- function(grp) {
    Cg <- counts(grp)
    row <- (lay$group - 1L) * S + lay$s
    hi <- Cg[cbind(row, lay$hi_col)]
    tot <- hi + Cg[cbind(row, lay$lo_col)]
    ifelse(tot > 0, hi / tot, NA_real_)
  }
  # Drop points empty in the observed data.
  keep <- !is.na(fun_all(pr$grp))
  lay <- lay[keep, ]
  lay$curve <- as.integer(factor(lay$curve, levels = unique(lay$curve)))
  list(lay = lay, fun = function(grp) {
    Cg <- counts(grp)
    row <- (lay$group - 1L) * S + lay$s
    hi <- Cg[cbind(row, lay$hi_col)]
    tot <- hi + Cg[cbind(row, lay$lo_col)]
    ifelse(tot > 0, hi / tot, NA_real_)
  })
}

#' Standardised category residuals per group, score group and category.
#' sg: integer score group per respondent (from the package's quantile
#' grouping into `class_intervals` groups).
residual_setup <- function(pr, class_intervals = 5L) {
  bin <- .cicc_class_bins(pr$r, "quantile", class_intervals, NULL)
  sg <- as.integer(droplevels(bin))
  Q <- max(sg)
  S <- length(pr$sc)
  sg_of_score <- vapply(seq_len(S), function(s) sg[pr$s_idx == s][1L],
                        integer(1))
  m <- lengths(pr$thr)
  K <- sum(m)
  P <- .cond_cat_probs(pr$thr)
  n0 <- agg_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n  # G x S, fixed
  lay <- list()
  EV <- list()
  for (i in seq_along(m)) {
    Pi <- P[[i]][pr$sc + 1L, , drop = FALSE]          # S x (m_i + 1)
    for (g in seq_len(pr$G)) {
      E <- rowsum(n0[g, ] * Pi, sg_of_score, reorder = TRUE)
      V <- rowsum(n0[g, ] * Pi * (1 - Pi), sg_of_score, reorder = TRUE)
      q <- as.integer(rownames(E))
      lay[[length(lay) + 1L]] <- data.frame(
        item = i, group = g, q = rep(q, times = m[i] + 1L),
        category = rep(0:m[i], each = length(q)),
        E = as.vector(E), V = as.vector(V))
    }
  }
  lay <- do.call(rbind, lay)
  lay <- lay[lay$V > 0, ]
  lay$curve <- as.integer(factor(paste(lay$item, lay$group),
                                 levels = unique(paste(lay$item, lay$group))))
  # Flat index of each layout row in the concatenated per-item count tables.
  sizes <- pr$G * Q * (m + 1L)
  offset <- c(0L, cumsum(sizes))[seq_along(m)]
  flat <- offset[lay$item] +
    ((lay$group - 1L) * Q + lay$q - 1L) * (m[lay$item] + 1L) +
    lay$category + 1L
  fun <- function(grp) {
    key <- (grp - 1L) * Q + sg                       # group x score group
    tabs <- unlist(lapply(seq_along(m), function(i) {
      tabulate((key - 1L) * (m[i] + 1L) + pr$X[, i] + 1L, nbins = sizes[i])
    }))
    O <- tabs[flat]
    (O - lay$E) / sqrt(lay$V)
  }
  list(lay = lay, fun = fun)
}
