# Category-by-score table under the PCM, conditional on the total score.
#
# Null 2 from the 2026-11-16 discussion: under the partial credit model,
# P(X_i = j | total score r) depends only on the item parameters, so the
# item-category by score table has model-expected counts without any theta.
# This is the table behind Glas's (1991, Cito report 91-2) CML test R1c for
# polytomous items and its scaled deviates d_rij. The new part here is the set
# of focused statistics (rest-score mean and SD within each category, and the
# category's count in the score tails) and a parametric bootstrap that samples
# response patterns exactly from P(x | r), so the observed score distribution
# is kept and no theta is drawn.

library(psychotools)

# Generating-polynomial coefficients for one item: exp(-cumsum(delta)),
# for j = 0..m. delta are the PCM thresholds (Andrich/"mode" type).
.eps <- function(delta) exp(-c(0, cumsum(delta)))

# Elementary symmetric functions for a set of items (convolution of the
# items' generating polynomials). Index 1 is score 0.
.gamma <- function(eps_list) {
  g <- 1
  for (e in eps_list) g <- stats::convolve(g, rev(e), type = "open")
  pmax(g, 0)
}

# P(X_i = j | r) for item i, as a (K+1) x (m_i+1) matrix, rows r = 0..K.
cond_cat_prob <- function(delta_list, i) {
  eps <- lapply(delta_list, .eps)
  g_all <- .gamma(eps)
  g_wo <- .gamma(eps[-i])
  K <- length(g_all) - 1L
  m <- length(eps[[i]]) - 1L
  P <- matrix(0, K + 1L, m + 1L)
  for (r in 0:K) for (j in 0:m) {
    s <- r - j
    if (s >= 0 && s < length(g_wo)) P[r + 1L, j + 1L] <- eps[[i]][j + 1L] * g_wo[s + 1L] / g_all[r + 1L]
  }
  P
}

# Sample response patterns exactly from P(x | r) given item thresholds.
# totals: integer vector of total scores. Returns an n x k matrix.
sample_given_totals <- function(delta_list, totals) {
  eps <- lapply(delta_list, .eps)
  k <- length(eps)
  suffix_g <- vector("list", k + 1L)
  suffix_g[[k + 1L]] <- 1
  for (i in k:1) suffix_g[[i]] <- .gamma(eps[i:k])
  n <- length(totals)
  out <- matrix(0L, n, k)
  s <- totals
  for (i in seq_len(k)) {
    m <- length(eps[[i]]) - 1L
    g_next <- suffix_g[[i + 1L]]
    g_here <- suffix_g[[i]]
    pr <- sapply(0:m, function(j) {
      idx <- s - j
      ok <- idx >= 0 & idx < length(g_next)
      v <- numeric(n)
      v[ok] <- eps[[i]][j + 1L] * g_next[idx[ok] + 1L]
      v
    })
    pr <- matrix(pr, nrow = n)
    pr <- pr / rowSums(pr)
    u <- stats::runif(n)
    x <- rowSums(u > t(apply(pr, 1, cumsum))) # 0..m
    out[, i] <- x
    s <- s - x
  }
  out
}

# CML fit, returning a list of threshold vectors (one per item).
fit_thresholds <- function(X) {
  m <- psychotools::pcmodel(X)
  lapply(psychotools::threshpar(m, type = "mode"), unname)
}

# Statistics for one focus item.
# Returns a data.frame with one row per category j:
#   n_obs, n_exp, mean_obs/mean_exp (rest score), sd_obs/sd_exp,
#   tail_obs/tail_exp (count in the score tails), Q (sum of squared
#   Glas deviates over score levels).
item_stats <- function(X, delta_list, i, tail_q = c(0.1, 0.9)) {
  K <- sum(vapply(delta_list, length, 1L))
  r <- rowSums(X)
  keep <- r > 0 & r < K
  X <- X[keep, , drop = FALSE]
  r <- r[keep]
  P <- cond_cat_prob(delta_list, i)
  m <- ncol(P) - 1L
  n_r <- tabulate(r + 1L, nbins = K + 1L)
  O <- table(factor(r, levels = 0:K), factor(X[, i], levels = 0:m))
  O <- matrix(O, K + 1L)
  E <- n_r * P
  V <- n_r * P * (1 - P)
  lo <- stats::quantile(r, tail_q[1], type = 1)
  hi <- stats::quantile(r, tail_q[2], type = 1)
  tail_rows <- (0:K) <= lo | (0:K) >= hi
  res <- lapply(0:m, function(j) {
    rest <- (0:K) - j
    ok <- rest >= 0 & n_r > 0
    o <- O[ok, j + 1L]; e <- E[ok, j + 1L]; v <- V[ok, j + 1L]; rr <- rest[ok]
    wmean <- function(w) sum(w * rr) / sum(w)
    wsd <- function(w) { mu <- wmean(w); sqrt(sum(w * (rr - mu)^2) / sum(w)) }
    d <- ifelse(v > 0, (o - e) / sqrt(v), 0)
    data.frame(
      category = j, n_obs = sum(o), n_exp = sum(e),
      mean_obs = wmean(o), mean_exp = wmean(e),
      sd_obs = wsd(o), sd_exp = wsd(e),
      tail_obs = sum(O[ok & tail_rows, j + 1L]),
      tail_exp = sum(E[ok & tail_rows, j + 1L]),
      Q = sum(d^2)
    )
  })
  do.call(rbind, res)
}

# Statistics compared across the bootstrap: observed minus expected.
.boot_vec <- function(s) {
  c(setNames(s$mean_obs - s$mean_exp, paste0("mean_", s$category)),
    setNames(s$sd_obs - s$sd_exp, paste0("sd_", s$category)),
    setNames(s$tail_obs - s$tail_exp, paste0("tail_", s$category)),
    setNames(s$Q, paste0("Q_", s$category)))
}

# Parametric bootstrap: sample patterns given the observed totals and the
# fitted thresholds, refit, recompute. Two-sided p for the signed
# differences, upper-tail p for Q.
category_test <- function(X, i, B = 200, tail_q = c(0.1, 0.9)) {
  delta <- fit_thresholds(X)
  K <- sum(vapply(delta, length, 1L))
  r <- rowSums(X)
  totals <- r[r > 0 & r < K]
  obs <- item_stats(X, delta, i, tail_q)
  t_obs <- .boot_vec(obs)
  t_boot <- t(replicate(B, {
    Xb <- sample_given_totals(delta, totals)
    db <- tryCatch(fit_thresholds(Xb), error = function(e) NULL)
    if (is.null(db) || any(vapply(db, length, 1L) != vapply(delta, length, 1L))) {
      rep(NA_real_, length(t_obs))
    } else .boot_vec(item_stats(Xb, db, i, tail_q))
  }))
  t_boot <- t_boot[stats::complete.cases(t_boot), , drop = FALSE]
  is_q <- startsWith(names(t_obs), "Q_")
  p <- vapply(seq_along(t_obs), function(k) {
    tb <- t_boot[, k]
    if (is_q[k]) (1 + sum(tb >= t_obs[k])) / (1 + length(tb))
    else min(1, 2 * min((1 + sum(tb >= t_obs[k])), (1 + sum(tb <= t_obs[k]))) / (1 + length(tb)))
  }, 0)
  list(stats = obs, thresholds = delta[[i]],
       p = setNames(p, names(t_obs)), B_used = nrow(t_boot))
}

# ── Data generation ─────────────────────────────────────────────

sim_pcm <- function(theta, delta_list) {
  sapply(delta_list, function(d) {
    m <- length(d)
    lp <- sapply(0:m, function(j) j * theta - sum(d[seq_len(j)]))
    lp <- matrix(lp, nrow = length(theta))
    p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
    rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
  })
}

# Condition B: categories `pair` of item i made interchangeable. Generate
# with the pair merged (one threshold removed), then split at random.
make_exchangeable <- function(theta, delta_list, i, pair = c(1, 2)) {
  d <- delta_list
  d[[i]] <- d[[i]][-pair[2]]
  X <- sim_pcm(theta, d)
  x <- X[, i]
  hi <- x > pair[1]
  x[hi] <- x[hi] + 1L
  merged <- X[, i] == pair[1]
  x[merged] <- pair[1] + stats::rbinom(sum(merged), 1, 0.5)
  X[, i] <- x
  X
}

# Condition C: with probability pi_dk a respondent answers item i in the
# middle category whatever their theta ("don't know" use).
make_dontknow <- function(theta, delta_list, i, pi_dk = 0.25) {
  X <- sim_pcm(theta, delta_list)
  mid <- length(delta_list[[i]]) %/% 2
  dk <- stats::runif(length(theta)) < pi_dk
  X[dk, i] <- mid
  X
}
