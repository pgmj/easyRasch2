# PAV calibration of the cumulative events X_i >= k with CML-conditional
# forecasts, and a test of it (2026-11-18). Built on R/functions.R.
#
# Forecast for person v, item i, event k: p = P(X_i >= k | r_v), from the CML
# item parameters, so no theta enters. Statistic: the miscalibration
# component of the CORP decomposition of the Brier score (Dimitriadis,
# Gneiting & Jordan 2021),
#   MCB = mean (p - y)^2 - mean (cep - y)^2,
# where cep is the PAV recalibration of y on p. MCB >= 0 and is 0 when the
# observed events are already isotonic in p and equal to it.
#
# Null distribution: the exact conditional bootstrap used for the category
# statistics (patterns sampled given the observed totals, CML thresholds
# refitted in every replicate), so the in-sample fit of the forecast is in
# the reference distribution. Within each item, the m events are combined
# with a single-step min-p adjustment (Westfall & Young 1993) over the same
# bootstrap replicates.

# Weighted PAV on y ordered by p, ties in p pooled first. Returns the
# conditional event probability (CEP) per unique p, with weights.
.pav <- function(p, y) {
  u <- sort(unique(p))
  idx <- match(p, u)
  w <- tabulate(idx, length(u))
  v <- as.numeric(rowsum(y, idx, reorder = TRUE)[, 1] / w)
  lv <- v; lw <- w; blocks <- as.list(seq_along(u)); i <- 1L
  while (i < length(lv)) {
    if (lv[i] > lv[i + 1L]) {
      lv[i] <- (lv[i] * lw[i] + lv[i + 1L] * lw[i + 1L]) / (lw[i] + lw[i + 1L])
      lw[i] <- lw[i] + lw[i + 1L]
      blocks[[i]] <- c(blocks[[i]], blocks[[i + 1L]])
      lv <- lv[-(i + 1L)]; lw <- lw[-(i + 1L)]; blocks <- blocks[-(i + 1L)]
      if (i > 1L) i <- i - 1L
    } else i <- i + 1L
  }
  cep <- numeric(length(u))
  for (b in seq_along(blocks)) cep[blocks[[b]]] <- lv[b]
  data.frame(p = u, cep = cep, w = w)
}

# Forecasts P(X_i >= k | r) for one item, as a (K+1) x m matrix (rows r).
event_probs <- function(delta_list, i) {
  P <- cond_cat_prob(delta_list, i)
  m <- ncol(P) - 1L
  vapply(1:m, function(k) rowSums(P[, (k + 1L):(m + 1L), drop = FALSE]),
         numeric(nrow(P)))
}

# MCB for every item and event: a k x m matrix (rows items).
# Persons with an extreme total (0 or K) are left out, as in item_stats().
mcb_matrix <- function(X, delta_list) {
  K <- sum(lengths(delta_list))
  r <- rowSums(X)
  keep <- r > 0 & r < K
  X <- X[keep, , drop = FALSE]; r <- r[keep]
  t(vapply(seq_along(delta_list), function(i) {
    E <- event_probs(delta_list, i)
    vapply(seq_len(ncol(E)), function(k) {
      p <- E[r + 1L, k]
      y <- as.numeric(X[, i] >= k)
      pv <- .pav(p, y)
      cep <- pv$cep[match(p, pv$p)]
      mean((p - y)^2) - mean((cep - y)^2)
    }, 0)
  }, numeric(length(delta_list[[1]]))))
}

# Single-step min-p over the events of one item, observed and replicates
# ranked in one pool (2026-11-18, second version).
#   T: m x (Bu + 1) matrix, column 1 observed, larger = worse fit.
# Each member's p for event j is the share of the pool at or above its value.
# The item p is the share of the pool whose smallest p is at or below the
# observed one. The first version divided the replicates' p by Bu and the
# observed p by Bu + 1, so replicates tied with the observed count were not
# counted and the item p was too small (.017 at .01 in a toy check).
.minp_pooled <- function(T) {
  N <- ncol(T)
  P <- t(apply(T, 1, function(t) (N - rank(t, ties.method = "min") + 1) / N))
  P <- matrix(P, nrow = nrow(T))
  mp <- apply(P, 2, min)
  list(p_raw = P[, 1], p_item = mean(mp <= mp[1]))
}

# One dataset: observed MCB, bootstrap with refit, raw p per item-event,
# item-level p (pooled single-step min-p over the events of the item), Holm
# over items. The bootstrap MCB array is returned (rounded to 7 significant
# digits) so other combinations can be computed without a rerun.
analyse_calibration <- function(X, B) {
  k <- ncol(X); m <- 4L
  if (any(vapply(seq_len(k), function(i) length(unique(X[, i])) < m + 1, TRUE)))
    return(list(skipped = "empty category in data"))
  delta <- tryCatch(fit_thresholds(X), error = function(e) NULL)
  if (is.null(delta)) return(list(skipped = "fit failed"))
  K <- k * m
  r <- rowSums(X)
  totals <- r[r > 0 & r < K]
  obs <- mcb_matrix(X, delta)
  boot <- replicate(B, {
    Xb <- sample_given_totals(delta, totals)
    db <- tryCatch(fit_thresholds(Xb), error = function(e) NULL)
    if (is.null(db) || any(lengths(db) != m)) matrix(NA_real_, k, m)
    else mcb_matrix(Xb, db)
  }, simplify = "array")                       # k x m x B
  ok <- apply(boot, 3, function(a) all(is.finite(a)))
  boot <- boot[, , ok, drop = FALSE]
  Bu <- dim(boot)[3]
  if (Bu == 0) return(list(skipped = NA_character_, B_used = 0L, mcb = obs))
  p_raw <- matrix(NA_real_, k, m)
  p_item <- numeric(k)
  for (i in seq_len(k)) {
    res <- .minp_pooled(cbind(obs[i, ], matrix(boot[i, , ], nrow = m)))
    p_raw[i, ] <- res$p_raw
    p_item[i] <- res$p_item
  }
  list(skipped = NA_character_, B_used = Bu, mcb = obs,
       p_raw = p_raw, p_item = p_item,
       p_item_holm = stats::p.adjust(p_item, "holm"),
       boot_mcb = signif(boot, 7))
}

# Studentised max-T over the events of each item (2026-11-19), the
# combination RMitemCatFit() uses. Each event's MCB is z-scored by the mean
# and SD of the pool (observed + Bu replicates), the item statistic is the
# largest z of its events, and the item p is the share of the pool at or
# above the observed statistic. Across items: Westfall-Young step-down on the
# item statistics, floored at the item p. obs: k x m, boot: k x m x Bu.
.stepdown_maxT <- function(s_obs, s_sim) {     # larger = more extreme
  B <- nrow(s_sim)
  ord <- order(s_obs, decreasing = TRUE)
  run <- rep(-Inf, B)
  p <- numeric(length(ord))
  for (j in rev(seq_along(ord))) {
    run <- pmax(run, s_sim[, ord[j]])
    p[j] <- (1 + sum(run >= s_obs[ord[j]])) / (B + 1)
  }
  out <- numeric(length(ord))
  out[ord] <- pmin(cummax(p), 1)
  out
}

.maxT_items <- function(obs, boot) {
  k <- nrow(obs); m <- ncol(obs); N <- dim(boot)[3] + 1L
  pool <- array(c(obs, boot), c(k, m, N))
  Tm <- matrix(NA_real_, N, k)
  for (i in seq_len(k)) {
    E <- matrix(pool[i, , ], nrow = m)
    Z <- t(apply(E, 1, function(t) {
      s <- stats::sd(t)
      if (!is.finite(s) || s < 1e-12) rep(0, length(t)) else (t - mean(t)) / s
    }))
    Tm[, i] <- apply(matrix(Z, nrow = m), 2, max)
  }
  p_item <- colMeans(sweep(Tm, 2L, Tm[1L, ], `>=`))
  padj <- pmax(.stepdown_maxT(Tm[1L, ], Tm[-1L, , drop = FALSE]), p_item)
  list(p_item = p_item, padj = padj, Tm = Tm)
}

# One null dataset for the large-B check: min-p and max-T item p, max-T
# step-down. Keeps the item statistics (N x k), not the full MCB array.
analyse_calibration_maxT <- function(X, B) {
  k <- ncol(X); m <- 4L
  if (any(vapply(seq_len(k), function(i) length(unique(X[, i])) < m + 1, TRUE)))
    return(list(skipped = "empty category in data"))
  delta <- tryCatch(fit_thresholds(X), error = function(e) NULL)
  if (is.null(delta)) return(list(skipped = "fit failed"))
  r <- rowSums(X)
  totals <- r[r > 0 & r < k * m]
  obs <- mcb_matrix(X, delta)
  boot <- replicate(B, {
    Xb <- sample_given_totals(delta, totals)
    db <- tryCatch(fit_thresholds(Xb), error = function(e) NULL)
    if (is.null(db) || any(lengths(db) != m)) matrix(NA_real_, k, m)
    else mcb_matrix(Xb, db)
  }, simplify = "array")
  ok <- apply(boot, 3, function(a) all(is.finite(a)))
  boot <- boot[, , ok, drop = FALSE]
  Bu <- dim(boot)[3]
  if (Bu == 0) return(list(skipped = NA_character_, B_used = 0L, mcb = obs))
  p_minp <- vapply(seq_len(k), function(i)
    .minp_pooled(cbind(obs[i, ], matrix(boot[i, , ], nrow = m)))$p_item, 0)
  mt <- .maxT_items(obs, boot)
  list(skipped = NA_character_, B_used = Bu, mcb = obs,
       p_item_minp = p_minp, p_item = mt$p_item, padj = mt$padj,
       Tm = signif(mt$Tm, 7))
}
