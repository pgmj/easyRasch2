# Prototype: threshold-level DIF tests for two groups (RMdifAnchor development).
# Development only, not package code.
#
# Partial credit model, group g, item i, threshold k:
#   delta_ik(g) = beta_i(g) + tau_ik(g),  beta_i = mean of the item's thresholds.
# Fixing a group's scale adds the same constant to every delta in that group.
#
# Two families of tests:
#
#   Spacing (DSF without anchoring). The adjacent distances
#   delta_i,k+1 - delta_ik do not change when the scale is shifted, so their
#   group differences need no anchor. One joint Wald test per item with
#   m_i - 1 df (m_i thresholds). Corrections across items: Holm, BH, and a
#   Monte-Carlo single-step min-p from the joint normal of all contrasts.
#
#   Anchored thresholds. Each group is shifted so that the mean location of
#   the anchor items is zero, then delta_ik(2) - delta_ik(1) is tested per
#   threshold. This is the location difference plus the spacing difference.
#   Corrections across all thresholds: Holm and single-step max-|z|
#   (Monte-Carlo by default, multcomp optional).
#   The covariance is singular (the anchor items' location differences sum
#   to zero), which both single-step methods handle.
#
# Also returned: anchored item locations (mean of the item's anchored
# threshold differences), with the first anchor item untested as in
# psychotools::anchortest().
#
# Anchors: "centred" (all items), "forward" (iterative forward, MTT) and
# "gini" (constant class, Gini, Strobl et al. 2021), selected at item level by
# psychotools::anchor().
#
# Requires psychotools, multcomp, mvtnorm.

dif_thr_fit <- function(data, group) {
  group <- factor(group)
  if (nlevels(group) != 2L) stop("Prototype handles two groups only.")
  lv <- levels(group)
  data <- as.matrix(data)
  m1 <- psychotools::pcmodel(data[group == lv[1L], , drop = FALSE])
  m2 <- psychotools::pcmodel(data[group == lv[2L], , drop = FALSE])
  t1 <- psychotools::threshpar(m1, type = "mode", relative = FALSE, vcov = TRUE)
  t2 <- psychotools::threshpar(m2, type = "mode", relative = FALSE, vcov = TRUE)
  th1 <- stats::coef(t1)
  th2 <- stats::coef(t2)
  if (!identical(names(th1), names(th2))) {
    stop("Threshold sets differ between groups (unused category?).")
  }
  items <- colnames(data)
  thr_item <- sub("-C[0-9]+$", "", names(th1))
  list(
    m1 = m1, m2 = m2, th1 = th1, th2 = th2,
    V = stats::vcov(t1) + stats::vcov(t2),
    items = items,
    thr_item = match(thr_item, items),
    thr_k = as.integer(ave(seq_along(thr_item), thr_item, FUN = seq_along)),
    levels = lv
  )
}

# Item-level anchor selection. Returns indices into fit$items, in the order
# psychotools gives them (the first is the one left untested).
dif_thr_anchor <- function(fit, anchor = c("centred", "forward", "gini")) {
  anchor <- match.arg(anchor)
  switch(anchor,
    centred = seq_along(fit$items),
    forward = psychotools::anchor(fit$m1, fit$m2, class = "forward",
                                  select = "MTT")$anchor_items,
    gini = psychotools::anchor(fit$m1, fit$m2, class = "constant",
                               select = "Gini")$anchor_items
  )
}

# Tie-break sensitivity. psychotools::anchor() draws random numbers only
# when it ranks tied items (in anchorselect()). The auxiliary tests and the
# anchor class step are deterministic, so the auxiliary tests are computed
# once and the ranking is redone for each seed. Returns one anchor vector per
# seed (same order as psychotools gives).
dif_thr_anchor_ties <- function(fit, anchor = c("forward", "gini"), seeds) {
  anchor <- match.arg(anchor)
  ns <- asNamespace("psychotools")
  select <- if (anchor == "forward") "MTT" else "Gini"
  class <- if (anchor == "forward") "forward" else "constant"
  len <- if (anchor == "forward") length(fit$items) else 1L
  aux <- ns$auxtests(obj1 = fit$m1, obj2 = fit$m2, select)$aux_tests
  cache <- list()
  lapply(seeds, function(sd) {
    set.seed(sd)
    ro <- ns$anchorselect(obj1 = fit$m1, obj2 = fit$m2, aux, select)$ranking_order
    key <- paste(ro, collapse = ",")
    if (is.null(cache[[key]])) {
      cache[[key]] <<- ns$anchorclass(obj1 = fit$m1, obj2 = fit$m2,
                                      ranking_order = ro, class = class,
                                      length = len, range = c(0.1, 0.8))$anchor_items
    }
    cache[[key]]
  })
}

# Spacing test per item. B = Monte-Carlo draws for the single-step min-p.
dif_thr_spacing <- function(fit, B = 5000L) {
  d <- fit$th2 - fit$th1
  n <- length(d)
  poly <- which(tabulate(fit$thr_item, length(fit$items)) >= 2L)
  C <- NULL
  blk <- integer(0)
  for (i in poly) {
    idx <- which(fit$thr_item == i)
    Ci <- matrix(0, length(idx) - 1L, n)
    for (r in seq_len(length(idx) - 1L)) {
      Ci[r, idx[r]] <- -1
      Ci[r, idx[r + 1L]] <- 1
    }
    C <- rbind(C, Ci)
    blk <- c(blk, rep(i, nrow(Ci)))
  }
  est <- as.numeric(C %*% d)
  S <- C %*% fit$V %*% t(C)
  Sinv <- lapply(poly, function(i) solve(S[blk == i, blk == i, drop = FALSE]))
  W <- vapply(seq_along(poly), function(j) {
    e <- est[blk == poly[j]]
    as.numeric(t(e) %*% Sinv[[j]] %*% e)
  }, numeric(1))
  df <- vapply(poly, function(i) sum(blk == i), integer(1))
  p <- stats::pchisq(W, df, lower.tail = FALSE)
  # Single-step min-p under the joint null of all spacing contrasts.
  Z <- mvtnorm::rmvnorm(B, sigma = S, method = "eigen")
  psim <- vapply(seq_along(poly), function(j) {
    Zi <- Z[, blk == poly[j], drop = FALSE]
    stats::pchisq(rowSums((Zi %*% Sinv[[j]]) * Zi), df[j], lower.tail = FALSE)
  }, numeric(B))
  minp <- apply(psim, 1L, min)
  padj_ss <- vapply(p, function(pp) (1 + sum(minp <= pp)) / (B + 1), numeric(1))
  data.frame(
    item_no = poly, df = df, W = W, p = p,
    padj_holm = stats::p.adjust(p, "holm"),
    padj_bh = stats::p.adjust(p, "BH"),
    padj_ss = padj_ss
  )
}

# Single-step (max-|z|) adjusted p-values. "multcomp" uses multcomp's
# pmvnorm integration. "mc" draws B vectors from the joint normal of the
# estimates and compares each |z| with the simulated max-|z|. It is the same
# adjustment, computed by simulation, and much faster with 30-40 correlated
# tests. Singular covariance matrices are fine for both.
.ss_pvalues <- function(est, V, method = c("mc", "multcomp"), B = 5000L) {
  method <- match.arg(method)
  if (method == "multcomp") {
    names(est) <- paste0("c", seq_along(est))
    dimnames(V) <- list(names(est), names(est))
    gl <- multcomp::glht(multcomp::parm(est, V), linfct = diag(length(est)))
    return(as.numeric(summary(gl, test = multcomp::adjusted("single-step"))$test$pvalues))
  }
  se <- sqrt(diag(V))
  V <- (V + t(V)) / 2
  Z <- mvtnorm::rmvnorm(B, sigma = V, method = "eigen")
  maxz <- apply(abs(sweep(Z, 2L, se, "/")), 1L, max)
  vapply(abs(est / se), function(a) (1 + sum(maxz >= a)) / (B + 1), numeric(1))
}

# Anchored threshold and location differences (group 2 minus group 1).
dif_thr_anchored <- function(fit, anchor_items, drop_first = TRUE,
                             ss_method = "mc", B = 5000L) {
  d <- fit$th2 - fit$th1
  n <- length(d)
  k <- length(fit$items)
  w <- numeric(n)
  for (j in anchor_items) {
    idx <- which(fit$thr_item == j)
    w[idx] <- 1 / (length(anchor_items) * length(idx))
  }
  L <- diag(n) - matrix(1, n, 1L) %*% matrix(w, 1L)
  D <- as.numeric(L %*% d)
  Cv <- L %*% fit$V %*% t(L)
  se <- sqrt(pmax(diag(Cv), 0))
  ok <- se > 1e-8
  z <- ifelse(ok, D / se, NA_real_)
  p <- 2 * stats::pnorm(-abs(z))
  thr <- data.frame(
    item_no = fit$thr_item, k = fit$thr_k, est = D, se = se, z = z, p = p,
    padj_holm = NA_real_, padj_ss = NA_real_
  )
  thr$padj_holm[ok] <- stats::p.adjust(p[ok], "holm")
  thr$padj_ss[ok] <- tryCatch(.ss_pvalues(D[ok], Cv[ok, ok, drop = FALSE], ss_method, B),
                              error = function(e) NA_real_)
  # Item locations: mean of each item's anchored threshold differences.
  M <- t(vapply(seq_len(k), function(i) as.numeric(fit$thr_item == i) /
                  sum(fit$thr_item == i), numeric(n)))
  Dl <- as.numeric(M %*% D)
  Cl <- M %*% Cv %*% t(M)
  tested <- if (drop_first) setdiff(seq_len(k), anchor_items[1L]) else seq_len(k)
  sel <- seq_len(k) %in% tested
  sel_l <- sel & sqrt(pmax(diag(Cl), 0)) > 1e-8
  loc <- data.frame(item_no = seq_len(k), est = Dl,
                    se = sqrt(pmax(diag(Cl), 0)), z = NA_real_, p = NA_real_,
                    padj_holm = NA_real_, padj_ss = NA_real_,
                    is_anchor = seq_len(k) %in% anchor_items)
  loc$z[sel_l] <- Dl[sel_l] / loc$se[sel_l]
  loc$p[sel_l] <- 2 * stats::pnorm(-abs(loc$z[sel_l]))
  loc$padj_holm[sel_l] <- stats::p.adjust(loc$p[sel_l], "holm")
  loc$padj_ss[sel_l] <- tryCatch(
    .ss_pvalues(Dl[sel_l], Cl[sel_l, sel_l, drop = FALSE], ss_method, B),
    error = function(e) NA_real_)
  list(thresholds = thr, locations = loc)
}
