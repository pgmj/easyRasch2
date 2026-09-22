## Exact conditional null for the PCA/t-test unidimensionality protocol
## (Smith 2002 / Tennant & Pallant 2006 / RUMM2030; see Hagell 2014).
##
## Companion to dev/rci_exact_null.qmd: the protocol statistic is the RCI
## with item subtests in place of occasions, so the same enumeration applies.
## Here the null is conditional on the TOTAL score, which under the Rasch
## model is free of theta (verified in the validation block below).
##
## Run: source() this file. Deterministic apart from the Monte Carlo anchor.


## ---- thresholds -------------------------------------------------------
mk <- function(k, ncat, spread = 1.2) {
  lapply(seq(-spread, spread, length.out = k),
         function(b) b + seq(-0.9, 0.9, length.out = ncat - 1L))
}

## ---- Lord-Wingersky: exact sum-score distribution given theta ----------
score_dist <- function(thr_list, theta) {
  d <- 1
  for (thr in thr_list) {
    p  <- easyRasch2:::.pcm_cat_probs(theta, thr)
    nd <- numeric(length(d) + length(p) - 1L)
    for (a in seq_along(d)) nd[a + seq_along(p) - 1L] <-
        nd[a + seq_along(p) - 1L] + d[a] * p
    d <- nd
  }
  d                                   # index = score + 1
}

## ---- score -> (theta, SE) lookup, WLE ---------------------------------
lut <- function(thr_list) {
  steps <- vapply(thr_list, length, integer(1L))
  m <- vapply(0:sum(steps), function(r)
    easyRasch2:::.theta_wle(easyRasch2:::.score_pattern(r, steps),
                            thr_list, c(-10, 10)), numeric(2L))
  list(theta = m[1L, ], se = m[2L, ])
}

## ---- the t grid for one split ----------------------------------------
## t[a+1, b+1] = (theta_A(a) - theta_B(b)) / sqrt(SE_A(a)^2 + SE_B(b)^2)
split_grid <- function(thr_all, idxA) {
  A <- thr_all[idxA]; B <- thr_all[-idxA]
  lA <- lut(A); lB <- lut(B)
  tg <- outer(seq_along(lA$theta), seq_along(lB$theta),
              function(i, j) (lA$theta[i] - lB$theta[j]) /
                sqrt(lA$se[i]^2 + lB$se[j]^2))
  dg <- outer(seq_along(lA$theta), seq_along(lB$theta),
              function(i, j) lA$theta[i] - lB$theta[j])
  list(A = A, B = B, lA = lA, lB = lB, t = tg, d = dg)
}

## ---- conditional-on-total-score null ---------------------------------
## Under the Rasch model P(r_A = a | r_A + r_B = r) is free of theta.
## Computed at theta_eval and normalised; theta-invariance is checked below.
cond_null <- function(g, theta_eval = 0) {
  dA <- score_dist(g$A, theta_eval)
  dB <- score_dist(g$B, theta_eval)
  maxA <- length(dA) - 1L; maxB <- length(dB) - 1L
  lapply(0:(maxA + maxB), function(r) {
    a  <- max(0L, r - maxB):min(maxA, r)
    w  <- dA[a + 1L] * dB[r - a + 1L]
    list(a = a, p = w / sum(w),
         t = g$t[cbind(a + 1L, r - a + 1L)],
         d = g$d[cbind(a + 1L, r - a + 1L)])
  })
}

## per-total-score exact type I error of a fixed critical value
alpha_by_r <- function(cn, crit = 1.96)
  vapply(cn, function(x) sum(x$p[abs(x$t) > crit]), numeric(1L))

## expected raw logit gap per total score (WLE-bias / centring check)
gap_by_r <- function(cn) vapply(cn, function(x) sum(x$p * x$d), numeric(1L))

## ---- total-score distribution under a latent density ------------------
score_weights <- function(thr_all, mu = 0, sigma = 1.4, nq = 121L) {
  q <- seq(mu - 5 * sigma, mu + 5 * sigma, length.out = nq)
  w <- dnorm(q, mu, sigma); w <- w / sum(w)
  Reduce(`+`, Map(function(t, wt) wt * score_dist(thr_all, t), q, w))
}

## ---- pooled exact critical values achieving nominal alpha -------------
crit_pooled <- function(cn, fr, alpha = 0.05) {
  tv <- unlist(lapply(cn, `[[`, "t"))
  pv <- unlist(Map(function(x, f) x$p * f, cn, fr))
  o  <- order(tv); tv <- tv[o]; pv <- pv[o] / sum(pv)
  cw <- cumsum(pv)
  c(lower = tv[which(cw >= alpha / 2)[1L]],
    upper = tv[which(cw >= 1 - alpha / 2)[1L]])
}

## ---- exact Poisson-binomial distribution of the flagged count ---------
pois_binom <- function(p) {
  d <- 1
  for (pi in p) d <- c(d * (1 - pi), 0) + c(0, d * pi)
  d                                    # index = count + 1
}

## ======================= validation anchors =========================

thr  <- mk(12, 4)
idxA <- seq(1, 12, by = 2)          # interleaved 6 / 6
g    <- split_grid(thr, idxA)

## --- anchor 1: theta-invariance of the conditional null ---------------
cn0 <- cond_null(g, 0); cnm <- cond_null(g, -2); cnp <- cond_null(g, 2.5)
dev <- max(vapply(seq_along(cn0), function(i)
  max(abs(cn0[[i]]$p - cnm[[i]]$p), abs(cn0[[i]]$p - cnp[[i]]$p)), numeric(1L)))
cat(sprintf("anchor 1  max |P(a|r) difference| across theta -2 / 0 / 2.5 : %.3e\n", dev))

## --- anchor 2: Monte Carlo cross-check of alpha_r ---------------------
set.seed(11)
N <- 400000
th <- rnorm(N, 0, 1.4)
dat <- easyRasch2:::sim_partial_score(thr, th)
rA <- rowSums(dat[, idxA, drop = FALSE]); rB <- rowSums(dat[, -idxA, drop = FALSE])
tt <- g$t[cbind(rA + 1L, rB + 1L)]
r  <- rA + rB
a_exact <- alpha_by_r(cn0, 1.96)
mc <- tapply(abs(tt) > 1.96, factor(r, levels = 0:(length(a_exact) - 1L)), mean)
nn <- as.vector(table(factor(r, levels = 0:(length(a_exact) - 1L))))
keep <- nn >= 2000
cmp <- data.frame(r = (0:(length(a_exact) - 1L))[keep], n = nn[keep],
                  exact = round(a_exact[keep], 4),
                  mc = round(as.vector(mc)[keep], 4))
cmp$diff <- round(cmp$mc - cmp$exact, 4)
cat("\nanchor 2  exact vs Monte Carlo alpha_r (scores with n >= 2000)\n")
print(cmp, row.names = FALSE)
cat(sprintf("max abs difference: %.4f   (MC se ~ %.4f)\n",
            max(abs(cmp$diff)), max(sqrt(0.05 * 0.95 / cmp$n))))

## --- anchor 3: pooled critical value recovers nominal alpha -----------
fr <- score_weights(thr, 0, 1.4)
cp <- crit_pooled(cn0, fr, 0.05)
ach <- sum(vapply(seq_along(cn0), function(i)
  fr[i] * sum(cn0[[i]]$p[cn0[[i]]$t < cp[["lower"]] | cn0[[i]]$t > cp[["upper"]]]),
  numeric(1L))) / sum(fr)
cat(sprintf("\nanchor 3  pooled exact crit = [%.3f, %.3f], achieved alpha = %.4f (nominal .05)\n",
            cp[["lower"]], cp[["upper"]], ach))

## ======================= configuration sweep ========================

run <- function(label, k, ncat, idxA, mu = 0, sigma = 1.4) {
  thr <- mk(k, ncat)
  g   <- split_grid(thr, idxA)
  cn  <- cond_null(g, 0)
  fr  <- score_weights(thr, mu, sigma); fr <- fr / sum(fr)
  a   <- alpha_by_r(cn, 1.96)
  gp  <- gap_by_r(cn)
  cp  <- crit_pooled(cn, fr, 0.05)
  nz  <- a > 0
  data.frame(
    config   = label,
    kA       = length(idxA), kB = k - length(idxA),
    alpha_bar= sum(fr * a),                       # true null rate of 1.96
    ratio    = 0.05 / sum(fr * a),                # how lenient the 5% mark is
    pct_zero = 100 * sum(fr[!nz]),                # respondents who cannot be flagged
    a_max    = max(a),
    crit_lo  = cp[["lower"]], crit_hi = cp[["upper"]],
    asym     = abs(cp[["lower"]]) - cp[["upper"]],
    max_gap  = max(abs(gp[nz])),                  # |E[thetaA - thetaB | r]| in logits
    stringsAsFactors = FALSE)
}

cfg <- list(
  run("12 items, 4 cat, 6/6",        12, 4, seq(1, 12, by = 2)),
  run("12 items, 4 cat, 4/8",        12, 4, seq(1, 12, by = 3)),
  run("20 items, 4 cat, 10/10",      20, 4, seq(1, 20, by = 2)),
  run("20 items, 4 cat, 7/13",       20, 4, seq(1, 20, by = 3)),
  run("20 items, 4 cat, 4/16",       20, 4, seq(1, 20, by = 5)),
  run("30 items, 4 cat, 15/15",      30, 4, seq(1, 30, by = 2)),
  run("40 items, 4 cat, 20/20",      40, 4, seq(1, 40, by = 2)),
  run("20 items, 2 cat, 10/10",      20, 2, seq(1, 20, by = 2)),
  run("20 items, 7 cat, 10/10",      20, 7, seq(1, 20, by = 2)),
  run("20 items, 4 cat, easy/hard",  20, 4, 1:10),
  run("20 items, 4 cat, 10/10, off target", 20, 4, seq(1, 20, by = 2), mu = -2)
)
res <- do.call(rbind, cfg)

cat("\n=== Exact null rate of the 1.96 rule, conditional on total score ===\n\n")
print(data.frame(
  config   = res$config,
  subtests = paste0(res$kA, "/", res$kB),
  `true alpha` = round(res$alpha_bar, 4),
  `5% is X times too lenient` = round(res$ratio, 2),
  `% cannot be flagged` = round(res$pct_zero, 1),
  check.names = FALSE), row.names = FALSE)

cat("\n=== Exact critical values that would give 5%, and centring ===\n\n")
print(data.frame(
  config = res$config,
  lower  = round(res$crit_lo, 3),
  upper  = round(res$crit_hi, 3),
  `asymmetry` = round(res$asym, 3),
  `max |E gap| logits` = round(res$max_gap, 3),
  check.names = FALSE), row.names = FALSE)


## ===================== count-level exact null =======================

count_null <- function(label, k, ncat, idxA, n, mu = 0, sigma = 1.4) {
  thr <- mk(k, ncat); g <- split_grid(thr, idxA); cn <- cond_null(g, 0)
  fr  <- score_weights(thr, mu, sigma); fr <- fr / sum(fr)
  a   <- alpha_by_r(cn, 1.96)
  ## allocate n respondents across total scores in expectation
  nr  <- round(n * fr); nr[nr < 0] <- 0
  p_i <- rep(a, times = nr)
  d   <- pois_binom(p_i)                       # exact, index = count + 1
  cnt <- seq_along(d) - 1L
  mean_p <- sum(d * cnt) / length(p_i)
  sd_pb  <- sqrt(sum(d * cnt^2) - sum(d * cnt)^2)
  sd_bin <- sqrt(length(p_i) * mean_p * (1 - mean_p))
  q95 <- cnt[which(cumsum(d) >= 0.95)[1L]] / length(p_i)
  p_over5 <- sum(d[cnt / length(p_i) > 0.05])  # P(observed proportion > 5%)
  data.frame(config = label, n = length(p_i),
             mean_prop = mean_p, crit_prop_95 = q95,
             P_obs_over_5pct = p_over5,
             sd_poisbinom = sd_pb, sd_binomial = sd_bin,
             stringsAsFactors = FALSE)
}

grid <- do.call(rbind, c(
  lapply(c(200, 500, 1000), function(n)
    count_null("20 items, 4 cat, 10/10", 20, 4, seq(1, 20, by = 2), n)),
  lapply(c(200, 500, 1000), function(n)
    count_null("12 items, 4 cat, 6/6", 12, 4, seq(1, 12, by = 2), n)),
  lapply(c(200, 500, 1000), function(n)
    count_null("20 items, 2 cat, 10/10", 20, 2, seq(1, 20, by = 2), n)),
  lapply(c(200, 500, 1000), function(n)
    count_null("40 items, 4 cat, 20/20", 40, 4, seq(1, 40, by = 2), n))))

cat("\n=== Exact null distribution of the observed proportion ===\n")
cat("    (1.96 rule, conditional on total scores; protocol benchmark is 5%)\n\n")
print(data.frame(
  config = grid$config, n = grid$n,
  `mean %`        = round(100 * grid$mean_prop, 2),
  `95th pct %`    = round(100 * grid$crit_prop_95, 2),
  `P(obs > 5%)`   = round(grid$P_obs_over_5pct, 4),
  `SD pois-binom` = round(grid$sd_poisbinom, 2),
  `SD binomial`   = round(grid$sd_binomial, 2),
  check.names = FALSE), row.names = FALSE)

## ============== extreme-score exclusion robustness ==================

## Marginal joint over (r_A, r_B) under a latent density, then alpha with and
## without respondents extreme on either subtest (what MLE software must drop).
joint_w <- function(thr_all, idxA, mu = 0, sigma = 1.4, nq = 121L) {
  A <- thr_all[idxA]; B <- thr_all[-idxA]
  q <- seq(mu - 5 * sigma, mu + 5 * sigma, length.out = nq)
  w <- dnorm(q, mu, sigma); w <- w / sum(w)
  Reduce(`+`, Map(function(t, wt) wt * outer(score_dist(A, t), score_dist(B, t)),
                  q, w))
}

chk <- function(label, k, ncat, idxA, mu = 0, sigma = 1.4) {
  thr <- mk(k, ncat); g <- split_grid(thr, idxA)
  W <- joint_w(thr, idxA, mu, sigma); W <- W / sum(W)
  sig <- abs(g$t) > 1.96
  nA <- nrow(W); nB <- ncol(W)
  keep <- matrix(TRUE, nA, nB)
  keep[c(1L, nA), ] <- FALSE; keep[, c(1L, nB)] <- FALSE   # extreme on either
  data.frame(config = label,
             alpha_all  = sum(W * sig),
             pct_extreme = 100 * (1 - sum(W[keep])),
             alpha_nonextreme = sum(W[keep] * sig[keep]) / sum(W[keep]),
             stringsAsFactors = FALSE)
}

out <- rbind(
  chk("12 items, 4 cat, 6/6",   12, 4, seq(1, 12, by = 2)),
  chk("20 items, 4 cat, 10/10", 20, 4, seq(1, 20, by = 2)),
  chk("20 items, 4 cat, 4/16",  20, 4, seq(1, 20, by = 5)),
  chk("30 items, 4 cat, 15/15", 30, 4, seq(1, 30, by = 2)),
  chk("40 items, 4 cat, 20/20", 40, 4, seq(1, 40, by = 2)),
  chk("20 items, 2 cat, 10/10", 20, 2, seq(1, 20, by = 2)),
  chk("20 items, 4 cat, easy/hard", 20, 4, 1:10),
  chk("20 items, 4 cat, off target", 20, 4, seq(1, 20, by = 2), mu = -2))

cat("\n=== Does dropping extreme-score respondents rescue the 5% benchmark? ===\n\n")
print(data.frame(
  config = out$config,
  `alpha, all retained`   = round(out$alpha_all, 4),
  `% extreme on a subtest`= round(out$pct_extreme, 1),
  `alpha, extremes dropped` = round(out$alpha_nonextreme, 4),
  check.names = FALSE), row.names = FALSE)
