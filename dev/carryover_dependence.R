## Response dependence between occasions and RMpersonChange()
##
## Companion to dev/rci_exact_null.qmd. Implements the Marais (2009) carry-over
## model, validates it against her published Table 4, and runs the studies
## reported in that note's "Response dependence between occasions" section.
##
## Nothing runs on source(). Call the study functions:
##   marais_anchor()   validation against Marais Table 4
##   drift_table()     mean observed change by true theta, no true change
##   study_A()         type I error and flag direction, no true change
##   study_B()         attenuation, spurious severity gradient, power
##   study_C()         detectability via a racked Q3 screen
## study_A() and study_B() take a few minutes on 10 cores.

suppressMessages(library(easyRasch2))
library(parallel)

## --- vectorised PCM sampler --------------------------------------------
## Shifting every threshold of an item by s is the same as shifting theta by
## -s, so a per-person shift is applied as an effective theta.
rpcm <- function(theta_eff, thr) {
  m   <- length(thr)
  num <- matrix(0, length(theta_eff), m + 1L)
  cs  <- cumsum(thr)
  for (c in seq_len(m)) num[, c + 1L] <- c * theta_eff - cs[c]
  num <- exp(num - apply(num, 1L, max))
  p   <- num / rowSums(num)
  cp  <- t(apply(p, 1L, cumsum))
  max.col(cp >= runif(length(theta_eff)), ties.method = "first") - 1L
}

## Marais (2009) Eq. 3: the time-2 difficulty of item j shifts by
## (1 - 2 * x_j1) * d, so a 1 at time 1 makes the item easier at time 2 and a 0
## makes it harder. The polytomous generalisation shifts every threshold by
## (1 - 2 * x_j1 / m_j) * d, which reduces to Marais when m_j = 1.
sim_pair <- function(thr_list, theta1, theta2, d = 0) {
  k  <- length(thr_list); n <- length(theta1)
  mj <- vapply(thr_list, length, integer(1L))
  x1 <- vapply(seq_len(k), function(j) rpcm(theta1, thr_list[[j]]), integer(n))
  sh <- if (d == 0) matrix(0, n, k) else (1 - 2 * sweep(x1, 2L, mj, "/")) * d
  x2 <- vapply(seq_len(k), function(j)
          rpcm(theta2 - sh[, j], thr_list[[j]]), integer(n))
  colnames(x1) <- colnames(x2) <- paste0("I", seq_len(k))
  list(t1 = as.data.frame(x1), t2 = as.data.frame(x2))
}

## Item locations over [-2.5, 2.5], following Marais. Note this is a wider
## spread than mk() in dev/rci_exact_null.qmd, which uses [-1.2, 1.2].
mk <- function(k, ncat, lo = -2.5, hi = 2.5) {
  b <- seq(lo, hi, length.out = k)
  if (ncat == 2L) lapply(b, function(bb) bb)
  else lapply(b, function(bb) bb + seq(-0.9, 0.9, length.out = ncat - 1L))
}

## --- fast RCI path (complete data, known item parameters) --------------
## Verified against RMpersonChange() at zero difference in RCI and zero
## classification disagreements over 1000 respondents.
lut <- function(thr_list) {
  steps <- vapply(thr_list, length, integer(1L))
  m <- vapply(0:sum(steps), function(r)
    easyRasch2:::.theta_wle(easyRasch2:::.score_pattern(r, steps),
                            thr_list, c(-10, 10)), numeric(2L))
  list(theta = m[1L, ], se = m[2L, ])
}
score_dist <- function(thr_list, theta) {
  dd <- 1
  for (thr in thr_list) {
    p  <- easyRasch2:::.pcm_cat_probs(theta, thr)
    nd <- numeric(length(dd) + length(p) - 1L)
    for (a in seq_along(dd)) nd[a + seq_along(p) - 1L] <-
        nd[a + seq_along(p) - 1L] + dd[a] * p
    dd <- nd
  }
  dd
}
## Exact two-sided critical values of the RCI under INDEPENDENCE, as in
## dev/rci_exact_null.qmd. This is what the package hands the user.
exact_crit <- function(thr_list, thetas, alpha = 0.05) {
  L <- lut(thr_list); th <- L$theta; se <- L$se
  rci <- outer(seq_along(th), seq_along(th),
               function(i, j) (th[j] - th[i]) / sqrt(se[i]^2 + se[j]^2))
  w <- matrix(0, length(th), length(th))
  for (t in thetas) { dd <- score_dist(thr_list, t); w <- w + outer(dd, dd) }
  w <- w / sum(w)
  o <- order(as.numeric(rci)); v <- as.numeric(rci)[o]
  cw <- cumsum(as.numeric(w)[o])
  c(lower = v[which(cw >= alpha / 2)[1L]],
    upper = v[which(cw >= 1 - alpha / 2)[1L]])
}
rci_fast <- function(pair, thr_list, crit) {
  L <- lut(thr_list)
  r1 <- rowSums(pair$t1); r2 <- rowSums(pair$t2)
  th1 <- L$theta[r1 + 1L]; s1 <- L$se[r1 + 1L]
  th2 <- L$theta[r2 + 1L]; s2 <- L$se[r2 + 1L]
  rci <- (th2 - th1) / sqrt(s1^2 + s2^2)
  data.frame(theta_t1 = th1, theta_t2 = th2, change = th2 - th1, rci = rci,
             flag = rci < crit[["lower"]] | rci > crit[["upper"]],
             up = rci > crit[["upper"]], down = rci < crit[["lower"]])
}

## --- validation against Marais Table 4 ---------------------------------
marais_anchor <- function(reps = 3, K = 25, N = 500, seed = 7) {
  set.seed(seed)
  thr <- mk(K, 2); pub <- c("-0.01", "0.14", "0.31", "0.47")
  dd <- c(0, 0.5, 1, 1.5)
  for (i in seq_along(dd)) {
    r <- replicate(reps, {
      th <- rnorm(N, 0.85, 1)
      p  <- sim_pair(thr, th, th, d = dd[i])
      rack <- cbind(p$t1, setNames(p$t2, paste0(names(p$t2), "_t2")))
      q <- as.matrix(RMlocdepQ3(rack, output = "dataframe"))
      mean(q[cbind(K + seq_len(K), seq_len(K))])
    })
    cat(sprintf("d = %.1f   same-item Q3 = %+.3f   Marais Table 4: %s\n",
                dd[i], mean(r), pub[i]))
  }
}

## --- the drift ---------------------------------------------------------
drift_table <- function(K = 20, ncat = 4, N = 4000, seed = 808) {
  set.seed(seed)
  thr <- mk(K, ncat); L <- lut(thr); tg <- seq(-2.5, 2.5, by = 0.5)
  tab <- sapply(c(0, 0.5, 1, 1.5), function(d)
    sapply(tg, function(t0) {
      th <- rep(t0, N); p <- sim_pair(thr, th, th, d = d)
      mean(L$theta[rowSums(p$t2) + 1L] - L$theta[rowSums(p$t1) + 1L])
    }))
  dimnames(tab) <- list(sprintf("% .1f", tg), paste0("d=", c(0, 0.5, 1, 1.5)))
  round(tab, 3)
}

## --- Study A: type I error and flag direction, no true change ----------
study_A <- function(N = 2000, REPS = 25, cores = 10, seed = 99) {
  scales <- list("20 items, 4 cat" = mk(20, 4), "6 items, 4 cat" = mk(6, 4),
                 "25 items, 2 cat" = mk(25, 2))
  targets <- c("on target" = 0, "items easy (mu +1.5)" = 1.5)
  set.seed(seed)
  grid <- expand.grid(scale = names(scales), targ = names(targets),
                      d = c(0, 0.5, 1, 1.5), stringsAsFactors = FALSE)
  one <- function(i) {
    g <- grid[i, ]; thr <- scales[[g$scale]]; mu <- targets[[g$targ]]
    crit <- exact_crit(thr, rnorm(4000, mu, 1))   # calibrated under independence
    res <- do.call(rbind, lapply(seq_len(REPS), function(rep) {
      th <- rnorm(N, mu, 1)
      f <- rci_fast(sim_pair(thr, th, th, d = g$d), thr, crit)
      f$theta_true <- th; f
    }))
    band <- cut(res$theta_true, c(-Inf, -1, 1, Inf),
                labels = c("below", "near", "above"))
    data.frame(scale = g$scale, targ = g$targ, d = g$d, rci_sd = sd(res$rci),
               flag = mean(res$flag),
               flag_below = mean(res$flag[band == "below"]),
               flag_near  = mean(res$flag[band == "near"]),
               flag_above = mean(res$flag[band == "above"]),
               up_above   = mean(res$up[band == "above"]),
               down_below = mean(res$down[band == "below"]),
               stringsAsFactors = FALSE)
  }
  do.call(rbind, mclapply(seq_len(nrow(grid)), one, mc.cores = cores))
}

## --- Study B: real change of +0.5 logits for everyone ------------------
study_B <- function(N = 2000, REPS = 25, ALPHA = 0.5, cores = 10, seed = 4242) {
  scales <- list("20 items, 4 cat" = mk(20, 4), "6 items, 4 cat" = mk(6, 4))
  set.seed(seed)
  grid <- expand.grid(scale = names(scales), d = c(0, 0.5, 1, 1.5),
                      stringsAsFactors = FALSE)
  one <- function(i) {
    g <- grid[i, ]; thr <- scales[[g$scale]]
    crit <- exact_crit(thr, rnorm(4000, 0, 1))
    res <- do.call(rbind, lapply(seq_len(REPS), function(rep) {
      th <- rnorm(N, 0, 1)
      f <- rci_fast(sim_pair(thr, th, th + ALPHA, d = g$d), thr, crit)
      f$theta_true <- th; f
    }))
    band <- cut(res$theta_true, c(-Inf, -1, 1, Inf),
                labels = c("below", "near", "above"))
    data.frame(scale = g$scale, d = g$d, mean_change = mean(res$change),
               ch_below = mean(res$change[band == "below"]),
               ch_near  = mean(res$change[band == "near"]),
               ch_above = mean(res$change[band == "above"]),
               power = mean(res$up),
               pw_below = mean(res$up[band == "below"]),
               pw_above = mean(res$up[band == "above"]),
               stringsAsFactors = FALSE)
  }
  do.call(rbind, mclapply(seq_len(nrow(grid)), one, mc.cores = cores))
}

## --- Study C: would a user detect it? ----------------------------------
## Rack both occasions into one 2k-item analysis and screen the same-item pairs
## with Yen-corrected Q3 against RMlocdepQ3Cutoff()'s simulated cutoff.
study_C <- function(K = 20, N = 500, REPS = 5, iterations = 150, cores = 10,
                    seed = 321) {
  set.seed(seed); thr <- mk(K, 4)
  do.call(rbind, lapply(c(0, 0.5, 1, 1.5), function(d) {
    r <- do.call(rbind, lapply(seq_len(REPS), function(rep) {
      th <- rnorm(N, 0, 1)
      p  <- sim_pair(thr, th, th, d = d)
      rack <- cbind(p$t1, setNames(p$t2, paste0(names(p$t2), "_t2")))
      cut <- RMlocdepQ3Cutoff(rack, iterations = iterations, n_cores = cores,
                              seed = rep)
      q   <- as.matrix(RMlocdepQ3(rack, output = "dataframe"))
      off <- q[lower.tri(q)]; off <- off[!is.na(off)]
      same <- (q - mean(off))[cbind(K + seq_len(K), seq_len(K))]
      ct <- as.numeric(cut$suggested_cutoff)
      data.frame(same = mean(same), cutoff = ct, nflag = sum(same > ct))
    }))
    data.frame(d = d, same_q3star = mean(r$same), cutoff = mean(r$cutoff),
               flagged = mean(r$nflag), pct_any = 100 * mean(r$nflag > 0))
  }))
}

## --- Study D: does heterogeneous true change look like carry-over? -----
## Carry-over is same-item specific. Heterogeneous change is a second person
## dimension and acts on the whole cross-occasion block. study_D_blocks()
## reports the Q3 means by block; study_D_screen() reports what the cutoff
## screen actually flags, under two centrings of the same-item statistic.
.q3_blocks <- function(t1, t2, K) {
  rack <- cbind(t1, setNames(t2, paste0(names(t2), "_t2")))
  q <- as.matrix(RMlocdepQ3(rack, output = "dataframe"))
  off <- q[lower.tri(q)]; off <- off[!is.na(off)]
  qs <- q - mean(off); i <- seq_len(K)
  dg <- qs[cbind(K + i, i)]
  cr <- qs[(K + 1):(2 * K), i]
  od <- cr[row(cr) != col(cr)]
  w1 <- qs[i, i][lower.tri(diag(K))]
  w2 <- qs[(K + 1):(2 * K), (K + 1):(2 * K)][lower.tri(diag(K))]
  c(same_item = mean(dg), other_cross = mean(od, na.rm = TRUE),
    contrast = mean(dg) - mean(od, na.rm = TRUE),
    within = mean(c(w1, w2), na.rm = TRUE))
}

.conds <- list(
  list(lab = "no change, no carry-over",   d = 0,   mu = 0,   sd = 0),
  list(lab = "uniform change +0.5",        d = 0,   mu = 0.5, sd = 0),
  list(lab = "change SD 0.6",              d = 0,   mu = 0.5, sd = 0.6),
  list(lab = "change SD 1.0",              d = 0,   mu = 0.5, sd = 1.0),
  list(lab = "carry-over d = 0.5",         d = 0.5, mu = 0,   sd = 0),
  list(lab = "carry-over d = 1.0",         d = 1,   mu = 0,   sd = 0),
  list(lab = "both, d = 1 and SD 0.6",     d = 1,   mu = 0.5, sd = 0.6))

study_D_blocks <- function(K = 20, N = 500, REPS = 20, seed = 2024) {
  set.seed(seed); thr <- mk(K, 4)
  do.call(rbind, lapply(.conds, function(cc) {
    r <- replicate(REPS, {
      th <- rnorm(N, 0, 1)
      p <- sim_pair(thr, th, th + cc$mu + rnorm(N, 0, cc$sd), d = cc$d)
      .q3_blocks(p$t1, p$t2, K)
    })
    data.frame(condition = cc$lab, t(round(rowMeans(r), 3)),
               stringsAsFactors = FALSE)
  }))
}

## Yen centring subtracts the mean of ALL off-diagonal Q3; cross centring
## subtracts the mean of the cross-occasion off-diagonal pairs only. The second
## has more power and is the one the help page recommends.
study_D_screen <- function(K = 20, N = 500, REPS = 5, iterations = 150,
                           cores = 10, seed = 555) {
  set.seed(seed); thr <- mk(K, 4)
  one <- function(t1, t2) {
    rack <- cbind(t1, setNames(t2, paste0(names(t2), "_t2")))
    ct <- as.numeric(RMlocdepQ3Cutoff(rack, iterations = iterations,
                                      n_cores = cores)$suggested_cutoff)
    q <- as.matrix(RMlocdepQ3(rack, output = "dataframe"))
    allo <- q[lower.tri(q)]; allo <- allo[!is.na(allo)]
    i <- seq_len(K); dg <- q[cbind(K + i, i)]
    cr <- q[(K + 1):(2 * K), i]
    cro <- mean(cr[row(cr) != col(cr)], na.rm = TRUE)
    c(yen = sum(dg - mean(allo) > ct), cross = sum(dg - cro > ct))
  }
  do.call(rbind, lapply(.conds, function(cc) {
    r <- replicate(REPS, {
      th <- rnorm(N, 0, 1)
      p <- sim_pair(thr, th, th + cc$mu + rnorm(N, 0, cc$sd), d = cc$d)
      one(p$t1, p$t2)
    })
    data.frame(condition = cc$lab,
               yen_flagged = mean(r["yen", ]), cross_flagged = mean(r["cross", ]),
               yen_any = 100 * mean(r["yen", ] > 0),
               cross_any = 100 * mean(r["cross", ] > 0),
               stringsAsFactors = FALSE)
  }))
}
