## Conditioning the RCI null on the two-occasion total
##
## Companion to dev/rci_exact_null.qmd, section "Can the sample dependence be
## removed?". Under the Rasch model the distribution of r1 given r1 + r2 is free
## of theta, so a per-respondent RCI null can be built with no latent
## distribution and no plug-in. This script works out whether that is worth
## doing. The answer is no; see the note.
##
## Reuses the generator and score-lookup helpers from carryover_dependence.R
## (mk, sim_pair, lut, score_dist, exact_crit).
source(file.path(dirname(sys.frame(1)$ofile %||% "."), "carryover_dependence.R"))

## Nothing runs on source(). Call:
##   theta_free_anchor()   invariance check across configurations
##   size_table()          achievable size, strict and mid-p
##   typeI_study()         type I error of five rules across latent densities
##   single_subject()      exact size given true theta, one respondent
##   power_table()         exact power at true theta = 0

rci_grid <- function(thr) {
  L <- lut(thr); R <- length(L$theta) - 1L
  list(R = R, theta = L$theta, se = L$se,
       rci = outer(seq_len(R + 1L), seq_len(R + 1L),
                   function(i, j) (L$theta[j] - L$theta[i]) /
                     sqrt(L$se[i]^2 + L$se[j]^2)))
}

## theta-free conditional null of the RCI given the two-occasion total.
## Evaluated at one convenient theta and normalised within each total, which
## gives gamma(a) gamma(s - a) / gamma_2(s).
cond_total_null <- function(thr, theta_eval = 0) {
  g <- rci_grid(thr); d <- score_dist(thr, theta_eval); R <- g$R
  lapply(0:(2 * R), function(s) {
    a <- max(0L, s - R):min(R, s)
    w <- d[a + 1L] * d[s - a + 1L]
    list(a = a, p = w / sum(w), rci = g$rci[cbind(a + 1L, s - a + 1L)])
  })
}

## conditional p-value for every attainable (r1, r2); midp halves the tie mass
ct_pmat <- function(thr, midp = FALSE) {
  cn <- cond_total_null(thr); R <- rci_grid(thr)$R
  P <- matrix(NA_real_, R + 1L, R + 1L)
  for (si in seq_along(cn)) {
    x <- cn[[si]]; ax <- abs(x$rci)
    for (u in seq_along(x$a)) {
      gt <- sum(x$p[ax > ax[u] + 1e-12])
      eq <- sum(x$p[abs(ax - ax[u]) <= 1e-12])
      P[x$a[u] + 1L, (si - 1L) - x$a[u] + 1L] <- if (midp) gt + 0.5 * eq else gt + eq
    }
  }
  P
}

## per-respondent critical values from a plug-in theta (conditional_crit)
plugin_crit <- function(thr, theta_null, alpha = 0.05) {
  g <- rci_grid(thr)
  t(vapply(theta_null, function(t0) {
    d <- score_dist(thr, t0); w <- outer(d, d)
    o <- order(as.numeric(g$rci)); v <- as.numeric(g$rci)[o]
    cw <- cumsum(as.numeric(w)[o] / sum(w))
    c(v[which(cw >= alpha / 2)[1L]], v[which(cw >= 1 - alpha / 2)[1L]])
  }, numeric(2L)))
}

theta_free_anchor <- function() {
  cfgs <- list("6 items, 4 cat" = mk(6, 4), "20 items, 4 cat" = mk(20, 4),
               "20 items, 2 cat" = mk(20, 2),
               "9 items, mixed 2/3/4 cat" = c(mk(3, 2), mk(3, 3), mk(3, 4)))
  for (nm in names(cfgs)) {
    thr <- cfgs[[nm]]
    a <- cond_total_null(thr, 0); b <- cond_total_null(thr, -2.5)
    d <- cond_total_null(thr, 2)
    dev <- max(vapply(seq_along(a), function(i)
      max(abs(a[[i]]$p - b[[i]]$p), abs(a[[i]]$p - d[[i]]$p)), numeric(1)))
    cat(sprintf("%-26s max |P(r1 | total) difference|: %.2e\n", nm, dev))
  }
}

.total_weights <- function(thr) {
  q <- seq(-5, 5, length.out = 121); w <- dnorm(q); w <- w / sum(w)
  ds <- Reduce(`+`, Map(function(t, wt) {
    d <- score_dist(thr, t)
    wt * as.numeric(stats::convolve(d, rev(d), type = "open"))
  }, q, w))
  ds / sum(ds)
}

size_table <- function(alpha = 0.05) {
  cfgs <- list("6 items, 4 cat" = mk(6, 4), "12 items, 4 cat" = mk(12, 4),
               "20 items, 4 cat" = mk(20, 4), "20 items, 2 cat" = mk(20, 2),
               "40 items, 4 cat" = mk(40, 4))
  do.call(rbind, lapply(names(cfgs), function(nm) {
    thr <- cfgs[[nm]]; ds <- .total_weights(thr)
    sz <- function(midp) {
      P <- ct_pmat(thr, midp); cn <- cond_total_null(thr)
      vapply(seq_along(cn), function(si) {
        x <- cn[[si]]
        sum(x$p[P[cbind(x$a + 1L, (si - 1L) - x$a + 1L)] <= alpha])
      }, numeric(1))
    }
    s0 <- sz(FALSE); s1 <- sz(TRUE)
    data.frame(scale = nm, strict = sum(ds * s0), midp = sum(ds * s1),
               cannot_reject_pct = 100 * sum(ds[s0 == 0]),
               stringsAsFactors = FALSE)
  }))
}

typeI_study <- function(N = 2000, REPS = 25, alpha = 0.05, cores = 8, seed = 606) {
  scales <- list("6 items, 4 cat" = mk(6, 4), "20 items, 4 cat" = mk(20, 4))
  dists <- list(c(0, 1), c(0, 1.5), c(-2, 1), c(1.5, 1))
  dlab <- c("N(0, 1)", "N(0, 1.5)", "N(-2, 1)", "N(1.5, 1)")
  set.seed(seed)
  do.call(rbind, lapply(names(scales), function(nm) {
    thr <- scales[[nm]]; gg <- rci_grid(thr)
    ref <- exact_crit(thr, rnorm(4000, 0, 1))        # calibrated once, reused
    Ps <- ct_pmat(thr, FALSE); Pm <- ct_pmat(thr, TRUE)
    do.call(rbind, parallel::mclapply(seq_along(dists), function(i) {
      mu <- dists[[i]][1]; sg <- dists[[i]][2]
      r <- do.call(rbind, lapply(seq_len(REPS), function(rep) {
        th <- rnorm(N, mu, sg); p <- sim_pair(thr, th, th, d = 0)
        r1 <- rowSums(p$t1); r2 <- rowSums(p$t2)
        th1 <- gg$theta[r1 + 1L]; s1 <- gg$se[r1 + 1L]
        th2 <- gg$theta[r2 + 1L]; s2 <- gg$se[r2 + 1L]
        rci <- (th2 - th1) / sqrt(s1^2 + s2^2)
        own <- exact_crit(thr, th)
        tn <- (th1 / s1^2 + th2 / s2^2) / (1 / s1^2 + 1 / s2^2)
        pc <- plugin_crit(thr, tn, alpha)
        c(pooled_own = mean(rci < own[["lower"]] | rci > own[["upper"]]),
          pooled_reused = mean(rci < ref[["lower"]] | rci > ref[["upper"]]),
          cond_crit = mean(rci < pc[, 1] | rci > pc[, 2]),
          ct_strict = mean(Ps[cbind(r1 + 1L, r2 + 1L)] <= alpha),
          ct_midp = mean(Pm[cbind(r1 + 1L, r2 + 1L)] <= alpha))
      }))
      own <- exact_crit(thr, rnorm(4000, mu, sg))
      data.frame(scale = nm, dist = dlab[i], crit_own = own[["upper"]],
                 t(colMeans(r)), stringsAsFactors = FALSE)
    }, mc.cores = cores))
  }))
}

## exact, no simulation: size or power for one respondent at a known location
.one_subject <- function(thr, t0, delta = 0, alpha = 0.05) {
  g <- rci_grid(thr); R <- g$R
  Ps <- ct_pmat(thr, FALSE); Pm <- ct_pmat(thr, TRUE)
  tn <- outer(seq_len(R + 1L), seq_len(R + 1L), function(i, j)
    (g$theta[i] / g$se[i]^2 + g$theta[j] / g$se[j]^2) /
      (1 / g$se[i]^2 + 1 / g$se[j]^2))
  pc <- plugin_crit(thr, as.numeric(tn), alpha)
  fpi <- matrix(as.numeric(g$rci) < pc[, 1] | as.numeric(g$rci) > pc[, 2],
                R + 1L, R + 1L)
  W <- outer(score_dist(thr, t0), score_dist(thr, t0 + delta))
  c(plug_in = sum(W[fpi]), strict = sum(W[Ps <= alpha]), midp = sum(W[Pm <= alpha]))
}

single_subject <- function(thr = mk(6, 4), thetas = seq(-3, 3, by = 1)) {
  out <- t(vapply(thetas, function(t0) .one_subject(thr, t0, 0), numeric(3)))
  rownames(out) <- sprintf("%+.0f", thetas); out
}

power_table <- function(thr = mk(6, 4), deltas = c(0, 0.5, 1, 1.5, 2)) {
  out <- t(vapply(deltas, function(dl) .one_subject(thr, 0, dl), numeric(3)))
  rownames(out) <- paste("change", deltas); out
}
