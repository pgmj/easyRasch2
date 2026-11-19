# Step 1 of RMdifImpact (dev/dif-impact-design.md, Part B, restarted 2026-11-17).
#
# Population (large-sample) values of the two-group comparison under three
# ways of handling one item with uniform DIF: keep it (pooled model), remove
# it, or split it (uniform-shift model, the item gets a focal-group shift).
# Three estimators of the group difference: mean WLE, latent mean (CML item
# parameters fixed, MML group moments) and raw sum score (under split, the
# focal group's sum score is true-score equated to the reference items).
#
# Nothing is simulated. Every fit is made to EXPECTED sufficient statistics
# (category counts and raw-score distribution per group, from quadrature over
# the true latent distribution), which gives the probability limit of each
# estimate. Standard errors at given group sizes combine
#   sampling: the within-group variance of the estimator over persons / n_g
#   param:    delta method on the CML parameters (vcov from the expected
#             information at n_ref, n_focal)
# assuming the two are independent (CML estimates condition on the raw
# scores). That assumption is checked by simulation in step 2.
#
# Run from the package root. Full run (216 general + 42 scenario settings)
# on the host, about 9 min on 10 workers:
#   cd R_packages/easyRasch2
#   DIF_THEORY_WORKERS=10 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
#     nohup Rscript dev/dif_impact_theory.R > dev/dif_impact_theory/run.log 2>&1 &
# Smoke test (7 scenario settings, container): DIF_THEORY_SMOKE=1, writes theory_smoke.rds.
# Results: dev/dif_impact_theory/theory_results.rds (+ dated copy).

suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE))
source("dev/dif_shift_cml.R")

out_dir <- "dev/dif_impact_theory"
dir.create(out_dir, showWarnings = FALSE)
n_workers <- as.integer(Sys.getenv("DIF_THEORY_WORKERS", "4"))

# ---------------------------------------------------------------------------
# Probability machinery on a theta grid
# ---------------------------------------------------------------------------

cat_probs_grid <- function(thr, grid) {
  m <- length(thr)
  eta <- outer(grid, 0:m) -
    matrix(c(0, cumsum(thr)), length(grid), m + 1L, byrow = TRUE)
  # No max-subtraction: with |theta| <= 15 and thresholds within +-3 the
  # exponents stay far below overflow, and this is the hot loop.
  p <- exp(eta)
  p / rowSums(p)
}

# P(raw score = r | theta) for r = 0..R, by convolution over items.
score_probs_grid <- function(thr_list, grid) {
  out <- matrix(1, length(grid), 1L)
  for (thr in thr_list) {
    pc <- cat_probs_grid(thr, grid)
    m <- ncol(pc) - 1L
    w <- ncol(out)
    new <- matrix(0, length(grid), w + m)
    for (c in 0:m) {
      idx <- (c + 1L):(c + w)
      new[, idx] <- new[, idx] + out * pc[, c + 1L]
    }
    out <- new
  }
  out
}

quad <- function(mu, sigma, n_nodes = 201L) {
  g <- seq(mu - 7 * sigma, mu + 7 * sigma, length.out = n_nodes)
  w <- stats::dnorm(g, mu, sigma)
  list(grid = g, w = w / sum(w))
}

# Sum over items of expected score, information and third central moment.
item_moments <- function(thr_list, theta) {
  E <- I <- J <- 0
  for (thr in thr_list) {
    p <- cat_probs_grid(thr, theta)
    cats <- 0:length(thr)
    e <- drop(p %*% cats)
    dv <- outer(-e, cats, "+")
    E <- E + e
    I <- I + rowSums(p * dv^2)
    J <- J + rowSums(p * dv^3)
  }
  list(E = E, I = I, J = J)
}

# Expected sufficient statistics per person in a group with N(mu, sigma).
pop_stats <- function(thr_list, mu, sigma) {
  q <- quad(mu, sigma)
  pr <- colSums(q$w * score_probs_grid(thr_list, q$grid))
  ctot <- unlist(lapply(thr_list, function(thr) {
    colSums(q$w * cat_probs_grid(thr, q$grid))[-1L]
  }))
  list(pr = pr, ctot = ctot)
}

# ---------------------------------------------------------------------------
# CML on expected counts (reuses the likelihood in dif_shift_cml.R)
# ---------------------------------------------------------------------------

make_prep <- function(m, dif_items, stats_list, n) {
  P <- sum(m)
  item_of <- rep.int(seq_along(m), m)
  cat_of <- unlist(lapply(m, seq_len))
  S <- matrix(0, P, length(dif_items))
  for (q in seq_along(dif_items)) {
    pos <- which(item_of == dif_items[q])
    S[pos, q] <- cat_of[pos]
  }
  stats <- Map(function(st, nn) list(ctot = nn * st$ctot, ptot = nn * st$pr),
               stats_list, n)
  list(m = m, P = P, item_of = item_of, cat_of = cat_of,
       dif_items = as.integer(dif_items), S = S, stats = stats)
}

fit_prep <- function(prep) {
  nd <- length(prep$dif_items)
  start <- rep(0, prep$P - 1L + nd)
  opt <- stats::optim(start, .shift_negll, .shift_grad, prep = prep,
                      method = "BFGS",
                      control = list(reltol = 1e-15, maxit = 10000L))
  H <- stats::optimHess(opt$par, .shift_negll, .shift_grad, prep = prep)
  list(par = opt$par, vcov = solve(H), conv = opt$convergence,
       grad_max = max(abs(.shift_grad(opt$par, prep))) /
         sum(vapply(prep$stats, function(s) sum(s$ptot), 0)))
}

# Andrich thresholds for both groups from the parameter vector, centred on the
# reference group's grand mean (as shift_cml() does).
par_to_thr <- function(par, prep, item_names) {
  u <- .shift_unpack(par, prep)
  thr <- lapply(split(u$psi, prep$item_of), function(p) diff(c(0, p)))
  shift <- mean(unlist(thr))
  thr <- lapply(thr, function(x) x - shift)
  names(thr) <- item_names
  thr2 <- thr
  for (q in seq_along(prep$dif_items)) {
    j <- prep$dif_items[q]
    thr2[[j]] <- thr2[[j]] + u$d[q]
  }
  list(ref = thr, focal = thr2)
}

# ---------------------------------------------------------------------------
# Estimators by raw score
# ---------------------------------------------------------------------------

# WLE for every raw score at once, by vectorised bisection on Warm's score
# function (same function as .theta_wle(), search range [-10, 10]).
wle_by_score <- function(thr_list, lo = -10, hi = 10, iter = 40L) {
  r <- 0:sum(lengths(thr_list))
  a <- rep(lo, length(r)); b <- rep(hi, length(r))
  for (it in seq_len(iter)) {
    mid <- (a + b) / 2
    mo <- item_moments(thr_list, mid)
    s <- r - mo$E + mo$J / (2 * mo$I)
    up <- s > 0
    a[up] <- mid[up]; b[!up] <- mid[!up]
  }
  (a + b) / 2
}

# theta with expected score r (r strictly inside 0..R).
inv_true_score <- function(thr_list, r, lo = -15, hi = 15, iter = 40L) {
  a <- rep(lo, length(r)); b <- rep(hi, length(r))
  for (it in seq_len(iter)) {
    mid <- (a + b) / 2
    up <- item_moments(thr_list, mid)$E < r
    a[up] <- mid[up]; b[!up] <- mid[!up]
  }
  (a + b) / 2
}

# True-score equating of focal raw scores onto the reference items:
# T_ref(T_focal^{-1}(r)); 0 and R map to themselves.
equate_focal <- function(thr_ref, thr_focal) {
  R <- sum(lengths(thr_ref))
  r <- 0:R
  out <- r
  inner <- r > 0 & r < R
  th <- inv_true_score(thr_focal, r[inner])
  out[inner] <- item_moments(thr_ref, th)$E
  out
}

MML_GRID <- seq(-10, 10, length.out = 401L)

# Log-likelihood of raw score r at each grid node, up to a pattern constant
# that does not involve theta: r * theta - sum_i log Z_i(theta).
mml_L <- function(thr_list, grid = MML_GRID) {
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

mml_negll <- function(par, L, pr, grid = MML_GRID) {
  lp <- stats::dnorm(grid, par[1L], exp(par[2L]), log = TRUE)
  lmx <- max(lp)
  lp <- lp - (lmx + log(sum(exp(lp - lmx))))
  M <- sweep(L, 2L, lp, "+")
  mx <- do.call(pmax, as.data.frame(M))  # row maxima (one per raw score)
  keep <- pr > 0
  -sum(pr[keep] * (mx + log(rowSums(exp(M - mx))))[keep])
}

# Latent mean and SD by MML with item parameters fixed; per-person
# information for the SE of the mean.
mml_fit <- function(thr_list, pr, start = c(0, 0), hessian = TRUE) {
  L <- mml_L(thr_list)
  opt <- stats::optim(start, mml_negll, L = L, pr = pr, method = "BFGS",
                      control = list(reltol = 1e-15, maxit = 2000L))
  var_mu <- NA_real_
  if (hessian) {
    H <- stats::optimHess(opt$par, mml_negll, L = L, pr = pr)
    var_mu <- solve(H)[1L, 1L]
  }
  list(mu = opt$par[1L], sd = exp(opt$par[2L]), par = opt$par,
       var_mu = var_mu)
}

# ---------------------------------------------------------------------------
# Group comparison under one strategy
# ---------------------------------------------------------------------------

# thr: list(ref, focal) used for scoring; pr: list(ref, focal) true raw-score
# distributions on the scored items; equate: equate the focal sum score.
group_comparison <- function(thr, pr, equate = FALSE, mml_start = NULL,
                             hessian = TRUE) {
  out <- list()
  # WLE
  w <- lapply(thr, wle_by_score)
  m <- mapply(function(wv, p) sum(p * wv), w, pr)
  v <- mapply(function(wv, p, mm) sum(p * (wv - mm)^2), w, pr, m)
  out$wle <- c(delta = m[[2]] - m[[1]], var_ref = v[[1]], var_focal = v[[2]],
               smd = (m[[2]] - m[[1]]) / sqrt(mean(v)),
               sd_pooled = sqrt(mean(v)))
  # latent
  if (is.null(mml_start)) mml_start <- list(c(0, 0), c(0, 0))
  lat <- Map(mml_fit, thr, pr, mml_start, MoreArgs = list(hessian = hessian))
  out$latent <- c(delta = lat[[2]]$mu - lat[[1]]$mu,
                  var_ref = lat[[1]]$var_mu, var_focal = lat[[2]]$var_mu,
                  smd = (lat[[2]]$mu - lat[[1]]$mu) /
                    sqrt((lat[[1]]$sd^2 + lat[[2]]$sd^2) / 2),
                  sd_pooled = sqrt((lat[[1]]$sd^2 + lat[[2]]$sd^2) / 2))
  # sum score
  R <- length(pr[[1]]) - 1L
  sv <- list(0:R, if (equate) equate_focal(thr$ref, thr$focal) else 0:R)
  ms <- mapply(function(s, p) sum(p * s), sv, pr)
  vs <- mapply(function(s, p, mm) sum(p * (s - mm)^2), sv, pr, ms)
  out$sum <- c(delta = ms[[2]] - ms[[1]], var_ref = vs[[1]],
               var_focal = vs[[2]], smd = (ms[[2]] - ms[[1]]) / sqrt(mean(vs)),
               sd_pooled = sqrt(mean(vs)))
  attr(out, "mml_par") <- lapply(lat, `[[`, "par")
  out
}

# Delta method: SE of Delta from item-parameter uncertainty. The estimator's
# Delta is recomputed at perturbed CML parameters with the raw-score
# distributions held fixed. Also returns the part due to the shift(s) alone.
param_se <- function(fit, prep, item_names, pr, equate, base, h = 1e-4) {
  np <- length(fit$par)
  f <- function(par) {
    gc <- group_comparison(par_to_thr(par, prep, item_names), pr, equate,
                           mml_start = attr(base, "mml_par"), hessian = FALSE)
    c(wle = gc$wle[["delta"]], latent = gc$latent[["delta"]],
      sum = gc$sum[["delta"]])
  }
  G <- matrix(0, 3L, np, dimnames = list(c("wle", "latent", "sum"), NULL))
  for (p in seq_len(np)) {
    e <- numeric(np); e[p] <- h
    G[, p] <- (f(fit$par + e) - f(fit$par - e)) / (2 * h)
  }
  total <- sqrt(pmax(diag(G %*% fit$vcov %*% t(G)), 0))
  nd <- length(prep$dif_items)
  shift_only <- if (nd) {
    idx <- prep$P - 1L + seq_len(nd)
    Gd <- G[, idx, drop = FALSE]
    sqrt(pmax(diag(Gd %*% fit$vcov[idx, idx, drop = FALSE] %*% t(Gd)), 0))
  } else c(wle = 0, latent = 0, sum = 0)
  # d Delta / d shift: how a unit error in the shift moves Delta.
  dshift <- if (nd) G[, prep$P, drop = TRUE] else c(wle = 0, latent = 0, sum = 0)
  list(total = total, shift_only = shift_only, dshift = dshift)
}

# ---------------------------------------------------------------------------
# One setting: all strategies at DIF d and at d = 0
# ---------------------------------------------------------------------------

make_scale <- function(type, k) {
  if (type == "dich") {
    loc <- seq(-2, 2, length.out = k)
    thr <- as.list(loc)
  } else {
    loc <- seq(-1.5, 1.5, length.out = k)
    thr <- lapply(loc, function(l) l + c(-1, 0, 1))
  }
  names(thr) <- sprintf("I%02d", seq_len(k))
  list(thr = thr, loc = loc)
}

pick_item <- function(loc, position) {
  if (position == "central") which.min(abs(loc)) else which.max(loc)
}

# DIF items for a scenario: item indices and the sign of each item's shift.
#   one:          one item, central (nearest 0), near_cut (nearest the cut) or
#                 offtarget (largest location)
#   two_same:     two items shifted the same way, central or near_cut
#   two_opposite: one item near the cut shifted +d, one near -cut shifted -d
#                 (cancels in the group mean, not locally)
scenario_items <- function(loc, scenario, position, cut) {
  near <- function(target, n = 1L) order(abs(loc - target))[seq_len(n)]
  switch(scenario,
    one = list(j = switch(position, central = near(0), near_cut = near(cut),
                          offtarget = which.max(loc)),
               sign = 1),
    two_same = list(j = if (position == "central") near(0, 2L) else near(cut, 2L),
                    sign = c(1, 1)),
    two_opposite = list(j = c(near(cut), near(-cut)), sign = c(1, -1)),
    stop("Unknown scenario: ", scenario, call. = FALSE)
  )
}

# Local between-group bias caused by the DIF items, from the true (or
# estimated split-model) parameters of the two groups.
#   b(theta)   = T_ref^{-1}(T_focal(theta)) - theta: where a focal person at
#                theta lands when scored with the reference parameters, i.e.
#                the horizontal distance between the groups' test
#                characteristic curves. Negative = the focal person is placed
#                too low. It is the between-group distance at equal theta and
#                does not depend on how the pooled fit splits the shift
#                between the groups (to first order).
#   gap(theta) = T_focal(theta) - T_ref(theta): expected raw-score difference
#                at equal theta (vertical distance).
# Maxima are taken over the central `q` range of the pooled latent
# distribution (mixture of the two groups), because the share of an
# off-target item keeps rising into empty tails.
local_bias <- function(thr_ref, thr_focal, mu, sd, w_focal, cut,
                       q = c(0.025, 0.975), n_grid = 401L) {
  pmix <- function(x) {
    (1 - w_focal) * stats::pnorm(x, mu[1], sd[1]) +
      w_focal * stats::pnorm(x, mu[2], sd[2])
  }
  rng <- vapply(q, function(p) {
    stats::uniroot(function(x) pmix(x) - p, c(-15, 15), tol = 1e-10)$root
  }, 0)
  grid <- seq(rng[1], rng[2], length.out = n_grid)
  at <- function(th) {
    Tf <- item_moments(thr_focal, th)$E
    Tr <- item_moments(thr_ref, th)$E
    b <- inv_true_score(thr_ref, Tf) - th
    sem <- 1 / sqrt(item_moments(thr_ref, th)$I)
    list(b = b, gap = Tf - Tr, b_sem = b / sem)
  }
  g <- at(grid)
  ib <- which.max(abs(g$b)); is <- which.max(abs(g$b_sem))
  ig <- which.max(abs(g$gap))
  ct <- at(cut)
  # Signed values at the argmax are kept for reading, but with DIF in both
  # directions |b_lo| and |b_hi| can be similar and the sign at the argmax
  # flips between samples; the unsigned *_maxabs are the summaries to use.
  c(range_lo = rng[1], range_hi = rng[2],
    b_lo = min(g$b), b_hi = max(g$b),
    b_maxabs = abs(g$b[ib]), b_sem_maxabs = abs(g$b_sem[is]),
    gap_maxabs = abs(g$gap[ig]),
    b_absmax = g$b[ib], theta_b_absmax = grid[ib],
    b_sem_absmax = g$b_sem[is], theta_b_sem_absmax = grid[is],
    gap_absmax = g$gap[ig], theta_gap_absmax = grid[ig],
    b_cut = ct$b, b_sem_cut = ct$b_sem, gap_cut = ct$gap)
}

# Classification at a cut on the theta scale (WLE >= cut), keep versus split.
# Raw scores are discrete, so a shift either moves the cut past a raw score
# or does not. Returns the share of each group whose classification differs,
# and the raw-score cut under each scoring.
cut_flip <- function(thr_keep, thr_split, pr, cut) {
  w_keep <- wle_by_score(thr_keep)
  out <- c()
  for (gi in 1:2) {
    w_split <- wle_by_score(thr_split[[gi]])
    flip <- (w_keep >= cut) != (w_split >= cut)
    lab <- c("ref", "focal")[gi]
    out[paste0("flip_prop_", lab)] <- sum(pr[[gi]][flip])
    out[paste0("raw_cut_split_", lab)] <- which(w_split >= cut)[1] - 1L
  }
  out["raw_cut_keep"] <- which(w_keep >= cut)[1] - 1L
  out
}

# j: DIF item indices; d: shift per item (same length; zeros give the no-DIF
# counterfactual with the same items removed or split).
strategies_at <- function(thr_true, j, d, mu2, n) {
  thr_focal <- thr_true
  for (q in seq_along(j)) thr_focal[[j[q]]] <- thr_focal[[j[q]]] + d[q]
  st <- list(pop_stats(thr_true, 0, 1), pop_stats(thr_focal, mu2, 1))
  red <- thr_true[-j]
  st_red <- list(pop_stats(red, 0, 1), pop_stats(thr_focal[-j], mu2, 1))
  m <- lengths(thr_true)
  nm <- names(thr_true)
  defs <- list(
    keep = list(prep = make_prep(m, integer(0), st, n), names = nm,
                pr = list(st[[1]]$pr, st[[2]]$pr), equate = FALSE),
    remove = list(prep = make_prep(m[-j], integer(0), st_red, n),
                  names = nm[-j],
                  pr = list(st_red[[1]]$pr, st_red[[2]]$pr), equate = FALSE),
    split = list(prep = make_prep(m, j, st, n), names = nm,
                 pr = list(st[[1]]$pr, st[[2]]$pr), equate = TRUE)
  )
  lapply(defs, function(df) {
    fit <- fit_prep(df$prep)
    thr <- par_to_thr(fit$par, df$prep, df$names)
    gc <- group_comparison(thr, df$pr, df$equate)
    ps <- param_se(fit, df$prep, df$names, df$pr, df$equate, gc)
    list(fit = fit, thr = thr, gc = gc, ps = ps,
         d_hat = if (length(df$prep$dif_items)) .shift_unpack(fit$par, df$prep)$d
                 else NA_real_)
  })
}

# First-order approximation of the keep bias for Delta (logits):
# reference persons shift by +c d I_j/I, focal by -(1-c) d I_j/I, where c is
# the share of d absorbed by the pooled location of item j.
# Several DIF items add up: the sum over items of each item's term. The
# returned c and shares are those of the first DIF item.
first_order_keep <- function(thr_true, j, d, mu2, keep_thr) {
  loc_true <- vapply(thr_true, mean, 0)
  loc_pool <- vapply(keep_thr, mean, 0)
  common <- setdiff(seq_along(thr_true), j)
  q0 <- quad(0, 1); q1 <- quad(mu2, 1)
  It0 <- item_moments(thr_true, q0$grid)$I
  It1 <- item_moments(thr_true, q1$grid)$I
  terms <- lapply(seq_along(j), function(i) {
    rel <- function(l) l[j[i]] - mean(l[common])
    c_abs <- unname((rel(loc_pool) - rel(loc_true)) / d[i])
    s_ref <- sum(q0$w * item_moments(thr_true[j[i]], q0$grid)$I / It0)
    s_foc <- sum(q1$w * item_moments(thr_true[j[i]], q1$grid)$I / It1)
    c(c_abs = c_abs, s_ref = s_ref, s_foc = s_foc,
      bias = -(1 - c_abs) * d[i] * s_foc - c_abs * d[i] * s_ref)
  })
  c(c_absorbed = terms[[1]][["c_abs"]], share_ref = terms[[1]][["s_ref"]],
    share_focal = terms[[1]][["s_foc"]],
    bias_first_order = sum(vapply(terms, function(t) t[["bias"]], 0)))
}

run_setting <- function(s) {
  if (is.null(s$scenario)) s$scenario <- "one"
  if (is.null(s$cut)) s$cut <- 1
  sc <- make_scale(s$type, s$k)
  si <- scenario_items(sc$loc, s$scenario, s$position, s$cut)
  j <- si$j
  dvec <- s$d * si$sign
  n <- c(s$n_ref, s$n_focal)
  at_d <- strategies_at(sc$thr, j, dvec, s$mu2, n)
  at_0 <- strategies_at(sc$thr, j, 0 * dvec, s$mu2, n)
  fo <- first_order_keep(sc$thr, j, dvec, s$mu2, at_d$keep$thr$ref)
  # Local bias from the true parameters, and classification at the cut
  # (population raw-score distributions on all items).
  thr_focal_true <- sc$thr
  for (q in seq_along(j)) thr_focal_true[[j[q]]] <- thr_focal_true[[j[q]]] + dvec[q]
  lb <- local_bias(sc$thr, thr_focal_true, mu = c(0, s$mu2), sd = c(1, 1),
                   w_focal = n[2] / sum(n), cut = s$cut)
  st_all <- list(pop_stats(sc$thr, 0, 1)$pr, pop_stats(thr_focal_true, s$mu2, 1)$pr)
  cf <- cut_flip(at_d$keep$thr$ref, at_d$split$thr, st_all, s$cut)
  fo <- c(fo, lb, cf)
  rows <- list()
  for (strat in names(at_d)) for (est in c("wle", "latent", "sum")) {
    a <- at_d[[strat]]; b <- at_0[[strat]]
    se_samp <- sqrt(a$gc[[est]][["var_ref"]] / n[1] +
                      a$gc[[est]][["var_focal"]] / n[2])
    se_par <- a$ps$total[[est]]
    rows[[length(rows) + 1L]] <- data.frame(
      strategy = strat, estimator = est,
      delta_d = a$gc[[est]][["delta"]], delta_0 = b$gc[[est]][["delta"]],
      impact = a$gc[[est]][["delta"]] - b$gc[[est]][["delta"]],
      smd_d = a$gc[[est]][["smd"]], smd_0 = b$gc[[est]][["smd"]],
      impact_smd = a$gc[[est]][["smd"]] - b$gc[[est]][["smd"]],
      se_samp = se_samp, se_param = se_par,
      se_shift = a$ps$shift_only[[est]], dshift = a$ps$dshift[[est]],
      se_total = sqrt(se_samp^2 + se_par^2),
      sd_pooled = a$gc[[est]][["sd_pooled"]],
      d_hat = a$d_hat[1],
      se_d_hat = if (strat == "split") sqrt(diag(a$fit$vcov))[length(a$fit$par) - length(j) + 1L] else NA_real_,
      conv = a$fit$conv, grad_max = a$fit$grad_max
    )
  }
  out <- do.call(rbind, rows)
  out$rmse_dif <- sqrt(out$impact^2 + out$se_total^2)
  # SE on the SMD scale (estimator's own pooled SD), so strategies with k and
  # k - 1 items and the sum score compare directly.
  out$se_total_smd <- out$se_total / out$sd_pooled
  # Consequence at these group sizes: bias in SE units and the coverage of a
  # nominal 95% Wald CI for the group difference.
  out$bias_in_se <- out$impact / out$se_total
  out$ci95_coverage <- stats::pnorm(1.96 - out$bias_in_se) -
    stats::pnorm(-1.96 - out$bias_in_se)
  cbind(as.data.frame(s)[rep(1L, nrow(out)), , drop = FALSE],
        item = paste(names(sc$thr)[j], collapse = ","),
        item_loc = paste(round(sc$loc[j], 3), collapse = ","),
        dif = paste(dvec, collapse = ","),
        out, t(fo)[rep(1L, nrow(out)), , drop = FALSE], row.names = NULL)
}

# ---------------------------------------------------------------------------
# Checks
# ---------------------------------------------------------------------------

run_checks <- function() {
  chk <- list()
  sc <- make_scale("poly", 6)
  # 1. Vectorised WLE equals the package .theta_wle().
  w_vec <- wle_by_score(sc$thr)
  R <- sum(lengths(sc$thr))
  w_pkg <- vapply(0:R, function(r) {
    resp <- c(r, rep(0, length(sc$thr) - 1L))
    .theta_wle(resp, sc$thr, c(-10, 10), tol = 1e-10)[["theta"]]
  }, 0)
  chk$wle_max_abs_diff <- max(abs(w_vec - w_pkg))
  # 2. Population shift CML recovers the truth (d = 0.7, focal mean 0.5).
  j <- 3L
  thr_f <- sc$thr; thr_f[[j]] <- thr_f[[j]] + 0.7
  st <- list(pop_stats(sc$thr, 0, 1), pop_stats(thr_f, 0.5, 1))
  prep <- make_prep(lengths(sc$thr), j, st, c(1000, 1000))
  fit <- fit_prep(prep)
  thr <- par_to_thr(fit$par, prep, names(sc$thr))
  chk$shift_d_hat <- .shift_unpack(fit$par, prep)$d
  chk$shift_thr_max_abs_err <- max(abs(unlist(thr$ref) - unlist(sc$thr)))
  # 3. MML with the true parameters recovers mean and SD.
  mm <- mml_fit(thr_f, st[[2]]$pr)
  chk$mml_mu <- mm$mu; chk$mml_sd <- mm$sd
  # 4. Equating with identical parameters is the identity.
  chk$equate_identity_max <- max(abs(equate_focal(sc$thr, sc$thr) - 0:R))
  # 5. One large simulated dataset agrees with the population keep values
  #    (dichotomous, k = 10, central item, d = 0.5, mu2 = 0.5, 100k + 100k).
  set.seed(20261117)
  sd10 <- make_scale("dich", 10)
  jj <- pick_item(sd10$loc, "central")
  sim <- sim_dif(1e5, 1e5, sd10$thr, dif = stats::setNames(0.5, names(sd10$thr)[jj]),
                 mu2 = 0.5)
  pooled <- .fit_cml_thresholds(sim$X)
  wtab <- wle_by_score(pooled)
  wle <- wtab[rowSums(sim$X) + 1L]
  chk$sim_keep_delta_wle <- mean(wle[sim$g == "g2"]) - mean(wle[sim$g == "g1"])
  pop <- strategies_at(sd10$thr, jj, 0.5, 0.5, c(1e5, 1e5))
  chk$pop_keep_delta_wle <- pop$keep$gc$wle[["delta"]]
  chk$pop_keep_se_samp <- sqrt(pop$keep$gc$wle[["var_ref"]] / 1e5 +
                                 pop$keep$gc$wle[["var_focal"]] / 1e5)
  chk$sim_split_d_hat <- shift_cml(sim$X, sim$g, jj, hessian = FALSE)$d
  chk
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

# Grid A: the general grid (one item, dichotomous and polytomous, k up to
# 20). First run 2026-11-17 with d = 0.25, 0.5, 1 and no local-bias columns;
# d changed to 0.5, 0.8, 1.2 the same day (user).
grid_general <- function() {
  grid <- expand.grid(
    type = c("dich", "poly"), k = c(5L, 10L, 20L),
    position = c("central", "offtarget"), d = c(0.5, 0.8, 1.2),
    mu2 = c(0, 0.5),
    nset = 1:3, stringsAsFactors = FALSE
  )
  nsets <- data.frame(nset = 1:3, n_ref = c(500, 1000, 2000),
                      n_focal = c(500, 200, 100))
  grid <- merge(grid, nsets, by = "nset")
  grid$nset <- NULL
  grid$scenario <- "one"
  grid$cut <- 1
  grid$cell_set <- "general"
  grid
}

# Grid B: the scenarios agreed for the simulation (2026-11-17). Polytomous,
# 4 categories, cut at theta = 1, focal mean 0.5. Shared with
# dif_impact_sim.R so every simulated cell has an exact population target.
grid_scenarios <- function() {
  configs <- data.frame(
    k = c(5L, 5L, 10L, 10L, 10L, 10L, 10L),
    scenario = c("one", "one", "one", "one", "two_same", "two_same",
                 "two_opposite"),
    position = c("central", "near_cut", "central", "near_cut", "central",
                 "near_cut", "near_cut"),
    stringsAsFactors = FALSE
  )
  # Sample sizes (user, 2026-11-17): 250/250 (small study), 500/500 and
  # 800/200 (same total, 4:1 imbalance). The 2500/2500 "detected but small"
  # cells were dropped: step 1 gives that message exactly.
  main <- merge(configs, expand.grid(d = c(0.5, 0.8), nset = 1:3))
  nsets <- data.frame(nset = 1:3, n_ref = c(250, 500, 800),
                      n_focal = c(250, 500, 200))
  main <- merge(main, nsets, by = "nset")
  main$nset <- NULL
  main$cell_set <- "main"
  out <- main
  out$type <- "poly"; out$mu2 <- 0.5; out$cut <- 1
  out[order(out$k, out$scenario, out$position, out$d, out$n_ref), ]
}

if (sys.nframe() == 0L || identical(Sys.getenv("DIF_THEORY_RUN"), "1")) {
  t0 <- Sys.time()
  checks <- run_checks()
  print(str(checks))

  cols <- c("type", "k", "scenario", "position", "d", "mu2", "n_ref",
            "n_focal", "cut", "cell_set")
  grid <- rbind(grid_general()[, cols], grid_scenarios()[, cols])
  smoke <- identical(Sys.getenv("DIF_THEORY_SMOKE"), "1")
  if (smoke) {
    grid <- grid[grid$cell_set == "main" & grid$d == 0.5 & grid$n_focal == 200, ]
  }
  rownames(grid) <- NULL
  settings <- split(grid, seq_len(nrow(grid)))

  res <- parallel::mclapply(settings, function(s) {
    tryCatch(run_setting(as.list(s)), error = function(e) {
      cbind(s, error = conditionMessage(e))
    })
  }, mc.cores = n_workers)
  ok <- vapply(res, function(x) !"error" %in% names(x), TRUE)
  if (any(!ok)) {
    message(sum(!ok), " settings failed")
    print(do.call(rbind, res[!ok]))
  }
  results <- do.call(rbind, res[ok])
  out <- list(results = results, checks = checks, grid = grid,
              run_time = Sys.time() - t0, date = Sys.Date(),
              session = utils::sessionInfo())
  stem <- if (smoke) "theory_smoke" else "theory_results"
  saveRDS(out, file.path(out_dir, paste0(stem, ".rds")))
  saveRDS(out, file.path(out_dir, sprintf("%s_%s.rds", stem,
                                          format(Sys.Date(), "%Y-%m-%d"))))
  message("Done in ", format(out$run_time))
}
