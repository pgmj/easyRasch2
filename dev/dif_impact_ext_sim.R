# Step 3 of RMdifImpact (dev/dif-impact-design.md, Part B, "Decisions after
# the simulation", 2026-11-17): validation of the three outputs added after
# the step 2 simulation, before they go into the package function.
#
#   1. Per-item impacts (two-item cells): the impact of keeping item j while
#      the other DIF items stay split, Delta(split S \ j) - Delta(split S),
#      with a delta-method SE on the shift of item j.
#   2. Classification at a cut, on the theta scale (WLE >= cut) and on the
#      sum-score scale (raw >= cut, focal raw scores equated under split).
#      Per group and strategy (keep, split): the share classified at or
#      above the cut (prev), false negatives (true status at or above the
#      cut, classified below) and false positives, as shares of the group,
#      and keep - split. Model-based: response model and latent N(mu, sd)
#      from the split model. Classification by raw score is a step function
#      of the parameters, so the delta method does not apply. Intervals are
#      percentiles of B draws of all split-model parameters (thresholds and
#      shifts, from the CML covariance matrix) and of each group's
#      (mu, log sd), linearised in the item parameters plus their own
#      sampling error. The keep WLE table follows the reference table and
#      is linearised in the shifts through population keep fits.
#   3. Projection to a target population (user-supplied n, means and SDs
#      per group): the keep impact on Delta (population CML fit to the
#      expected statistics under the estimated split model), the SE of Delta
#      in a study of the target size, bias / SE, and the rates at the cut.
#      Impact SE by the delta method on the shifts, rates by the same draws
#      with the target distributions held fixed.
#
# Oracle DIF items only (the flagged-set behaviour is known from step 2),
# DIF data only. Cells: grid_scenarios() (42 cells), cut at theta = 1, the
# sum-score cut is the raw score nearest the reference expected score at
# theta = 1. Projection target: n = 40 000 / 10 000, reference N(0, 1),
# focal N(0.25, 1.25^2), different from every simulated cell.
#
# Run from the package root. Smoke test (container, 4 workers):
#   DIF_EXT_SMOKE=1 LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript dev/dif_impact_ext_sim.R
# Full run on the host (10 workers, about 3-5 h; restarts skip finished
# cells, which are saved one file per cell in ext_results_cells/):
#   DIF_EXT_WORKERS=10 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
#     nohup Rscript dev/dif_impact_ext_sim.R > dev/dif_impact_ext_sim/run.log 2>&1 &

source("dev/dif_impact_sim.R")  # theory + step 2 helpers; run blocks do not fire

ext_dir <- "dev/dif_impact_ext_sim"
dir.create(ext_dir, showWarnings = FALSE)

N_DRAWS <- 300L
TARGET <- list(n = c(40000, 10000), mu = c(0, 0.25), sd = c(1, 1.25))

# Theta grid for the classification integrals and for reading WLEs and the
# true-score equating off the item moments (linear interpolation, checked
# against wle_by_score() and equate_focal() in the smoke test). The
# true-status indicator is interpolated across the grid cell that holds the
# cut, so the rates are smooth in mu and sd.
TH_H <- 0.02
TH_GRID <- seq(-10, 10, by = TH_H)

# ---------------------------------------------------------------------------
# Classification at a cut
# ---------------------------------------------------------------------------

sp_conv <- function(SP, pc) {
  m <- ncol(pc) - 1L
  w <- ncol(SP)
  new <- matrix(0, nrow(SP), w + m)
  for (c in 0:m) {
    idx <- (c + 1L):(c + w)
    new[, idx] <- new[, idx] + SP * pc[, c + 1L]
  }
  new
}

madd <- function(a, b) list(E = a$E + b$E, I = a$I + b$I, J = a$J + b$J)

wle_from_moments <- function(M, R) {
  stats::approx(M$E - M$J / (2 * M$I), TH_GRID, xout = 0:R, rule = 2)$y
}

# Everything the rates need from one set of split-model thresholds: WLE
# tables per group, equated focal raw scores, score probabilities on the
# grid, and the true-status point for the sum-score cut. The common items
# are shared by both groups, so their moments and score distribution are
# computed once.
rate_parts <- function(thr_split, items, cut_raw) {
  common <- setdiff(seq_along(thr_split$ref), items)
  R <- sum(lengths(thr_split$ref))
  Mc <- item_moments(thr_split$ref[common], TH_GRID)
  Mr <- madd(Mc, item_moments(thr_split$ref[items], TH_GRID))
  Mf <- madd(Mc, item_moments(thr_split$focal[items], TH_GRID))
  SPc <- score_probs_grid(thr_split$ref[common], TH_GRID)
  SPr <- SPc; SPf <- SPc
  for (j in items) {
    SPr <- sp_conv(SPr, cat_probs_grid(thr_split$ref[[j]], TH_GRID))
    SPf <- sp_conv(SPf, cat_probs_grid(thr_split$focal[[j]], TH_GRID))
  }
  r <- 0:R
  inner <- r > 0 & r < R
  eq <- r
  eq[inner] <- stats::approx(TH_GRID, Mr$E,
                             xout = stats::approx(Mf$E, TH_GRID, xout = r[inner])$y)$y
  list(R = R, wle_ref = wle_from_moments(Mr, R), wle_focal = wle_from_moments(Mf, R),
       eq = eq, SP = list(ref = SPr, focal = SPf),
       theta_c_sum = stats::approx(Mr$E, TH_GRID, xout = cut_raw - 0.5)$y)
}

rates_one <- function(SP, cls, mu, sd, theta_c) {
  w <- stats::dnorm(TH_GRID, mu, sd) * TH_H
  above <- pmin(pmax((TH_GRID - theta_c) / TH_H + 0.5, 0), 1)
  pp <- drop(SP %*% as.numeric(cls))
  c(prev = sum(w * pp), fn = sum(w * above * (1 - pp)),
    fp = sum(w * (1 - above) * pp))
}

# wk: keep WLE table on the split metric (offset through the common items
# already added). True status on the sum-score scale: theta at or above
# T_ref^{-1}(cut_raw - 0.5).
cut_rates <- function(pt, wk, mu, sd, cut, cut_raw) {
  th_c <- c(theta = cut, sum = pt$theta_c_sum)
  r <- 0:pt$R
  cls <- list(
    theta = list(keep = list(ref = wk >= cut, focal = wk >= cut),
                 split = list(ref = pt$wle_ref >= cut, focal = pt$wle_focal >= cut)),
    sum = list(keep = list(ref = r >= cut_raw, focal = r >= cut_raw),
               split = list(ref = r >= cut_raw, focal = pt$eq >= cut_raw - 1e-9)))
  out <- c()
  for (sc in names(cls)) for (st in c("keep", "split")) for (gi in 1:2) {
    g <- c("ref", "focal")[gi]
    v <- rates_one(pt$SP[[g]], cls[[sc]][[st]][[g]], mu[gi], sd[gi], th_c[[sc]])
    out[paste(sc, st, g, names(v), sep = ".")] <- v
  }
  for (sc in names(cls)) for (g in c("ref", "focal")) {
    for (m in c("prev", "fn", "fp")) {
      out[paste(sc, "diff", g, m, sep = ".")] <-
        out[[paste(sc, "keep", g, m, sep = ".")]] -
        out[[paste(sc, "split", g, m, sep = ".")]]
    }
  }
  out
}

with_shifts <- function(thr_split, items, d) {
  for (q in seq_along(items)) {
    thr_split$focal[[items[q]]] <- thr_split$ref[[items[q]]] + d[q]
  }
  thr_split
}

# Upper Cholesky root: draws = z %*% chol(V) are reproducible across
# machines (unique root).
mvn_draws <- function(B, mean, V) {
  z <- matrix(stats::rnorm(B * length(mean)), B)
  sweep(z %*% chol(V), 2L, mean, "+")
}

# ---------------------------------------------------------------------------
# Population fits (keep model) and projection to a target population
# ---------------------------------------------------------------------------

fit_prep_nohess <- function(prep) {
  nd <- length(prep$dif_items)
  opt <- stats::optim(rep(0, prep$P - 1L + nd), .shift_negll, .shift_grad,
                      prep = prep, method = "BFGS",
                      control = list(reltol = 1e-15, maxit = 10000L))
  list(par = opt$par, conv = opt$convergence)
}

# Pooled (keep) thresholds that CML converges to when the data follow
# thr_split with latent N(mu, sd) per group and group sizes n.
pop_keep <- function(thr_split, design) {
  st <- list(pop_stats(thr_split$ref, design$mu[1], design$sd[1]),
             pop_stats(thr_split$focal, design$mu[2], design$sd[2]))
  prep <- make_prep(lengths(thr_split$ref), integer(0), st, design$n)
  fk <- fit_prep_nohess(prep)
  list(thr = par_to_thr(fk$par, prep, names(thr_split$ref))$ref,
       pr = list(st[[1]]$pr, st[[2]]$pr), conv = fk$conv)
}

# Keep impact on Delta in the target population, and the SE of Delta
# (split) in a study of the target size.
project <- function(thr_split, target) {
  pk <- pop_keep(thr_split, target)
  gk <- group_comparison(list(ref = pk$thr, focal = pk$thr), pk$pr, FALSE,
                         hessian = FALSE)
  gs <- group_comparison(thr_split, pk$pr, TRUE, hessian = TRUE)
  est <- c("wle", "latent", "sum")
  imp <- vapply(est, function(e) gk[[e]][["delta"]] - gs[[e]][["delta"]], 0)
  se_t <- vapply(est, function(e) {
    sqrt(gs[[e]][["var_ref"]] / target$n[1] + gs[[e]][["var_focal"]] / target$n[2])
  }, 0)
  list(impact = imp, se_target = se_t, thr_keep = pk$thr, conv = pk$conv)
}

# ---------------------------------------------------------------------------
# Estimates from one dataset
# ---------------------------------------------------------------------------

lat_fit <- function(thr, pr, n) {
  f <- mml_fit(thr, pr, hessian = FALSE)
  H <- stats::optimHess(f$par, mml_negll, L = mml_L(thr), pr = pr)
  list(par = f$par, vcov = solve(H) / n)
}

estimate_ext <- function(X, g, items, cut, cut_raw, h = 1e-3) {
  nm <- colnames(X)
  foc <- g == "focal"
  n <- c(sum(!foc), sum(foc))

  # keep and split fits. shift_cml() sorts the DIF items, so items and
  # d_hat are taken from the fit (column order) from here on.
  thr_keep <- .fit_cml_thresholds(X)
  names(thr_keep) <- nm
  keep_sc <- score_strategy(X, g, list(ref = thr_keep, focal = thr_keep), 0,
                            FALSE, NULL)
  sf <- shift_cml(X, g, items, hessian = TRUE)
  if (is.null(sf$vcov)) stop("Singular Hessian in the split fit.", call. = FALSE)
  items <- sf$prep$dif_items
  nd <- length(items)
  common <- setdiff(seq_len(ncol(X)), items)
  P <- sf$prep$P
  idx <- P - 1L + seq_len(nd)
  Vd <- sf$vcov[idx, idx, drop = FALSE]
  d_hat <- unname(sf$d)
  par_hat <- c(sf$psi[-1L], d_hat)
  thr <- list(ref = sf$thresholds$ref, focal = sf$thresholds$focal)
  R <- sum(lengths(thr$ref))
  sp <- score_strategy(X, g, thr, 0, TRUE, NULL)
  unit <- function(q) h * (seq_len(nd) == q)

  # joint impact and the gradient of Delta(split) in the shifts (step 2)
  imp <- keep_sc$group[, "delta"] - sp$group[, "delta"]
  Gi <- matrix(0, 3L, nd)
  for (q in seq_len(nd)) {
    sp_p <- score_strategy(X, g, with_shifts(thr, items, d_hat + unit(q)), 0,
                           TRUE, NULL)
    sp_m <- score_strategy(X, g, with_shifts(thr, items, d_hat - unit(q)), 0,
                           TRUE, NULL)
    Gi[, q] <- -(sp_p$group[, "delta"] - sp_m$group[, "delta"]) / (2 * h)
  }
  se_imp <- sqrt(pmax(diag(Gi %*% Vd %*% t(Gi)), 0))
  rows <- list(data.frame(context = "sample",
                          quantity = paste0("impact_", rownames(keep_sc$group)),
                          estimate = unname(imp), se = se_imp,
                          lo = NA_real_, hi = NA_real_))

  # 1. per-item impacts: keep item q common, the other DIF items split
  if (nd >= 2L) {
    for (q in seq_len(nd)) {
      sq <- shift_cml(X, g, items[-q], hessian = FALSE)
      sc_q <- score_strategy(X, g, list(ref = sq$thresholds$ref,
                                        focal = sq$thresholds$focal), 0, TRUE, NULL)
      rows[[length(rows) + 1L]] <- data.frame(
        context = "sample",
        quantity = paste0("item_", nm[items[q]], "_impact_", rownames(sp$group)),
        estimate = unname(sc_q$group[, "delta"] - sp$group[, "delta"]),
        se = abs(Gi[, q]) * sqrt(Vd[q, q]), lo = NA_real_, hi = NA_real_)
    }
  }

  # 2. classification at the cut, in the sample's own distributions.
  # Latent moments per group, and their Jacobian in the item parameters
  # (forward differences).
  pr <- lapply(c("ref", "focal"), function(lev) {
    tabulate(rowSums(X)[g == lev] + 1L, R + 1L) / sum(g == lev)
  })
  lat <- lapply(1:2, function(gi) lat_fit(thr[[gi]], pr[[gi]], n[gi]))
  ms_hat <- c(lat[[1]]$par, lat[[2]]$par)  # mu_ref, log sd_ref, mu_f, log sd_f
  np <- length(par_hat)
  Jms <- matrix(0, 4L, np)
  for (p in seq_len(np)) {
    tp <- par_to_thr(par_hat + h * (seq_len(np) == p), sf$prep, nm)
    Jms[, p] <- (c(mml_fit(tp$ref, pr[[1]], lat[[1]]$par, hessian = FALSE)$par,
                   mml_fit(tp$focal, pr[[2]], lat[[2]]$par, hessian = FALSE)$par) -
                   ms_hat) / h
  }
  pt_hat <- rate_parts(thr, items, cut_raw)
  wk <- wle_by_score(thr_keep) + link_offset(thr_keep, thr$ref, common)
  rate_est <- cut_rates(pt_hat, wk, ms_hat[c(1, 3)], exp(ms_hat[c(2, 4)]),
                        cut, cut_raw)

  # The keep WLE table moves with the reference table (shared threshold
  # error) and with the shifts (population keep fits at the sample design).
  keep_wle <- function(thr_k) wle_by_score(thr_k) + link_offset(thr_k, thr$ref, common)
  design_s <- list(n = n, mu = ms_hat[c(1, 3)], sd = exp(ms_hat[c(2, 4)]))
  Jw_s <- matrix(0, R + 1L, nd)
  for (q in seq_len(nd)) {
    Jw_s[, q] <- (keep_wle(pop_keep(with_shifts(thr, items, d_hat + unit(q)), design_s)$thr) -
                    keep_wle(pop_keep(with_shifts(thr, items, d_hat - unit(q)), design_s)$thr)) /
      (2 * h)
  }

  D <- mvn_draws(N_DRAWS, par_hat, sf$vcov)
  Ems <- cbind(mvn_draws(N_DRAWS, c(0, 0), lat[[1]]$vcov),
               mvn_draws(N_DRAWS, c(0, 0), lat[[2]]$vcov))
  pts <- lapply(seq_len(N_DRAWS), function(b) {
    rate_parts(par_to_thr(D[b, ], sf$prep, nm), items, cut_raw)
  })
  shift_wk <- function(b) {
    pts[[b]]$wle_ref - pt_hat$wle_ref + drop(Jw_s %*% (D[b, idx] - d_hat))
  }
  draws <- t(vapply(seq_len(N_DRAWS), function(b) {
    ms <- ms_hat + drop(Jms %*% (D[b, ] - par_hat)) + Ems[b, ]
    cut_rates(pts[[b]], wk + shift_wk(b), ms[c(1, 3)], exp(ms[c(2, 4)]),
              cut, cut_raw)
  }, rate_est))
  rows[[length(rows) + 1L]] <- data.frame(
    context = "sample", quantity = names(rate_est), estimate = unname(rate_est),
    se = apply(draws, 2, stats::sd),
    lo = apply(draws, 2, stats::quantile, .025),
    hi = apply(draws, 2, stats::quantile, .975))

  # 3. projection to the target population
  pj <- project(thr, TARGET)
  wk_p <- keep_wle(pj$thr_keep)
  Gp <- matrix(0, 3L, nd)
  Jw_p <- matrix(0, R + 1L, nd)
  for (q in seq_len(nd)) {
    pp <- project(with_shifts(thr, items, d_hat + unit(q)), TARGET)
    pm <- project(with_shifts(thr, items, d_hat - unit(q)), TARGET)
    Gp[, q] <- (pp$impact - pm$impact) / (2 * h)
    Jw_p[, q] <- (keep_wle(pp$thr_keep) - keep_wle(pm$thr_keep)) / (2 * h)
  }
  se_p <- sqrt(pmax(diag(Gp %*% Vd %*% t(Gp)), 0))
  est <- names(pj$impact)
  rows[[length(rows) + 1L]] <- data.frame(
    context = "projection",
    quantity = c(paste0("impact_", est), paste0("bias_in_se_", est),
                 paste0("se_target_", est)),
    estimate = c(pj$impact, pj$impact / pj$se_target, pj$se_target),
    se = c(se_p, se_p / pj$se_target, rep(NA_real_, 3L)),
    lo = NA_real_, hi = NA_real_)
  rate_p <- cut_rates(pt_hat, wk_p, TARGET$mu, TARGET$sd, cut, cut_raw)
  draws_p <- t(vapply(seq_len(N_DRAWS), function(b) {
    wkb <- wk_p + pts[[b]]$wle_ref - pt_hat$wle_ref +
      drop(Jw_p %*% (D[b, idx] - d_hat))
    cut_rates(pts[[b]], wkb, TARGET$mu, TARGET$sd, cut, cut_raw)
  }, rate_p))
  rows[[length(rows) + 1L]] <- data.frame(
    context = "projection", quantity = names(rate_p),
    estimate = unname(rate_p), se = apply(draws_p, 2, stats::sd),
    lo = apply(draws_p, 2, stats::quantile, .025),
    hi = apply(draws_p, 2, stats::quantile, .975))
  out <- do.call(rbind, rows)
  out$d_hat <- paste(round(d_hat, 4), collapse = ",")
  out
}

# ---------------------------------------------------------------------------
# Population targets
# ---------------------------------------------------------------------------

pop_split_delta <- function(thr_true, thr_focal, split_items, mu2, n) {
  st <- list(pop_stats(thr_true, 0, 1), pop_stats(thr_focal, mu2, 1))
  prep <- make_prep(lengths(thr_true), split_items, st, n)
  f <- fit_prep_nohess(prep)
  th <- par_to_thr(f$par, prep, names(thr_true))
  gc <- group_comparison(th, list(st[[1]]$pr, st[[2]]$pr),
                         equate = length(split_items) > 0L, hessian = FALSE)
  list(delta = vapply(c("wle", "latent", "sum"), function(e) gc[[e]][["delta"]], 0),
       thr = th)
}

cut_raw_for <- function(cell) {
  sc <- make_scale(cell$type, cell$k)
  round(item_moments(sc$thr, cell$cut)$E)
}

ext_truth <- function(cell) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  j <- si$j
  dvec <- cell$d * si$sign
  n <- c(cell$n_ref, cell$n_focal)
  common <- setdiff(seq_along(sc$thr), j)
  thr_split <- with_shifts(list(ref = sc$thr, focal = sc$thr), j, dvec)
  cut_raw <- cut_raw_for(cell)
  keep <- pop_split_delta(sc$thr, thr_split$focal, integer(0), cell$mu2, n)
  split <- pop_split_delta(sc$thr, thr_split$focal, j, cell$mu2, n)
  tg <- c()
  tg[paste0("sample.impact_", names(keep$delta))] <- keep$delta - split$delta
  if (length(j) >= 2L) {
    for (q in seq_along(j)) {
      pq <- pop_split_delta(sc$thr, thr_split$focal, j[-q], cell$mu2, n)
      tg[paste0("sample.item_", names(sc$thr)[j[q]], "_impact_",
                names(pq$delta))] <-
        pq$delta - split$delta
    }
  }
  pt <- rate_parts(thr_split, j, cut_raw)
  wk <- function(thr_k) wle_by_score(thr_k) + link_offset(thr_k, sc$thr, common)
  r_s <- cut_rates(pt, wk(keep$thr$ref), c(0, cell$mu2), c(1, 1), cell$cut,
                   cut_raw)
  tg[paste0("sample.", names(r_s))] <- r_s
  pj <- project(thr_split, TARGET)
  tg[paste0("projection.impact_", names(pj$impact))] <- pj$impact
  tg[paste0("projection.bias_in_se_", names(pj$impact))] <- pj$impact / pj$se_target
  tg[paste0("projection.se_target_", names(pj$impact))] <- pj$se_target
  r_p <- cut_rates(pt, wk(pj$thr_keep), TARGET$mu, TARGET$sd, cell$cut, cut_raw)
  tg[paste0("projection.", names(r_p))] <- r_p
  list(target = tg, cut_raw = cut_raw)
}

# ---------------------------------------------------------------------------
# Replications and summaries
# ---------------------------------------------------------------------------

ext_rep <- function(cell, cut_raw, rep_id) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  dif <- stats::setNames(cell$d * si$sign, names(sc$thr)[si$j])
  tw <- sim_twin(sc$thr, dif, cell$n_ref, cell$n_focal, mu_focal = cell$mu2)
  e <- estimate_ext(tw$X_dif, tw$g, si$j, cell$cut, cut_raw)
  cbind(rep = rep_id, e)
}

ext_cell <- function(cell, reps, seed, workers) {
  tr <- ext_truth(cell)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  out <- parallel::mclapply(seq_len(reps), function(r) {
    tryCatch(ext_rep(cell, tr$cut_raw, r), error = function(e) {
      list(error = conditionMessage(e), rep = r)
    })
  }, mc.cores = workers, mc.set.seed = TRUE)
  err <- vapply(out, function(x) is.list(x) && !is.data.frame(x), TRUE)
  list(cell = cell, truth = tr,
       estimates = do.call(rbind, out[!err]),
       errors = if (any(err)) vapply(out[err], `[[`, "", "error") else character(0))
}

summarise_ext <- function(res) {
  e <- res$estimates
  key <- paste(e$context, e$quantity, sep = ".")
  tg <- res$truth$target
  do.call(rbind, lapply(split(e, key), function(x) {
    k <- paste(x$context[1], x$quantity[1], sep = ".")
    t0 <- if (k %in% names(tg)) tg[[k]] else NA_real_
    has_ci <- all(is.finite(x$lo))
    ok <- is.finite(x$se)
    data.frame(
      context = x$context[1], quantity = x$quantity[1], reps = nrow(x),
      target = t0, mean_est = mean(x$estimate), bias = mean(x$estimate) - t0,
      emp_sd = stats::sd(x$estimate),
      mean_se = if (any(ok)) mean(x$se[ok]) else NA_real_,
      se_ratio = if (any(ok)) mean(x$se[ok]) / stats::sd(x$estimate) else NA_real_,
      coverage = if (has_ci) mean(x$lo <= t0 & t0 <= x$hi)
                 else if (any(ok)) mean(abs(x$estimate[ok] - t0) <= 1.96 * x$se[ok])
                 else NA_real_,
      interval = if (has_ci) "draws" else if (any(ok)) "delta" else "none",
      share_zero = mean(x$estimate == 0),
      stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  smoke <- identical(Sys.getenv("DIF_EXT_SMOKE"), "1")
  workers <- as.integer(Sys.getenv("DIF_EXT_WORKERS", "4"))
  cells <- grid_scenarios()
  reps <- 1000L
  if (smoke) {
    pick <- (cells$k == 5 & cells$scenario == "one" & cells$position == "near_cut" &
               cells$d == 0.5 & cells$n_focal == 200) |
      (cells$scenario == "two_opposite" & cells$d == 0.8 & cells$n_focal == 500) |
      (cells$scenario == "two_same" & cells$position == "central" &
         cells$d == 0.5 & cells$n_focal == 250)
    cells <- cells[pick, ]
    reps <- as.integer(Sys.getenv("DIF_EXT_SMOKE_REPS", "8"))
  }
  rownames(cells) <- NULL
  stem <- if (smoke) "ext_smoke" else "ext_results"
  # One file per cell, so a stopped run restarts where it left off (finished
  # cells are read back, not rerun).
  cell_dir <- file.path(ext_dir, paste0(stem, "_cells"))
  dir.create(cell_dir, showWarnings = FALSE)
  t0 <- Sys.time()
  results <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    f_cell <- file.path(cell_dir, sprintf("cell_%02d.rds", i))
    if (file.exists(f_cell)) {
      results[[i]] <- readRDS(f_cell)
      message(sprintf("cell %d/%d read from %s", i, nrow(cells), f_cell))
      next
    }
    cell <- as.list(cells[i, ])
    res <- ext_cell(cell, reps, seed = 20261120L + i, workers = workers)
    res$summary <- if (!is.null(res$estimates)) summarise_ext(res) else NULL
    saveRDS(res, f_cell)
    results[[i]] <- res
    message(sprintf("cell %d/%d done (%s), errors: %d", i, nrow(cells),
                    format(Sys.time() - t0, digits = 3), length(res$errors)))
    if (length(res$errors)) print(utils::head(unique(res$errors)))
  }
  out <- list(results = results, cells = cells, reps = reps,
              n_draws = N_DRAWS, target_population = TARGET,
              run_time = Sys.time() - t0, date = Sys.Date(),
              session = utils::sessionInfo())
  saveRDS(out, file.path(ext_dir, paste0(stem, ".rds")))
  saveRDS(out, file.path(ext_dir, sprintf("%s_%s.rds", stem,
                                          format(Sys.Date(), "%Y-%m-%d"))))
  message("Done in ", format(out$run_time))
}
