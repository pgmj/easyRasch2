# Step 3b of RMdifImpact (dev/dif-impact-design.md, Part B, "Step 3
# results", 2026-11-18). Reruns the two outputs that step 3 did not
# validate.
#
#   1. Rates at a cut, conditional on the fitted classification rule. The
#      rule (which raw scores are classified at or above the cut, per scale,
#      strategy and group) is taken from the fit and held fixed:
#        theta cut, keep:  keep WLE (linked to the split metric) >= cut
#        theta cut, split: group's split WLE >= cut
#        sum cut, keep:    raw >= cut_raw
#        sum cut, split:   raw >= cut_raw (ref), equated raw >= cut_raw (focal)
#      Its rates (prev, FN, FP as shares of the group, and keep - split) are
#      estimated under the split model with each group's latent N(mu, sd),
#      in the sample (estimated mu, sd) and in the target population of
#      step 3. With the rule fixed, the rates are smooth in the parameters.
#      Target per replication: the same rule's rates under the true item
#      parameters and true latent distributions. Two intervals:
#        delta: central differences in all split-model parameters, with
#               (mu, log sd) moving through their Jacobian, plus the
#               sampling error of (mu, log sd) themselves (sample context)
#        draws: percentiles of N_DRAWS draws, as in step 3
#   2. Per-item impacts (two-item cells) with the ratio SE
#      |impact_j / dif_hat_j| * SE(dif_hat_j), SE from the CML covariance.
#
# Also recorded: the fitted rule as a raw-score cut against the rule from
# the population tables (how often the fitted rule is the population one).
#
# Oracle DIF items, DIF data, grid_scenarios() (42 cells), seeds
# 20261130 + cell. Default 500 replications (Monte Carlo SE of a .95
# coverage .010).
#
# Run from the package root. Smoke test (container, 4 workers):
#   DIF_3B_SMOKE=1 LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript dev/dif_impact_ext_sim_3b.R
# Full run on the host (about 1.5 h on 10 workers, 2.5 h on 6; restarts
# skip finished cells, saved one file per cell in ext3b_results_cells/):
#   DIF_3B_WORKERS=10 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
#     nohup Rscript dev/dif_impact_ext_sim_3b.R > dev/dif_impact_ext_sim/run_3b.log 2>&1 &

source("dev/dif_impact_ext_sim.R")  # step 3 helpers; its run block does not fire

# ---------------------------------------------------------------------------
# Fixed classification rules and their rates
# ---------------------------------------------------------------------------

make_rule <- function(pt, wk, cut, cut_raw) {
  r <- 0:pt$R
  list(theta = list(keep = list(ref = wk >= cut, focal = wk >= cut),
                    split = list(ref = pt$wle_ref >= cut,
                                 focal = pt$wle_focal >= cut)),
       sum = list(keep = list(ref = r >= cut_raw, focal = r >= cut_raw),
                  split = list(ref = r >= cut_raw,
                               focal = pt$eq >= cut_raw - 1e-9)))
}

# Lowest raw score classified at or above the cut (R + 1 if none).
rule_cuts <- function(rule) {
  out <- c()
  for (sc in names(rule)) for (st in c("keep", "split")) for (g in c("ref", "focal")) {
    cl <- rule[[sc]][[st]][[g]]
    stopifnot(!is.unsorted(cl))  # monotone in the raw score
    out[paste(sc, st, g, sep = ".")] <- if (any(cl)) which(cl)[1] - 1L else length(cl)
  }
  out
}

# Rates of a fixed rule. pt supplies the score distributions given theta
# and the sum-score true-status cut; both follow the item parameters.
fixed_rates <- function(pt, rule, mu, sd, cut) {
  th_c <- c(theta = cut, sum = pt$theta_c_sum)
  out <- c()
  for (sc in names(rule)) for (st in c("keep", "split")) for (gi in 1:2) {
    g <- c("ref", "focal")[gi]
    v <- rates_one(pt$SP[[g]], rule[[sc]][[st]][[g]], mu[gi], sd[gi], th_c[[sc]])
    out[paste(sc, st, g, names(v), sep = ".")] <- v
  }
  for (sc in names(rule)) for (g in c("ref", "focal")) {
    for (m in c("prev", "fn", "fp")) {
      out[paste(sc, "diff", g, m, sep = ".")] <-
        out[[paste(sc, "keep", g, m, sep = ".")]] -
        out[[paste(sc, "split", g, m, sep = ".")]]
    }
  }
  out
}

# ---------------------------------------------------------------------------
# Estimates from one dataset
# ---------------------------------------------------------------------------

estimate_3b <- function(X, g, items, cut, cut_raw, h = 1e-3) {
  nm <- colnames(X)
  foc <- g == "focal"
  n <- c(sum(!foc), sum(foc))

  thr_keep <- .fit_cml_thresholds(X)
  names(thr_keep) <- nm
  sf <- shift_cml(X, g, items, hessian = TRUE)
  if (is.null(sf$vcov)) stop("Singular Hessian in the split fit.", call. = FALSE)
  items <- sf$prep$dif_items  # sorted by shift_cml()
  nd <- length(items)
  common <- setdiff(seq_len(ncol(X)), items)
  P <- sf$prep$P
  idx <- P - 1L + seq_len(nd)
  V <- sf$vcov
  dif_hat <- unname(sf$d)
  par_hat <- c(sf$psi[-1L], dif_hat)
  np <- length(par_hat)
  thr <- list(ref = sf$thresholds$ref, focal = sf$thresholds$focal)
  R <- sum(lengths(thr$ref))
  rows <- list()

  # 1. per-item impacts, ratio SE
  if (nd >= 2L) {
    sp <- score_strategy(X, g, thr, 0, TRUE, NULL)
    for (q in seq_len(nd)) {
      sq <- shift_cml(X, g, items[-q], hessian = FALSE)
      sc_q <- score_strategy(X, g, list(ref = sq$thresholds$ref,
                                        focal = sq$thresholds$focal), 0, TRUE, NULL)
      imp_q <- unname(sc_q$group[, "delta"] - sp$group[, "delta"])
      rows[[length(rows) + 1L]] <- data.frame(
        context = "sample",
        quantity = paste0("item_", nm[items[q]], "_impact_", rownames(sp$group)),
        estimate = imp_q,
        se_delta = abs(imp_q / dif_hat[q]) * sqrt(V[idx[q], idx[q]]))
    }
  }

  # 2. the fitted rule
  pr <- lapply(c("ref", "focal"), function(lev) {
    tabulate(rowSums(X)[g == lev] + 1L, R + 1L) / sum(g == lev)
  })
  lat <- lapply(1:2, function(gi) lat_fit(thr[[gi]], pr[[gi]], n[gi]))
  ms_hat <- c(lat[[1]]$par, lat[[2]]$par)  # mu_ref, log sd_ref, mu_f, log sd_f
  Vms <- matrix(0, 4L, 4L)
  Vms[1:2, 1:2] <- lat[[1]]$vcov
  Vms[3:4, 3:4] <- lat[[2]]$vcov
  pt_hat <- rate_parts(thr, items, cut_raw)
  wk <- wle_by_score(thr_keep) + link_offset(thr_keep, thr$ref, common)
  rule <- make_rule(pt_hat, wk, cut, cut_raw)
  rs <- function(pt, ms) fixed_rates(pt, rule, ms[c(1, 3)], exp(ms[c(2, 4)]), cut)
  rp <- function(pt) fixed_rates(pt, rule, TARGET$mu, TARGET$sd, cut)
  est_s <- rs(pt_hat, ms_hat)
  est_p <- rp(pt_hat)

  # (mu, log sd) Jacobian in the item parameters (forward differences)
  Jms <- matrix(0, 4L, np)
  for (p in seq_len(np)) {
    tp <- par_to_thr(par_hat + h * (seq_len(np) == p), sf$prep, nm)
    Jms[, p] <- (c(mml_fit(tp$ref, pr[[1]], lat[[1]]$par, hessian = FALSE)$par,
                   mml_fit(tp$focal, pr[[2]], lat[[2]]$par, hessian = FALSE)$par) -
                   ms_hat) / h
  }

  # delta: total derivative in each item parameter (mu, sd follow through
  # Jms), then the own error of (mu, log sd)
  Gs <- matrix(0, length(est_s), np)
  Gp <- matrix(0, length(est_p), np)
  for (p in seq_len(np)) {
    e <- h * (seq_len(np) == p)
    pt_u <- rate_parts(par_to_thr(par_hat + e, sf$prep, nm), items, cut_raw)
    pt_l <- rate_parts(par_to_thr(par_hat - e, sf$prep, nm), items, cut_raw)
    Gs[, p] <- (rs(pt_u, ms_hat + h * Jms[, p]) - rs(pt_l, ms_hat - h * Jms[, p])) / (2 * h)
    Gp[, p] <- (rp(pt_u) - rp(pt_l)) / (2 * h)
  }
  Gm <- matrix(0, length(est_s), 4L)
  for (m in 1:4) {
    e <- h * (1:4 == m)
    Gm[, m] <- (rs(pt_hat, ms_hat + e) - rs(pt_hat, ms_hat - e)) / (2 * h)
  }
  se_s <- sqrt(pmax(rowSums((Gs %*% V) * Gs) + rowSums((Gm %*% Vms) * Gm), 0))
  se_p <- sqrt(pmax(rowSums((Gp %*% V) * Gp), 0))

  # draws
  D <- mvn_draws(N_DRAWS, par_hat, V)
  Ems <- cbind(mvn_draws(N_DRAWS, c(0, 0), lat[[1]]$vcov),
               mvn_draws(N_DRAWS, c(0, 0), lat[[2]]$vcov))
  dr_s <- matrix(0, N_DRAWS, length(est_s))
  dr_p <- matrix(0, N_DRAWS, length(est_p))
  for (b in seq_len(N_DRAWS)) {
    pt_b <- rate_parts(par_to_thr(D[b, ], sf$prep, nm), items, cut_raw)
    dr_s[b, ] <- rs(pt_b, ms_hat + drop(Jms %*% (D[b, ] - par_hat)) + Ems[b, ])
    dr_p[b, ] <- rp(pt_b)
  }
  q_lo <- function(M) apply(M, 2, stats::quantile, .025, names = FALSE)
  q_hi <- function(M) apply(M, 2, stats::quantile, .975, names = FALSE)

  rows[[length(rows) + 1L]] <- data.frame(
    context = "sample", quantity = names(est_s), estimate = unname(est_s),
    se_delta = se_s, se_draws = apply(dr_s, 2, stats::sd),
    lo_draws = q_lo(dr_s), hi_draws = q_hi(dr_s))
  rows[[length(rows) + 1L]] <- data.frame(
    context = "projection", quantity = names(est_p), estimate = unname(est_p),
    se_delta = se_p, se_draws = apply(dr_p, 2, stats::sd),
    lo_draws = q_lo(dr_p), hi_draws = q_hi(dr_p))
  rc <- rule_cuts(rule)
  rows[[length(rows) + 1L]] <- data.frame(context = "rule", quantity = names(rc),
                                          estimate = unname(rc))
  out <- do.call(rbind, lapply(rows, function(x) {
    for (col in c("se_delta", "se_draws", "lo_draws", "hi_draws")) {
      if (is.null(x[[col]])) x[[col]] <- NA_real_
    }
    x[c("context", "quantity", "estimate", "se_delta", "se_draws",
        "lo_draws", "hi_draws")]
  }))
  list(est = out, rule = rule, dif_hat = dif_hat)
}

# ---------------------------------------------------------------------------
# Truth per cell, and per replication for the fitted rule
# ---------------------------------------------------------------------------

truth_3b <- function(cell) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  j <- si$j
  thr_split <- with_shifts(list(ref = sc$thr, focal = sc$thr), j,
                           cell$d * si$sign)
  cut_raw <- cut_raw_for(cell)
  n <- c(cell$n_ref, cell$n_focal)
  common <- setdiff(seq_along(sc$thr), j)
  pt <- rate_parts(thr_split, j, cut_raw)
  keep <- pop_split_delta(sc$thr, thr_split$focal, integer(0), cell$mu2, n)
  wk <- wle_by_score(keep$thr$ref) + link_offset(keep$thr$ref, sc$thr, common)
  tg_item <- c()
  if (length(j) >= 2L) {
    split <- pop_split_delta(sc$thr, thr_split$focal, j, cell$mu2, n)
    for (q in seq_along(j)) {
      pq <- pop_split_delta(sc$thr, thr_split$focal, j[-q], cell$mu2, n)
      tg_item[paste0("item_", names(sc$thr)[j[q]], "_impact_",
                     names(pq$delta))] <- pq$delta - split$delta
    }
  }
  list(pt = pt, cut_raw = cut_raw, mu = c(0, cell$mu2), sd = c(1, 1),
       rule_cuts = rule_cuts(make_rule(pt, wk, cell$cut, cut_raw)),
       item = tg_item)
}

rep_3b <- function(cell, tr, rep_id) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  dif <- stats::setNames(cell$d * si$sign, names(sc$thr)[si$j])
  tw <- sim_twin(sc$thr, dif, cell$n_ref, cell$n_focal, mu_focal = cell$mu2)
  e <- estimate_3b(tw$X_dif, tw$g, si$j, cell$cut, tr$cut_raw)
  x <- e$est
  # targets: the fitted rule under the true parameters and distributions
  tg_s <- fixed_rates(tr$pt, e$rule, tr$mu, tr$sd, cell$cut)
  tg_p <- fixed_rates(tr$pt, e$rule, TARGET$mu, TARGET$sd, cell$cut)
  x$target <- NA_real_
  i <- x$context == "sample" & x$quantity %in% names(tg_s)
  x$target[i] <- tg_s[x$quantity[i]]
  i <- x$context == "projection"
  x$target[i] <- tg_p[x$quantity[i]]
  i <- x$context == "rule"
  x$target[i] <- tr$rule_cuts[x$quantity[i]]
  i <- x$quantity %in% names(tr$item)
  x$target[i] <- tr$item[x$quantity[i]]
  x$dif_hat <- paste(round(e$dif_hat, 4), collapse = ",")
  cbind(rep = rep_id, x)
}

run_cell_3b <- function(cell, reps, seed, workers) {
  tr <- truth_3b(cell)
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  out <- parallel::mclapply(seq_len(reps), function(r) {
    tryCatch(rep_3b(cell, tr, r), error = function(e) {
      list(error = conditionMessage(e), rep = r)
    })
  }, mc.cores = workers, mc.set.seed = TRUE)
  err <- vapply(out, function(x) is.list(x) && !is.data.frame(x), TRUE)
  list(cell = cell, truth = tr[c("cut_raw", "rule_cuts", "item")],
       estimates = do.call(rbind, out[!err]),
       errors = if (any(err)) vapply(out[err], `[[`, "", "error") else character(0))
}

# Error = estimate - target (the target varies by replication for the
# rates). SE ratios are against the SD of the error.
summarise_3b <- function(res) {
  e <- res$estimates
  key <- paste(e$context, e$quantity, sep = ".")
  do.call(rbind, lapply(split(e, key), function(x) {
    err <- x$estimate - x$target
    sd_err <- stats::sd(err)
    rule <- x$context[1] == "rule"
    data.frame(
      context = x$context[1], quantity = x$quantity[1], reps = nrow(x),
      mean_target = mean(x$target), mean_est = mean(x$estimate),
      bias = mean(err), sd_err = sd_err,
      se_ratio_delta = if (rule) NA_real_ else mean(x$se_delta) / sd_err,
      cov_delta = if (rule) NA_real_
                  else mean(abs(err) <= 1.96 * x$se_delta + 1e-12),
      se_ratio_draws = if (rule || all(is.na(x$se_draws))) NA_real_
                       else mean(x$se_draws) / sd_err,
      cov_draws = if (rule || all(is.na(x$lo_draws))) NA_real_
                  else mean(x$lo_draws - 1e-12 <= x$target &
                              x$target <= x$hi_draws + 1e-12),
      rule_match = if (rule) mean(err == 0) else NA_real_,
      share_zero = mean(x$estimate == 0 & x$target == 0),
      stringsAsFactors = FALSE)
  }))
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  smoke <- identical(Sys.getenv("DIF_3B_SMOKE"), "1")
  workers <- as.integer(Sys.getenv("DIF_3B_WORKERS", "4"))
  reps <- as.integer(Sys.getenv("DIF_3B_REPS", "500"))
  cells <- grid_scenarios()
  rownames(cells) <- NULL
  ids <- seq_len(nrow(cells))
  if (smoke) {
    # a k5 cell, the cell-14 keep knife-edge, the cell-34 focal knife-edge
    ids <- c(4L, 14L, 34L)
    reps <- as.integer(Sys.getenv("DIF_3B_SMOKE_REPS", "8"))
  }
  stem <- if (smoke) "ext3b_smoke" else "ext3b_results"
  cell_dir <- file.path(ext_dir, paste0(stem, "_cells"))
  dir.create(cell_dir, showWarnings = FALSE)
  message(sprintf("%d cells, %d replications, %d workers", length(ids), reps,
                  workers))
  t0 <- Sys.time()
  results <- vector("list", length(ids))
  for (ii in seq_along(ids)) {
    i <- ids[ii]
    f_cell <- file.path(cell_dir, sprintf("cell_%02d.rds", i))
    if (file.exists(f_cell)) {
      results[[ii]] <- readRDS(f_cell)
      message(sprintf("cell %d/%d read from %s", i, nrow(cells), f_cell))
      next
    }
    cell <- as.list(cells[i, ])
    res <- run_cell_3b(cell, reps, seed = 20261130L + i, workers = workers)
    res$summary <- if (!is.null(res$estimates)) summarise_3b(res) else NULL
    saveRDS(res, f_cell)
    results[[ii]] <- res
    message(sprintf("cell %d/%d done (%s), errors: %d", i, nrow(cells),
                    format(Sys.time() - t0, digits = 3), length(res$errors)))
    if (length(res$errors)) print(utils::head(unique(res$errors)))
  }
  out <- list(results = results, cells = cells[ids, ], reps = reps,
              n_draws = N_DRAWS, target_population = TARGET, workers = workers,
              run_time = Sys.time() - t0, date = Sys.Date(),
              session = utils::sessionInfo())
  saveRDS(out, file.path(ext_dir, paste0(stem, ".rds")))
  saveRDS(out, file.path(ext_dir, sprintf("%s_%s.rds", stem,
                                          format(Sys.Date(), "%Y-%m-%d"))))
  message("Done in ", format(out$run_time))
}
