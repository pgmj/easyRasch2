# Steps 2 and 3 of RMdifImpact (dev/dif-impact-design.md, Part B,
# 2026-11-17): simulation of keep / remove / split with uniform DIF, two
# groups, and of the estimates a real-data function would report, with the
# oracle DIF items and with the items RMdifAnchor() flags.
#
# Each replication draws one dataset WITH the DIF and its twin WITHOUT it
# from the same persons and the same uniform draws (common random numbers),
# so only the DIF items' focal-group responses differ.
#
# Per replication and dataset (DIF and no-DIF twin):
#   group:     keep / remove / split (oracle items), Delta for WLE mean,
#              latent mean and sum score, the paired impact and naive SEs.
#   person:    WLE on the true metric per group: RMSE, bias, signed shift
#              from the twin, classification at the cut.
#   estimates: what the function would report, for the oracle item set and
#              for the RMdifAnchor() flags: impact = Delta(keep) -
#              Delta(split) per estimator, and the local between-group bias
#              (max |b|, range, value at the cut, raw-point gap), each with a
#              delta-method SE on the estimated shifts.
#
# Cells: grid_scenarios() in dif_impact_theory.R (42 polytomous cells), each
# with exact population targets from run_setting().
#
# Complete data only (the latent mean is fitted to the raw-score
# distribution).
#
# Run from the package root. Smoke test (container, 4 workers):
#   DIF_SIM_SMOKE=1 LANG=C.UTF-8 LC_ALL=C.UTF-8 Rscript dev/dif_impact_sim.R
# Full run on the host (10 workers, about 4-4.5 h; 2 s per replication at
# k = 5, about 5 s at k = 10):
#   DIF_SIM_WORKERS=10 LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 \
#     nohup Rscript dev/dif_impact_sim.R > dev/dif_impact_sim/run.log 2>&1 &

source("dev/dif_impact_theory.R")  # helpers only; its run block does not fire

sim_dir <- "dev/dif_impact_sim"
dir.create(sim_dir, showWarnings = FALSE)

# Local-bias quantities reported: the range of b (signed), the unsigned
# maxima, and the values at the cut.
LB_QUANT <- c("b_lo", "b_hi", "b_maxabs", "b_sem_maxabs", "gap_maxabs",
              "b_cut", "gap_cut")

# ---------------------------------------------------------------------------
# Data generation with common random numbers
# ---------------------------------------------------------------------------

# Responses from uniforms U (n x k): category = number of cumulative
# probabilities below the uniform.
sim_from_u <- function(theta, thr_list, U) {
  X <- vapply(seq_along(thr_list), function(i) {
    p <- cat_probs_grid(thr_list[[i]], theta)
    cp <- t(apply(p, 1L, cumsum))
    as.integer(rowSums(U[, i] > cp[, -ncol(cp), drop = FALSE]))
  }, integer(length(theta)))
  if (!is.matrix(X)) X <- matrix(X, nrow = length(theta))
  colnames(X) <- names(thr_list)
  X
}

# thr: true reference thresholds (named list); dif: named vector of shifts
# for the focal group; theta from N(mu, sd) per group.
sim_twin <- function(thr, dif, n_ref, n_focal, mu_ref = 0, mu_focal = 0,
                     sd_ref = 1, sd_focal = 1) {
  g <- factor(rep(c("ref", "focal"), c(n_ref, n_focal)),
              levels = c("ref", "focal"))
  theta <- c(stats::rnorm(n_ref, mu_ref, sd_ref),
             stats::rnorm(n_focal, mu_focal, sd_focal))
  U <- matrix(stats::runif(length(theta) * length(thr)), length(theta))
  thr_f <- thr
  for (nm in names(dif)) thr_f[[nm]] <- thr_f[[nm]] + dif[[nm]]
  foc <- g == "focal"
  X0 <- sim_from_u(theta, thr, U)
  X1 <- X0
  X1[foc, ] <- sim_from_u(theta[foc], thr_f, U[foc, , drop = FALSE])
  list(X_dif = X1, X_null = X0, g = g, theta = theta)
}

# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------

item_loc <- function(thr) vapply(thr, mean, 0)

# Offset that moves a fit onto the true metric through the DIF-free items.
link_offset <- function(thr_est, thr_true, common) {
  mean(item_loc(thr_true)[common]) - mean(item_loc(thr_est)[common])
}

score_strategy <- function(X, g, thr_by_group, offset, equate, cut) {
  R <- sum(lengths(thr_by_group$ref))
  r <- rowSums(X)
  # A fit that downcoded an empty category has fewer thresholds than the data
  # has categories; scoring raw sums with it would be wrong.
  if (length(thr_by_group$ref) != ncol(X) ||
      any(apply(X, 2, max) > lengths(thr_by_group$ref))) {
    stop("Threshold lengths do not match the data (empty category).",
         call. = FALSE)
  }
  foc <- g == "focal"
  n <- c(sum(!foc), sum(foc))
  wtab <- lapply(thr_by_group, wle_by_score)
  semtab <- Map(function(thr, w) 1 / sqrt(item_moments(thr, w)$I),
                thr_by_group, wtab)
  wle <- ifelse(foc, wtab$focal[r + 1L], wtab$ref[r + 1L]) + offset
  sem <- ifelse(foc, semtab$focal[r + 1L], semtab$ref[r + 1L])
  lat <- lapply(c("ref", "focal"), function(lev) {
    pr <- tabulate(r[g == lev] + 1L, R + 1L) / sum(g == lev)
    mml_fit(thr_by_group[[lev]], pr)
  })
  s <- r
  if (equate) {
    s[foc] <- equate_focal(thr_by_group$ref, thr_by_group$focal)[r[foc] + 1L]
  }
  grp <- function(v) {
    m <- c(mean(v[!foc]), mean(v[foc]))
    vv <- c(stats::var(v[!foc]), stats::var(v[foc]))
    c(delta = m[2] - m[1], se = sqrt(vv[1] / n[1] + vv[2] / n[2]),
      smd = (m[2] - m[1]) / sqrt(mean(vv)))
  }
  group <- rbind(
    wle = grp(wle),
    latent = c(delta = lat[[2]]$mu - lat[[1]]$mu,
               se = sqrt(lat[[1]]$var_mu / n[1] + lat[[2]]$var_mu / n[2]),
               smd = (lat[[2]]$mu - lat[[1]]$mu) /
                 sqrt((lat[[1]]$sd^2 + lat[[2]]$sd^2) / 2)),
    sum = grp(s)
  )
  list(group = group, wle = wle, sem = sem,
       lat_mu = c(lat[[1]]$mu, lat[[2]]$mu),
       lat_sd = c(lat[[1]]$sd, lat[[2]]$sd),
       class = if (!is.null(cut)) wle >= cut else NULL)
}

# Oracle strategies on one dataset (the true DIF items are handled).
fit_and_score <- function(X, g, thr_true, dif_items, cut) {
  common <- setdiff(seq_len(ncol(X)), dif_items)
  thr_keep <- .fit_cml_thresholds(X)
  names(thr_keep) <- colnames(X)
  keep <- score_strategy(X, g, list(ref = thr_keep, focal = thr_keep),
                         link_offset(thr_keep, thr_true, common), FALSE, cut)
  Xr <- X[, common, drop = FALSE]
  thr_rem <- .fit_cml_thresholds(Xr)
  names(thr_rem) <- colnames(Xr)
  remove <- score_strategy(Xr, g, list(ref = thr_rem, focal = thr_rem),
                           link_offset(thr_rem, thr_true[common],
                                       seq_along(common)), FALSE, cut)
  sf <- shift_cml(X, g, dif_items, hessian = FALSE)
  split <- score_strategy(X, g, list(ref = sf$thresholds$ref,
                                     focal = sf$thresholds$focal),
                          link_offset(sf$thresholds$ref, thr_true, common),
                          TRUE, cut)
  list(keep = keep, remove = remove, split = split)
}

# ---------------------------------------------------------------------------
# The estimates a real-data function would report
# ---------------------------------------------------------------------------

# Split the given items (all others common), estimate the impact of keeping
# them as Delta(keep) - Delta(split) for each estimator, and the local
# between-group bias from the split parameters, over the central 95% of the
# estimated latent mixture. SEs by the delta method on the estimated shifts
# only (step 1: the shift error carries over 90% of the parameter variance).
estimate_impact <- function(X, g, items, keep_sc, cut, h = 1e-3) {
  sf <- shift_cml(X, g, items, hessian = TRUE)
  nd <- length(items)
  idx <- sf$prep$P - 1L + seq_len(nd)
  vc_d <- if (!is.null(sf$vcov)) sf$vcov[idx, idx, drop = FALSE] else
    matrix(NA_real_, nd, nd)
  thr <- list(ref = sf$thresholds$ref, focal = sf$thresholds$focal)
  sp <- score_strategy(X, g, thr, 0, TRUE, NULL)
  w_focal <- mean(g == "focal")
  lb_fun <- function(th) {
    local_bias(th$ref, th$focal, mu = sp$lat_mu, sd = sp$lat_sd,
               w_focal = w_focal, cut = cut)[LB_QUANT]
  }
  imp <- keep_sc$group[, "delta"] - sp$group[, "delta"]
  lb <- lb_fun(thr)
  Gi <- matrix(0, 3L, nd); Gl <- matrix(0, length(lb), nd)
  for (q in seq_len(nd)) {
    tp <- thr; tm <- thr
    tp$focal[[items[q]]] <- tp$focal[[items[q]]] + h
    tm$focal[[items[q]]] <- tm$focal[[items[q]]] - h
    dp <- score_strategy(X, g, tp, 0, TRUE, NULL)$group[, "delta"]
    dm <- score_strategy(X, g, tm, 0, TRUE, NULL)$group[, "delta"]
    Gi[, q] <- -(dp - dm) / (2 * h)
    Gl[, q] <- (lb_fun(tp) - lb_fun(tm)) / (2 * h)
  }
  se_imp <- sqrt(pmax(diag(Gi %*% vc_d %*% t(Gi)), 0))
  se_lb <- sqrt(pmax(diag(Gl %*% vc_d %*% t(Gl)), 0))
  data.frame(
    quantity = c(paste0("impact_", rownames(keep_sc$group)), LB_QUANT),
    estimate = c(unname(imp), unname(lb)),
    se = c(se_imp, se_lb),
    d_hat = paste(round(sf$d, 4), collapse = ","),
    stringsAsFactors = FALSE
  )
}

# Items RMdifAnchor() flags (default settings: forward consensus anchor,
# item level, single-step FWER, alpha .05).
# RMdifAnchor() stops when a category is unused in one group (common at
# k = 5 with a small focal group); that replication then has no flagged set,
# and the reason is recorded.
anchor_flags <- function(X, g) {
  a <- tryCatch(suppressMessages(RMdifAnchor(as.data.frame(X), g,
                                             output = "list")),
                error = function(e) e)
  if (inherits(a, "error")) {
    return(list(items = integer(0), error = conditionMessage(a)))
  }
  list(items = which(a$table$Flagged %in% TRUE), error = NA_character_)
}

# ---------------------------------------------------------------------------
# One replication
# ---------------------------------------------------------------------------

one_rep <- function(cell, rep_id) {
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  j <- si$j
  dvec <- cell$d * si$sign
  dif <- stats::setNames(dvec, names(sc$thr)[j])
  tw <- sim_twin(sc$thr, dif, cell$n_ref, cell$n_focal, mu_focal = cell$mu2)
  foc <- tw$g == "focal"
  fits <- list(dif = fit_and_score(tw$X_dif, tw$g, sc$thr, j, cell$cut),
               null = fit_and_score(tw$X_null, tw$g, sc$thr, j, cell$cut))
  a <- fits$dif; b <- fits$null

  # group level, oracle strategies
  rows <- list()
  for (strat in c("keep", "remove", "split")) {
    for (est in c("wle", "latent", "sum")) {
      ga <- a[[strat]]$group[est, ]; gb <- b[[strat]]$group[est, ]
      rows[[length(rows) + 1L]] <- data.frame(
        rep = rep_id, strategy = strat, estimator = est,
        delta_dif = ga[["delta"]], delta_null = gb[["delta"]],
        impact = ga[["delta"]] - gb[["delta"]],
        smd_dif = ga[["smd"]], smd_null = gb[["smd"]],
        impact_smd = ga[["smd"]] - gb[["smd"]],
        se_naive = ga[["se"]], stringsAsFactors = FALSE
      )
    }
  }
  group <- do.call(rbind, rows)

  # person level, WLE on the true metric
  prow <- list()
  for (strat in c("keep", "remove", "split")) {
    sa <- a[[strat]]; sb <- b[[strat]]
    shift_sem <- (sa$wle - sb$wle) / sb$sem
    for (lev in c("ref", "focal")) {
      idx <- if (lev == "focal") foc else !foc
      err <- sa$wle[idx] - tw$theta[idx]
      prow[[length(prow) + 1L]] <- data.frame(
        rep = rep_id, strategy = strat, group = lev,
        rmse_theta = sqrt(mean(err^2)), bias_theta = mean(err),
        mean_sem = mean(sa$sem[idx]),
        mean_shift_sem = mean(shift_sem[idx]),
        class_disagree = mean(sa$class[idx] != sb$class[idx]),
        class_error = mean(sa$class[idx] != (tw$theta[idx] >= cell$cut)),
        stringsAsFactors = FALSE
      )
    }
  }
  person <- do.call(rbind, prow)

  # estimates the function would report
  erows <- list()
  for (dat in c("dif", "null")) {
    X <- if (dat == "dif") tw$X_dif else tw$X_null
    keep_sc <- fits[[dat]]$keep
    af <- anchor_flags(X, tw$g)
    flagged <- af$items
    sets <- list(oracle = j, flagged = flagged)
    for (set in names(sets)) {
      items <- sets[[set]]
      meta <- data.frame(rep = rep_id, data = dat, set = set,
                         n_flagged = length(flagged),
                         tp = sum(flagged %in% j), fp = sum(!flagged %in% j),
                         items = paste(names(sc$thr)[items], collapse = ","),
                         anchor_error = !is.na(af$error),
                         stringsAsFactors = FALSE)
      if (!length(items)) {
        erows[[length(erows) + 1L]] <- cbind(meta, quantity = NA_character_,
                                             estimate = NA_real_,
                                             se = NA_real_, d_hat = NA_character_)
        next
      }
      e <- estimate_impact(X, tw$g, items, keep_sc, cell$cut)
      erows[[length(erows) + 1L]] <- cbind(meta[rep(1L, nrow(e)), ], e,
                                           row.names = NULL)
    }
  }
  list(group = group, person = person, estimates = do.call(rbind, erows))
}

run_cell <- function(cell, reps, seed, workers) {
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  out <- parallel::mclapply(seq_len(reps), function(r) {
    tryCatch(one_rep(cell, r), error = function(e) {
      list(error = conditionMessage(e), rep = r)
    })
  }, mc.cores = workers, mc.set.seed = TRUE)
  err <- vapply(out, function(x) !is.null(x$error), TRUE)
  list(cell = cell,
       group = do.call(rbind, lapply(out[!err], `[[`, "group")),
       person = do.call(rbind, lapply(out[!err], `[[`, "person")),
       estimates = do.call(rbind, lapply(out[!err], `[[`, "estimates")),
       errors = if (any(err)) vapply(out[err], `[[`, "", "error") else character(0))
}

# ---------------------------------------------------------------------------
# Population targets and summaries
# ---------------------------------------------------------------------------

# Exact targets for one cell from run_setting(): the impact a keep - split
# estimate converges to (per estimator), and the local-bias quantities.
cell_truth <- function(cell) {
  r <- run_setting(cell)
  imp <- vapply(c("wle", "latent", "sum"), function(est) {
    r$delta_d[r$strategy == "keep" & r$estimator == est] -
      r$delta_d[r$strategy == "split" & r$estimator == est]
  }, 0)
  lb <- unlist(r[1L, LB_QUANT])
  list(target = c(stats::setNames(imp, paste0("impact_", names(imp))), lb),
       setting = r)
}

# Oracle strategies: impact, empirical vs naive SE, naive CI coverage.
summarise_cell <- function(res) {
  gr <- res$group
  key <- interaction(gr$strategy, gr$estimator, drop = TRUE, sep = ":")
  s <- do.call(rbind, lapply(split(gr, key), function(x) {
    data.frame(
      strategy = x$strategy[1], estimator = x$estimator[1], reps = nrow(x),
      mean_impact = mean(x$impact),
      mcse_impact = stats::sd(x$impact) / sqrt(nrow(x)),
      mean_impact_smd = mean(x$impact_smd),
      emp_se = stats::sd(x$delta_dif), mean_se_naive = mean(x$se_naive),
      cover_naive = mean(abs(x$delta_dif - mean(x$delta_null)) <=
                           1.96 * x$se_naive),
      stringsAsFactors = FALSE
    )
  }))
  rownames(s) <- NULL
  s
}

# The reported estimates against the population targets (DIF data), and
# their spread on the no-DIF twin (noise floor).
summarise_estimates <- function(res, target) {
  e <- res$estimates
  out <- list()
  for (dat in c("dif", "null")) for (set in c("oracle", "flagged")) {
    x <- e[e$data == dat & e$set == set, ]
    reps_total <- length(unique(x$rep))
    # no flags, or RMdifAnchor() could not run
    reps_none <- length(unique(x$rep[is.na(x$quantity)]))
    x <- x[!is.na(x$quantity), ]
    for (qn in names(target)) {
      y <- x[x$quantity == qn, ]
      tg <- if (dat == "dif") target[[qn]] else 0
      ok <- is.finite(y$estimate) & is.finite(y$se)
      out[[length(out) + 1L]] <- data.frame(
        data = dat, set = set, quantity = qn,
        reps = reps_total, reps_no_flags = reps_none, reps_used = nrow(y),
        target = tg, mean_est = mean(y$estimate),
        bias = mean(y$estimate) - tg,
        emp_sd = stats::sd(y$estimate), mean_se = mean(y$se[ok]),
        se_ratio = mean(y$se[ok]) / stats::sd(y$estimate[ok]),
        coverage = mean(abs(y$estimate[ok] - tg) <= 1.96 * y$se[ok]),
        q95_abs = unname(stats::quantile(abs(y$estimate), .95)),
        stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

summarise_detection <- function(res, n_dif) {
  e <- res$estimates
  e <- e[!duplicated(e[, c("rep", "data")]), ]
  do.call(rbind, lapply(c("dif", "null"), function(dat) {
    x <- e[e$data == dat, ]
    data.frame(data = dat, reps = nrow(x),
               anchor_error = mean(x$anchor_error),
               any_flag = mean(x$n_flagged > 0),
               all_true_flagged = if (dat == "dif") mean(x$tp == n_dif) else NA,
               mean_tp = mean(x$tp), any_fp = mean(x$fp > 0),
               mean_fp = mean(x$fp))
  }))
}

# ---------------------------------------------------------------------------
# Run
# ---------------------------------------------------------------------------

if (sys.nframe() == 0L) {
  smoke <- identical(Sys.getenv("DIF_SIM_SMOKE"), "1")
  workers <- as.integer(Sys.getenv("DIF_SIM_WORKERS", "4"))
  cells <- grid_scenarios()
  reps <- 1000L
  if (smoke) {
    pick <- (cells$k == 5 & cells$scenario == "one" & cells$position == "near_cut" &
               cells$d == 0.5 & cells$n_focal == 200) |
      (cells$scenario == "two_opposite" & cells$d == 0.5 & cells$n_focal == 500) |
      (cells$k == 10 & cells$scenario == "one" & cells$position == "central" &
         cells$d == 0.5 & cells$n_focal == 250)
    cells <- cells[pick, ]
    reps <- 40L
  }
  rownames(cells) <- NULL
  t0 <- Sys.time()
  results <- vector("list", nrow(cells))
  for (i in seq_len(nrow(cells))) {
    cell <- as.list(cells[i, ])
    tr <- cell_truth(cell)
    res <- run_cell(cell, reps, seed = 20261118L + i, workers = workers)
    res$truth <- tr
    res$summary <- summarise_cell(res)
    res$estimate_summary <- summarise_estimates(res, tr$target)
    sc <- make_scale(cell$type, cell$k)
    res$detection <- summarise_detection(
      res, length(scenario_items(sc$loc, cell$scenario, cell$position,
                                 cell$cut)$j))
    results[[i]] <- res
    message(sprintf("cell %d/%d done (%s), errors: %d", i, nrow(cells),
                    format(Sys.time() - t0, digits = 3), length(res$errors)))
    if (length(res$errors)) print(utils::head(unique(res$errors)))
  }
  out <- list(results = results, cells = cells, reps = reps,
              run_time = Sys.time() - t0, date = Sys.Date(),
              session = utils::sessionInfo())
  stem <- if (smoke) "sim_smoke" else "sim_results"
  saveRDS(out, file.path(sim_dir, paste0(stem, ".rds")))
  saveRDS(out, file.path(sim_dir, sprintf("%s_%s.rds", stem,
                                          format(Sys.Date(), "%Y-%m-%d"))))
  message("Done in ", format(out$run_time))
}
