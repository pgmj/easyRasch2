## Prototype of RMitemImpact() internals (not package code)
##
## Design: dev/item-impact-interval-design.md (function) and
## dev/item-impact-validation-design.md (study, section 3).
## Decisions 2026-10-12: infit only, all five metrics, complete cases.
## 2026-10-13: overfit_mechanism = "auto" (response dependence when the item's
## largest residual correlation exceeds the Q3 cutoff, steep slope otherwise).
##
## Steps: conditional infit MSQ per item, FWB intervals (per item and
## simultaneous), item selection, impact curve per item side by simulation
## from the fitted model (same-item counterfactual with common random numbers,
## as in infit_impact_validation.R), mapping of the MSQ interval through the
## curve, keep or remove, and the group check.
##
## The response generator and the per-replicate analysis follow one_rep() in
## infit_impact_validation.R line by line, with the same order of random
## draws, so that sim_rep() with the same seed reproduces one_rep() exactly.
## Checks: dev/item_impact_prototype_check.R.
##
## Requires devtools::load_all("..") for .fit_cml_thresholds(),
## .estimate_thetas(), .wle_theta_pool(), .cond_moments(),
## .q3_residual_matrix() and RMlocdepQ3Cutoff().

IMP_RHOS    <- c(0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1)
IMP_SLOPES  <- c(1, 1.15, 1.3, 1.5, 1.75, 2, 2.5, 3)
IMP_LAMBDAS <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2)
IMP_D_GROUP <- 0.5
IMP_METRICS <- c("kl", "shift", "group", "cbias", "rmse")
## Grid extensions, calibrated only when the MSQ interval reaches past the
## end of the main grid and MSQ is still moving away from 1 there.
IMP_EXT <- list(xdim = c(0.15, 0), joint = c(0.15, 0), slope = c(4, 5),
                ld = c(3, 4))

## ---- response generation (identical to infit_impact_validation.R) ---------
imp_gen_item <- function(theta, thr, u, a = 1) {
  m <- length(thr)
  eta <- a * (outer(theta, 0:m) - matrix(c(0, cumsum(thr)), length(theta),
                                         m + 1L, byrow = TRUE))
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta) / rowSums(exp(eta))
  cum <- t(apply(P, 1, cumsum))
  rowSums(cum[, -(m + 1L), drop = FALSE] < u)
}

imp_analyse <- function(X, gen_loc, common) {
  thr <- .fit_cml_thresholds(X)
  est_loc <- vapply(thr, mean, numeric(1))
  shift <- mean(gen_loc[common]) - mean(est_loc[common])
  th <- .estimate_thetas(X, thr, method = "WLE")
  list(thr = thr, theta = th$theta + shift, sem = th$sem, sum = rowSums(X))
}

## ---- conditional infit with case weights ------------------------------------
## Unit weights give the ordinary conditional infit (equals iarm::out_infit()).
imp_infit_weighted <- function(X, w) {
  m <- psychotools::pcmodel(as.data.frame(X), weights = w)
  thr <- lapply(psychotools::threshpar(m), as.numeric)
  mom <- .cond_moments(thr)
  s <- rowSums(X)
  keep <- s > 0 & s < mom$max_score
  E <- mom$E[s[keep] + 1, , drop = FALSE]
  V <- mom$V[s[keep] + 1, , drop = FALSE]
  wk <- w[keep]
  colSums(wk * (X[keep, , drop = FALSE] - E)^2) / colSums(wk * V)
}

imp_cond_infit <- function(X, thr, j) {
  mom <- .cond_moments(thr)
  s <- rowSums(X)
  keep <- s > 0 & s < mom$max_score
  E <- mom$E[s[keep] + 1, j]
  V <- mom$V[s[keep] + 1, j]
  sum((X[keep, j] - E)^2) / sum(V)
}

## ---- FWB intervals for the MSQ ---------------------------------------------
imp_fwb <- function(X, iterations, conf_level) {
  n <- nrow(X)
  obs <- imp_infit_weighted(X, rep(1, n))
  draws <- matrix(NA_real_, iterations, ncol(X))
  for (b in seq_len(iterations)) {
    w <- stats::rexp(n)
    w <- w / mean(w)
    draws[b, ] <- tryCatch(imp_infit_weighted(X, w),
                           error = function(e) rep(NA_real_, ncol(X)))
  }
  draws <- draws[stats::complete.cases(draws), , drop = FALSE]
  a <- (1 - conf_level) / 2
  sdv <- apply(draws, 2, stats::sd)
  tmax <- apply(abs(sweep(draws, 2, obs)) / rep(sdv, each = nrow(draws)), 1, max)
  cmax <- unname(stats::quantile(tmax, conf_level))
  list(obs = obs, draws = draws, sd = sdv,
       lo = apply(draws, 2, stats::quantile, a),
       hi = apply(draws, 2, stats::quantile, 1 - a),
       slo = obs - cmax * sdv, shi = obs + cmax * sdv, cmax = cmax)
}

## ---- one calibration replicate ----------------------------------------------
## model: list(thr, thetas, loc) from the fitted data (or the truth).
## focal: the item(s) whose responses are made to misfit.
## mech = "xdim": focal responds to theta2 with correlation par (underfit).
## mech = "slope": focal has discrimination par.
## mech = "joint": every focal item responds to one shared theta2.
## mech = "ld": focal responds to theta1 + par * x_partner (uniform response
##   dependence), the partner unchanged.
## The order of random draws matches one_rep() in infit_impact_validation.R:
## for "ld", one_rep()'s target i is the partner and its j is the focal item.
sim_rep <- function(model, focal, mech, par, n, seed, partner = NULL) {
  set.seed(seed)
  k <- length(model$thr)
  mu <- mean(model$thetas); sg <- sd(model$thetas)
  th1 <- sample(model$thetas, n, replace = TRUE)
  ee  <- sample(model$thetas, n, replace = TRUE)
  uses_th2 <- mech %in% c("xdim", "joint")
  th2 <- if (uses_th2) mu + par * (th1 - mu) + sqrt(1 - par^2) * (ee - mu) else th1
  G   <- rbinom(n, 1, 0.5)
  th2g <- if (uses_th2) th2 + G * IMP_D_GROUP * sg * sqrt(1 - par^2) else th1
  U <- matrix(runif(n * k), n, k)

  Xb <- vapply(seq_len(k), function(j) imp_gen_item(th1, model$thr[[j]], U[, j]),
               numeric(n))
  Xa <- Xb
  if (uses_th2) {
    for (ii in focal) Xa[, ii] <- imp_gen_item(th2, model$thr[[ii]], U[, ii])
  } else if (mech == "slope") {
    Xa[, focal] <- imp_gen_item(th1, model$thr[[focal]], U[, focal], a = par)
  } else {
    Xa[, focal] <- imp_gen_item(th1 + par * Xb[, partner], model$thr[[focal]],
                                U[, focal])
  }
  if (uses_th2) {
    Xag <- Xb
    for (ii in focal) Xag[, ii] <- imp_gen_item(th2g, model$thr[[ii]], U[, ii])
  }

  changed <- if (mech == "ld") c(partner, focal) else focal
  common <- setdiff(seq_len(k), changed)
  keep_c <- setdiff(seq_len(k), focal)
  A  <- imp_analyse(Xa, model$loc, common)
  Bf <- imp_analyse(Xb, model$loc, common)
  C  <- imp_analyse(Xb[, keep_c, drop = FALSE], model$loc[keep_c],
                    match(common, keep_c))
  if (uses_th2) AG <- imp_analyse(Xag, model$loc, common)

  msq_each <- vapply(focal, function(q) imp_cond_infit(Xa, A$thr, q), numeric(1))
  rmse <- function(x) sqrt(mean((x - th1)^2))
  shift_sem <- abs(A$theta - Bf$theta) / Bf$sem
  r_se <- (A$sem / Bf$sem)^2
  kl <- 0.5 * ((A$theta - Bf$theta) / Bf$sem)^2 + 0.5 * (r_se - 1 - log(r_se))
  gdiff <- function(x) mean(x[G == 1]) - mean(x[G == 0])
  dsem <- function(x) sd(x$theta - th1) / sqrt(mean(x$sem^2))
  key <- factor(th1, levels = sort(unique(model$thetas)))
  list(
    summary = data.frame(
      msq = mean(msq_each),
      rmse_a = rmse(A$theta), rmse_b = rmse(Bf$theta), rmse_c = rmse(C$theta),
      shift_mean = mean(shift_sem), kl_mean = mean(kl),
      E_wle = if (uses_th2) (gdiff(AG$theta) - gdiff(Bf$theta)) / sg / IMP_D_GROUP
              else NA_real_,
      D_a = dsem(A), D_b = dsem(Bf)
    ),
    fsum = data.frame(
      n = as.numeric(table(key)),
      bias_a = as.numeric(tapply(A$theta - th1, key, sum, default = 0)),
      bias_b = as.numeric(tapply(Bf$theta - th1, key, sum, default = 0)),
      sem_b  = as.numeric(tapply(Bf$sem, key, sum, default = 0))
    )
  )
}

## ---- impact curve for one item side -----------------------------------------
## One row per grid point: mean MSQ and every metric over the replicates.
## Metrics are stored as departures from no impact, so the null grid point
## (rho = 1, a = 1, lambda = 0) is 0 on every metric and its MSQ is set to 1.
calibrate_curve <- function(model, focal, mech, n_sim, reps, seed,
                            partner = NULL, n_cores = 1, f_min_n = 200,
                            grid = NULL) {
  if (is.null(grid))
    grid <- switch(mech, xdim = IMP_RHOS, joint = IMP_RHOS, slope = IMP_SLOPES,
                   ld = IMP_LAMBDAS)
  null_par <- switch(mech, xdim = 1, joint = 1, slope = 1, ld = 0)
  jobs <- expand.grid(g = seq_along(grid), r = seq_len(reps))
  res <- parallel::mclapply(seq_len(nrow(jobs)), function(z) {
    g <- jobs$g[z]; r <- jobs$r[z]
    tryCatch(sim_rep(model, focal, mech, grid[g], n_sim,
                     seed = seed + 1000L * g + r, partner = partner),
             error = function(e) conditionMessage(e))
  }, mc.cores = n_cores)
  ok <- vapply(res, is.list, logical(1))
  rows <- lapply(seq_along(grid), function(g) {
    idx <- which(ok & jobs$g == g)
    if (!length(idx)) return(NULL)
    s <- do.call(rbind, lapply(res[idx], `[[`, "summary"))
    f <- Reduce(`+`, lapply(res[idx], `[[`, "fsum"))
    f <- f[f$n >= f_min_n, , drop = FALSE]
    d <- (f$bias_a - f$bias_b) / f$n
    data.frame(par = grid[g], reps = length(idx), msq = mean(s$msq),
               kl = mean(s$kl_mean), shift = mean(s$shift_mean),
               group = mean(s$E_wle), cbias = max(abs(d) / (f$sem_b / f$n)),
               rmse = mean(s$rmse_a / s$rmse_b) - 1,
               dsem = mean(s$D_a / s$D_b) - 1,
               keep_remove = mean(s$rmse_a / s$rmse_c))
  })
  cv <- do.call(rbind, rows)
  z <- cv$par == null_par
  cv$msq[z] <- 1
  cv[z, c("kl", "shift", "cbias", "rmse", "dsem")] <- 0
  if (mech %in% c("xdim", "joint")) cv$group[z] <- 0
  cv <- cv[order(abs(cv$par - null_par)), ]   # path from the null outwards
  attr(cv, "mech") <- mech
  attr(cv, "failed") <- sum(!ok)
  cv
}

## Main grid plus the extension when the interval reaches past the grid end
## and MSQ has not yet turned back toward 1.
calibrate_side <- function(model, focal, mech, lo, hi, n_sim, reps, seed,
                           partner = NULL, n_cores = 1) {
  cv <- calibrate_curve(model, focal, mech, n_sim, reps, seed, partner, n_cores)
  up <- max(cv$msq) > 1
  far <- if (up) which.max(cv$msq) else which.min(cv$msq)
  past <- if (up) hi > max(cv$msq) else lo < min(cv$msq)
  if (past && far == nrow(cv)) {
    ext <- calibrate_curve(model, focal, mech, n_sim, reps, seed + 5e5, partner,
                           n_cores, grid = IMP_EXT[[mech]])
    cv <- rbind(cv, ext)
    attr(cv, "mech") <- mech
  }
  cv
}

## ---- mapping an MSQ interval through a curve --------------------------------
## The curve is a path in (MSQ, impact) indexed by the mechanism parameter. It
## is interpolated finely along the parameter, and the impact interval is the
## range of impact over the path points whose MSQ lies in [lo, hi]. This
## equals the monotone mapping when MSQ is monotone in the parameter, and
## stays valid when it is not (response dependence in phq9, where MSQ turns
## back toward 1 at large lambda while impact keeps rising): every parameter
## value consistent with the interval contributes. beyond_grid flags an
## interval reaching past the most extreme MSQ on the path.
map_interval <- function(cv, metric, lo, hi, n_fine = 400) {
  t <- seq(1, nrow(cv), length.out = n_fine)
  m <- stats::approx(seq_len(nrow(cv)), cv$msq, xout = t)$y
  v <- stats::approx(seq_len(nrow(cv)), cv[[metric]], xout = t)$y
  up <- max(cv$msq) > 1
  # MSQ turns back toward 1 before the grid end: the parameter is not
  # identified from the MSQ, and the upper bound is set by the grid end.
  turns_back <- (if (up) which.max(cv$msq) else which.min(cv$msq)) < nrow(cv)
  beyond <- if (up) hi > max(m) else lo < min(m)
  inside <- m >= lo & m <= hi
  if (!any(inside)) {
    if (!beyond) return(c(lo = 0, hi = 0, beyond_grid = 0, turns_back = turns_back))
    far <- if (up) which.max(m) else which.min(m)
    return(c(lo = v[far], hi = max(v), beyond_grid = 1, turns_back = turns_back))
  }
  # past the grid end the upper bound is at least the largest value on the path
  c(lo = min(v[inside]), hi = if (beyond) max(v) else max(v[inside]),
    beyond_grid = as.numeric(beyond), turns_back = as.numeric(turns_back))
}

## Impact at a point MSQ (for keep/remove and curve-error checks): the value
## on the monotone branch nearest the null.
map_point <- function(cv, metric, msq) {
  cut <- if (max(cv$msq) > 1) which.max(cv$msq) else which.min(cv$msq)
  b <- cv[seq_len(cut), ]
  stats::approx(b$msq, b[[metric]], xout = msq, rule = 2, ties = mean)$y
}

## ---- the prototype function --------------------------------------------------
item_impact <- function(X, items = NULL, iterations = 200, simultaneous = TRUE,
                        conf_level = 0.95,
                        overfit_mechanism = c("auto", "dependence", "slope"),
                        items_extra = NULL, partner_override = NULL,
                        n_sim = 5000, reps = 20, group_check = TRUE,
                        q3_iterations = 200, n_cores = 1, seed = 1,
                        verbose = FALSE) {
  overfit_mechanism <- match.arg(overfit_mechanism)
  X <- as.matrix(X)
  X <- X[stats::complete.cases(X), , drop = FALSE]
  k <- ncol(X)
  timing <- c(fwb = 0, q3_cutoff = 0, calibration = 0, group = 0)

  t0 <- Sys.time()
  set.seed(seed)
  cal_seed <- sample.int(1e8, 1)   # base for the calibration seeds
  fwb <- imp_fwb(X, iterations, conf_level)
  timing["fwb"] <- as.numeric(Sys.time() - t0, units = "secs")
  lo <- if (simultaneous) fwb$slo else fwb$lo
  hi <- if (simultaneous) fwb$shi else fwb$hi
  flagged <- which(lo > 1 | hi < 1)
  cal_items <- if (is.null(items)) flagged else items
  cal_items <- sort(union(cal_items, items_extra))

  pool <- .wle_theta_pool(X)
  model <- list(thr = pool$thr_list, thetas = pool$thetas,
                loc = vapply(pool$thr_list, mean, numeric(1)))
  q3 <- .q3_residual_matrix(X, "CML")
  q3s <- q3 - mean(q3[upper.tri(q3)])
  ## Q3 cutoff (99th percentile of max Q3 - mean Q3 under the fitted model),
  ## needed by "auto" for overfit items and by the group check. Computed once.
  need_over <- overfit_mechanism == "auto" && any(lo[cal_items] < 1)
  need_group <- group_check && sum(cal_items %in% which(lo > 1)) >= 2
  q3_cut <- NA_real_
  t0 <- Sys.time()
  if (need_over || need_group)
    q3_cut <- RMlocdepQ3Cutoff(as.data.frame(X), iterations = q3_iterations,
                               parallel = FALSE, cutoff_method = "quantile",
                               seed = seed)$suggested_cutoff
  timing["q3_cutoff"] <- as.numeric(Sys.time() - t0, units = "secs")

  t0 <- Sys.time()
  curves <- list()
  partner <- rep(NA_integer_, k)
  partner_q3 <- rep(NA_real_, k)
  over_mech <- rep(NA_character_, k)
  for (i in cal_items) {
    sides <- c(if (hi[i] > 1) "under", if (lo[i] < 1) "over")
    for (sd_ in sides) {
      p <- NULL
      if (sd_ == "under") {
        mech <- "xdim"
      } else {
        qi <- q3s[i, ]; qi[i] <- NA
        p <- if (!is.null(partner_override) && !is.na(partner_override[i]))
               partner_override[i]
             else if (all(is.na(qi))) setdiff(order(abs(model$loc - model$loc[i])), i)[1]
             else which.max(qi)
        partner_q3[i] <- qi[p]
        # auto: dependence only when the pair's residual correlation exceeds
        # the Q3 cutoff, otherwise a steep slope (the generic overfit)
        mech <- switch(overfit_mechanism, slope = "slope", dependence = "ld",
                       auto = if (isTRUE(qi[p] > q3_cut)) "ld" else "slope")
        if (mech == "ld") partner[i] <- p else p <- NULL
        over_mech[i] <- mech
      }
      if (verbose) message(sprintf("calibrating item %d, %s side (%s)", i, sd_, mech))
      curves[[sprintf("%d_%s", i, sd_)]] <-
        calibrate_side(model, i, mech, lo[i], hi[i], n_sim, reps,
                       seed = cal_seed + i * 1e4 + (sd_ == "over") * 5e3,
                       partner = p, n_cores = n_cores)
    }
  }
  timing["calibration"] <- as.numeric(Sys.time() - t0, units = "secs")

  metrics <- c(IMP_METRICS, "dsem")
  tab <- data.frame(item = seq_len(k), name = colnames(X), msq = fwb$obs,
                    lo = lo, hi = hi, flagged = seq_len(k) %in% flagged,
                    calibrated = seq_len(k) %in% cal_items,
                    over_mech = over_mech, partner = partner,
                    partner_q3 = partner_q3, q3_cutoff = q3_cut)
  for (m in metrics) tab[[paste0(m, "_lo")]] <- tab[[paste0(m, "_hi")]] <- NA_real_
  tab$beyond_grid <- tab$turns_back <- NA
  tab$keep_remove <- NA_real_
  for (i in cal_items) {
    cvs <- curves[grep(sprintf("^%d_", i), names(curves))]
    bg <- tb <- FALSE
    for (m in metrics) {
      iv <- lapply(cvs, function(cv) {
        if (all(is.na(cv[[m]]))) return(NULL)
        map_interval(cv, m, lo[i], hi[i])
      })
      iv <- do.call(rbind, iv)
      if (is.null(iv)) next
      tab[i, paste0(m, "_lo")] <- min(iv[, "lo"])
      tab[i, paste0(m, "_hi")] <- max(iv[, "hi"])
      bg <- bg || any(iv[, "beyond_grid"] == 1)
      tb <- tb || any(iv[, "turns_back"] == 1)
    }
    tab$beyond_grid[i] <- bg
    tab$turns_back[i] <- tb
    side_cv <- cvs[[if (fwb$obs[i] > 1) 1 else length(cvs)]]
    tab$keep_remove[i] <- map_point(side_cv, "keep_remove", fwb$obs[i])
  }

  ## group check: underfit-side calibrated items clustered on Q3
  group <- NULL
  t0 <- Sys.time()
  under <- intersect(cal_items, which(fwb$obs > 1 & lo > 1))
  if (group_check && length(under) >= 2) {
    cut <- q3_cut
    adj <- q3s[under, under, drop = FALSE] > cut
    adj[is.na(adj)] <- FALSE
    comp <- igraph_free_components(adj)
    sets <- lapply(unique(comp), function(z) under[comp == z])
    sets <- sets[lengths(sets) >= 2]
    if (length(sets)) {
      group <- lapply(sets, function(s) {
        gm <- rowMeans(fwb$draws[, s, drop = FALSE])
        g_obs <- mean(fwb$obs[s])
        a <- (1 - conf_level) / 2
        glo <- unname(stats::quantile(gm, a)); ghi <- unname(stats::quantile(gm, 1 - a))
        cv <- calibrate_side(model, s, "joint", glo, ghi, n_sim, reps,
                             seed = cal_seed + 9e5 + sum(s), n_cores = n_cores)
        out <- data.frame(items = paste(s, collapse = ","), msq = g_obs,
                          lo = glo, hi = ghi)
        for (m in IMP_METRICS) {
          iv <- map_interval(cv, m, glo, ghi)
          out[[paste0(m, "_lo")]] <- iv["lo"]; out[[paste0(m, "_hi")]] <- iv["hi"]
        }
        out$keep_remove <- map_point(cv, "keep_remove", g_obs)
        attr(out, "curve") <- cv
        out
      })
    }
    if (is.null(group)) group <- list()   # checked, no group found
    attr(group, "q3_cutoff") <- cut
  }
  timing["group"] <- as.numeric(Sys.time() - t0, units = "secs")

  list(table = tab, group = group, curves = curves, fwb = fwb, model = model,
       q3 = q3, timing = timing,
       settings = list(iterations = iterations, simultaneous = simultaneous,
                       conf_level = conf_level, overfit_mechanism = overfit_mechanism,
                       n_sim = n_sim, reps = reps, n = nrow(X)))
}

## Connected components of a logical adjacency matrix, without igraph.
igraph_free_components <- function(adj) {
  n <- nrow(adj); comp <- rep(NA_integer_, n); cc <- 0L
  for (s in seq_len(n)) {
    if (!is.na(comp[s])) next
    cc <- cc + 1L; stack <- s
    while (length(stack)) {
      v <- stack[1]; stack <- stack[-1]
      if (!is.na(comp[v])) next
      comp[v] <- cc
      stack <- c(stack, which(adj[v, ] & is.na(comp)))
    }
  }
  comp
}
