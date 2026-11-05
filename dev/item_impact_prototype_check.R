## Checks of item_impact_prototype.R (validation design, stage 1)
##
## 1. sim_rep() reproduces one_rep() from infit_impact_validation.R exactly
##    (same seed, same inputs), for all four mechanisms, and one cached
##    replicate at n = 20,000.
## 2. Unit-weight FWB infit equals the conditional infit used in the
##    calibration.
## 3. Curves at the proposed defaults (n_sim = 5000, reps = 20) with the true
##    model agree with the cached curves (n = 20,000, 40 replicates).
## 4. End-to-end runs on single datasets: null, underfit, overfit by response
##    dependence, a joint group. Output, timing, and the group check.
##
## Run from easyRasch2/dev: source("item_impact_prototype_check.R")
## Results saved to item_impact_validation/prototype_check.rds (+ dated).

Sys.setenv(IMPACT_SMOKE = "1")
source("infit_impact_validation.R")   # BASES, one_rep(), N_CAL (runs nothing)
source("item_impact_prototype.R")
OUT <- "item_impact_validation"
dir.create(OUT, showWarnings = FALSE)
NC <- 10
checks <- list()

seed_of <- function(base, tpos, rho, r, n, mech)
  1e7 * match(base$design, c("d20", "phq9")) + 1e5 * tpos +
    round(1000 * rho) * 10 + r + 1e9 * (n != N_CAL) + n +
    5e8 * (mech == "slope") + 6e8 * (mech == "joint") + 7e8 * (mech == "ld")

## ---- 1. exact reproduction ---------------------------------------------------
cols <- c("msq", "rmse_a", "rmse_b", "rmse_c", "shift_mean", "kl_mean", "E_wle",
          "D_a", "D_b")
repro <- do.call(rbind, lapply(list(
  list("d20", 1, 0.7, "xdim"), list("phq9", 3, 0.5, "xdim"),
  list("d20", 2, 2, "slope"), list("phq9", 1, 1, "ld"),
  list("d20", 1, 0.6, "joint"), list("phq9", 1, 0.6, "joint")), function(a) {
  b <- BASES[[a[[1]]]]; tpos <- a[[2]]; par <- a[[3]]; mech <- a[[4]]
  n <- 2000; r <- 1
  ref <- one_rep(b, tpos, par, r, n, mech)
  i <- b$targets[tpos]
  focal <- if (mech == "joint") b$targets else if (mech == "ld")
    setdiff(order(abs(b$loc - b$loc[i])), i)[1] else i
  partner <- if (mech == "ld") i else NULL
  new <- sim_rep(list(thr = b$thr, thetas = b$thetas, loc = b$loc), focal, mech,
                 par, n, seed_of(b, tpos, par, r, n, mech), partner)
  d_sum <- max(abs(unlist(ref$summary[cols]) - unlist(new$summary[cols])), na.rm = TRUE)
  d_f <- max(abs(as.matrix(ref$fsum[, c("n", "bias_a", "bias_b", "sem_b")]) -
                   as.matrix(new$fsum)))
  data.frame(design = b$design, mech = mech, par = par, max_diff_summary = d_sum,
             max_diff_fsum = d_f)
}))
print(repro)

# one cached replicate at n = 20,000 (environment unchanged since the run)
cached <- readRDS(file.path("infit_impact_validation", "d20_on_rho070.rds"))$summary
b <- BASES$d20
new <- sim_rep(list(thr = b$thr, thetas = b$thetas, loc = b$loc), b$targets[1],
               "xdim", 0.7, N_CAL, seed_of(b, 1, 0.7, 1, N_CAL, "xdim"))
cc <- intersect(cols, names(cached))
cache_diff <- max(abs(unlist(cached[cached$rep == 1, cc]) -
                        unlist(new$summary[cc])))
cat(sprintf("cached replicate, max abs difference: %.2g\n", cache_diff))
checks$repro <- repro
checks$cache_diff <- cache_diff

## ---- 2. unit-weight infit -----------------------------------------------------
X <- BASES$phq9$X
thr <- .fit_cml_thresholds(X)
w1 <- imp_infit_weighted(X, rep(1, nrow(X)))
direct <- vapply(seq_len(ncol(X)), function(j) imp_cond_infit(X, thr, j), numeric(1))
checks$infit_diff <- max(abs(w1 - direct))
cat(sprintf("unit-weight vs direct conditional infit: %.2g\n", checks$infit_diff))

# (saved at the end)

## ---- 3. curves at the proposed defaults vs the cache ------------------------
## Cached curve: cell means at n = 20,000 and 40 replicates, F from the pooled
## sums with the same 200-person floor.
cached_curve <- function(design, target, mech) {
  pat <- switch(mech, xdim = "rho", slope = "slope", ld = "ldlambda")
  grid <- switch(mech, xdim = IMP_RHOS, slope = IMP_SLOPES, ld = IMP_LAMBDAS)
  do.call(rbind, lapply(grid, function(p) {
    x <- readRDS(file.path("infit_impact_validation",
                           sprintf("%s_%s_%s%03d.rds", design, target, pat, round(100 * p))))
    s <- x$summary; f <- x$fsum; f <- f[f$n >= 200, ]
    d <- (f$bias_a - f$bias_b) / f$n
    data.frame(par = p, msq = mean(s$msq), kl = mean(s$kl_mean),
               shift = mean(s$shift_mean), cbias = max(abs(d) / (f$sem_b / f$n)),
               rmse = mean(s$rmse_a / s$rmse_b) - 1,
               keep_remove = mean(s$rmse_a / s$rmse_c))
  }))
}
curve_cmp <- do.call(rbind, lapply(list(
  list("d20", 1, "xdim"), list("phq9", 1, "xdim"),
  list("d20", 1, "slope"), list("phq9", 1, "ld")), function(a) {
  b <- BASES[[a[[1]]]]; tpos <- a[[2]]; mech <- a[[3]]
  i <- b$targets[tpos]
  focal <- if (mech == "ld") setdiff(order(abs(b$loc - b$loc[i])), i)[1] else i
  partner <- if (mech == "ld") i else NULL
  t0 <- Sys.time()
  cv <- calibrate_curve(list(thr = b$thr, thetas = b$thetas, loc = b$loc), focal,
                        mech, n_sim = 5000, reps = 20, seed = 424242,
                        partner = partner, n_cores = NC)
  secs <- as.numeric(Sys.time() - t0, units = "secs")
  ref <- cached_curve(b$design, b$target_label[tpos], mech)
  m <- merge(cv[, c("par", "msq", "kl", "shift", "cbias", "rmse", "keep_remove")],
             ref, by = "par", suffixes = c("", "_cache"))
  m$design <- b$design; m$mech <- mech; m$seconds <- secs
  m
}))
print(curve_cmp[, c("design", "mech", "par", "msq", "msq_cache", "kl", "kl_cache",
                    "shift", "shift_cache", "cbias", "cbias_cache",
                    "keep_remove", "keep_remove_cache", "seconds")], digits = 3)
checks$curve_cmp <- curve_cmp

## ---- 4. end-to-end on single datasets ---------------------------------------
## Data generated from the base truth (resampled WLE persons, base thresholds).
gen_data <- function(b, mech, par, focal, n, seed, partner = NULL) {
  set.seed(seed)
  k <- length(b$thr); mu <- mean(b$thetas)
  th1 <- sample(b$thetas, n, replace = TRUE)
  ee <- sample(b$thetas, n, replace = TRUE)
  th2 <- mu + par * (th1 - mu) + sqrt(1 - par^2) * (ee - mu)
  U <- matrix(runif(n * k), n, k)
  X <- vapply(seq_len(k), function(j) imp_gen_item(th1, b$thr[[j]], U[, j]), numeric(n))
  if (mech %in% c("xdim", "joint")) for (ii in focal) X[, ii] <- imp_gen_item(th2, b$thr[[ii]], U[, ii])
  if (mech == "ld") X[, focal] <- imp_gen_item(th1 + par * X[, partner], b$thr[[focal]], U[, focal])
  colnames(X) <- colnames(b$X)
  X
}
e2e_cases <- list(
  null_d20  = list("d20", "null", 1, NULL, NULL),
  xdim_d20  = list("d20", "xdim", 0.5, BASES$d20$targets[1], NULL),
  ld_phq9   = list("phq9", "ld", 1,
                   setdiff(order(abs(BASES$phq9$loc - BASES$phq9$loc[BASES$phq9$targets[1]])),
                           BASES$phq9$targets[1])[1], BASES$phq9$targets[1]),
  joint_d20 = list("d20", "joint", 0.6, BASES$d20$targets, NULL)
)
truth_file <- list(xdim_d20 = "d20_on_rho050.rds", ld_phq9 = "phq9_on_ldlambda100.rds",
                   joint_d20 = "d20_joint_jointrho060.rds")
e2e <- lapply(names(e2e_cases), function(nm) {
  a <- e2e_cases[[nm]]
  X <- gen_data(BASES[[a[[1]]]], a[[2]], a[[3]], a[[4]], 2000, seed = 777, partner = a[[5]])
  t0 <- Sys.time()
  res <- item_impact(X, n_cores = NC, seed = 11)
  res$total_seconds <- as.numeric(Sys.time() - t0, units = "secs")
  res$case <- nm; res$focal <- a[[4]]; res$partner_true <- a[[5]]
  cat(sprintf("\n== %s: %.0f s (fwb %.0f, calibration %.0f, group %.0f)\n", nm,
              res$total_seconds, res$timing["fwb"], res$timing["calibration"],
              res$timing["group"]))
  tb <- res$table[res$table$calibrated | res$table$flagged, ]
  print(tb[, c("item", "msq", "lo", "hi", "partner", "kl_lo", "kl_hi", "shift_lo",
               "shift_hi", "group_lo", "group_hi", "cbias_lo", "cbias_hi",
               "rmse_lo", "rmse_hi", "dsem_lo", "dsem_hi", "beyond_grid",
               "turns_back", "keep_remove")], digits = 3)
  if (!is.null(truth_file[[nm]])) {
    x <- readRDS(file.path("infit_impact_validation", truth_file[[nm]]))
    s <- x$summary; f <- x$fsum; f <- f[f$n >= 200, ]
    d <- (f$bias_a - f$bias_b) / f$n
    cat("truth (cached cell, n = 20,000, 40 replicates):\n")
    print(data.frame(msq = mean(s$msq), kl = mean(s$kl_mean), shift = mean(s$shift_mean),
                     group = mean(s$E_wle), cbias = max(abs(d) / (f$sem_b / f$n)),
                     rmse = mean(s$rmse_a / s$rmse_b) - 1,
                     keep_remove = mean(s$rmse_a / s$rmse_c)), digits = 3)
  }
  if (!is.null(res$group)) { cat("group rows:\n"); print(do.call(rbind, res$group), digits = 3)
    cat("Q3 cutoff:", attr(res$group, "q3_cutoff"), "\n") }
  res$X <- X
  res
})
names(e2e) <- names(e2e_cases)
checks$e2e <- e2e

## ---- 5. overfit side: dependence vs slope on the same datasets ---------------
## Cascade items (xdim_d20, joint_d20) have no true partner. Records the
## partner's Q3 (relative to the mean) against the Q3 cutoff.
slope_cmp <- lapply(c("xdim_d20", "ld_phq9", "joint_d20"), function(nm) {
  r0 <- checks$e2e[[nm]]
  r1 <- item_impact(r0$X, n_cores = NC, seed = 11, overfit_mechanism = "slope",
                    group_check = FALSE)
  ov <- which(r0$table$calibrated & r0$table$msq < 1)
  q3s <- r0$q3 - mean(r0$q3[upper.tri(r0$q3)])
  cut <- RMlocdepQ3Cutoff(as.data.frame(r0$X), iterations = 200, parallel = FALSE,
                          cutoff_method = "quantile", seed = 11)$suggested_cutoff
  out <- data.frame(case = nm, item = ov, msq = r0$table$msq[ov],
                    partner = r0$table$partner[ov],
                    partner_q3 = q3s[cbind(ov, r0$table$partner[ov])], q3_cutoff = cut,
                    kl_dep = sprintf("%.4f-%.4f", r0$table$kl_lo[ov], r0$table$kl_hi[ov]),
                    kl_slope = sprintf("%.4f-%.4f", r1$table$kl_lo[ov], r1$table$kl_hi[ov]),
                    shift_dep = sprintf("%.3f-%.3f", r0$table$shift_lo[ov], r0$table$shift_hi[ov]),
                    shift_slope = sprintf("%.3f-%.3f", r1$table$shift_lo[ov], r1$table$shift_hi[ov]),
                    beyond_slope = r1$table$beyond_grid[ov])
  out
})
slope_cmp <- do.call(rbind, slope_cmp)
cat("\n== overfit side, dependence vs slope\n")
print(slope_cmp, digits = 3)
checks$slope_cmp <- slope_cmp

saveRDS(checks, file.path(OUT, "prototype_check.rds"))
saveRDS(checks, file.path(OUT, sprintf("prototype_check_%s.rds", Sys.Date())))
cat("done\n")
