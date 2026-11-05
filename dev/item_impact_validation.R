## Validation study for item impact intervals (RMitemImpact prototype)
##
## Design: dev/item-impact-validation-design.md. Prototype:
## dev/item_impact_prototype.R (overfit_mechanism = "auto" default, decided
## 2026-10-13).
##
## Arms (section 6 of the design note):
##   A  coverage with matched mechanism: A1 null, A2 xdim .7, A3 xdim .5,
##      A4 xdim .7 off-target, A5 response dependence lambda = 1;
##      n = 500, 2000, 5000
##   B  B1 steep slope a = 2 (analysed auto, dependence, slope), B2 flat slope
##      (MSQ about 1.25, analysed with the default), B3 = A5 data analysed with
##      the true partner forced; n = 2000
##   C  C1 three items on one shared theta2 (rho .6), C2 three items each on
##      its own theta2 (rho .6); n = 2000
##   D  FWB intervals only, B = 1000: null and A2 on both bases, n = 500,
##      2000, 5000 (d20 uses this study's generator, not the Part A data)
## Data come from the base truth of infit_impact_validation.R (resampled WLE
## persons and base thresholds). Truth cells are that study's cache
## (n = 20,000, 40 replicates); B2's cells are added to it here.
##
## Every dataset is analysed with the function's default selection (items
## whose simultaneous interval excludes 1) plus the target item(s), so the
## target has an interval at every n and cascade items are recorded too.
##
## Run from easyRasch2/dev:
##   Sys.setenv(IMPACT_VAL_PILOT = "1"); source("item_impact_validation.R")
## runs the timing pilot (5 datasets per cell, into item_impact_validation/pilot/).
## IMPACT_VAL_SMOKE = "1" defines everything and runs nothing.
## Cells are cached undated and dated.

Sys.setenv(IMPACT_SMOKE = "1")
source("infit_impact_validation.R")    # BASES, run_cell() for truth, N_CAL
source("item_impact_prototype.R")

NC_OUTER <- 10
PILOT <- identical(Sys.getenv("IMPACT_VAL_PILOT"), "1")
R_DATA <- if (PILOT) 5 else 100
R_D    <- if (PILOT) 5 else 200
VAL_DIR <- if (PILOT) "item_impact_validation/pilot" else "item_impact_validation"
dir.create(VAL_DIR, showWarnings = FALSE, recursive = TRUE)

FLAT_A <- c(d20 = 0.55, phq9 = 0.60)   # pilot 2026-10-13: MSQ 1.25 and 1.27

## ---- conditions -------------------------------------------------------------
## focal: item(s) made to misfit. partner: for "ld", the item it depends on.
cond_spec <- function(cond, b) {
  on <- b$targets[1]; off <- b$targets[3]
  near <- setdiff(order(abs(b$loc - b$loc[on])), on)[1]
  switch(cond,
    A1 = list(mech = "none", par = NA, focal = integer(0)),
    A2 = list(mech = "xdim", par = 0.7, focal = on),
    A3 = list(mech = "xdim", par = 0.5, focal = on),
    A4 = list(mech = "xdim", par = 0.7, focal = off),
    A5 = list(mech = "ld", par = 1, focal = near, partner = on),
    B1 = list(mech = "slope", par = 2, focal = on),
    B2 = list(mech = "slope", par = FLAT_A[[b$design]], focal = on),
    C1 = list(mech = "joint", par = 0.6, focal = b$targets),
    C2 = list(mech = "indep", par = 0.6, focal = b$targets))
}
CONDS <- c("A1", "A2", "A3", "A4", "A5", "B1", "B2", "C1", "C2")

## Analyses per condition: list of item_impact() argument sets.
analyses_for <- function(cond, sp) {
  base <- list(default = list())
  if (cond == "B1") base <- list(auto = list(), dependence = list(overfit_mechanism = "dependence"),
                                 slope = list(overfit_mechanism = "slope"))
  base
}

## ---- data generation ----------------------------------------------------------
gen_val_data <- function(b, sp, n, seed) {
  set.seed(seed)
  k <- length(b$thr); mu <- mean(b$thetas)
  th1 <- sample(b$thetas, n, replace = TRUE)
  ee <- sample(b$thetas, n, replace = TRUE)
  U <- matrix(runif(n * k), n, k)
  X <- vapply(seq_len(k), function(j) imp_gen_item(th1, b$thr[[j]], U[, j]), numeric(n))
  if (sp$mech %in% c("xdim", "joint")) {
    th2 <- mu + sp$par * (th1 - mu) + sqrt(1 - sp$par^2) * (ee - mu)
    for (ii in sp$focal) X[, ii] <- imp_gen_item(th2, b$thr[[ii]], U[, ii])
  } else if (sp$mech == "indep") {
    for (ii in sp$focal) {
      e_i <- sample(b$thetas, n, replace = TRUE)
      th2 <- mu + sp$par * (th1 - mu) + sqrt(1 - sp$par^2) * (e_i - mu)
      X[, ii] <- imp_gen_item(th2, b$thr[[ii]], U[, ii])
    }
  } else if (sp$mech == "slope") {
    X[, sp$focal] <- imp_gen_item(th1, b$thr[[sp$focal]], U[, sp$focal], a = sp$par)
  } else if (sp$mech == "ld") {
    X[, sp$focal] <- imp_gen_item(th1 + sp$par * X[, sp$partner], b$thr[[sp$focal]],
                                  U[, sp$focal])
  }
  colnames(X) <- colnames(b$X)
  X
}

data_seed <- function(cond, b, n, r)
  1e8 * (2 * (match(cond, CONDS) - 1) + match(b$design, c("d20", "phq9"))) +
    100 * n + r

## Keep what the report needs, drop the FWB draws.
slim <- function(res) {
  res$fwb <- res$fwb[c("obs", "sd", "lo", "hi", "slo", "shi", "cmax")]
  res$model <- NULL
  res
}

analyse_dataset <- function(cond, b, n, r) {
  sp <- cond_spec(cond, b)
  seed <- data_seed(cond, b, n, r)
  X <- gen_val_data(b, sp, n, seed)
  out <- list(cond = cond, design = b$design, n = n, rep = r, focal = sp$focal,
              partner = sp$partner, analyses = list())
  for (an in names(analyses_for(cond, sp))) {
    args <- c(list(X = X, items_extra = sp$focal, n_cores = 1, seed = seed + 1),
              analyses_for(cond, sp)[[an]])
    t0 <- Sys.time()
    res <- do.call(item_impact, args)
    res$seconds <- as.numeric(Sys.time() - t0, units = "secs")
    out$analyses[[an]] <- slim(res)
  }
  if (cond == "A5" && n == 2000) {           # B3: true partner forced
    po <- rep(NA_integer_, ncol(X)); po[sp$focal] <- sp$partner
    t0 <- Sys.time()
    res <- item_impact(X, items_extra = sp$focal, n_cores = 1, seed = seed + 1,
                       overfit_mechanism = "dependence", partner_override = po,
                       group_check = FALSE)
    res$seconds <- as.numeric(Sys.time() - t0, units = "secs")
    out$analyses$true_partner <- slim(res)
  }
  out
}

fwb_only <- function(cond, b, n, r, B = 1000) {
  sp <- cond_spec(cond, b)
  seed <- data_seed(cond, b, n, r)
  X <- gen_val_data(b, sp, n, seed)
  set.seed(seed + 2)
  t0 <- Sys.time()
  f <- imp_fwb(X, B, 0.95)
  sec <- as.numeric(Sys.time() - t0, units = "secs")
  # B = 200 from the first 200 draws of the same run
  d200 <- f$draws[seq_len(min(200, nrow(f$draws))), , drop = FALSE]
  sd200 <- apply(d200, 2, sd)
  t200 <- apply(abs(sweep(d200, 2, f$obs)) / rep(sd200, each = nrow(d200)), 1, max)
  c200 <- unname(quantile(t200, 0.95))
  data.frame(cond = cond, design = b$design, n = n, rep = r, item = seq_len(ncol(X)),
             obs = f$obs, sd = f$sd, lo = f$lo, hi = f$hi, slo = f$slo, shi = f$shi,
             cmax = f$cmax,
             sd200 = sd200, lo200 = apply(d200, 2, quantile, 0.025),
             hi200 = apply(d200, 2, quantile, 0.975),
             slo200 = f$obs - c200 * sd200, shi200 = f$obs + c200 * sd200, cmax200 = c200,
             draws_ok = nrow(f$draws), seconds = sec)
}

## ---- cells --------------------------------------------------------------------
run_val_cell <- function(arm, cond, b, n) {
  f <- file.path(VAL_DIR, sprintf("%s_%s_%s_n%04d.rds", arm, cond, b$design, n))
  if (file.exists(f)) return(invisible(readRDS(f)))
  t0 <- Sys.time()
  R <- if (arm == "D") R_D else R_DATA
  res <- parallel::mclapply(seq_len(R), function(r) {
    tryCatch(if (arm == "D") fwb_only(cond, b, n, r) else analyse_dataset(cond, b, n, r),
             error = function(e) conditionMessage(e))
  }, mc.cores = NC_OUTER, mc.preschedule = FALSE)
  ok <- vapply(res, function(x) is.list(x), logical(1))
  out <- list(results = if (arm == "D") do.call(rbind, res[ok]) else res[ok],
              failed = unlist(res[!ok]),
              minutes = as.numeric(Sys.time() - t0, units = "mins"),
              arm = arm, cond = cond, design = b$design, n = n, R = R)
  saveRDS(out, f)
  saveRDS(out, sub("[.]rds$", paste0("_", Sys.Date(), ".rds"), f))
  cat(sprintf("%s %s %s n=%d: %d/%d ok, %.1f min\n", arm, cond, b$design, n,
              sum(ok), R, out$minutes))
  invisible(out)
}

CELLS <- rbind(
  expand.grid(arm = "A", cond = c("A1", "A2", "A3", "A4", "A5"),
              design = c("d20", "phq9"), n = c(500, 2000, 5000), stringsAsFactors = FALSE),
  expand.grid(arm = "B", cond = c("B1", "B2"), design = c("d20", "phq9"), n = 2000,
              stringsAsFactors = FALSE),
  expand.grid(arm = "C", cond = c("C1", "C2"), design = c("d20", "phq9"), n = 2000,
              stringsAsFactors = FALSE),
  expand.grid(arm = "D", cond = c("A1", "A2"), design = c("d20", "phq9"),
              n = c(500, 2000, 5000), stringsAsFactors = FALSE))

if (!identical(Sys.getenv("IMPACT_VAL_SMOKE"), "1")) {
  # truth cells for B2 (added to the impact-validation cache)
  for (b in BASES) run_cell(b, 1, FLAT_A[[b$design]], mech = "slope")
  for (z in seq_len(nrow(CELLS)))
    run_val_cell(CELLS$arm[z], CELLS$cond[z], BASES[[CELLS$design[z]]], CELLS$n[z])
  cat("done\n")
}
