## Validation of the impact-calibrated infit band (upper side, M2)
##
## Design: dev/infit-impact-delta-design.md (decisions in section 9).
## For a target item i and a second-dimension correlation rho, three datasets
## share every person and every uniform draw:
##   (a) item i responds to theta2 (misfit, correlation rho with theta1)
##   (b) item i responds to theta1 (the same item, fitting)
##   (c) item i removed
## plus a pair (ag, bg) in which a random half of the sample is shifted on the
## unique part of theta2 by D_GROUP SDs, with theta1 unchanged (metric E).
## Each dataset is analysed as a user would: CML refit, WLE with SEM, linked to
## the generating scale through the items other than i.
##
## Designs: d20 (20 dichotomous items, the Part A design, targets at 0, -1 and
## -2 logits relative to the person mean) and phq9 (9 items with 4
## categories, the package data, targets on, mid and off target).
## Persons: resampled WLE from a base dataset (dgp = "resample" convention).
##
## Run from easyRasch2/dev: source("infit_impact_validation.R")
## Results are cached per cell in infit_impact_validation/ (undated and dated
## copies). Report: infit_impact_validation.qmd, which only reads the cache.
## Smoke test: Sys.setenv(IMPACT_SMOKE = "1") before sourcing defines
## everything and runs nothing.

suppressMessages({
  devtools::load_all("..", quiet = TRUE)   # .cond_moments, .estimate_thetas, ...
  library(parallel)
})

NCORES   <- 10
N_CAL    <- 20000      # population proxy
R_CAL    <- 40         # replicates per cell
N_STAB   <- c(2000, 5000)  # stability check, on-target items only
SLOPES   <- c(1, 1.15, 1.3, 1.5, 1.75, 2, 2.5, 3)  # overfit arm
LAMBDAS  <- c(0, 0.25, 0.5, 0.75, 1, 1.5, 2)          # response dependence arm
RHOS     <- c(0.3, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1)
D_GROUP  <- 0.5        # group shift on the unique part of theta2, in SD units
OUT_DIR  <- "infit_impact_validation"
dir.create(OUT_DIR, showWarnings = FALSE)

## ---- base designs ----------------------------------------------------------
make_base <- function(design) {
  if (design == "d20") {
    set.seed(20261012)
    locs <- seq(-2, 2, length.out = 20)
    X <- (matrix(runif(2000 * 20), 2000) <
            plogis(outer(rnorm(2000, 0, 1.5), locs, "-"))) * 1L
    colnames(X) <- paste0("V", 1:20)
  } else {
    e <- new.env()
    load("../data/phq9.rda", envir = e)
    X <- as.matrix(e$phq9[, paste0("q", 1:9)])
    X <- X[stats::complete.cases(X), ]
  }
  pool <- .wle_theta_pool(X)
  loc <- vapply(pool$thr_list, mean, numeric(1))
  rel <- loc - mean(pool$thetas)
  targets <- if (design == "d20") {
    out <- integer(0)
    for (tg in c(0, -1, -2)) out <- c(out, setdiff(order(abs(rel - tg)), out)[1])
    out
  } else {
    o <- order(abs(rel))
    c(o[1], o[ceiling(length(o) / 2)], o[length(o)])
  }
  list(design = design, X = X, thr = pool$thr_list, thetas = pool$thetas,
       loc = loc, rel = rel, targets = targets,
       target_label = c("on", "mid", "off"))
}

## ---- response generation with common uniforms ------------------------------
## Category = number of cumulative category probabilities below U.
gen_item <- function(theta, thr, u, a = 1) {
  m <- length(thr)
  # PCM category logits: c * theta - sum of the first c thresholds, times the
  # discrimination a (GPCM). a = 1 is the PCM, identical to the earlier runs.
  eta <- a * (outer(theta, 0:m) - matrix(c(0, cumsum(thr)), length(theta),
                                         m + 1L, byrow = TRUE))
  eta <- eta - apply(eta, 1, max)
  P <- exp(eta) / rowSums(exp(eta))
  cum <- t(apply(P, 1, cumsum))
  rowSums(cum[, -(m + 1L), drop = FALSE] < u)
}

## ---- analysis of one dataset -----------------------------------------------
## CML refit, WLE + SEM, linked to the generating scale via the common items.
analyse <- function(X, gen_loc, common) {
  thr <- .fit_cml_thresholds(X)
  est_loc <- vapply(thr, mean, numeric(1))
  shift <- mean(gen_loc[common]) - mean(est_loc[common])
  th <- .estimate_thetas(X, thr, method = "WLE")
  list(thr = thr, theta = th$theta + shift, sem = th$sem, sum = rowSums(X))
}

cond_infit <- function(X, thr, j) {
  mom <- .cond_moments(thr)
  s <- rowSums(X)
  keep <- s > 0 & s < mom$max_score
  E <- mom$E[s[keep] + 1, j]
  V <- mom$V[s[keep] + 1, j]
  sum((X[keep, j] - E)^2) / sum(V)
}

## ---- one replicate ----------------------------------------------------------
## mech = "xdim": `rho` is the correlation of theta2 with theta1 (underfit).
## mech = "slope": `rho` is the discrimination a of item i (a > 1 overfit),
## item i still responds to theta1, and metric E is not defined.
## mech = "joint": all three target items respond to one shared theta2 with
## correlation `rho` (Part A's design). `tpos` is ignored. Counterfactual: all
## three fitting. Dataset (c) removes all three.
## mech = "ld": uniform response dependence (log-linear Rasch model term,
## Kreiner). Item j, the item nearest item i in location, responds to
## theta1 + lambda * x_i, with `rho` = lambda. Counterfactual: the same pair
## independent. Dataset (c) removes item j. Metric E is not defined.
one_rep <- function(base, tpos, rho, r, n = N_CAL, mech = "xdim") {
  set.seed(1e7 * match(base$design, c("d20", "phq9")) + 1e5 * tpos +
             round(1000 * rho) * 10 + r + 1e9 * (n != N_CAL) + n +
             5e8 * (mech == "slope") + 6e8 * (mech == "joint") +
             7e8 * (mech == "ld"))
  i <- base$targets[tpos]
  if (mech == "joint") i <- base$targets
  if (mech == "ld") {
    j <- setdiff(order(abs(base$loc - base$loc[i])), i)[1]
  }
  k <- length(base$thr)
  mu <- mean(base$thetas); sg <- sd(base$thetas)
  th1 <- sample(base$thetas, n, replace = TRUE)
  ee  <- sample(base$thetas, n, replace = TRUE)
  uses_th2 <- mech %in% c("xdim", "joint")
  th2 <- if (uses_th2) mu + rho * (th1 - mu) + sqrt(1 - rho^2) * (ee - mu) else th1
  # Groups differ by D_GROUP SDs on the item-specific factor e only, so
  # theta1 is equal across groups and theta2 differs by
  # sqrt(1 - rho^2) * D_GROUP SDs (zero at rho = 1).
  G   <- rbinom(n, 1, 0.5)
  th2g <- if (uses_th2) th2 + G * D_GROUP * sg * sqrt(1 - rho^2) else th1
  U <- matrix(runif(n * k), n, k)

  Xb <- vapply(seq_len(k), function(j) gen_item(th1, base$thr[[j]], U[, j]),
               numeric(n))
  Xa <- Xb
  if (mech %in% c("xdim", "joint")) {
    for (ii in i) Xa[, ii] <- gen_item(th2, base$thr[[ii]], U[, ii])
  } else if (mech == "slope") {
    Xa[, i] <- gen_item(th1, base$thr[[i]], U[, i], a = rho)
  } else {
    Xa[, j] <- gen_item(th1 + rho * Xb[, i], base$thr[[j]], U[, j])
  }
  colnames(Xa) <- colnames(Xb) <- colnames(base$X)
  if (uses_th2) {
    Xag <- Xb
    for (ii in i) Xag[, ii] <- gen_item(th2g, base$thr[[ii]], U[, ii])
    colnames(Xag) <- colnames(base$X)
  }

  # Items whose responses change between (a) and (b), and the item removed in (c)
  changed <- if (mech == "ld") c(i, j) else i
  removed <- if (mech == "ld") j else i
  common <- setdiff(seq_len(k), changed)
  keep_c <- setdiff(seq_len(k), removed)
  A  <- analyse(Xa, base$loc, common)
  Bf <- analyse(Xb, base$loc, common)
  C  <- analyse(Xb[, keep_c, drop = FALSE], base$loc[keep_c],
                match(common, keep_c))
  if (uses_th2) AG <- analyse(Xag, base$loc, common)

  focal <- if (mech == "ld") j else i
  msq_each <- vapply(focal, function(q) cond_infit(Xa, A$thr, q), numeric(1))
  msq <- mean(msq_each)
  msq_b <- mean(vapply(focal, function(q) cond_infit(Xb, Bf$thr, q), numeric(1)))
  msq_partner <- if (mech == "ld") cond_infit(Xa, A$thr, i) else NA_real_
  other_msq <- vapply(common, function(q) cond_infit(Xa, A$thr, q), numeric(1))

  rmse <- function(x) sqrt(mean((x - th1)^2))
  shift_sem <- abs(A$theta - Bf$theta) / Bf$sem
  r_se <- (A$sem / Bf$sem)^2
  kl <- 0.5 * ((A$theta - Bf$theta) / Bf$sem)^2 + 0.5 * (r_se - 1 - log(r_se))
  gdiff <- function(x) mean(x[G == 1]) - mean(x[G == 0])
  # D: true error SD over the reported SEM (1 = reported precision honest)
  dsem <- function(x) sd(x$theta - th1) / sqrt(mean(x$sem^2))

  # F: conditional bias by true theta1 (the pool is finite, so condition on
  # its distinct values), summed for pooling over replicates.
  key <- factor(th1, levels = sort(unique(base$thetas)))
  fsum <- data.frame(
    theta = as.numeric(levels(key)),
    n = as.numeric(table(key)),
    bias_a = as.numeric(tapply(A$theta - th1, key, sum, default = 0)),
    bias_b = as.numeric(tapply(Bf$theta - th1, key, sum, default = 0)),
    sem_b  = as.numeric(tapply(Bf$sem, key, sum, default = 0))
  )

  list(
    summary = data.frame(
      design = base$design, target = if (mech == "joint") "joint" else base$target_label[tpos],
      item = i[1], n = n,
      rel_loc = base$rel[i[1]], mech = mech, rho = rho, rep = r,
      msq = msq, msq_b = msq_b,
      msq_1 = msq_each[1],
      msq_2 = if (length(msq_each) >= 2) msq_each[2] else NA_real_,
      msq_3 = if (length(msq_each) >= 3) msq_each[3] else NA_real_,
      msq_partner = msq_partner,
      ld_partner = if (mech == "ld") j else NA_integer_,
      other_msq_min = min(other_msq), other_msq_max = max(other_msq),
      rmse_a = rmse(A$theta), rmse_b = rmse(Bf$theta), rmse_c = rmse(C$theta),
      spear_a = cor(A$theta, th1, method = "spearman"),
      spear_b = cor(Bf$theta, th1, method = "spearman"),
      spear_c = cor(C$theta, th1, method = "spearman"),
      spear_sum_a = cor(A$sum, th1, method = "spearman"),
      spear_sum_b = cor(Bf$sum, th1, method = "spearman"),
      spear_sum_c = cor(C$sum, th1, method = "spearman"),
      shift_mean = mean(shift_sem), shift_gt1 = mean(shift_sem > 1),
      kl_mean = mean(kl),
      kl_shift = mean(0.5 * ((A$theta - Bf$theta) / Bf$sem)^2),
      E_wle = if (uses_th2) (gdiff(AG$theta) - gdiff(Bf$theta)) / sg / D_GROUP else NA_real_,
      E_sum = if (uses_th2) (gdiff(AG$sum) - gdiff(Bf$sum)) / sd(Bf$sum) / D_GROUP else NA_real_,
      D_a = dsem(A), D_b = dsem(Bf)
    ),
    fsum = fsum
  )
}

run_cell <- function(base, tpos, rho, n = N_CAL, mech = "xdim") {
  f <- file.path(OUT_DIR, sprintf("%s_%s_%s%03d%s.rds", base$design,
                                  if (mech == "joint") "joint" else base$target_label[tpos],
                                  switch(mech, xdim = "rho", slope = "slope",
                                         joint = "jointrho", ld = "ldlambda"),
                                  round(100 * rho),
                                  if (n == N_CAL) "" else sprintf("_n%05d", n)))
  if (file.exists(f)) return(invisible(readRDS(f)))
  t0 <- Sys.time()
  res <- mclapply(seq_len(R_CAL), function(r) {
    tryCatch(one_rep(base, tpos, rho, r, n, mech), error = function(e) conditionMessage(e))
  }, mc.cores = NCORES)
  ok <- vapply(res, is.list, logical(1))
  out <- list(
    summary = do.call(rbind, lapply(res[ok], `[[`, "summary")),
    fsum = Reduce(function(x, y) {
      x[, c("n", "bias_a", "bias_b", "sem_b")] <-
        x[, c("n", "bias_a", "bias_b", "sem_b")] +
        y[, c("n", "bias_a", "bias_b", "sem_b")]
      x
    }, lapply(res[ok], `[[`, "fsum")),
    failed = unlist(res[!ok]),
    minutes = as.numeric(Sys.time() - t0, units = "mins"),
    base_info = list(design = base$design, target = base$target_label[tpos], n = n,
                     mech = mech,
                     item = base$targets[tpos], rel_loc = base$rel[base$targets[tpos]],
                     n_items = length(base$thr), theta_sd = sd(base$thetas))
  )
  saveRDS(out, f)
  saveRDS(out, sub("[.]rds$", paste0("_", Sys.Date(), ".rds"), f))
  cat(sprintf("%s %s %s=%.2f n=%d: %d reps, %.1f min, msq %.3f\n", base$design,
              base$target_label[tpos], mech, rho, n, sum(ok), out$minutes,
              mean(out$summary$msq)))
  invisible(out)
}

BASES <- lapply(c("d20", "phq9"), make_base)
names(BASES) <- c("d20", "phq9")
saveRDS(lapply(BASES, function(b) b[c("design", "thr", "loc", "rel", "targets",
                                       "target_label")]),
        file.path(OUT_DIR, "bases.rds"))

if (!identical(Sys.getenv("IMPACT_SMOKE"), "1")) {
  for (b in BASES) for (tpos in 1:3) for (rho in RHOS) run_cell(b, tpos, rho)
  for (n in N_STAB) for (b in BASES) for (rho in RHOS) run_cell(b, 1, rho, n)
  # Overfit arm (added 2026-10-12): item i with discrimination a >= 1.
  for (b in BASES) for (tpos in 1:3) for (a in SLOPES) run_cell(b, tpos, a, mech = "slope")
  # Joint arm (added 2026-10-12): three items misfit together on one theta2.
  for (b in BASES) for (rho in RHOS) run_cell(b, 1, rho, mech = "joint")
  # Local dependence arm (added 2026-10-12): uniform response dependence.
  for (b in BASES) for (tpos in 1:3) for (lam in LAMBDAS) run_cell(b, tpos, lam, mech = "ld")
  cat("done\n")
}
