## S-X2 against conditional infit, revised
##
## Design: dev/sx2-comparison-v2-design.md. Replaces the comparison in
## the July 2026 S-X2 study (removed, see git history). All methods are run on the same freshly generated
## dataset in every replication:
##   infit       RMitemInfitCutoff() B = 400 + RMitemInfit(p_value = TRUE),
##               WY (padj_infit) and BH on the marginal bootstrap p
##   sx2         mirt S-X2, p on mirt's df and on df + 1
##   sx2_boot    parametric bootstrap of S-X2, B = 400, statistic
##               z = qnorm(1 - p_df1), upper-tail WY via .bootstrap_pvalues()
##   rmsea_boot  same draws, RMSEA.S_X2 as the statistic
##
## Run from easyRasch2/dev: source("sx2_comparison_v2.R"). About 6 to 7 hours
## on 10 cores. Each cell is cached in sx2_comparison_v2/ as soon as it
## finishes, so an interrupted run resumes at the next unfinished cell.
## Progress is written to sx2_comparison_v2/progress.log.
##
## Smoke test: Sys.setenv(SX2V2_SMOKE = "1") before sourcing runs 2
## replications of two cells with B = 20 into sx2_comparison_v2_smoke/.

suppressMessages({
  library(easyRasch2)
  library(mirt)
  library(parallel)
})

SMOKE    <- identical(Sys.getenv("SX2V2_SMOKE"), "1")
B        <- if (SMOKE) 20L else 400L
NCORES   <- 10L
THETA_SD <- 1.5
OUT_DIR  <- if (SMOKE) "sx2_comparison_v2_smoke" else "sx2_comparison_v2"

## ---- items -------------------------------------------------------------
N_ITEMS <- 20L
LOCS    <- seq(-2, 2, length.out = N_ITEMS)
# misfitting items nearest 0, -1 and -2 logits, as in the original study
MISFIT_ORDER <- integer(0)
for (tg in c(0, -1, -2)) {
  MISFIT_ORDER <- c(MISFIT_ORDER,
                    setdiff(order(abs(LOCS - tg)), MISFIT_ORDER)[1L])
}
ITEM_NAMES <- paste0("I", seq_len(N_ITEMS))

## ---- cells -------------------------------------------------------------
NS <- c(150, 500, 1000, 2000)
cells <- rbind(
  expand.grid(k = 0L, rho = NA_real_, n = NS, reps = 500L),
  expand.grid(k = 1L, rho = c(0.50, 0.75), n = NS, reps = 200L),
  expand.grid(k = 3L, rho = c(0.15, 0.50, 0.75), n = NS, reps = 200L)
)
cells$id <- with(cells, sprintf("k%d_rho%s_n%d", k,
                                ifelse(is.na(rho), "NA", sprintf("%.2f", rho)), n))
if (SMOKE) {
  cells <- cells[cells$id %in% c("k0_rhoNA_n150", "k3_rho0.50_n500"), ]
  cells$reps <- 2L
}

## ---- data generation -----------------------------------------------------
# Fitting items respond to theta1, misfitting items to theta2. Both SD 1.5,
# correlation rho. Dichotomous Rasch responses.
gen_data <- function(n, k, rho) {
  th1 <- rnorm(n, 0, THETA_SD)
  th <- matrix(th1, n, N_ITEMS)
  if (k > 0L) {
    th2 <- rho * th1 + sqrt(1 - rho^2) * rnorm(n, 0, THETA_SD)
    th[, MISFIT_ORDER[seq_len(k)]] <- th2
  }
  p <- plogis(th - matrix(LOCS, n, N_ITEMS, byrow = TRUE))
  d <- as.data.frame((matrix(runif(n * N_ITEMS), n) < p) * 1L)
  names(d) <- ITEM_NAMES
  d
}

both_categories <- function(d) all(vapply(d, function(x) {
  s <- sum(x)
  s > 0 && s < length(x)
}, logical(1)))

## ---- S-X2 helpers --------------------------------------------------------
fit_rasch <- function(d) mirt(d, 1, itemtype = "Rasch", verbose = FALSE)

sx2_stats <- function(m) {
  f <- itemfit(m, fit_stats = "S_X2")
  S <- f$S_X2
  df <- f$df.S_X2
  p_df1 <- pchisq(S, df + 1, lower.tail = FALSE)
  data.frame(sx2 = S, df_mirt = df, p_mirt = f$p.S_X2, p_df1 = p_df1,
             rmsea = f$RMSEA.S_X2)
}

# z on a common scale whatever the df after collapsing
p_to_z <- function(p) qnorm(1 - pmin(pmax(p, 1e-12), 1 - 1e-12))

sx2_bootstrap <- function(m0, n) {
  cf <- coef(m0, simplify = TRUE)
  dpar <- cf$items[, "d"]                # P(x = 1) = plogis(theta + d)
  sd_hat <- sqrt(cf$cov[1, 1])
  z_sim <- rmsea_sim <- matrix(NA_real_, B, N_ITEMS,
                               dimnames = list(NULL, ITEM_NAMES))
  for (b in seq_len(B)) {
    th <- rnorm(n, 0, sd_hat)
    p <- plogis(outer(th, dpar, "+"))
    db <- as.data.frame((matrix(runif(n * N_ITEMS), n) < p) * 1L)
    names(db) <- ITEM_NAMES
    if (!both_categories(db)) next
    st <- tryCatch(sx2_stats(fit_rasch(db)), error = function(e) NULL)
    if (is.null(st)) next
    z_sim[b, ] <- p_to_z(st$p_df1)
    rmsea_sim[b, ] <- st$rmsea
  }
  list(z = z_sim, rmsea = rmsea_sim,
       dropped = sum(rowSums(is.finite(z_sim)) == 0))
}

## ---- one replication ---------------------------------------------------
one_rep <- function(rep, cell) {
  seed <- 1e6 * match(cell$id, cells$id) + rep
  set.seed(seed)
  d <- gen_data(cell$n, cell$k, cell$rho)
  if (!both_categories(d)) return(list(res = NULL, reason = "constant item"))

  t0 <- Sys.time()
  out <- data.frame(item = seq_len(N_ITEMS), location = LOCS,
                    misfit = seq_len(N_ITEMS) %in% MISFIT_ORDER[seq_len(cell$k)])

  # conditional infit, package defaults
  fi <- tryCatch({
    co <- RMitemInfitCutoff(d, iterations = B, parallel = FALSE,
                            seed = seed, verbose = FALSE)
    RMitemInfit(d, cutoff = co, p_value = TRUE, correction = "fwer",
                output = "dataframe")
  }, error = function(e) NULL)
  if (is.null(fi)) {
    out$infit_msq <- out$p_infit <- out$padj_infit <- NA_real_
  } else {
    i <- match(ITEM_NAMES, fi$Item)
    out$infit_msq  <- fi$Infit_MSQ[i]
    out$p_infit    <- fi$p_infit[i]
    out$padj_infit <- fi$padj_infit[i]
  }

  # S-X2, asymptotic and bootstrap
  sx <- tryCatch({
    m0 <- fit_rasch(d)
    st <- sx2_stats(m0)
    bt <- sx2_bootstrap(m0, cell$n)
    z_obs <- setNames(p_to_z(st$p_df1), ITEM_NAMES)
    r_obs <- setNames(st$rmsea, ITEM_NAMES)
    pz <- easyRasch2:::.bootstrap_pvalues(z_obs, bt$z, "fwer", "upper")
    pr <- easyRasch2:::.bootstrap_pvalues(r_obs, bt$rmsea, "fwer", "upper")
    cbind(st, p_sx2boot = pz$p, padj_sx2boot = pz$padj,
          p_rmseaboot = pr$p, padj_rmseaboot = pr$padj,
          boot_dropped = bt$dropped)
  }, error = function(e) NULL)
  if (is.null(sx)) {
    for (v in c("sx2", "df_mirt", "p_mirt", "p_df1", "rmsea", "p_sx2boot",
                "padj_sx2boot", "p_rmseaboot", "padj_rmseaboot",
                "boot_dropped")) out[[v]] <- NA_real_
  } else {
    out <- cbind(out, sx)
  }

  list(res = out, secs = as.numeric(Sys.time() - t0, units = "secs"),
       infit_failed = is.null(fi), sx2_failed = is.null(sx))
}

## ---- run, one cached file per cell ---------------------------------------
dir.create(OUT_DIR, showWarnings = FALSE)
log_file <- file.path(OUT_DIR, "progress.log")
logline <- function(...) {
  msg <- paste0(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "  ", ...)
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

logline("start: ", nrow(cells), " cells, B = ", B,
        ", easyRasch2 ", as.character(packageVersion("easyRasch2")),
        ", mirt ", as.character(packageVersion("mirt")))

for (ci in seq_len(nrow(cells))) {
  cell <- cells[ci, ]
  f <- file.path(OUT_DIR, paste0(cell$id, ".rds"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  res <- mclapply(seq_len(cell$reps), one_rep, cell = cell,
                  mc.cores = NCORES, mc.preschedule = FALSE)
  saveRDS(list(cell = cell, res = res, B = B, locs = LOCS,
               misfit_order = MISFIT_ORDER,
               easyRasch2 = as.character(packageVersion("easyRasch2")),
               mirt = as.character(packageVersion("mirt")),
               elapsed_min = as.numeric(Sys.time() - t0, units = "mins")), f)
  n_err <- sum(vapply(res, function(r) inherits(r, "try-error") ||
                        is.null(r$res), logical(1)))
  logline(sprintf("cell %d/%d %s done in %.1f min, %d reps without results",
                  ci, nrow(cells), cell$id,
                  as.numeric(Sys.time() - t0, units = "mins"), n_err))
}
logline("all cells done")
