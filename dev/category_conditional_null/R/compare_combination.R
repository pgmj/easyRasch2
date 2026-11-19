# Compare two ways of combining the events of an item (2026-11-19), on the
# stored MCB values of the calibration rerun (results/cal_cells, boot_mcb).
#   minp: single-step min-p over the events, as in analyse_calibration().
#   maxT: each event's MCB studentised by the mean and SD of the pool
#         (observed + B simulated), item statistic = max over events, item
#         p = share of the pool at or above it. As .bootstrap_pvalues() in
#         easyRasch2 with moments = "joint".
# Across items both use the Westfall-Young step-down: on -item p for minp,
# on the item statistic for maxT. No new simulation.
#   Rscript R/compare_combination.R   (CORES=4 default)
# Writes results/compare_combination.rds (+ dated copy).

library(parallel)
cores <- as.integer(Sys.getenv("CORES", "4"))
alpha <- 0.05

stepdown <- function(s_obs, s_sim) {     # larger = more extreme; s_sim B x k
  B <- nrow(s_sim)
  ord <- order(s_obs, decreasing = TRUE)
  run <- rep(-Inf, B)
  p <- numeric(length(ord))
  for (j in rev(seq_along(ord))) {
    run <- pmax(run, s_sim[, ord[j]])
    p[j] <- (1 + sum(run >= s_obs[ord[j]])) / (B + 1)
  }
  out <- numeric(length(ord))
  out[ord] <- pmin(cummax(p), 1)
  out
}

combine <- function(obs, boot) {
  k <- nrow(obs); m <- ncol(obs); N <- dim(boot)[3] + 1L
  pool <- array(c(obs, boot), c(k, m, N))
  PI_min <- PI_max <- T_max <- matrix(NA_real_, N, k)
  for (i in seq_len(k)) {
    E <- matrix(pool[i, , ], nrow = m)                      # m x N
    Pev <- t(apply(E, 1, function(t) (N - rank(t, ties.method = "min") + 1) / N))
    Pev <- matrix(Pev, nrow = m)
    mp <- apply(Pev, 2, min)
    PI_min[, i] <- rank(mp, ties.method = "max") / N
    Z <- t(apply(E, 1, function(t) {
      s <- stats::sd(t)
      if (!is.finite(s) || s < 1e-12) rep(0, length(t)) else (t - mean(t)) / s
    }))
    Z <- matrix(Z, nrow = m)
    T_max[, i] <- apply(Z, 2, max)
    PI_max[, i] <- (N - rank(T_max[, i], ties.method = "min") + 1) / N
  }
  list(p_min = PI_min[1, ], p_max = PI_max[1, ],
       padj_min = pmax(stepdown(-PI_min[1, ], -PI_min[-1, , drop = FALSE]), PI_min[1, ]),
       padj_max = pmax(stepdown(T_max[1, ], T_max[-1, , drop = FALSE]), PI_max[1, ]),
       N = N)
}

files <- list.files("results/cal_cells", pattern = "_k(5|10)\\.rds$", full.names = TRUE)
rows <- mclapply(files, function(f) {
  x <- readRDS(f); cl <- x$cell
  do.call(rbind, lapply(x$results, function(r) {
    if (is.null(r) || !is.na(r$skipped) || is.null(r$boot_mcb) || r$B_used == 0) return(NULL)
    cb <- combine(r$mcb, r$boot_mcb)
    data.frame(cond = cl$cond, null = cl$null, n = cl$n, k = cl$k, rep = r$rep,
               B_used = r$B_used, item = seq_along(cb$p_min),
               p_min = cb$p_min, p_max = cb$p_max,
               padj_min = cb$padj_min, padj_max = cb$padj_max,
               p_min_stored = r$p_item)
  }))
}, mc.cores = cores)
res <- do.call(rbind, rows)
stopifnot(nrow(res) > 0)
# boot_mcb is stored to 7 significant digits, so ties can differ from the
# unrounded values the run used. Report the agreement instead of asserting it.
d <- abs(res$p_min - res$p_min_stored)
message("min-p vs stored: identical in ", round(100 * mean(d < 1e-12), 2), "% of items, max diff ", signif(max(d), 3),
        ", rejections at .05 differ in ", sum((res$p_min < alpha) != (res$p_min_stored < alpha)), " items")
out <- list(items = res, alpha = alpha, date = Sys.Date())
saveRDS(out, "results/compare_combination.rds")
saveRDS(out, sprintf("results/compare_combination_%s.rds", Sys.Date()))
message("item rows: ", nrow(res), ", datasets: ", nrow(unique(res[, c("cond", "n", "k", "rep")])))
