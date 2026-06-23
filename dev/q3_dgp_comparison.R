# =====================================================================
# Comparative calibration study: Q3 cut-off DGPs ("resample" vs "conditional")
#
# Question: which data-generating process for RMlocdepQ3Cutoff() yields a
# correctly calibrated null cut-off for the global Q3 statistic (max - mean)?
#
# Design: define a KNOWN well-fitting model (PCM fitted to phq9, normal latent
# with MML-estimated SD, NO local dependence). Data simulated from it are null
# by construction, so the TRUE null distribution of the statistic is obtained
# by Monte-Carlo. We then ask, for each DGP:
#   (1) does its simulated null distribution match the true null? (Part 1)
#   (2) is its 99th-percentile cut-off unbiased for the true 99th percentile,
#       and does it control the false-positive rate? (Part 2)
# =====================================================================
suppressMessages(pkgload::load_all(".", quiet = TRUE))

set.seed(2024)

# ---- knobs (modest defaults; raise for a publication-grade run) ------------
N        <- 300L     # sample size of each simulated dataset
M_truth  <- 1500L    # datasets to estimate the TRUE null distribution
K        <- 20L      # datasets to estimate each DGP's cut-off bias/variance
B        <- 250L     # inner bootstrap iterations per cut-off

# ---- define the true (null) model from phq9 --------------------------------
load("data/phq9.rda")
obs <- as.matrix(phq9[, paste0("q", 1:9)])
thr_list <- .fit_cml_thresholds(obs)              # true item thresholds (CML)

# true latent SD via marginal ML
grid   <- seq(-6, 6, length.out = 81L)
loglik <- .grid_loglik(obs, .logp_tables(thr_list, grid), grid)
sigma  <- .estimate_prior_sd(loglik, grid, 0)
cat(sprintf("True model: %d items, latent SD = %.3f, N = %d per dataset\n\n",
            length(thr_list), sigma, N))

# simulate one null dataset from the true model (normal latent, no LD)
sim_null <- function() {
  theta <- stats::rnorm(N, 0, sigma)
  as.data.frame(sim_partial_score(thr_list, theta))
}
# global Q3 statistic, computed exactly as RMlocdepQ3Cutoff() does internally
stat_maxmean <- function(D) {
  m <- .q3_residual_matrix(D, estimator = "CML")  # NA diagonal
  max(m, na.rm = TRUE) - mean(m, na.rm = TRUE)
}

# ---- Part 1: TRUE null distribution of the statistic -----------------------
cat("Part 1: estimating the true null distribution (", M_truth, "datasets)...\n")
truth <- vapply(seq_len(M_truth), function(i) stat_maxmean(sim_null()), numeric(1))
P_true <- stats::quantile(truth, c(0.90, 0.95, 0.99), names = TRUE)
cat(sprintf("  true null max-mean Q3: mean=%.3f  P95=%.3f  P99=%.3f\n\n",
            mean(truth), P_true[["95%"]], P_true[["99%"]]))

# one representative dataset -> each DGP's simulated null distribution
D0 <- sim_null()
cut_r0 <- RMlocdepQ3Cutoff(D0, iterations = 1000, parallel = FALSE,
                           seed = 1, dgp = "resample")
cut_c0 <- RMlocdepQ3Cutoff(D0, iterations = 1000, parallel = FALSE,
                           seed = 1, dgp = "conditional")
cat("One dataset, 1000 inner iterations -- simulated null vs truth:\n")
cmp <- rbind(
  truth       = P_true,
  resample    = stats::quantile(cut_r0$results$diff, c(.90, .95, .99)),
  conditional = stats::quantile(cut_c0$results$diff, c(.90, .95, .99))
)
print(round(cmp, 3))

# ---- Part 2: cut-off bias + realized false-positive rate -------------------
cat(sprintf("\nPart 2: cut-off bias over %d datasets (B = %d inner)...\n", K, B))
res <- data.frame(stat = numeric(K), cut_r = numeric(K), cut_c = numeric(K))
for (k in seq_len(K)) {
  Dk <- sim_null()
  res$stat[k]  <- stat_maxmean(Dk)
  res$cut_r[k] <- RMlocdepQ3Cutoff(Dk, iterations = B, parallel = FALSE,
                                   seed = k, dgp = "resample")$suggested_cutoff
  res$cut_c[k] <- RMlocdepQ3Cutoff(Dk, iterations = B, parallel = FALSE,
                                   seed = k, dgp = "conditional")$suggested_cutoff
  cat(sprintf("  %2d/%d\r", k, K))
}
cat("\n")

p99 <- P_true[["99%"]]
summ <- data.frame(
  DGP        = c("resample", "conditional"),
  mean_cutoff = c(mean(res$cut_r), mean(res$cut_c)),
  sd_cutoff   = c(stats::sd(res$cut_r), stats::sd(res$cut_c)),
  bias_vs_P99 = c(mean(res$cut_r) - p99, mean(res$cut_c) - p99),
  # realized false-positive rate: fresh-null stat exceeding the DGP cut-off
  realized_fpr = c(mean(res$stat > res$cut_r), mean(res$stat > res$cut_c))
)
cat(sprintf("\nTrue P99 = %.3f. Nominal FPR at the 99th-pctl cut-off = 0.01.\n", p99))
print(round.data.frame <- within(summ, {
  mean_cutoff  <- round(mean_cutoff, 3)
  sd_cutoff    <- round(sd_cutoff, 3)
  bias_vs_P99  <- round(bias_vs_P99, 3)
  realized_fpr <- round(realized_fpr, 3)
}))

saveRDS(list(truth = truth, P_true = P_true, cut_r0 = cut_r0$results$diff,
             cut_c0 = cut_c0$results$diff, per_dataset = res, summary = summ),
        "dev/q3_dgp_comparison_results.rds")
cat("\nSaved dev/q3_dgp_comparison_results.rds\n")
