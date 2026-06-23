# =====================================================================
# Q3 DGP comparison via the per-pair bootstrap p-value + family-wise
# correction (the powerful, localizing way to call LD -- not the global
# max-mean cutoff).
#
# For each DGP ("resample" vs "conditional") this measures, under a KNOWN
# well-fitting model:
#   * FWER  -- P(any item pair flagged) under the null; target = alpha.
# and, under a model with ONE injected dependent pair:
#   * power -- P(the true dependent pair flagged);
#   * off-target FWER -- P(any OTHER pair flagged).
#
# Run from the package root:  Rscript dev/q3_dgp_fwer.R
# (uses pkgload::load_all so the unreleased conditional DGP is available).
# =====================================================================

# ------------------------------- SETTINGS -----------------------------------
n_cores     <- 8L          # parallel workers (outer loop over datasets)
K           <- 200L        # datasets per scenario (FWER precision ~ sqrt(a(1-a)/K))
B           <- 1000L       # inner bootstrap iterations per cut-off
alpha       <- 0.05        # family-wise significance level
correction  <- "fwer"      # "fwer" (Westfall-Young), "fdr_bh", "fdr_by", "none"
N           <- 300L        # respondents per simulated dataset
ld_tau      <- 0.8         # shared-nuisance SD on the dependent pair (effect size)
run_power   <- TRUE        # also run the injected-LD (power) scenario
master_seed <- 2024L
out_file    <- "dev/q3_dgp_fwer_results.rds"
# ----------------------------------------------------------------------------

suppressMessages(pkgload::load_all(".", quiet = TRUE))
library(parallel)

# ---- true (null) model from phq9 -------------------------------------------
load("data/phq9.rda")
obs        <- as.matrix(phq9[, paste0("q", 1:9)])
thr_list   <- .fit_cml_thresholds(obs)
grid       <- seq(-6, 6, length.out = 81L)
sigma      <- .estimate_prior_sd(.grid_loglik(obs, .logp_tables(thr_list, grid), grid),
                                 grid, 0)
k_items    <- length(thr_list)
item_names <- paste0("I", seq_len(k_items))
ld_pair    <- c(1L, 2L)                       # injected dependent pair (items I1,I2)
cat(sprintf("True model: %d items, latent SD = %.3f, N = %d | K=%d B=%d cores=%d\n",
            k_items, sigma, N, K, B, n_cores))

# ---- data generators -------------------------------------------------------
sim_null <- function() {
  D <- sim_partial_score(thr_list, rnorm(N, 0, sigma))
  colnames(D) <- item_names; as.data.frame(D)
}
# inject LD: a shared nuisance u added to the linear predictor of the pair only
sim_ld <- function() {
  th <- rnorm(N, 0, sigma)
  D  <- sim_partial_score(thr_list, th)
  u  <- rnorm(N, 0, ld_tau)
  for (i in ld_pair) D[, i] <- sim_poly_item(thr_list[[i]], th + u)
  colnames(D) <- item_names; as.data.frame(D)
}

# ---- one dataset, one DGP --> flags ----------------------------------------
run_one <- function(D, dgp) {
  cut   <- RMlocdepQ3Cutoff(D, iterations = B, parallel = FALSE, dgp = dgp)
  pairs <- suppressWarnings(
    RMlocdepQ3(D, cutoff = cut, p_value = TRUE, correction = correction,
               alpha = alpha, output = "dataframe")$pairs)
  flagged <- pairs$padj_q3 < alpha & !is.na(pairs$padj_q3)
  is_ld   <- (pairs$Item1 == item_names[ld_pair[1]] &
                pairs$Item2 == item_names[ld_pair[2]]) |
             (pairs$Item1 == item_names[ld_pair[2]] &
                pairs$Item2 == item_names[ld_pair[1]])
  c(any_flag    = as.integer(any(flagged)),
    n_flag      = sum(flagged),
    ld_flag     = as.integer(any(flagged &  is_ld)),   # power (under LD)
    other_flag  = as.integer(any(flagged & !is_ld)))   # off-target false positive
}

# ---- run a scenario across K datasets (parallel over datasets) -------------
run_scenario <- function(gen, seeds) {
  mclapply(seq_along(seeds), function(k) {
    set.seed(seeds[k]); D <- gen()
    lapply(c("resample", "conditional"), function(dgp)
      tryCatch(run_one(D, dgp),
               error = function(e) c(any_flag = NA, n_flag = NA,
                                     ld_flag = NA, other_flag = NA)))
  }, mc.cores = n_cores)
}
summarise <- function(scen) {
  grab <- function(j, field) vapply(scen, function(x) x[[j]][[field]], numeric(1))
  data.frame(
    DGP          = c("resample", "conditional"),
    FWER         = c(mean(grab(1, "any_flag"),   na.rm = TRUE),
                     mean(grab(2, "any_flag"),   na.rm = TRUE)),
    mean_n_flag  = c(mean(grab(1, "n_flag"),     na.rm = TRUE),
                     mean(grab(2, "n_flag"),     na.rm = TRUE)),
    power_ld     = c(mean(grab(1, "ld_flag"),    na.rm = TRUE),
                     mean(grab(2, "ld_flag"),    na.rm = TRUE)),
    offtarget    = c(mean(grab(1, "other_flag"), na.rm = TRUE),
                     mean(grab(2, "other_flag"), na.rm = TRUE))
  )
}

set.seed(master_seed)
null_seeds <- sample.int(.Machine$integer.max, K)
ld_seeds   <- sample.int(.Machine$integer.max, K)

cat("\n== NULL scenario (target FWER = alpha) ==\n")
t_null <- system.time(null_scen <- run_scenario(sim_null, null_seeds))["elapsed"]
null_summary <- summarise(null_scen)
print(null_summary[, c("DGP", "FWER", "mean_n_flag")])
cat(sprintf("(%.0f s; FWER binomial SE ~ %.3f at alpha=%.2f)\n",
            t_null, sqrt(alpha * (1 - alpha) / K), alpha))

ld_scen <- ld_summary <- NULL
if (run_power) {
  cat("\n== INJECTED-LD scenario (pair I1-I2, tau =", ld_tau, ") ==\n")
  t_ld <- system.time(ld_scen <- run_scenario(sim_ld, ld_seeds))["elapsed"]
  ld_summary <- summarise(ld_scen)
  print(ld_summary[, c("DGP", "power_ld", "offtarget", "mean_n_flag")])
  cat(sprintf("(%.0f s)\n", t_ld))
}

saveRDS(list(
  settings = list(n_cores = n_cores, K = K, B = B, alpha = alpha,
                  correction = correction, N = N, ld_tau = ld_tau,
                  sigma = sigma, k_items = k_items, ld_pair = ld_pair,
                  master_seed = master_seed),
  null_scenario = null_scen, null_summary = null_summary,
  ld_scenario   = ld_scen,   ld_summary   = ld_summary
), out_file)
cat("\nSaved", out_file, "\n")
