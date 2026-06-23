# =====================================================================
# Infit DGP comparison via the per-ITEM bootstrap p-value + family-wise
# correction (infit is judged per item -- expected infit depends on item
# location -- never by a global cutoff).
#
# For each DGP ("resample" vs "conditional") this measures, under a KNOWN
# well-fitting model:
#   * FWER  -- P(any item flagged) under the null; target = alpha.
# and, under a model with ONE injected misfitting (noisy/underfit) item:
#   * power -- P(the true misfitting item flagged);
#   * off-target FWER -- P(any OTHER item flagged).
#
# Run from the package root:  Rscript dev/infit_dgp_fwer.R
# =====================================================================

# ------------------------------- SETTINGS -----------------------------------
n_cores     <- 8L          # parallel workers (outer loop over datasets)
K           <- 200L        # datasets per scenario
B           <- 1000L       # inner bootstrap iterations per cut-off
alpha       <- 0.05        # family-wise significance level
correction  <- "fwer"      # "fwer" (Westfall-Young), "fdr_bh", "fdr_by", "none"
N           <- 300L        # respondents per simulated dataset
noise_p     <- 0.30        # fraction of responses randomised on the misfit item
run_power   <- TRUE        # also run the injected-misfit (power) scenario
master_seed <- 2024L
out_file    <- "dev/infit_dgp_fwer_results.rds"
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
k_items     <- length(thr_list)
item_names  <- paste0("I", seq_len(k_items))
misfit_item <- 5L                              # injected noisy/underfit item (I5)
cat(sprintf("True model: %d items, latent SD = %.3f, N = %d | K=%d B=%d cores=%d\n",
            k_items, sigma, N, K, B, n_cores))

# ---- data generators -------------------------------------------------------
sim_null <- function() {
  D <- sim_partial_score(thr_list, rnorm(N, 0, sigma))
  colnames(D) <- item_names; as.data.frame(D)
}
# inject underfit: randomise a fraction of the target item's responses (noise)
sim_misfit <- function() {
  D    <- sim_partial_score(thr_list, rnorm(N, 0, sigma))
  m    <- length(thr_list[[misfit_item]])      # max category index
  flip <- stats::runif(N) < noise_p
  D[flip, misfit_item] <- sample(0:m, sum(flip), replace = TRUE)
  colnames(D) <- item_names; as.data.frame(D)
}

# ---- one dataset, one DGP --> flags ----------------------------------------
run_one <- function(D, dgp) {
  cut <- RMitemInfitCutoff(D, iterations = B, parallel = FALSE, dgp = dgp)
  tab <- suppressWarnings(
    RMitemInfit(D, cutoff = cut, p_value = TRUE, correction = correction,
                alpha = alpha, output = "dataframe"))
  flagged <- tab$padj_infit < alpha & !is.na(tab$padj_infit)
  is_tgt  <- tab$Item == item_names[misfit_item]
  c(any_flag   = as.integer(any(flagged)),
    n_flag     = sum(flagged),
    tgt_flag   = as.integer(any(flagged &  is_tgt)),   # power (under misfit)
    other_flag = as.integer(any(flagged & !is_tgt)))   # off-target false positive
}

run_scenario <- function(gen, seeds) {
  mclapply(seq_along(seeds), function(k) {
    set.seed(seeds[k]); D <- gen()
    lapply(c("resample", "conditional"), function(dgp)
      tryCatch(run_one(D, dgp),
               error = function(e) c(any_flag = NA, n_flag = NA,
                                     tgt_flag = NA, other_flag = NA)))
  }, mc.cores = n_cores)
}
summarise <- function(scen) {
  grab <- function(j, field) vapply(scen, function(x) x[[j]][[field]], numeric(1))
  data.frame(
    DGP         = c("resample", "conditional"),
    FWER        = c(mean(grab(1, "any_flag"),   na.rm = TRUE),
                    mean(grab(2, "any_flag"),   na.rm = TRUE)),
    mean_n_flag = c(mean(grab(1, "n_flag"),     na.rm = TRUE),
                    mean(grab(2, "n_flag"),     na.rm = TRUE)),
    power_item  = c(mean(grab(1, "tgt_flag"),   na.rm = TRUE),
                    mean(grab(2, "tgt_flag"),   na.rm = TRUE)),
    offtarget   = c(mean(grab(1, "other_flag"), na.rm = TRUE),
                    mean(grab(2, "other_flag"), na.rm = TRUE))
  )
}

set.seed(master_seed)
null_seeds   <- sample.int(.Machine$integer.max, K)
misfit_seeds <- sample.int(.Machine$integer.max, K)

cat("\n== NULL scenario (target FWER = alpha) ==\n")
t_null <- system.time(null_scen <- run_scenario(sim_null, null_seeds))["elapsed"]
null_summary <- summarise(null_scen)
print(null_summary[, c("DGP", "FWER", "mean_n_flag")])
cat(sprintf("(%.0f s; FWER binomial SE ~ %.3f at alpha=%.2f)\n",
            t_null, sqrt(alpha * (1 - alpha) / K), alpha))

misfit_scen <- misfit_summary <- NULL
if (run_power) {
  cat("\n== INJECTED-MISFIT scenario (item I5, noise_p =", noise_p, ") ==\n")
  t_m <- system.time(misfit_scen <- run_scenario(sim_misfit, misfit_seeds))["elapsed"]
  misfit_summary <- summarise(misfit_scen)
  print(misfit_summary[, c("DGP", "power_item", "offtarget", "mean_n_flag")])
  cat(sprintf("(%.0f s)\n", t_m))
}

saveRDS(list(
  settings = list(n_cores = n_cores, K = K, B = B, alpha = alpha,
                  correction = correction, N = N, noise_p = noise_p,
                  sigma = sigma, k_items = k_items, misfit_item = misfit_item,
                  master_seed = master_seed),
  null_scenario   = null_scen,   null_summary   = null_summary,
  misfit_scenario = misfit_scen, misfit_summary = misfit_summary
), out_file)
cat("\nSaved", out_file, "\n")
