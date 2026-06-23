# =====================================================================
# At what N do asymptotic CML threshold SEs become unreliable, so that the
# fractional weighted bootstrap (FWB) is worth recommending?
#
# Truth (Monte-Carlo): simulate M datasets from a fixed PCM at sample size N,
# refit by CML (psychotools), and take the SD of the centred thresholds across
# datasets -> the true finite-sample SE. Then compare, per threshold:
#   * asymptotic SE  (psychotools threshpar vcov, centred)  + Wald CI coverage
#   * FWB SE/CI      (exp / beta / power weights)           + percentile coverage
# against that truth, on-target and off-target (a shifted sample -> sparse
# extreme categories, where FWB should help most).
#
# Decision rule: recommend FWB where asymptotic 95% coverage drops materially
# below nominal (say < 0.93) while FWB coverage stays near 0.95.
#
# Run from the package root (8 cores):  Rscript dev/fwb_threshold_se_simN.R
# =====================================================================

# ------------------------------- SETTINGS -----------------------------------
n_cores   <- 8L
N_grid    <- c(50L, 80L, 120L, 200L, 350L, 600L, 1000L)
M_truth   <- 1500L                 # datasets for the true SD + asymptotic coverage
m_fwb     <- 200L                  # subset of those also given an FWB run
R         <- 400L                  # FWB weight draws per dataset
fwb_types <- c("exp", "beta", "power")
mu_off    <- 1.0                   # off-target latent-mean shift (logits)
conf      <- 0.95
seed0     <- 2024L
out_file  <- "dev/fwb_threshold_se_simN_results.rds"
# ----------------------------------------------------------------------------

suppressMessages({ pkgload::load_all(".", quiet = TRUE); library(psychotools) })
library(parallel)

# ---- true model from phq9 (9 items x 4 categories) -------------------------
load("data/phq9.rda")
obs       <- as.matrix(phq9[, paste0("q", 1:9)])
thr_list  <- .fit_cml_thresholds(obs)           # centred Andrich thresholds (truth)
true_thr  <- unlist(thr_list, use.names = FALSE)
grid      <- seq(-6, 6, length.out = 81L)
sigma     <- .estimate_prior_sd(.grid_loglik(obs, .logp_tables(thr_list, grid), grid),
                                grid, 0)
P <- length(true_thr)
z <- stats::qnorm(1 - (1 - conf) / 2)
cat(sprintf("Truth: %d items, %d thresholds, latent SD=%.2f, mu_off=%.1f\n\n",
            length(thr_list), P, sigma, mu_off))

# ---- FWB weights (mean 1), centred-threshold helpers -----------------------
fwb_w <- function(n, type) {
  w <- switch(type,
    exp   = stats::rexp(n),
    beta  = 4 * stats::rbeta(n, 0.5, 1.5),
    power = (2 + sqrt(2)) * stats::rbeta(n, sqrt(2) - 1, 1))
  w / mean(w)
}
centred_thr  <- function(fit) {                  # centred Andrich threshold vector
  t <- unlist(lapply(psychotools::threshpar(fit), as.numeric), use.names = FALSE)
  t - mean(t)
}
centred_se   <- function(fit) {                  # asymptotic SE of centred thresholds
  V <- attr(psychotools::threshpar(fit, vcov = TRUE), "vcov")
  C <- diag(P) - matrix(1 / P, P, P)
  sqrt(diag(C %*% V %*% t(C)))
}
sim_data <- function(N, mu) {
  D <- sim_partial_score(thr_list, stats::rnorm(N, mu, sigma))
  as.data.frame(D)
}

# ---- one Monte-Carlo dataset: estimate, asymptotic SE/coverage, (FWB) ------
# Fully fault-tolerant: returns NULL if the dataset is degenerate (e.g. a null
# category makes threshpar return fewer than P thresholds), so such datasets
# are simply discarded. A safe quantile guards all-NA FWB columns.
sq <- function(x, p) if (all(is.na(x))) NA_real_ else stats::quantile(x, p, na.rm = TRUE)
one_dataset <- function(seed, N, mu, do_fwb) {
  tryCatch({
    set.seed(seed)
    D   <- sim_data(N, mu)
    fit <- psychotools::pcmodel(D)
    est <- centred_thr(fit)
    if (length(est) != P) return(NULL)            # null category -> discard
    ase <- tryCatch(centred_se(fit), error = function(e) rep(NA_real_, P))
    res <- list(est = est, ase = ase,
                acov = as.integer(abs(true_thr - est) <= z * ase))
    if (do_fwb) {
      res$fwb_cov <- tryCatch({
        fc <- lapply(fwb_types, function(ty) {
          draws <- vapply(seq_len(R), function(r) {
            f <- tryCatch(psychotools::pcmodel(D, weights = fwb_w(nrow(D), ty)),
                          error = function(e) NULL)
            if (is.null(f)) return(rep(NA_real_, P))
            ct <- centred_thr(f); if (length(ct) != P) rep(NA_real_, P) else ct
          }, numeric(P))                            # P x R
          lo <- apply(draws, 1L, sq, (1 - conf) / 2)
          hi <- apply(draws, 1L, sq, 1 - (1 - conf) / 2)
          list(cov = as.integer(true_thr >= lo & true_thr <= hi),
               sd  = apply(draws, 1L, stats::sd, na.rm = TRUE))
        })
        names(fc) <- fwb_types; fc
      }, error = function(e) NULL)
    }
    res
  }, error = function(e) NULL)
}

run_cell <- function(N, mu) {
  seeds <- sample.int(.Machine$integer.max, M_truth)
  rows  <- mclapply(seq_len(M_truth), function(k)
    one_dataset(seeds[k], N, mu, do_fwb = (k <= m_fwb)),
    mc.cores = n_cores)
  rows <- rows[vapply(rows, function(x) is.list(x) && !is.null(x$est), logical(1))]
  est  <- do.call(rbind, lapply(rows, `[[`, "est"))     # n x P
  true_sd   <- apply(est, 2L, stats::sd)
  asy_se    <- colMeans(do.call(rbind, lapply(rows, `[[`, "ase")), na.rm = TRUE)
  asy_cov   <- colMeans(do.call(rbind, lapply(rows, `[[`, "acov")), na.rm = TRUE)
  fwb_rows  <- Filter(function(x) !is.null(x$fwb_cov), rows)
  fwb_cov   <- sapply(fwb_types, function(ty)
    mean(colMeans(do.call(rbind, lapply(fwb_rows, function(x) x$fwb_cov[[ty]]$cov)))))
  fwb_se    <- sapply(fwb_types, function(ty)
    mean(colMeans(do.call(rbind, lapply(fwb_rows, function(x) x$fwb_cov[[ty]]$sd)))))
  list(N = N, mu = mu, n_ok = length(rows), n_fwb = length(fwb_rows),
       true_sd = mean(true_sd), asy_se = mean(asy_se),
       asy_cov = mean(asy_cov), fwb_cov = fwb_cov, fwb_se = fwb_se)
}

# ---- run the grid ----------------------------------------------------------
set.seed(seed0)
cells <- list()
for (mu in c(0, mu_off)) for (N in N_grid) {
  t <- system.time(cell <- run_cell(N, mu))["elapsed"]
  cells[[length(cells) + 1L]] <- cell
  cat(sprintf("mu=%.1f N=%4d | true_sd=%.3f asy_se=%.3f | cov: asy=%.3f exp=%.3f beta=%.3f power=%.3f | %.0fs\n",
              mu, N, cell$true_sd, cell$asy_se, cell$asy_cov,
              cell$fwb_cov[["exp"]], cell$fwb_cov[["beta"]], cell$fwb_cov[["power"]], t))
}

saveRDS(list(settings = list(N_grid = N_grid, M_truth = M_truth, m_fwb = m_fwb,
                             R = R, fwb_types = fwb_types, mu_off = mu_off,
                             conf = conf, sigma = sigma, P = P),
             cells = cells), out_file)
cat("\nSaved", out_file, "\n")
