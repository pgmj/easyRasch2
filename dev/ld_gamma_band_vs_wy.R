# =====================================================================
# RMlocdepGamma(): flagging against the interval versus the Westfall-Young
# corrected p-value.
#
# The partial gamma counterpart of dev/ld_q3_band_vs_wy.R, written the
# same day and for the same question. It matters separately because the
# companion study rasch_q3fwer measured partial gamma only through its
# ASYMPTOTIC rule (iarm's two-sided test over a family of 72), never
# through the interval that RMlocdepGammaCutoff() produces, so neither
# half of the package's own simulated path had been measured.
#
# Both rules are read from the SAME cutoff object per replication, so the
# comparison is paired. Conditions and pair bookkeeping are identical to
# the Q3 script: complete null, and one weak dependent pair (q1 + q2)
# giving 1 dependent pair, 14 sharing an item and 21 disjoint.
#
# Note that the tested statistic is `gamma_pair`, the larger of the pair's
# two rest-score directions, on both rules. The two direction tables carry
# the same flag for a pair, so only direction 1 is scored here.
#
# `hdci_width` and `iterations` are passed explicitly rather than left to
# the defaults, so the script keeps reproducing the pre-1.2.0 settings it
# was written to measure.
#
# RESULTS (n = 400, B = 250, hdci_width = 0.99, 500 replications each)
#
#   complete null                            interval      WY
#     at least one pair flagged                 .430      .044
#       of which above the upper bound          .268         -
#       of which below the lower bound          .236         -
#     mean pairs flagged per analysis            0.55      0.04
#   (1 - .99^36 = .304 predicted for the interval)
#
#   one weak dependent pair
#     detects the dependent pair                 .908      .788
#     >=1 false flag among the 21 disjoint       .278      .038
#     >=1 flag among the 14 sharing an item      .332      .004
#
# The interval overshoots its own prediction by more than the Q3 one does
# (.430 against .344), which is the .99 interval being estimated from 250
# draws rather than 500. The .332 for pairs sharing an item is the
# spurious negative halo of Kreiner & Christensen (2011, Example 3),
# reproduced here in the package. Re-run at the 1.2.0 defaults (B = 400,
# hdci_width = 0.95, p-values on) the familywise error rate under the
# complete null is .048 (SE .014).
#
# Run from the package root:  Rscript dev/ld_gamma_band_vs_wy.R
# =====================================================================

# ------------------------------- SETTINGS -----------------------------------
N            <- 400L    # respondents per simulated dataset
REPS         <- 500L    # replications per condition
ITER         <- 250L    # bootstrap iterations (the pre-1.2.0 gamma default)
HDCI         <- 0.99    # interval width (the pre-1.2.0 default)
RHO          <- 0.282   # copula correlation, calibrated to Q3 = .08
DEP          <- c(1L, 2L)   # the dependent pair, q1 + q2
N_CORES      <- 10L
MASTER_SEED  <- 91260822L
OUT_FILE     <- "dev/ld_gamma_band_vs_wy_results.rds"
# ----------------------------------------------------------------------------

suppressMessages(pkgload::load_all(".", quiet = TRUE))
source("dev/ld_flagging_helpers.R")
library(parallel)

arm <- phq9_arm()
dep_names <- names(arm$thr)[DEP]

one <- function(task) {
  set.seed(MASTER_SEED + task$r + if (task$dep) 500000L else 0L)
  X <- gen_data(N, arm$thr, arm$sd_theta,
                pair = if (task$dep) DEP else NULL, rho = RHO)
  quiet <- function(e) suppressWarnings(suppressMessages(e))
  cut <- try(quiet(RMlocdepGammaCutoff(X, iterations = ITER, parallel = FALSE,
                                       hdci_width = HDCI)), silent = TRUE)
  if (inherits(cut, "try-error")) return(NULL)
  get <- function(pv) try(quiet(RMlocdepGamma(
    X, cutoff = cut, p_value = pv, correction = "fwer",
    output = "dataframe")$direction1), silent = TRUE)
  band <- get(FALSE); wy <- get(TRUE)
  if (inherits(band, "try-error") || inherits(wy, "try-error")) return(NULL)

  is_dep <- function(d) (d$Item1 %in% dep_names) & (d$Item2 %in% dep_names)
  shares <- function(d) xor(d$Item1 %in% dep_names, d$Item2 %in% dep_names)
  null_of <- function(d) if (task$dep) !is_dep(d) & !shares(d) else rep(TRUE, nrow(d))
  # The interval is two-sided, so the direction of a flag is read off the
  # tested statistic against the bound it crossed.
  above <- function(d) d$flagged & d$gamma_pair > d$gamma_high
  below <- function(d) d$flagged & d$gamma_pair < d$gamma_low

  c(dep        = task$dep,
    npairs     = nrow(band),
    band_det   = if (task$dep) any(above(band)[is_dep(band)]) else NA,
    wy_det     = if (task$dep) any(wy$flagged[is_dep(wy)]) else NA,
    band_fwer  = any(band$flagged[null_of(band)]),
    band_up    = any(above(band)[null_of(band)]),
    band_down  = any(below(band)[null_of(band)]),
    wy_fwer    = any(wy$flagged[null_of(wy)]),
    band_nflag = sum(band$flagged[null_of(band)]),
    wy_nflag   = sum(wy$flagged[null_of(wy)]),
    band_share = if (task$dep) any(band$flagged[shares(band)]) else NA,
    wy_share   = if (task$dep) any(wy$flagged[shares(wy)]) else NA)
}

tasks <- c(lapply(seq_len(REPS), function(r) list(r = r, dep = FALSE)),
           lapply(seq_len(REPS), function(r) list(r = r, dep = TRUE)))
t0 <- Sys.time()
res <- mclapply(tasks, one, mc.cores = N_CORES)
M <- as.data.frame(do.call(rbind, Filter(Negate(is.null), res)))
elapsed <- difftime(Sys.time(), t0, units = "min")
saveRDS(list(results = M, settings = list(
  n = N, reps = REPS, iterations = ITER, hdci_width = HDCI, rho = RHO,
  dep_pair = dep_names, seed = MASTER_SEED, elapsed_min = as.numeric(elapsed)
)), OUT_FILE)

line <- function(label, v) {
  p <- mean(v, na.rm = TRUE); k <- sum(!is.na(v))
  cat(sprintf("%-46s %.3f  (SE %.3f, %d reps)\n", label, p, rate_se(p, k), k))
}
A <- M[M$dep == 0, ]; B <- M[M$dep == 1, ]
cat(sprintf("\n==== RMlocdepGamma, n = %d, %d pairs, B = %d, %.1f min ====\n",
            N, M$npairs[1], ITER, as.numeric(elapsed)))
cat("\n-- complete null --\n")
line("interval, at least one pair flagged", A$band_fwer)
line("interval, above the upper bound", A$band_up)
line("interval, below the lower bound", A$band_down)
line("Westfall-Young, at least one flagged", A$wy_fwer)
cat(sprintf("%-46s %.2f vs %.2f\n", "mean pairs flagged, interval vs WY",
            mean(A$band_nflag), mean(A$wy_nflag)))
cat(sprintf("%-46s %.3f\n", "1 - hdci_width^pairs, predicted",
            1 - HDCI^M$npairs[1]))
cat("\n-- one weak dependent pair --\n")
line("interval, detects the dependent pair", B$band_det)
line("Westfall-Young, detects it", B$wy_det)
line("interval, >=1 false flag, disjoint pairs", B$band_fwer)
line("Westfall-Young, >=1 false flag, disjoint", B$wy_fwer)
line("interval, >=1 flag, pairs sharing an item", B$band_share)
line("Westfall-Young, >=1 flag, sharing an item", B$wy_share)
