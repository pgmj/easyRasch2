# Null check of the standardised residuals in RMitemCatFitPlot(type =
# "residuals"). Question: under the PCM, how often is |z| > 1.96, and is z
# close to N(0, 1)? The plot's z uses the CML estimates and the binomial
# variance sum(n_r p (1 - p)), which treats the thresholds as known. The
# category counts are sufficient for the thresholds, so at the estimates
# each category's residuals sum to zero over score groups; we expect
# under-dispersion. z_true uses the generating thresholds and isolates that
# effect from discreteness.
#
# Run from easyRasch2/:  Rscript dev/category_conditional_null/R/residual_null_check.R [reps] [cores]
# Container: reps = 200, cores = 4. Host: reps = 1000, cores = 10.

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 200L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L

suppressMessages(devtools::load_all(".", quiet = TRUE))
out_dir   <- "dev/category_conditional_null/results"
cells_dir <- file.path(out_dir, "residual_null_cells")
dir.create(cells_dir, showWarnings = FALSE, recursive = TRUE)

GROUPS <- c(3L, 5L, 10L)
cells <- expand.grid(n = c(200L, 500L, 1000L), k = c(5L, 10L), m = c(2L, 4L))

# Thresholds: item locations spread over [-1, 1], steps evenly spaced over
# [-1.5, 1.5] with a fixed jitter. Fixed per (k, m).
true_thr <- function(k, m) {
  withr::with_seed(20261120 + 10 * k + m, .rng_kind = "Mersenne-Twister", {
    loc <- seq(-1, 1, length.out = k)
    lapply(loc, function(l) l + seq(-1.5, 1.5, length.out = m) +
             stats::rnorm(m, 0, 0.2))
  })
}

sim_pcm <- function(theta, thr) {
  vapply(thr, function(t) {
    lp <- outer(theta, 0:length(t)) -
      matrix(c(0, cumsum(t)), length(theta), length(t) + 1L, byrow = TRUE)
    p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
    rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
  }, numeric(length(theta)))
}

one_rep <- function(n, k, m, thr) {
  X <- sim_pcm(stats::rnorm(n), thr)
  colnames(X) <- paste0("I", seq_len(k))
  prep <- tryCatch(.catfit_prepare(as.data.frame(X)), error = function(e) NULL)
  # Skip a dataset with an empty category (fewer thresholds estimated than
  # generated), as the main study does.
  if (is.null(prep) || any(lengths(prep$thr) != m)) return(NULL)
  P_true <- .cond_cat_probs(thr)
  do.call(rbind, lapply(GROUPS, function(g) {
    # Equal-count groups as RMitemCatFitPlot() forms them.
    sg <- .score_group_labels(prep$r, .cicc_class_bins(prep$r, "quantile", g))
    est <- .catplot_residuals(prep$X, prep$r, prep$P, sg)
    tru <- .catplot_residuals(prep$X, prep$r, P_true, sg)
    stopifnot(identical(est$observed, tru$observed))
    data.frame(groups = g, item = est$item, score_group = as.integer(est$score_group),
               category = est$category, extreme = est$category %in% c(0L, m),
               expected = est$expected, z = est$z, z_true = tru$z)
  }))
}

RNGkind("L'Ecuyer-CMRG")
t0 <- Sys.time()
for (ci in seq_len(nrow(cells))) {
  cl <- cells[ci, ]
  f <- file.path(cells_dir, sprintf("n%d_k%d_m%d.rds", cl$n, cl$k, cl$m))
  if (file.exists(f)) next
  thr <- true_thr(cl$k, cl$m)
  set.seed(20261120 + ci)  # per cell, so one cell can be rerun alone
  res <- parallel::mclapply(seq_len(REPS), function(rp) {
    d <- one_rep(cl$n, cl$k, cl$m, thr)
    if (!is.null(d)) d$rep <- rp
    d
  }, mc.cores = CORES, mc.set.seed = TRUE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE)
  if (any(bad)) stop("errors in cell ", ci, ": ", res[bad][[1]])
  n_skip <- sum(vapply(res, is.null, TRUE))
  d <- do.call(rbind, res)
  d$n <- cl$n; d$k <- cl$k; d$m <- cl$m
  attr(d, "n_skipped") <- n_skip
  saveRDS(d, f)
  message(sprintf("cell %d/%d n=%d k=%d m=%d skipped=%d  %.1f min", ci, nrow(cells),
                  cl$n, cl$k, cl$m, n_skip, as.numeric(Sys.time() - t0, units = "mins")))
}

# ── Summary ─────────────────────────────────────────────────────
all <- do.call(rbind, lapply(list.files(cells_dir, full.names = TRUE), readRDS))
all <- all[!is.na(all$z) & !is.na(all$z_true), ]
all$e_bin <- cut(all$expected, c(0, 5, 20, Inf), labels = c("E < 5", "5-20", "> 20"),
                 right = FALSE)
summ <- function(d) data.frame(
  cells = nrow(d),
  beyond_est = mean(abs(d$z) > 1.96), beyond_true = mean(abs(d$z_true) > 1.96),
  sd_est = stats::sd(d$z), sd_true = stats::sd(d$z_true),
  mean_est = mean(d$z),
  upper_est = mean(d$z > 1.96), lower_est = mean(d$z < -1.96)
)
by_vars <- function(d, v) {
  sp <- split(d, d[v], drop = TRUE)
  out <- do.call(rbind, lapply(sp, summ))
  cbind(do.call(rbind, lapply(sp, function(x) x[1, v, drop = FALSE])), out)
}
summary_list <- list(
  reps = REPS,
  by_cell   = by_vars(all, c("groups", "n", "k", "m")),
  by_bin    = by_vars(all, c("groups", "k", "e_bin")),
  by_extreme = by_vars(all, c("groups", "k", "m", "extreme")),
  skipped = vapply(list.files(cells_dir, full.names = TRUE),
                   function(f) attr(readRDS(f), "n_skipped"), 0L)
)
saveRDS(summary_list, file.path(out_dir, "residual_null_summary.rds"))
saveRDS(summary_list, file.path(out_dir, paste0("residual_null_summary_", Sys.Date(), ".rds")))
message("done, ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
