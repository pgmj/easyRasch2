# Study B: permutation bands under the Rasch model without DIF.
# Pointwise vs simultaneous exceedance. See ../PLAN.md.
#
# Run from the package root.
#   Smoke (container, 2 cells x 200 datasets, 4 cores):
#     SMOKE=1 NCORES=4 LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/run_B.R
#   Full (host, 54 cells x 1000 datasets, about 1.5 h on 10 cores):
#     NCORES=10 LANG=en_US.UTF-8 nohup Rscript dev/catfit_dif/R/run_B.R \
#       > dev/catfit_dif/results/run_B.log 2>&1 &
# Cells already on disk are skipped, so an interrupted run can be restarted.
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/catfit_dif/R/perm_functions.R")

smoke <- identical(Sys.getenv("SMOKE"), "1")
cores <- as.integer(Sys.getenv("NCORES", "4"))
B <- 999L
R <- if (smoke) 200L else 1000L
alphas <- c(0.05, 0.10)
n_bins <- c(0, 4, 9, 19, 49, Inf)  # n per point: 1-4, 5-9, 10-19, 20-49, 50+

groups <- list(two_equal = c(0.5, 0.5), two_80_20 = c(0.8, 0.2),
               three_equal = rep(1 / 3, 3))
grid <- expand.grid(impact = c(0, 0.5, 1), groups = names(groups),
                    n = c(250L, 500L, 1000L), k = c(5L, 10L),
                    stringsAsFactors = FALSE)
grid$cell <- sprintf("c%02d_imp%.1f_%s_n%d_k%d", seq_len(nrow(grid)),
                     grid$impact, grid$groups, grid$n, grid$k)
grid$id <- seq_len(nrow(grid))
if (smoke) {
  grid <- grid[(grid$impact == 1 & grid$groups == "two_80_20" &
                  grid$n == 250L & grid$k == 5L) |
                 (grid$impact == 0.5 & grid$groups == "three_equal" &
                    grid$n == 500L & grid$k == 10L), ]
}
out_dir <- if (smoke) "dev/catfit_dif/results/permband_B_smoke_cells" else
  "dev/catfit_dif/results/permband_B_cells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

one_rep <- function(row, rep) {
  # Seed set per replication, never inside a generator.
  set.seed(20261125L + 10000L * row$id + rep)
  d <- gen_null(row$n, row$k, row$impact, groups[[row$groups]])
  hash <- digest::digest(d$X)
  pr <- prep_dif(d$X, d$grp)
  if (is.null(pr)) return(list(failed = TRUE, hash = hash))
  b <- perm_bands(pr, B = B)
  pt <- b$points
  # Certain events (P = 0 or 1) cannot vary: left out of all rates.
  pt <- pt[pt$predicted > 1e-12 & pt$predicted < 1 - 1e-12, ]
  pt$nbin <- cut(pt$n, n_bins, labels = c("1-4", "5-9", "10-19", "20-49",
                                           "50+"))
  pts <- do.call(rbind, lapply(alphas, function(a) {
    s <- stats::aggregate(cbind(points = 1, outside = pt$p_point <= a) ~
                            k + group + nbin, data = pt, FUN = sum)
    s$alpha <- a
    s
  }))
  cv <- b$curves
  cv$n_points_var <- vapply(cv$curve, function(c) {
    sum(pt$curve == c)
  }, numeric(1))
  list(failed = FALSE, hash = hash, points = pts,
       curves = cv[, c("item", "k", "group", "p_sim", "min_p_point",
                       "n_points", "n_points_var")])
}

t0 <- Sys.time()
for (j in seq_len(nrow(grid))) {
  row <- grid[j, ]
  f <- file.path(out_dir, paste0(row$cell, ".rds"))
  if (file.exists(f)) { cat("skip", row$cell, "\n"); next }
  tc <- Sys.time()
  res <- parallel::mclapply(seq_len(R), function(rep) {
    tryCatch(one_rep(row, rep),
             error = function(e) list(failed = TRUE, error = conditionMessage(e)))
  }, mc.cores = cores)
  saveRDS(list(cell = row, B = B, R = R, reps = res), f)
  cat(sprintf("%s done in %.1f min (failed: %d)\n", row$cell,
              as.numeric(difftime(Sys.time(), tc, units = "mins")),
              sum(vapply(res, `[[`, logical(1), "failed"))))
}
cat(sprintf("All cells done in %.1f min\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins"))))
