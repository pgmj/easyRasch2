# Full grid of the band study for the adjacent and residual plot types (the
# calibration type has its full grid in run_B.R). Same 54 cells, seeds and
# datasets as run_B.R. Cells on disk are skipped, so a restart resumes.
# Container (user request 2026-11-26, overnight, 4 cores):
#   NCORES=4 LANG=en_US.UTF-8 nohup Rscript dev/catfit_dif/R/run_types_full.R \
#     > dev/catfit_dif/results/run_types_full.log 2>&1 &
# Then: LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/summarise_types_full.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/catfit_dif/R/perm_functions.R")

cores <- as.integer(Sys.getenv("NCORES", "4"))
B <- 999L
R <- 1000L
alphas <- c(0.05, 0.10)
groups <- list(two_equal = c(0.5, 0.5), two_80_20 = c(0.8, 0.2),
               three_equal = rep(1 / 3, 3))
grid <- expand.grid(impact = c(0, 0.5, 1), groups = names(groups),
                    n = c(250L, 500L, 1000L), k = c(5L, 10L),
                    stringsAsFactors = FALSE)
grid$cell <- sprintf("c%02d_imp%.1f_%s_n%d_k%d", seq_len(nrow(grid)),
                     grid$impact, grid$groups, grid$n, grid$k)
grid$id <- seq_len(nrow(grid))
out_dir <- "dev/catfit_dif/results/permband_types_cells"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

one_rep <- function(row, rep) {
  set.seed(20261125L + 10000L * row$id + rep)
  d <- gen_null(row$n, row$k, row$impact, groups[[row$groups]])
  hash <- digest::digest(d$X)
  pr <- prep_dif(d$X, d$grp)
  if (is.null(pr)) return(list(failed = TRUE, hash = hash))
  setups <- list(adjacent = adjacent_setup(pr),
                 residuals = residual_setup(pr, 5L))
  res <- lapply(names(setups), function(nm) {
    st <- setups[[nm]]
    b <- perm_bands_generic(pr, st$lay, st$fun, B = B)
    pv <- b$points$p_point[b$points$varies]
    data.frame(type = nm, alpha = alphas,
               points = length(pv),
               points_out = vapply(alphas, function(a) sum(pv <= a), 0),
               curves = nrow(b$curves),
               any_point = vapply(alphas, function(a) {
                 sum(b$curves$min_p_point <= a)
               }, 0),
               simult = vapply(alphas, function(a) {
                 sum(b$curves$p_sim <= a)
               }, 0),
               points_per_curve = mean(b$curves$n_points_var))
  })
  list(failed = FALSE, hash = hash, res = do.call(rbind, res))
}

t0 <- Sys.time()
for (j in seq_len(nrow(grid))) {
  row <- grid[j, ]
  f <- file.path(out_dir, paste0(row$cell, ".rds"))
  if (file.exists(f)) { cat("skip", row$cell, "\n"); next }
  tc <- Sys.time()
  reps <- parallel::mclapply(seq_len(R), function(rep) {
    tryCatch(one_rep(row, rep), error = function(e) {
      list(failed = TRUE, error = conditionMessage(e))
    })
  }, mc.cores = cores)
  saveRDS(list(cell = row, B = B, R = R, reps = reps), f)
  cat(sprintf("%s done in %.1f min (failed: %d) at %s\n", row$cell,
              as.numeric(difftime(Sys.time(), tc, units = "mins")),
              sum(vapply(reps, `[[`, logical(1), "failed")),
              format(Sys.time(), "%H:%M")))
}
cat(sprintf("All cells done in %.1f h\n",
            as.numeric(difftime(Sys.time(), t0, units = "hours"))))
