# Calibration test (PAV miscalibration of the cumulative events, CML forecasts,
# exact conditional bootstrap with refit). Same 84 cells, generators,
# replications and B as R/run_study.R, new datasets. Run on the host:
#   cd R_packages/easyRasch2/dev/category_conditional_null
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript R/run_calibration_study.R > results/run_calibration_study.log 2>&1 &
# Resumable: each cell is saved to results/cal_cells/<cell>.rds and skipped
# if it exists. Expected time: about the same as run_study.R (2 h 16 min on
# 10 cores), since the refits dominate.
#
# Second run (2026-11-18 evening): item p from the pooled min-p (see
# .minp_pooled() in R/calibration_functions.R), and the bootstrap MCB array is
# kept per replication (boot_mcb). Results of the first run, whose item p was
# too small in the lower tail, are in results/old_minp_bug/. Same seeds, so
# the datasets are the same as in the first run (compare data_id).
#
# Cores: CORES=8 Rscript R/run_calibration_study.R (default 10, smoke 4).
# A different core count gives different datasets than the first run.
#
# Smoke mode (container): SMOKE=1 Rscript R/run_calibration_study.R, then
# SMOKE=1 Rscript R/summarise_calibration.R.

source("R/study_functions.R")
source("R/calibration_functions.R")
library(parallel)

smoke <- identical(Sys.getenv("SMOKE"), "1")
# CORES=<n> overrides the default (4 in smoke mode, 10 otherwise).
cores_env <- Sys.getenv("CORES")
cores <- if (nzchar(cores_env)) as.integer(cores_env) else if (smoke) 4L else 10L
if (is.na(cores) || cores < 1L) stop("CORES must be a positive integer, got '", cores_env, "'.")
message("Running on ", cores, " cores.")

ns <- if (smoke) 250 else c(250, 500, 1000)
ks <- if (smoke) 5 else c(5, 10)
# At least two replicates per worker in smoke mode (see run_study.R).
reps_null <- if (smoke) 12 else 500
reps_alt  <- if (smoke) 8 else 200
B_null    <- if (smoke) 19 else 199
B_alt     <- if (smoke) 19 else 99

conds <- conditions()
grid <- merge(conds, expand.grid(n = ns, k = ks))
grid$cell <- sprintf("%s_n%d_k%d", grid$cond, grid$n, grid$k)
grid <- grid[order(!grid$null, grid$n, grid$k), ]

outdir <- if (smoke) "results/cal_smoke_cells" else "results/cal_cells"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

RNGkind("L'Ecuyer-CMRG")
for (g in seq_len(nrow(grid))) {
  row <- grid[g, ]
  f <- file.path(outdir, paste0(row$cell, ".rds"))
  if (file.exists(f)) next
  reps <- if (row$null) reps_null else reps_alt
  B <- if (row$null) B_null else B_alt
  t0 <- Sys.time()
  set.seed(as.integer(g) * 1000L + 20261118L %% 1000L)
  res <- mclapply(seq_len(reps), function(rp) {
    X <- gen_data(row$cond, row$n, row$k)
    out <- tryCatch(analyse_calibration(X, B),
                    error = function(e) list(skipped = conditionMessage(e)))
    out$rep <- rp
    out$data_id <- sum(X * seq_along(X)) %% 1e9   # repeated-dataset guard
    out
  }, mc.cores = cores, mc.set.seed = TRUE)
  obj <- list(cell = row, reps = reps, B = B, results = res,
              date = Sys.Date(), session = utils::sessionInfo())
  saveRDS(obj, f)
  saveRDS(obj, sub("\\.rds$", paste0("_", Sys.Date(), ".rds"), f))
  message(sprintf("[%s] %s done: %d reps in %.1f min", format(Sys.time(), "%H:%M"),
                  row$cell, reps, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
message("all cells done")
