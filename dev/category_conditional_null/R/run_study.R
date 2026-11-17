# Main study. Run on the host, from this folder:
#   cd R_packages/easyRasch2/dev/category_conditional_null
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 nohup Rscript R/run_study.R > results/run_study.log 2>&1 &
# Resumable: each cell is saved to results/cells/<cell>.rds and skipped if it
# exists. Delete a cell file to rerun it.
#
# Smoke mode (container): SMOKE=1 Rscript R/run_study.R, then
# SMOKE=1 Rscript R/summarise.R to check for repeated datasets.

source("R/study_functions.R")
library(parallel)

smoke <- identical(Sys.getenv("SMOKE"), "1")
cores <- if (smoke) 4L else 10L

ns <- if (smoke) 250 else c(250, 500, 1000)
ks <- if (smoke) 5 else c(5, 10)
# Smoke mode runs at least two replicates per worker (mclapply splits reps
# evenly over cores). A bug that only shows from a worker's second replicate
# (the 2026-11-16 RNG bug) is invisible with one replicate per worker.
reps_null <- if (smoke) 12 else 500
reps_alt  <- if (smoke) 8 else 200
B_null    <- if (smoke) 19 else 199
B_alt     <- if (smoke) 19 else 99

conds <- conditions()
grid <- merge(conds, expand.grid(n = ns, k = ks))
grid$cell <- sprintf("%s_n%d_k%d", grid$cond, grid$n, grid$k)
# Null cells first: they decide whether the rest is interpretable.
grid <- grid[order(!grid$null, grid$n, grid$k), ]

outdir <- if (smoke) "results/smoke_cells" else "results/cells"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

RNGkind("L'Ecuyer-CMRG")
for (g in seq_len(nrow(grid))) {
  row <- grid[g, ]
  f <- file.path(outdir, paste0(row$cell, ".rds"))
  if (file.exists(f)) next
  reps <- if (row$null) reps_null else reps_alt
  B <- if (row$null) B_null else B_alt
  t0 <- Sys.time()
  set.seed(as.integer(g) * 1000L + 7L)
  res <- mclapply(seq_len(reps), function(rp) {
    X <- gen_data(row$cond, row$n, row$k)
    out <- tryCatch(analyse_dataset(X, B), error = function(e) list(skipped = conditionMessage(e)))
    out$rep <- rp
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
