# Large-B check of the studentised max-T calibration test under the null
# (2026-11-19). The rerun (B = 199) gave a family-wise rate of .074 in
# A0_pcm_ordered n1000 k5 and .048 or less in every other null cell. This run
# repeats the two k = 5 null cells at n = 1000 (A0 ordered and A1 disordered)
# with B = 999 and 1000 datasets per cell, new datasets. With 1000 datasets
# the Monte Carlo SE of a .05 rate is .007. n = 250 and 500 were left out
# (all at .048 or less at B = 199), decided with the user 2026-11-19.
#
# Run on the host:
#   cd R_packages/easyRasch2/dev/category_conditional_null
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 CORES=10 nohup Rscript R/run_calibration_B999.R > results/run_calibration_B999.log 2>&1 &
# Then: Rscript R/summarise_calibration_B999.R
#
# Expected time on 10 cores: about 3.1 h (the B = 199 n1000 k5 cell took
# 9.4 min at 500 datasets, x5 for B, x2 for datasets, two conditions). Cells
# are resumable (results/cal_b999_cells/<cell>.rds, skipped if it
# exists), so the run can be stopped and restarted.
# Environment: CORES (default 10, smoke 4), REPS (default 1000), B (default
# 999), NS (sample sizes, comma-separated, default 1000, e.g. NS=250,500,1000
# adds the smaller cells back, about 7 h in all), SMOKE=1 for a container
# smoke test (n = 250 only, 8 reps, B = 19).
# Stored per dataset: item statistics (B_used + 1) x 5, item p (min-p and
# max-T), max-T step-down p. About 40 MB per cell.

source("R/study_functions.R")
source("R/calibration_functions.R")
library(parallel)

smoke <- identical(Sys.getenv("SMOKE"), "1")
env_int <- function(name, default) {
  v <- Sys.getenv(name)
  out <- if (nzchar(v)) suppressWarnings(as.integer(v)) else as.integer(default)
  if (is.na(out) || out < 1L) stop(name, " must be a positive integer, got '", v, "'.")
  out
}
cores <- env_int("CORES", if (smoke) 4L else 10L)
reps  <- env_int("REPS", if (smoke) 8L else 1000L)
B     <- env_int("B", if (smoke) 19L else 999L)
message("Running on ", cores, " cores, ", reps, " datasets per cell, B = ", B, ".")

ns_env <- Sys.getenv("NS")
ns <- if (smoke) 250L else if (nzchar(ns_env)) suppressWarnings(as.integer(strsplit(ns_env, ",")[[1]])) else 1000L
if (anyNA(ns) || any(!ns %in% c(250L, 500L, 1000L)))
  stop("NS must be a comma-separated subset of 250, 500, 1000, got '", ns_env, "'.")
message("Sample sizes: ", paste(ns, collapse = ", "), ".")
grid <- expand.grid(cond = c("A0_pcm_ordered", "A1_pcm_disordered"), n = ns, k = 5,
                    stringsAsFactors = FALSE)
grid$cell <- sprintf("%s_n%d_k%d", grid$cond, grid$n, grid$k)
grid <- grid[order(grid$n), ]

outdir <- if (smoke) "results/cal_b999_smoke_cells" else "results/cal_b999_cells"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

RNGkind("L'Ecuyer-CMRG")
for (g in seq_len(nrow(grid))) {
  row <- grid[g, ]
  f <- file.path(outdir, paste0(row$cell, ".rds"))
  if (file.exists(f)) next
  t0 <- Sys.time()
  # Seed base differs from run_calibration_study.R, so the datasets are new.
  set.seed(match(row$cell, sort(grid$cell)) * 1000L + 20261119L %% 1000L + 7L)
  res <- mclapply(seq_len(reps), function(rp) {
    X <- gen_data(row$cond, row$n, row$k)
    out <- tryCatch(analyse_calibration_maxT(X, B),
                    error = function(e) list(skipped = conditionMessage(e)))
    out$rep <- rp
    out$data_id <- sum(X * seq_along(X)) %% 1e9
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
