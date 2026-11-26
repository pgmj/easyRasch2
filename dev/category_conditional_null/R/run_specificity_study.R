# Specificity study run (2026-11-23). Every condition of the main study plus
# the four nuisance conditions (specificity_functions.R), at n = 250 and 500
# and k = 5 and 10. Each dataset gets both tests:
#   RMitemCatFit()'s calibration test (main-study implementation,
#     analyse_calibration(), combined by max-T with the Westfall-Young
#     step-down as in RMitemCatFit());
#   conditional infit (iarm via easyRasch2::RMitemInfit()), whose p-values
#     come from an oracle reference set per cell (summarise_specificity.R).
# Replication as the main calibration run: PCM 500 datasets with B = 199,
# other conditions 200 with B = 99.
#
# Run from easyRasch2/dev/category_conditional_null/:
#   CORES=10 Rscript R/run_specificity_study.R            # full run
#   SMOKE=1 CORES=4 Rscript R/run_specificity_study.R     # smoke test
# CORES defaults to 4. Completed cells are skipped, so a stopped run resumes.
# Then: CORES=10 Rscript R/summarise_specificity.R

source("R/specificity_functions.R")
library(parallel)
stopifnot(requireNamespace("easyRasch2", quietly = TRUE),
          requireNamespace("iarm", quietly = TRUE))
message("easyRasch2 ", utils::packageVersion("easyRasch2"))

smoke <- identical(Sys.getenv("SMOKE"), "1")
cores <- as.integer(Sys.getenv("CORES", "4"))
if (is.na(cores) || cores < 1L) stop("CORES must be a positive integer.")
message("Running on ", cores, " cores", if (smoke) " (smoke test)", ".")

ns <- if (smoke) 250L else c(250L, 500L)
ks <- if (smoke) 5L else c(5L, 10L)
reps_null <- if (smoke) 8L else 500L
reps_alt  <- if (smoke) 4L else 200L
B_null    <- if (smoke) 19L else 199L
B_alt     <- if (smoke) 19L else 99L
R_REF     <- if (smoke) 50L else 1000L
conds <- spec_conditions()
if (smoke) conds <- conds[conds$cond %in% c("A0_pcm_ordered", "SWAP_1_2", NUISANCE$cond), ]

outdir <- if (smoke) "results/spec_smoke_cells" else "results/spec_cells"
refdir <- if (smoke) "results/spec_smoke_cells/ref" else "results/spec_infit_ref"
dir.create(refdir, recursive = TRUE, showWarnings = FALSE)

RNGkind("L'Ecuyer-CMRG")

# ── Infit reference sets, one per n, k and base thresholds ─────────
refs <- expand.grid(n = ns, k = ks, base = c("WIDE", "NARROW"),
                    stringsAsFactors = FALSE)
for (j in seq_len(nrow(refs))) {
  rf <- refs[j, ]
  f <- file.path(refdir, sprintf("ref_%s_n%d_k%d.rds", rf$base, rf$n, rf$k))
  if (file.exists(f)) next
  t0 <- Sys.time()
  set.seed(20261123L + 100L * j)
  ref <- infit_reference(rf$n, rf$k, rf$base, R_REF, cores)
  saveRDS(list(n = rf$n, k = rf$k, base = rf$base, ref = ref), f)
  message(sprintf("[%s] infit reference %s n%d k%d: %d datasets in %.1f min",
                  format(Sys.time(), "%H:%M"), rf$base, rf$n, rf$k, nrow(ref),
                  as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}

# ── Study cells ───────────────────────────────────────────────────
grid <- merge(conds, expand.grid(n = ns, k = ks))
grid$cell <- sprintf("%s_n%d_k%d", grid$cond, grid$n, grid$k)
grid <- grid[order(!grid$null, grid$n, grid$k, grid$cond), ]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

for (g in seq_len(nrow(grid))) {
  row <- grid[g, ]
  f <- file.path(outdir, paste0(row$cell, ".rds"))
  if (file.exists(f)) next
  reps <- if (row$null) reps_null else reps_alt
  B <- if (row$null) B_null else B_alt
  t0 <- Sys.time()
  # Seed from the cell name, so a cell's datasets do not depend on the grid.
  set.seed(20261123L + sum(utf8ToInt(row$cell) * seq_along(utf8ToInt(row$cell))))
  res <- mclapply(seq_len(reps), function(rp) {
    X <- gen_spec(row$cond, row$n, row$k)
    inf <- infit_vec(X)
    cal <- tryCatch(analyse_calibration(X, B),
                    error = function(e) list(skipped = conditionMessage(e)))
    out <- list(rep = rep_len(rp, 1L), infit = inf, skipped = cal$skipped,
                B_used = cal$B_used)
    if (is.na(cal$skipped) && !is.null(cal$boot_mcb) && cal$B_used > 0) {
      cp <- catfit_p(cal$mcb, cal$boot_mcb)
      out$catfit_p <- cp$p
      out$catfit_padj <- cp$padj
      out$p_raw <- cal$p_raw
    }
    out
  }, mc.cores = cores, mc.set.seed = TRUE)
  bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE)
  if (any(bad)) stop("errors in ", row$cell, ": ", res[bad][[1]])
  obj <- list(cell = row, reps = reps, B = B, base = base_of(row$cond),
              results = res, date = Sys.Date(), session = utils::sessionInfo())
  saveRDS(obj, f)
  message(sprintf("[%s] %s done: %d reps in %.1f min", format(Sys.time(), "%H:%M"),
                  row$cell, reps, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
message("all cells done")
