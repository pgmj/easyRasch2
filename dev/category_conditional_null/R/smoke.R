# Smoke test (container: 4 cores, small B). Not a study run.
# Rscript R/smoke.R
source("R/functions.R")
library(parallel)
cores <- 4

n <- 500
others <- list(c(-2, -1, 0, 1), c(-1, -.5, .5, 1), c(-1.5, -.5, .5, 1.5),
               c(-1, 0, 1, 2), c(-2, -.5, 0, 1))
dl_disord <- c(list(c(-1.5, 0.5, -0.5, 1.5)), others)  # A: PCM, disordered
dl_order  <- c(list(c(-1.5, -0.5, 0.5, 1.5)), others)  # B, C base

run_one <- function(cond, seed, B) {
  set.seed(seed)
  th <- rnorm(n)
  X <- switch(cond,
    A = sim_pcm(th, dl_disord),
    B = make_exchangeable(th, dl_order, 1, pair = c(1, 2)),
    C = make_dontknow(th, dl_order, 1, pi_dk = 0.25))
  category_test(X, 1, B = B)
}

# One dataset per condition, B = 200
one <- mclapply(c("A", "B", "C"), run_one, seed = 11, B = 200, mc.cores = 3)
names(one) <- c("A", "B", "C")
for (k in names(one)) {
  cat("\n== Condition", k, "  item-1 thresholds:",
      paste(round(one[[k]]$thresholds, 2), collapse = " "), "\n")
  s <- one[[k]]$stats
  s[, -1] <- round(s[, -1], 2)
  print(s, row.names = FALSE)
  print(round(one[[k]]$p, 3))
}

# Null calibration under A: 20 datasets, B = 100
nullp <- mclapply(1:20, function(s) run_one("A", 1000 + s, 100)$p, mc.cores = cores)
nullp <- do.call(rbind, nullp)
cat("\nCondition A, share of p < .05 over 20 datasets (expect ~.05, MC sd ~.05):\n")
print(round(colMeans(nullp < 0.05), 2))

dir.create("results", showWarnings = FALSE)
saveRDS(list(one = one, nullp = nullp), "results/smoke.rds")
saveRDS(list(one = one, nullp = nullp), paste0("results/smoke_", Sys.Date(), ".rds"))
