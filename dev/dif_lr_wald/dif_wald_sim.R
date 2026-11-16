# Calibration and power of anchored Wald DIF tests (two groups, item level).
#
# Questions:
#   1. Null: is the per-item Wald test calibrated (rejection near .05 with no
#      adjustment, z with SD near 1), and do Holm and single-step keep the
#      family-wise error rate at .05?
#   2. DIF: how often is a real DIF item found (power), and how often are
#      DIF-free items flagged because of the DIF in others (artificial DIF)?
#      The centred anchor (all items, as RMdifLR() uses today) is expected to
#      flag DIF-free items when DIF is unbalanced. The constant (Gini) and
#      forward (MTT) anchors from psychotools are meant to avoid that.
#
# The "all-other" anchor is left out: its z and p equal the centred ones,
# since item i against the mean of all items is (k - 1) / k times item i
# against the mean of the others.
#
# Design (48 cells):
#   n per group       150, 300, 600
#   structure         10 items x 4 categories, 20 dichotomous items
#   impact            group 2 latent mean 0 or 0.5 (SD 1 in both groups)
#   DIF scenario      null, one item +0.6, three items +0.6 (unbalanced),
#                     two items +0.6 / -0.6 (balanced)
# DIF is uniform: 0.6 logits added to every threshold of the item in group 2.
#
# Usage (from the easyRasch2 package root, on the host):
#   nohup Rscript dev/dif_lr_wald/dif_wald_sim.R 1000 10 > dev/dif_lr_wald/run.log 2>&1 &
# Arguments: R (replications per cell, default 1000), cores (default 10),
#   out_dir (default dev/dif_lr_wald).
# Needs psychotools and multcomp. Runs with parallel::mclapply (macOS/Linux).
# Each replication has its own seed, so results do not depend on `cores`.

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 10L
out_dir <- if (length(args) >= 3L) args[3L] else "dev/dif_lr_wald"

suppressMessages({
  library(psychotools)
  library(multcomp)
})
source(file.path("dev", "dif_lr_wald", "dif_wald.R"))

anchors <- c("centred", "constant", "forward")
dif_size <- 0.6

# --- Data generation ---------------------------------------------------------
# Partial credit model responses. `tau` is a list of threshold vectors.
sim_pcm <- function(theta, tau) {
  sapply(tau, function(t) {
    m <- length(t)
    eta <- outer(theta, 0:m) - matrix(c(0, cumsum(t)), length(theta), m + 1L,
                                      byrow = TRUE)
    eta <- eta - apply(eta, 1L, max)
    p <- exp(eta)
    p <- p / rowSums(p)
    u <- stats::runif(length(theta))
    rowSums(u > t(apply(p, 1L, cumsum))[, -(m + 1L), drop = FALSE])
  })
}

structures <- list(
  "10x4" = list(
    tau = lapply(seq(-1.5, 1.5, length.out = 10), function(l) l + c(-1, 0, 1)),
    dif_items = c(3L, 6L, 8L)
  ),
  "20x2" = list(
    tau = as.list(seq(-1.5, 1.5, length.out = 20)),
    dif_items = c(5L, 11L, 16L)
  )
)

# Shift per DIF item for each scenario (positions index `dif_items`).
scenario_shift <- function(scenario) {
  switch(scenario,
    null = numeric(0),
    one = dif_size,
    unbal3 = rep(dif_size, 3L),
    bal2 = c(dif_size, -dif_size)
  )
}

design <- expand.grid(
  n = c(150L, 300L, 600L),
  structure = names(structures),
  impact = c(0, 0.5),
  scenario = c("null", "one", "unbal3", "bal2"),
  stringsAsFactors = FALSE
)
design$cell <- seq_len(nrow(design))

# --- One replication -----------------------------------------------------------
one_rep <- function(cell, rep) {
  d <- design[cell, ]
  st <- structures[[d$structure]]
  set.seed(cell * 100000L + rep)
  shift <- scenario_shift(d$scenario)
  dif_items <- st$dif_items[seq_along(shift)]
  tau2 <- st$tau
  for (j in seq_along(dif_items)) tau2[[dif_items[j]]] <- tau2[[dif_items[j]]] + shift[j]
  x1 <- sim_pcm(stats::rnorm(d$n, 0, 1), st$tau)
  x2 <- sim_pcm(stats::rnorm(d$n, d$impact, 1), tau2)
  k <- length(st$tau)
  colnames(x1) <- colnames(x2) <- paste0("I", seq_len(k))
  # Every category must be used in both groups, as RMdifLR() requires.
  ncat <- vapply(st$tau, length, integer(1)) + 1L
  ok <- all(vapply(seq_len(k), function(i) {
    length(unique(x1[, i])) == ncat[i] && length(unique(x2[, i])) == ncat[i]
  }, logical(1)))
  if (!ok) {
    return(data.frame(cell = cell, rep = rep, status = "empty_category"))
  }
  res <- tryCatch(
    dif_wald(dif_wald_fit(rbind(x1, x2), rep(1:2, each = d$n)), anchors = anchors),
    error = function(e) conditionMessage(e)
  )
  if (is.character(res)) {
    return(data.frame(cell = cell, rep = rep, status = paste("error:", res)))
  }
  res$cell <- cell
  res$rep <- rep
  res$status <- "ok"
  res$item_no <- match(res$Item, paste0("I", seq_len(k)))
  res$true_shift <- 0
  res$true_shift[res$item_no %in% dif_items] <-
    shift[match(res$item_no[res$item_no %in% dif_items], dif_items)]
  res$Item <- NULL
  res
}

# --- Run ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Design: %d cells x %d replications on %d cores\n", nrow(design), R, cores))
t_start <- Sys.time()
all_res <- vector("list", nrow(design))
for (cell in design$cell) {
  t0 <- Sys.time()
  reps <- parallel::mclapply(seq_len(R), function(r) one_rep(cell, r),
                             mc.cores = cores, mc.preschedule = TRUE)
  bad <- vapply(reps, function(x) inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("Worker failure in cell ", cell, ": ", reps[[which(bad)[1L]]])
  ok_rows <- Filter(function(x) "Anchor" %in% names(x), reps)
  st_rows <- do.call(rbind, lapply(reps, function(x) x[1L, c("cell", "rep", "status")]))
  all_res[[cell]] <- list(results = do.call(rbind, ok_rows), status = st_rows)
  d <- design[cell, ]
  cat(sprintf(
    "[%s] cell %2d/%d  n=%d %s impact=%.1f %-6s  ok %d/%d  %.1f min\n",
    format(Sys.time(), "%H:%M"), cell, nrow(design), d$n, d$structure,
    d$impact, d$scenario, sum(st_rows$status == "ok"), R,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}
results <- do.call(rbind, lapply(all_res, `[[`, "results"))
status <- do.call(rbind, lapply(all_res, `[[`, "status"))
status <- merge(status, design, by = "cell")

# --- Summary -----------------------------------------------------------------
# Flags at alpha = .05. An untested anchor item (NA) counts as not flagged.
alpha <- 0.05
flag <- function(p) !is.na(p) & p < alpha
results$f_none <- flag(results$p)
results$f_holm <- flag(results$padj_holm)
results$f_ss <- flag(results$padj_ss)
results$dif_item <- results$true_shift != 0

# One row per replication and anchor (aggregate is fast on ~150,000 groups,
# where rbind over split data frames is not).
results$ff_none <- results$f_none & !results$dif_item
results$ff_holm <- results$f_holm & !results$dif_item
results$ff_ss <- results$f_ss & !results$dif_item
results$found_ss <- results$f_ss & results$dif_item
results$found_holm <- results$f_holm & results$dif_item
results$dif_tested <- results$dif_item & !is.na(results$p)
per_rep <- stats::aggregate(
  cbind(ff_none, ff_holm, ff_ss, found_ss, found_holm, dif_item, dif_tested) ~
    cell + rep + Anchor,
  data = results, FUN = sum, na.action = stats::na.pass
)
names(per_rep)[names(per_rep) == "dif_item"] <- "n_dif"
names(per_rep)[names(per_rep) == "ff_ss"] <- "n_false_ss"
per_rep$any_false_none <- per_rep$ff_none > 0
per_rep$any_false_holm <- per_rep$ff_holm > 0
per_rep$any_false_ss <- per_rep$n_false_ss > 0
results[c("ff_none", "ff_holm", "ff_ss", "found_ss", "found_holm", "dif_tested")] <- NULL
summ <- do.call(rbind, lapply(
  split(per_rep, list(per_rep$cell, per_rep$Anchor), drop = TRUE),
  function(x) {
    z <- results$z[results$cell == x$cell[1L] & results$Anchor == x$Anchor[1L] &
                     !results$dif_item & !is.na(results$z)]
    data.frame(
      cell = x$cell[1L], Anchor = x$Anchor[1L], reps = nrow(x),
      fwer_none = mean(x$any_false_none),
      fwer_holm = mean(x$any_false_holm),
      fwer_ss = mean(x$any_false_ss),
      mean_false_ss = mean(x$n_false_ss),
      power_ss = if (x$n_dif[1L] > 0) mean(x$found_ss / x$n_dif) else NA_real_,
      power_holm = if (x$n_dif[1L] > 0) mean(x$found_holm / x$n_dif) else NA_real_,
      dif_untested = if (x$n_dif[1L] > 0) mean(x$dif_tested < x$n_dif) else NA_real_,
      sd_z_free = stats::sd(z)
    )
  }
))
summ <- merge(design, summ, by = "cell")
summ <- summ[order(summ$scenario, summ$structure, summ$impact, summ$n, summ$Anchor), ]
fails <- as.data.frame(table(cell = status$cell, status = status$status))
fails <- fails[fails$Freq > 0 & fails$status != "ok", ]

out <- list(
  design = design, results = results, status = status, per_rep = per_rep,
  summary = summ, failures = fails, R = R, alpha = alpha, dif_size = dif_size,
  anchors = anchors,
  session = utils::sessionInfo(),
  run_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))
)
stamp <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(out, file.path(out_dir, paste0("dif_wald_sim_", stamp, ".rds")))
saveRDS(out, file.path(out_dir, "dif_wald_sim.rds"))
utils::write.csv(summ, file.path(out_dir, "dif_wald_sim_summary.csv"), row.names = FALSE)

options(width = 200)
cat("\nFailures (cell, status, count):\n")
print(fails, row.names = FALSE)
cat("\nSummary (fwer_* = any DIF-free item flagged; power_* = share of DIF items flagged):\n")
print(summ[, c("n", "structure", "impact", "scenario", "Anchor", "reps",
               "fwer_none", "fwer_holm", "fwer_ss", "mean_false_ss",
               "power_holm", "power_ss", "dif_untested", "sd_z_free")],
      digits = 3, row.names = FALSE)
cat(sprintf("\nTotal run time: %.1f min\n", out$run_minutes))
