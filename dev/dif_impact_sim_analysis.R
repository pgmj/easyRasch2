# Analysis of the RMdifImpact() step 2 simulation (dev/dif_impact_sim.R).
#
# Reads dev/dif_impact_sim/sim_results.rds (about 40 MB, not tracked) and
# writes dev/dif_impact_sim/sim_results_slim.rds (tracked) with every table
# the design log and a later qmd need, so neither has to load the full file.
#
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript dev/dif_impact_sim_analysis.R
#
# Tables in the slim file:
#   cells        the 42 cells
#   targets      population targets per cell (cell_truth())
#   strategies   oracle keep / remove / split per estimator: impact, empirical
#                SE, naive CI coverage, and the split / remove SE ratios on
#                three scales (native units, SMD, coefficient of variation)
#   person       person-level WLE results per strategy and group
#   estimates    reported estimates vs targets (oracle and flagged sets,
#                DIF and no-DIF data), from summarise_estimates()
#   flag_class   flagged-set estimates split by what RMdifAnchor() flagged:
#                exact (all true items, nothing else), missed (a true item
#                missing, no false flag), fp (at least one false flag)
#   null_flagged estimates on the no-DIF twin when RMdifAnchor() flagged
#                something (all of them false flags)
#   detection    summarise_detection() per cell

sim_dir <- file.path("dev", "dif_impact_sim")
o <- readRDS(file.path(sim_dir, "sim_results.rds"))

cell_id <- function(c) {
  sprintf("k%d %s %s d%.1f %d/%d", c$k, c$scenario, c$position, c$d,
          c$n_ref, c$n_focal)
}
with_cell <- function(x, df) {
  c0 <- x$cell
  cbind(cell = cell_id(c0), k = c0$k, scenario = c0$scenario,
        position = c0$position, d = c0$d, n_ref = c0$n_ref,
        n_focal = c0$n_focal, df, stringsAsFactors = FALSE)
}
n_dif <- function(cell) if (cell$scenario == "one") 1L else 2L

stopifnot(all(vapply(o$results, function(x) length(x$errors), 0L) == 0L))

# ---------------------------------------------------------------------------
# Oracle strategies
# ---------------------------------------------------------------------------

strategies <- do.call(rbind, lapply(o$results, function(x) {
  g <- x$group
  s <- x$summary
  key <- paste(s$strategy, s$estimator)
  sd_smd <- tapply(g$smd_dif, paste(g$strategy, g$estimator), stats::sd)
  cv <- tapply(g$delta_dif, paste(g$strategy, g$estimator),
               function(v) stats::sd(v) / abs(mean(v)))
  s$emp_se_smd <- unname(sd_smd[key])
  s$cv <- unname(cv[key])
  with_cell(x, s)
}))

ratio <- function(col) {
  sp <- strategies[strategies$strategy == "split", ]
  rm <- strategies[strategies$strategy == "remove", ]
  stopifnot(identical(paste(sp$cell, sp$estimator),
                      paste(rm$cell, rm$estimator)))
  sp[[col]] / rm[[col]]
}
split_remove <- strategies[strategies$strategy == "split",
                           c("cell", "k", "scenario", "position", "d",
                             "n_ref", "n_focal", "estimator")]
split_remove$se_native <- ratio("emp_se")
split_remove$se_smd <- ratio("emp_se_smd")
split_remove$cv <- ratio("cv")

# ---------------------------------------------------------------------------
# Person level
# ---------------------------------------------------------------------------

person <- do.call(rbind, lapply(o$results, function(x) {
  a <- stats::aggregate(cbind(rmse_theta, bias_theta, mean_sem,
                              mean_shift_sem, class_disagree, class_error) ~
                          strategy + group, x$person, mean)
  with_cell(x, a)
}))

# ---------------------------------------------------------------------------
# Reported estimates
# ---------------------------------------------------------------------------

estimates <- do.call(rbind, lapply(o$results, function(x) {
  with_cell(x, x$estimate_summary)
}))
targets <- do.call(rbind, lapply(o$results, function(x) {
  with_cell(x, data.frame(quantity = names(x$truth$target),
                          target = unname(x$truth$target)))
}))
detection <- do.call(rbind, lapply(o$results, function(x) {
  with_cell(x, x$detection)
}))

flag_quant <- c("impact_latent", "impact_wle", "impact_sum", "b_maxabs",
                "b_sem_maxabs", "b_cut")

flag_class <- do.call(rbind, lapply(o$results, function(x) {
  nd <- n_dif(x$cell)
  tg <- x$truth$target
  e <- x$estimates[x$estimates$set == "flagged" & x$estimates$data == "dif" &
                     !is.na(x$estimates$quantity), ]
  e$class <- ifelse(e$fp > 0, "fp", ifelse(e$tp < nd, "missed", "exact"))
  reps <- length(unique(x$estimates$rep))
  rows <- list()
  for (q in flag_quant) for (cl in c("exact", "missed", "fp")) {
    z <- e[e$quantity == q & e$class == cl, ]
    ok <- is.finite(z$estimate) & is.finite(z$se)
    rows[[length(rows) + 1L]] <- data.frame(
      quantity = q, class = cl, reps = reps, n = nrow(z),
      share = nrow(z) / reps, target = tg[[q]],
      mean_est = if (nrow(z)) mean(z$estimate) else NA_real_,
      bias = if (nrow(z)) mean(z$estimate) - tg[[q]] else NA_real_,
      rel_bias = if (nrow(z)) (mean(z$estimate) - tg[[q]]) / abs(tg[[q]])
                 else NA_real_,
      coverage = if (any(ok)) mean(abs(z$estimate[ok] - tg[[q]]) <=
                                     1.96 * z$se[ok]) else NA_real_)
  }
  with_cell(x, do.call(rbind, rows))
}))

null_flagged <- do.call(rbind, lapply(o$results, function(x) {
  e <- x$estimates[x$estimates$set == "flagged" & x$estimates$data == "null" &
                     !is.na(x$estimates$quantity), ]
  es <- x$estimate_summary
  reps <- length(unique(x$estimates$rep))
  rows <- lapply(flag_quant, function(q) {
    z <- e[e$quantity == q, ]
    ok <- is.finite(z$estimate) & is.finite(z$se)
    data.frame(
      quantity = q, reps = reps, n = nrow(z), share = nrow(z) / reps,
      mean_abs = mean(abs(z$estimate)),
      q95_abs = unname(stats::quantile(abs(z$estimate), .95)),
      ci_excludes_0 = mean(abs(z$estimate[ok]) > 1.96 * z$se[ok]),
      oracle_null_q95 = es$q95_abs[es$data == "null" & es$set == "oracle" &
                                     es$quantity == q],
      target_dif = x$truth$target[[q]])
  })
  with_cell(x, do.call(rbind, rows))
}))

slim <- list(cells = o$cells, reps = o$reps, run_time = o$run_time,
             date = o$date, targets = targets, strategies = strategies,
             split_remove = split_remove, person = person,
             estimates = estimates, flag_class = flag_class,
             null_flagged = null_flagged, detection = detection)
saveRDS(slim, file.path(sim_dir, "sim_results_slim.rds"))
saveRDS(slim, file.path(sim_dir, sprintf("sim_results_slim_%s.rds",
                                         format(o$date, "%Y-%m-%d"))))

# ---------------------------------------------------------------------------
# Printed overview (the numbers quoted in dev/dif-impact-design.md)
# ---------------------------------------------------------------------------

rng <- function(v) sprintf("%.3f to %.3f", min(v), max(v))
cat("Split / remove SE ratio, by estimator\n")
for (est in c("latent", "wle", "sum")) {
  s <- split_remove[split_remove$estimator == est, ]
  cat(sprintf("  %-6s native %s | SMD %s | CV %s\n", est, rng(s$se_native),
              rng(s$se_smd), rng(s$cv)))
}

or <- estimates[estimates$data == "dif" & estimates$set == "oracle", ]
cat("\nOracle set, DIF data: |bias|, SE ratio, coverage\n")
for (q in c("impact_latent", "impact_wle", "impact_sum", "b_maxabs",
            "b_sem_maxabs", "b_lo", "b_hi", "b_cut", "gap_maxabs", "gap_cut")) {
  for (sc in c("not two_opposite", "two_opposite")) {
    y <- or[or$quantity == q & ((or$scenario == "two_opposite") ==
                                  (sc == "two_opposite")), ]
    cat(sprintf("  %-13s %-16s max|bias| %.4f  SE ratio %s  coverage %s\n",
                q, sc, max(abs(y$bias)), rng(y$se_ratio), rng(y$coverage)))
  }
}

cat("\nFlagged set by class (impact_latent, b_maxabs)\n")
fc <- flag_class[flag_class$quantity %in% c("impact_latent", "b_maxabs") &
                   flag_class$n >= 10, ]
print(stats::aggregate(cbind(share, rel_bias, coverage) ~ quantity + class +
                         scenario, fc, function(v) round(range(v), 3)))

cat("\nNo-DIF data, flagged reps (false flags)\n")
nf <- null_flagged[null_flagged$quantity %in% c("impact_latent", "b_maxabs"), ]
print(stats::aggregate(cbind(share, mean_abs, q95_abs, ci_excludes_0) ~
                         quantity + k, nf, function(v) round(range(v), 3)))

cat("\nDetection, DIF data\n")
dd <- detection[detection$data == "dif", ]
print(dd[, c("cell", "anchor_error", "any_flag", "all_true_flagged",
             "any_fp")], row.names = FALSE, digits = 3)
