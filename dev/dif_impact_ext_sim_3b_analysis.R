# Analysis of the RMdifImpact() step 3b simulation (dev/dif_impact_ext_sim_3b.R).
#
# Reads dev/dif_impact_ext_sim/ext3b_results.rds (about 30 MB, not tracked)
# and writes dev/dif_impact_ext_sim/ext3b_results_slim.rds (tracked).
#
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript dev/dif_impact_ext_sim_3b_analysis.R
#
# Tables in the slim file:
#   summary   summarise_3b() per cell (rates, per-item impacts, rule rows),
#             with scale, strategy and a degenerate flag. Degenerate rows
#             have an error SD of zero or estimate = target = 0 in every
#             replication (keep - split for the reference group on the sum
#             scale, where both strategies use the same raw cut by
#             construction).
#   diff_nz   keep - split rates, restricted to replications where the
#             fitted keep and split rules differ (estimate or target not
#             zero). Where the rules agree, the difference is exactly zero
#             and so is its SE, which inflates the unconditional coverage.

source("dev/dif_impact_ext_sim_3b.R")
o <- readRDS(file.path(ext_dir, "ext3b_results.rds"))
stopifnot(all(vapply(o$results, function(x) length(x$errors), 0L) == 0L),
          all(vapply(o$results, function(x) length(unique(x$estimates$rep)),
                     0L) == o$reps))

cells <- o$cells
cell_cols <- c("k", "type", "scenario", "position", "d", "n_ref", "n_focal")
summary <- do.call(rbind, lapply(seq_along(o$results), function(i) {
  cbind(cell = i, cells[i, cell_cols], o$results[[i]]$summary,
        row.names = NULL)
}))
summary$scale <- ifelse(grepl("^item_", summary$quantity), "item",
                        sub("\\..*", "", summary$quantity))
summary$strategy <- ifelse(summary$scale == "item", "item",
                           sub("^[a-z]+\\.([a-z]+).*", "\\1", summary$quantity))
summary$degenerate <- summary$context != "rule" &
  (summary$sd_err < 1e-10 | summary$share_zero == 1)

diff_nz <- do.call(rbind, lapply(seq_along(o$results), function(i) {
  e <- o$results[[i]]$estimates
  e <- e[grepl("\\.diff\\.", e$quantity) & e$context != "rule", ]
  e <- e[!(e$estimate == 0 & e$target == 0), ]
  if (!nrow(e)) return(NULL)
  do.call(rbind, lapply(split(e, paste(e$context, e$quantity)), function(x) {
    err <- x$estimate - x$target
    data.frame(cell = i, cells[i, cell_cols], context = x$context[1],
               quantity = x$quantity[1], scale = sub("\\..*", "", x$quantity[1]),
               reps_nonzero = nrow(x),
               cov_delta = mean(abs(err) <= 1.96 * x$se_delta + 1e-12),
               cov_draws = mean(x$lo_draws - 1e-12 <= x$target &
                                  x$target <= x$hi_draws + 1e-12),
               se_ratio_delta = if (nrow(x) > 2) mean(x$se_delta) / stats::sd(err)
                                else NA_real_,
               row.names = NULL)
  }))
}))

slim <- list(cells = cells, reps = o$reps, n_draws = o$n_draws,
             target_population = o$target_population, run_time = o$run_time,
             date = o$date, summary = summary, diff_nz = diff_nz)
saveRDS(slim, file.path(ext_dir, "ext3b_results_slim.rds"))
saveRDS(slim, file.path(ext_dir, sprintf("ext3b_results_slim_%s.rds",
                                         format(o$date, "%Y-%m-%d"))))

# ---------------------------------------------------------------------------
# Printed overview (the numbers quoted in dev/dif-impact-design.md)
# ---------------------------------------------------------------------------

rng <- function(v) {
  v <- v[!is.na(v)]
  sprintf("%.3f to %.3f (median %.3f, n = %d)", min(v), max(v),
          stats::median(v), length(v))
}
s <- summary[summary$context != "rule" & !summary$degenerate, ]
cat(sprintf("%d cells x %d replications, %d draws; Monte Carlo SE of a .95 coverage %.4f\n",
            nrow(cells), o$reps, o$n_draws, sqrt(.95 * .05 / o$reps)))
cat(sprintf("Degenerate rows dropped: %d of %d\n\n", sum(summary$degenerate),
            sum(summary$context != "rule")))

cat("Rates per strategy (keep, split), fixed rule\n")
for (ctx in c("sample", "projection")) for (sc_ in c("sum", "theta")) {
  y <- s[s$context == ctx & s$scale == sc_ & s$strategy %in% c("keep", "split"), ]
  cat(sprintf("  %-10s %-5s delta coverage %s\n", ctx, sc_, rng(y$cov_delta)))
  cat(sprintf("  %-10s %-5s draws coverage %s\n", "", "", rng(y$cov_draws)))
  cat(sprintf("  %-10s %-5s delta SE ratio %s\n", "", "", rng(y$se_ratio_delta)))
}
y <- s[s$strategy %in% c("keep", "split"), ]
lo <- .95 - 2 * sqrt(.95 * .05 / o$reps)
cat(sprintf("  rows at or below %.3f (2 Monte Carlo SE): sample %d of %d, projection %d of %d\n",
            lo, sum(y$cov_delta[y$context == "sample"] <= lo), sum(y$context == "sample"),
            sum(y$cov_delta[y$context == "projection"] <= lo),
            sum(y$context == "projection")))

cat("\nDelta SE / draws SE (mean over replications), all non-degenerate rate rows\n")
y <- s[s$scale != "item", ]
cat(sprintf("  %s\n", rng(y$se_ratio_delta / y$se_ratio_draws)))

cat("\nKeep - split, replications where the fitted rules differ (>= 100 such reps)\n")
for (ctx in c("sample", "projection")) for (sc_ in c("sum", "theta")) {
  y <- diff_nz[diff_nz$context == ctx & diff_nz$scale == sc_ &
                 diff_nz$reps_nonzero >= 100, ]
  p <- diff_nz[diff_nz$context == ctx & diff_nz$scale == sc_, ]
  cat(sprintf("  %-10s %-5s delta coverage %s; pooled %.3f over %d reps\n", ctx,
              sc_, rng(y$cov_delta),
              sum(p$cov_delta * p$reps_nonzero) / sum(p$reps_nonzero),
              sum(p$reps_nonzero)))
}

cat("\nPer-item impacts, ratio SE |impact_j / dif_j| * SE(dif_j)\n")
y <- s[s$scale == "item", ]
cat(sprintf("  coverage %s\n  SE ratio %s\n", rng(y$cov_delta),
            rng(y$se_ratio_delta)))

cat("\nFitted rule equals the population rule (share of replications, per cell)\n")
y <- summary[summary$context == "rule", ]
for (q in sort(unique(y$quantity))) {
  cat(sprintf("  %-18s %s\n", q, rng(y$rule_match[y$quantity == q])))
}

cat("\nStep 3 knife-edge cells 14 and 34, theta keep/split, delta coverage\n")
for (i in c(14, 34)) for (ctx in c("sample", "projection")) {
  y <- s[s$cell == i & s$context == ctx & s$scale == "theta" &
           s$strategy %in% c("keep", "split"), ]
  cat(sprintf("  cell %d %-10s %s\n", i, ctx, rng(y$cov_delta)))
}
