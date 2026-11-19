# Analysis of the RMdifImpact() step 3c simulation (dev/dif_impact_ext_sim_3c.R):
# rates at a sum-score cut with the package's rule, every estimate from the
# package's .dimp_estimate().
#
# Reads dev/dif_impact_ext_sim/ext3c_results.rds (not tracked) and writes
# dev/dif_impact_ext_sim/ext3c_results_slim.rds (tracked), plus a dated copy.
#
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript dev/dif_impact_ext_sim_3c_analysis.R
#
# Tables in the slim file:
#   summary   summarise_3c() per cell: coverage and SE ratio per rate, over
#             all replications (subset "all") and, for keep - split, over
#             the replications where that group's fitted rules differ
#             (subset "rule_differs"). With strategy and a degenerate flag
#             (error SD zero, or estimate = target = 0 in every replication).
#   rule      per cell: share of replications where the focal group's
#             fitted keep and split rules differ, share where each fitted
#             rule equals the package rule at the true parameters, and the
#             number of replications with a warning.
#   diff_nz   keep - split over the replications where the group's rules
#             differ, per cell, context and quantity (reps_nonzero, coverage,
#             SE ratio), as in the 3b slim file.

ext_dir <- "dev/dif_impact_ext_sim"
o <- readRDS(file.path(ext_dir, "ext3c_results.rds"))
stopifnot(identical(o$scales, "sum"),
          all(vapply(o$results, function(x) length(x$errors), 0L) == 0L),
          all(vapply(o$results, function(x) length(unique(x$estimates$rep)),
                     0L) == o$reps),
          all(vapply(o$results, function(x) {
            e <- x$estimates
            !anyNA(e$target[e$context != "rule"]) &&
              !anyNA(e$se[e$context != "rule"])
          }, TRUE)))

cells <- o$cells
cell_cols <- c("k", "type", "scenario", "position", "d", "n_ref", "n_focal")

summary <- do.call(rbind, lapply(seq_along(o$results), function(i) {
  cbind(cell = i, cells[i, cell_cols], o$results[[i]]$summary,
        row.names = NULL)
}))
summary$strategy <- sub("^([a-z]+)\\..*", "\\1", summary$quantity)
summary$group <- sub("^[a-z]+\\.([a-z]+)\\..*", "\\1", summary$quantity)
summary$rate <- sub(".*\\.", "", summary$quantity)
summary$degenerate <- summary$sd_err < 1e-10 | summary$share_zero == 1

rule <- do.call(rbind, lapply(seq_along(o$results), function(i) {
  e <- o$results[[i]]$estimates
  per_rep <- e[!duplicated(e$rep), ]
  r <- e[e$context == "rule", ]
  match <- vapply(split(r, r$quantity), function(x) mean(x$estimate == x$target),
                  numeric(1L))
  data.frame(cell = i, cells[i, cell_cols],
             cut_raw = o$results[[i]]$truth$cut_raw,
             true_split_focal = o$results[[i]]$truth$rule_cuts_sum[["split.focal"]],
             focal_rule_differs = mean(per_rep$rule_differs_focal),
             ref_rule_differs = mean(per_rep$rule_differs_ref),
             match_keep_focal = match[["keep.focal"]],
             match_split_focal = match[["split.focal"]],
             match_split_ref = match[["split.ref"]],
             reps_with_warning = sum(per_rep$warnings > 0),
             row.names = NULL)
}))

diff_nz <- summary[summary$subset == "rule_differs", ]
names(diff_nz)[names(diff_nz) == "reps"] <- "reps_nonzero"

slim <- list(cells = cells, reps = o$reps, target_population = o$target_population,
             run_time = o$run_time, date = o$date,
             package_version = o$package_version,
             summary = summary[summary$subset == "all", ], rule = rule,
             diff_nz = diff_nz)
saveRDS(slim, file.path(ext_dir, "ext3c_results_slim.rds"))
saveRDS(slim, file.path(ext_dir, sprintf("ext3c_results_slim_%s.rds",
                                         format(o$date, "%Y-%m-%d"))))

# ---------------------------------------------------------------------------
# Printed overview
# ---------------------------------------------------------------------------

rng <- function(v) {
  v <- v[!is.na(v)]
  sprintf("%.3f to %.3f (median %.3f, n = %d)", min(v), max(v),
          stats::median(v), length(v))
}
mcse <- sqrt(.95 * .05 / o$reps)
lo <- .95 - 2 * mcse
s <- slim$summary[!slim$summary$degenerate, ]
cat(sprintf("%d cells x %d replications (easyRasch2 %s), run time %s\n",
            nrow(cells), o$reps, o$package_version, format(o$run_time)))
cat(sprintf("Monte Carlo SE of a .95 coverage %.4f; 2 MCSE band from %.3f\n", mcse, lo))
cat(sprintf("Degenerate rows dropped: %d of %d; replications with a warning: %d\n\n",
            sum(slim$summary$degenerate), nrow(slim$summary),
            sum(rule$reps_with_warning)))

cat("Rates of keep and split, fixed fitted rule, sum-score cut\n")
for (ctx in c("sample", "projection")) {
  y <- s[s$context == ctx & s$strategy %in% c("keep", "split"), ]
  cat(sprintf("  %-10s coverage %s\n", ctx, rng(y$coverage)))
  cat(sprintf("  %-10s SE ratio %s\n", "", rng(y$se_ratio)))
  cat(sprintf("  %-10s rows at or below %.3f: %d of %d\n", "", lo,
              sum(y$coverage <= lo), nrow(y)))
}

cat("\nKeep - split, replications where the group's rules differ\n")
for (ctx in c("sample", "projection")) {
  y <- diff_nz[diff_nz$context == ctx, ]
  y100 <- y[y$reps_nonzero >= 100, ]
  cat(sprintf("  %-10s pooled coverage %.3f over %d reps; per row (>= 100 reps) %s\n",
              ctx, sum(y$coverage * y$reps_nonzero) / sum(y$reps_nonzero),
              sum(y$reps_nonzero), rng(y100$coverage)))
  cat(sprintf("  %-10s rows (>= 100 reps) below .93: %d of %d\n", "",
              sum(y100$coverage < .93), nrow(y100)))
}

cat("\nFocal rule: share of replications where keep and split differ, per cell\n")
for (sn in unique(rule$scenario)) for (dd in sort(unique(rule$d))) {
  y <- rule[rule$scenario == sn & rule$d == dd, ]
  if (!nrow(y)) next
  cat(sprintf("  %-13s d = %.1f  %s\n", sn, dd, rng(y$focal_rule_differs)))
}
cat(sprintf("  reference rule ever differs: %s\n", any(rule$ref_rule_differs > 0)))

cat("\nFitted split focal rule = package rule at the true parameters\n")
cat(sprintf("  %s\n", rng(rule$match_split_focal)))

# Step 3b (dev rule, eq >= cut) on the same cells, for the dev write-up
f3b <- file.path(ext_dir, "ext3b_results_slim.rds")
if (file.exists(f3b)) {
  b <- readRDS(f3b)
  bs <- b$summary[b$summary$context != "rule" & !b$summary$degenerate &
                    b$summary$scale == "sum" &
                    b$summary$strategy %in% c("keep", "split"), ]
  cat("\nStep 3b, same cells, dev rule: keep/split sum-score coverage\n")
  for (ctx in c("sample", "projection")) {
    cat(sprintf("  %-10s %s\n", ctx, rng(bs$cov_delta[bs$context == ctx])))
  }
  sa <- b$summary[b$summary$context == "sample" &
                    b$summary$quantity == "sum.diff.focal.prev", ]
  cat(sprintf("  focal rule differs (share of reps per cell): %s\n",
              rng(1 - sa$share_zero)))
}
