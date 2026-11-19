# Analysis of the RMdifImpact() step 3 simulation (dev/dif_impact_ext_sim.R).
#
# Reads dev/dif_impact_ext_sim/ext_results.rds (about 75 MB, not tracked) and
# writes dev/dif_impact_ext_sim/ext_results_slim.rds (tracked).
#
#   LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8 Rscript dev/dif_impact_ext_sim_analysis.R
#
# Tables in the slim file:
#   cells       the 42 cells
#   summary     summarise_ext() per cell: target, bias, SE ratio, coverage
#   item_check  per-item impacts with the ratio SE |impact_j / d_j| * SE(d_j),
#               SE(d_j) stood in for by the empirical SD of d_j over
#               replications (the stored SE used |Gi_j| * SE(d_j))
#   knife       per cell, strategy and group: distance from the theta cut to
#               the nearest point of the population WLE table, next to the
#               coverage of the theta-scale rates

source("dev/dif_impact_ext_sim.R")
o <- readRDS(file.path(ext_dir, "ext_results.rds"))
stopifnot(all(vapply(o$results, function(x) length(x$errors), 0L) == 0L))

cells <- o$cells
cell_cols <- c("k", "scenario", "position", "d", "n_ref", "n_focal")
summary <- do.call(rbind, lapply(seq_along(o$results), function(i) {
  cbind(cell = i, cells[i, cell_cols], o$results[[i]]$summary,
        row.names = NULL)
}))
summary$q <- sub("^item_[^_]+_", "item_", summary$quantity)

# ---------------------------------------------------------------------------
# Per-item impacts: ratio SE, checked post hoc
# ---------------------------------------------------------------------------

item_check <- do.call(rbind, lapply(which(cells$scenario != "one"), function(i) {
  r <- o$results[[i]]
  e <- r$estimates[grepl("^item_", r$estimates$quantity), ]
  dh <- do.call(rbind, lapply(strsplit(e$d_hat, ","), as.numeric))
  it <- sub("^item_([^_]+)_.*", "\\1", e$quantity)
  dj <- dh[cbind(seq_len(nrow(e)), match(it, sort(unique(it))))]
  sd_d <- tapply(dj, e$quantity, stats::sd)[e$quantity]
  se_new <- abs(e$estimate / dj) * sd_d
  tg <- r$truth$target[paste0("sample.", e$quantity)]
  do.call(rbind, lapply(split(seq_len(nrow(e)), e$quantity), function(z) {
    data.frame(cell = i, cells[i, cell_cols], quantity = e$quantity[z[1]],
               emp_sd = stats::sd(e$estimate[z]),
               se_ratio_old = mean(e$se[z]) / stats::sd(e$estimate[z]),
               se_ratio_new = mean(se_new[z]) / stats::sd(e$estimate[z]),
               cov_old = mean(abs(e$estimate[z] - tg[z]) <= 1.96 * e$se[z]),
               cov_new = mean(abs(e$estimate[z] - tg[z]) <= 1.96 * se_new[z]),
               row.names = NULL)
  }))
}))

# ---------------------------------------------------------------------------
# Knife-edge: population WLE tables against the theta cut
# ---------------------------------------------------------------------------

knife <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  cell <- as.list(cells[i, ])
  sc <- make_scale(cell$type, cell$k)
  si <- scenario_items(sc$loc, cell$scenario, cell$position, cell$cut)
  j <- si$j
  thr_split <- with_shifts(list(ref = sc$thr, focal = sc$thr), j,
                           cell$d * si$sign)
  pt <- rate_parts(thr_split, j, cut_raw_for(cell))
  keep <- pop_split_delta(sc$thr, thr_split$focal, integer(0), cell$mu2,
                          c(cell$n_ref, cell$n_focal))
  common <- setdiff(seq_along(sc$thr), j)
  wk <- wle_by_score(keep$thr$ref) + link_offset(keep$thr$ref, sc$thr, common)
  tabs <- list(keep = list(ref = wk, focal = wk),
               split = list(ref = pt$wle_ref, focal = pt$wle_focal))
  s <- summary[summary$cell == i, ]
  do.call(rbind, lapply(c("keep", "split"), function(st) {
    do.call(rbind, lapply(c("ref", "focal"), function(g) {
      w <- tabs[[st]][[g]]
      dist <- w - cell$cut
      near <- which.min(abs(dist))
      cov <- s$coverage[s$context == "sample" &
                          s$quantity %in% paste("theta", st, g,
                                                c("prev", "fn", "fp"), sep = ".")]
      cov_p <- s$coverage[s$context == "projection" &
                            s$quantity %in% paste("theta", st, g,
                                                  c("prev", "fn", "fp"), sep = ".")]
      data.frame(cell = i, cells[i, cell_cols], strategy = st, group = g,
                 nearest_score = near - 1L, dist = dist[near],
                 cov_min_sample = min(cov), cov_min_projection = min(cov_p),
                 row.names = NULL)
    }))
  }))
}))

slim <- list(cells = cells, reps = o$reps, n_draws = o$n_draws,
             target_population = o$target_population, run_time = o$run_time,
             date = o$date, summary = summary, item_check = item_check,
             knife = knife)
saveRDS(slim, file.path(ext_dir, "ext_results_slim.rds"))
saveRDS(slim, file.path(ext_dir, sprintf("ext_results_slim_%s.rds",
                                         format(o$date, "%Y-%m-%d"))))

# ---------------------------------------------------------------------------
# Printed overview (the numbers quoted in dev/dif-impact-design.md)
# ---------------------------------------------------------------------------

rng <- function(v) sprintf("%.3f to %.3f (median %.3f)", min(v), max(v),
                           stats::median(v))
cat("Impacts (joint, projected) and projected bias/SE\n")
for (ctx in c("sample", "projection")) {
  for (q in c("impact_latent", "impact_wle", "impact_sum", "bias_in_se_wle")) {
    y <- summary[summary$context == ctx & summary$q == q, ]
    if (!nrow(y)) next
    cat(sprintf("  %-10s %-15s coverage %s | SE ratio %s\n", ctx, q,
                rng(y$coverage), rng(y$se_ratio)))
  }
}

cat("\nPer-item impacts, stored SE |Gi_j| SE(d_j) vs ratio SE\n")
cat(sprintf("  stored: coverage %s | SE ratio %s\n", rng(item_check$cov_old),
            rng(item_check$se_ratio_old)))
cat(sprintf("  ratio:  coverage %s | SE ratio %s\n", rng(item_check$cov_new),
            rng(item_check$se_ratio_new)))

cat("\nRates at the cut, coverage of the draws intervals (keep / split only)\n")
for (sc_ in c("sum", "theta")) for (ctx in c("sample", "projection")) {
  y <- summary[summary$context == ctx & grepl(paste0("^", sc_, "\\.(keep|split)"),
                                              summary$q), ]
  for (kk in c(5, 10)) {
    z <- y[y$k == kk, ]
    cat(sprintf("  %-5s %-10s k%-2d coverage %s\n", sc_, ctx, kk,
                rng(z$coverage)))
  }
}

cat("\nKnife-edge: theta-rate coverage by distance of the nearest WLE to the cut\n")
knife$band <- cut(abs(knife$dist), c(0, 0.01, 0.05, 0.1, Inf),
                  include.lowest = TRUE)
print(stats::aggregate(cbind(cov_min_sample, cov_min_projection) ~ band,
                       knife, function(v) round(c(n = length(v), min = min(v),
                                                  median = stats::median(v)), 3)))
