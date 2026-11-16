# Calibration and power of threshold-level DIF tests (two groups).
#
# Questions:
#   1. Spacing test (no anchor): is the per-item Wald chi-square calibrated,
#      and do Holm and the single-step min-p keep the family-wise error rate
#      at .05 over items without spacing DSF? Uniform DIF must not trigger it.
#   2. Anchored threshold differences: with the forward (MTT) and Gini
#      anchors, is the family-wise error rate over thresholds without DIF
#      kept at .05, and what is the power? The centred anchor (the old
#      RMdifLR() threshold scale) is included as the reference that fails.
#   3. The anchored item-location test is computed too, to compare what it
#      finds for DSF items with what the threshold tests find. It also extends
#      the item-level study in ../dif_lr_wald/ to unequal group sizes.
#   4. Tie-break sensitivity: psychotools breaks ties in the anchor ranking at
#      random. In the first `tie_reps` replications of each cell, the forward
#      and Gini rankings are redone under 20 seeds, and the anchors and
#      single-step flags are compared with the main run.
#
# Design (144 cells):
#   n               mean group size 150, 300, 600
#   ratio           equal (n, n) or unequal (1.5 n, 0.5 n; group 2 is the
#                   smaller, focal group with the impact and the DIF)
#   structure       10 items x 4 categories (3 thresholds),
#                   20 items x 3 categories (2 thresholds)
#   impact          group 2 latent mean 0 or 0.5 (SD 1 in both groups)
#   scenario        shifts of group 2's thresholds, per DIF item:
#     null          none
#     unif1         one item, every threshold +0.6 (location DIF only)
#     unif3         three items, every threshold +0.6 (unbalanced)
#     sym1          one item, lowest -0.6 and highest +0.6
#                   (symmetric DSF, location unchanged)
#     nonsym1       one item, highest threshold +0.6
#                   (non-symmetric DSF, location +0.2 or +0.3)
#     nonsym3       three items, highest threshold +0.6 (unbalanced)
#
# Null sets for the error rates:
#   spacing     items whose shift vector is constant (no spacing change)
#   thresholds  thresholds with shift 0
#   locations   items with mean shift 0
#
# Usage (from the easyRasch2 package root, on the host):
#   nohup Rscript dev/dif_anchor_threshold/dif_thresh_sim.R 1000 10 > dev/dif_anchor_threshold/run.log 2>&1 &
# Arguments: R (replications per cell, default 1000), cores (default 10),
#   out_dir (default dev/dif_anchor_threshold), B (Monte-Carlo draws for the
#   single-step and min-p adjustments, default 5000), tie_reps (replications
#   per cell with the tie-break check, default 250).
# The Monte-Carlo adjustments of the anchored tests use a fixed seed per
# replication, so in the tie-break check flags differ only through the anchor.
# Needs psychotools, multcomp and mvtnorm. Runs with parallel::mclapply.
# Each replication has its own seed, so results do not depend on `cores`.

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 10L
out_dir <- if (length(args) >= 3L) args[3L] else "dev/dif_anchor_threshold"
B <- if (length(args) >= 4L) as.integer(args[4L]) else 5000L
tie_reps <- if (length(args) >= 5L) as.integer(args[5L]) else 250L
n_tie_seeds <- 20L

suppressMessages({
  library(psychotools)
  library(multcomp)
})
source(file.path("dev", "dif_anchor_threshold", "dif_thresh.R"))

anchors <- c("centred", "forward", "gini")
alpha <- 0.05
size <- 0.6

# --- Data generation ---------------------------------------------------------
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
  "20x3" = list(
    tau = lapply(seq(-1.5, 1.5, length.out = 20), function(l) l + c(-0.75, 0.75)),
    dif_items = c(5L, 11L, 16L)
  )
)

# Shift vector for one DIF item with m thresholds.
item_shift <- function(type, m) {
  switch(type,
    unif = rep(size, m),
    sym = c(-size, rep(0, m - 2L), size),
    nonsym = c(rep(0, m - 1L), size)
  )
}
scenario_items <- function(scenario) {
  switch(scenario,
    null = character(0),
    unif1 = "unif", unif3 = rep("unif", 3L),
    sym1 = "sym",
    nonsym1 = "nonsym", nonsym3 = rep("nonsym", 3L)
  )
}

design <- expand.grid(
  n = c(150L, 300L, 600L),
  ratio = c("equal", "unequal"),
  structure = names(structures),
  impact = c(0, 0.5),
  scenario = c("null", "unif1", "unif3", "sym1", "nonsym1", "nonsym3"),
  stringsAsFactors = FALSE
)
design$cell <- seq_len(nrow(design))

# One summary row per family, anchor and correction.
summ_rows <- function(family, anchor, padj_list, null_set, z = NULL) {
  do.call(rbind, lapply(names(padj_list), function(cr) {
    pa <- padj_list[[cr]]
    tested <- !is.na(pa)
    flag <- tested & pa < alpha
    data.frame(
      family = family, anchor = anchor, correction = cr,
      n_null = sum(null_set & tested), n_false = sum(flag & null_set),
      n_alt = sum(!null_set), n_alt_tested = sum(!null_set & tested),
      n_found = sum(flag & !null_set),
      z_n = if (is.null(z)) NA_real_ else sum(null_set & !is.na(z)),
      z_sum = if (is.null(z)) NA_real_ else sum(z[null_set], na.rm = TRUE),
      z_ss = if (is.null(z)) NA_real_ else sum(z[null_set]^2, na.rm = TRUE)
    )
  }))
}

# --- One replication -----------------------------------------------------------
one_rep <- function(cell, rep) {
  d <- design[cell, ]
  st <- structures[[d$structure]]
  set.seed(cell * 100000L + rep)
  types <- scenario_items(d$scenario)
  dif_items <- st$dif_items[seq_along(types)]
  k <- length(st$tau)
  shift <- lapply(st$tau, function(t) rep(0, length(t)))
  for (j in seq_along(dif_items)) {
    shift[[dif_items[j]]] <- item_shift(types[j], length(st$tau[[dif_items[j]]]))
  }
  tau2 <- Map(`+`, st$tau, shift)
  n1 <- if (d$ratio == "equal") d$n else round(1.5 * d$n)
  n2 <- 2L * d$n - n1
  x1 <- sim_pcm(stats::rnorm(n1, 0, 1), st$tau)
  x2 <- sim_pcm(stats::rnorm(n2, d$impact, 1), tau2)
  colnames(x1) <- colnames(x2) <- paste0("I", seq_len(k))
  ncat <- vapply(st$tau, length, integer(1)) + 1L
  ok <- all(vapply(seq_len(k), function(i) {
    length(unique(x1[, i])) == ncat[i] && length(unique(x2[, i])) == ncat[i]
  }, logical(1)))
  if (!ok) return(list(status = "empty_category"))
  mc_seed <- 1000000000L + cell * 100000L + rep

  res <- tryCatch({
    fit <- dif_thr_fit(rbind(x1, x2), rep(1:2, c(n1, n2)))
    s_thr <- unlist(shift)
    sp_null <- vapply(shift, function(s) all(s == s[1L]), logical(1))
    loc_null <- vapply(shift, function(s) mean(s) == 0, logical(1))
    dif_any <- vapply(shift, function(s) any(s != 0), logical(1))
    sp <- dif_thr_spacing(fit, B = B)
    out <- list(summ_rows(
      "spacing", "none",
      list(none = sp$p, holm = sp$padj_holm, bh = sp$padj_bh, ss = sp$padj_ss),
      sp_null[sp$item_no],
      z = stats::qnorm(sp$p, lower.tail = FALSE)
    ))
    ties <- list()
    flags_of <- function(r) list(
      loc = !is.na(r$locations$padj_ss) & r$locations$padj_ss < alpha,
      thr = !is.na(r$thresholds$padj_ss) & r$thresholds$padj_ss < alpha
    )
    for (an in anchors) {
      A <- dif_thr_anchor(fit, an)
      set.seed(mc_seed)
      r <- dif_thr_anchored(fit, A, drop_first = an != "centred", B = B)
      th <- r$thresholds
      lo <- r$locations
      if (an != "centred" && rep <= tie_reps) {
        f0 <- flags_of(r)
        key0 <- paste(A[1L], paste(sort(A), collapse = ","))
        As <- dif_thr_anchor_ties(fit, an, mc_seed + seq_len(n_tie_seeds))
        cache <- list()
        tf <- lapply(As, function(At) {
          key <- paste(At[1L], paste(sort(At), collapse = ","))
          if (key == key0) return(f0)
          if (is.null(cache[[key]])) {
            set.seed(mc_seed)
            cache[[key]] <<- flags_of(dif_thr_anchored(fit, At, B = B))
          }
          cache[[key]]
        })
        loc_alt <- !loc_null
        thr_alt <- s_thr != 0
        loc_found <- vapply(tf, function(f) sum(f$loc & loc_alt), numeric(1))
        thr_found <- vapply(tf, function(f) sum(f$thr & thr_alt), numeric(1))
        ties[[an]] <- data.frame(
          anchor = an, n_seeds = n_tie_seeds,
          n_unique_set = length(unique(lapply(As, sort))),
          p_set_diff = mean(vapply(As, function(At) !identical(sort(At), sort(A)), logical(1))),
          p_first_diff = mean(vapply(As, function(At) At[1L] != A[1L], logical(1))),
          p_loc_flagdiff = mean(vapply(tf, function(f) any(f$loc != f0$loc), logical(1))),
          p_thr_flagdiff = mean(vapply(tf, function(f) any(f$thr != f0$thr), logical(1))),
          loc_found_main = sum(f0$loc & loc_alt),
          loc_found_min = min(loc_found), loc_found_max = max(loc_found),
          thr_found_main = sum(f0$thr & thr_alt),
          thr_found_min = min(thr_found), thr_found_max = max(thr_found),
          loc_false_max = max(vapply(tf, function(f) sum(f$loc & !loc_alt), numeric(1))),
          thr_false_max = max(vapply(tf, function(f) sum(f$thr & !thr_alt), numeric(1)))
        )
      }
      rows <- rbind(
        summ_rows("thresholds", an,
                  list(none = th$p, holm = th$padj_holm, ss = th$padj_ss),
                  s_thr[seq_len(nrow(th))] == 0, z = th$z),
        summ_rows("locations", an,
                  list(none = lo$p, holm = lo$padj_holm, ss = lo$padj_ss),
                  loc_null, z = lo$z)
      )
      rows$anchor_len <- length(A)
      rows$anchor_locdif <- sum(!loc_null[A])
      rows$anchor_anyshift <- sum(dif_any[A])
      rows$ss_failed <- c(all(is.na(th$padj_ss)), all(is.na(lo$padj_ss)))[
        match(rows$family, c("thresholds", "locations"))]
      out[[length(out) + 1L]] <- rows
    }
    out[[1L]]$anchor_len <- NA_integer_
    out[[1L]]$anchor_locdif <- NA_integer_
    out[[1L]]$anchor_anyshift <- NA_integer_
    out[[1L]]$ss_failed <- FALSE
    list(rows = do.call(rbind, out),
         ties = if (length(ties)) do.call(rbind, ties) else NULL)
  }, error = function(e) conditionMessage(e))
  if (is.character(res)) return(list(status = paste("error:", res)))
  res$rows$cell <- cell
  res$rows$rep <- rep
  if (!is.null(res$ties)) {
    res$ties$cell <- cell
    res$ties$rep <- rep
  }
  res$status <- "ok"
  res
}

# --- Run ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Design: %d cells x %d replications on %d cores, B = %d, tie check in %d replications per cell\n",
            nrow(design), R, cores, B, tie_reps))
t_start <- Sys.time()
all_res <- vector("list", nrow(design))
for (cell in design$cell) {
  t0 <- Sys.time()
  reps <- parallel::mclapply(seq_len(R), function(r) one_rep(cell, r),
                             mc.cores = cores, mc.preschedule = TRUE)
  bad <- vapply(reps, function(x) inherits(x, "try-error"), logical(1))
  if (any(bad)) stop("Worker failure in cell ", cell, ": ", reps[[which(bad)[1L]]])
  st_rows <- data.frame(cell = cell, rep = seq_len(R),
                        status = vapply(reps, `[[`, character(1), "status"))
  all_res[[cell]] <- list(
    per_rep = do.call(rbind, lapply(reps, `[[`, "rows")),
    ties = do.call(rbind, lapply(reps, `[[`, "ties")),
    status = st_rows
  )
  d <- design[cell, ]
  cat(sprintf(
    "[%s] cell %3d/%d  n=%d %-7s %s impact=%.1f %-7s  ok %d/%d  %.1f min\n",
    format(Sys.time(), "%H:%M"), cell, nrow(design), d$n, d$ratio, d$structure,
    d$impact, d$scenario, sum(st_rows$status == "ok"), R,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}
per_rep <- do.call(rbind, lapply(all_res, `[[`, "per_rep"))
ties <- do.call(rbind, lapply(all_res, `[[`, "ties"))
status <- do.call(rbind, lapply(all_res, `[[`, "status"))
status <- merge(status, design, by = "cell")

# --- Summary -----------------------------------------------------------------
# fwer    share of replications with any flag in the null set
# power   mean share of non-null units flagged (NA when there are none)
# untested  share of replications where a non-null unit got no test
# sd_z    SD of the null-set z (spacing: qnorm of p, one-sided), pooled
# anchor_locdif    share of replications whose anchor holds an item with
#                  location DIF (mean shift not 0), which contaminates it
# anchor_anyshift  the same for any shifted threshold, including symmetric
#                  DSF items, whose location is DIF-free
key <- c("cell", "family", "anchor", "correction")
grp <- split(per_rep, per_rep[key], drop = TRUE)
summ <- do.call(rbind, lapply(grp, function(x) {
  zn <- sum(x$z_n)
  data.frame(
    x[1L, key], reps = nrow(x),
    fwer = mean(x$n_false > 0),
    mean_false = mean(x$n_false),
    power = if (x$n_alt[1L] > 0) mean(x$n_found / x$n_alt) else NA_real_,
    untested = if (x$n_alt[1L] > 0) mean(x$n_alt_tested < x$n_alt) else NA_real_,
    sd_z = if (zn > 1) sqrt((sum(x$z_ss) - sum(x$z_sum)^2 / zn) / (zn - 1)) else NA_real_,
    anchor_len = mean(x$anchor_len),
    anchor_locdif = mean(x$anchor_locdif > 0),
    anchor_anyshift = mean(x$anchor_anyshift > 0),
    ss_failed = sum(x$ss_failed)
  )
}))
summ <- merge(design, summ, by = "cell")
summ <- summ[order(summ$family, summ$scenario, summ$ratio, summ$structure,
                   summ$impact, summ$n, summ$anchor, summ$correction), ]

# Tie-break summary per cell and anchor:
# set_diff, first_diff    mean share of seeds giving another anchor set, or
#                         another untested first anchor item
# loc/thr_flagdiff        mean share of seeds that change any single-step flag
# loc/thr_any_change      share of replications where some seed changes a flag
# loc/thr_found_varies    share of replications where the number of true DIF
#                         units found differs between seeds
tsum <- do.call(rbind, lapply(split(ties, ties[c("cell", "anchor")], drop = TRUE), function(x) {
  data.frame(
    cell = x$cell[1L], anchor = x$anchor[1L], reps = nrow(x),
    unique_sets = mean(x$n_unique_set),
    set_diff = mean(x$p_set_diff), first_diff = mean(x$p_first_diff),
    loc_flagdiff = mean(x$p_loc_flagdiff), thr_flagdiff = mean(x$p_thr_flagdiff),
    loc_any_change = mean(x$p_loc_flagdiff > 0),
    thr_any_change = mean(x$p_thr_flagdiff > 0),
    loc_found_varies = mean(x$loc_found_min < x$loc_found_max),
    thr_found_varies = mean(x$thr_found_min < x$thr_found_max)
  )
}))
tsum <- merge(design, tsum, by = "cell")
tsum <- tsum[order(tsum$anchor, tsum$scenario, tsum$ratio, tsum$structure,
                   tsum$impact, tsum$n), ]
fails <- as.data.frame(table(cell = status$cell, status = status$status))
fails <- fails[fails$Freq > 0 & fails$status != "ok", ]

out <- list(
  design = design, per_rep = per_rep, status = status, summary = summ,
  ties = ties, tie_summary = tsum, tie_reps = tie_reps, n_tie_seeds = n_tie_seeds,
  failures = fails, R = R, B = B, alpha = alpha, size = size,
  anchors = anchors, session = utils::sessionInfo(),
  run_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))
)
stamp <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(out, file.path(out_dir, paste0("dif_thresh_sim_", stamp, ".rds")))
saveRDS(out, file.path(out_dir, "dif_thresh_sim.rds"))
utils::write.csv(summ, file.path(out_dir, "dif_thresh_sim_summary.csv"), row.names = FALSE)
utils::write.csv(tsum, file.path(out_dir, "dif_thresh_sim_ties.csv"), row.names = FALSE)

options(width = 200)
cat("\nFailures (cell, status, count):\n")
print(fails, row.names = FALSE)
cat("\nSummary, single-step correction (fwer = any null unit flagged; power = share of non-null units flagged):\n")
print(summ[summ$correction %in% c("ss", "none"),
           c("family", "scenario", "ratio", "structure", "impact", "n", "anchor",
             "correction", "reps", "fwer", "power", "untested", "sd_z",
             "anchor_locdif", "ss_failed")],
      digits = 3, row.names = FALSE)
cat("\nTie-break sensitivity (forward and Gini):\n")
print(tsum[, c("anchor", "scenario", "ratio", "structure", "impact", "n", "reps",
               "unique_sets", "set_diff", "first_diff", "loc_any_change",
               "thr_any_change", "loc_found_varies", "thr_found_varies")],
      digits = 3, row.names = FALSE)
cat(sprintf("\nTotal run time: %.1f min\n", out$run_minutes))
