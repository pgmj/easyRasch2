# Consensus forward anchor over random tie-breaks (two groups).
#
# Study 2 (../dif_anchor_threshold) showed that the forward (MTT) anchor
# depends on psychotools' random tie-break: the same data give different
# anchors and different flags under different seeds. This study checks a
# consensus anchor built from S seeds (dif_consensus.R).
#
# Questions:
#   1. Does the consensus anchor keep the family-wise error rate of the
#      anchored location and threshold tests (single-step) near .05, and how
#      does its power compare with the single-seed forward anchor and Gini?
#   2. Is it reproducible? Two disjoint seed sets (A and B, S seeds each) are
#      compared: same anchor set, same reference item, same flags. S = 1 is
#      the single-seed forward anchor, for comparison.
#   3. How many seeds are needed (S = 50, 100, 200, nested within each set)?
#
# Methods, per replication:
#   single      forward anchor from one seed (first seed of set A)
#   cons50, cons100, cons200
#               consensus from the first 50, 100, 200 seeds of set A; the
#               reference (first) anchor item is untested, as in psychotools
#   cons200_all consensus from 200 seeds, every item tested when the anchor
#               has two or more items (no untested reference item)
#   gini        Gini anchor point (Strobl et al., 2021), for comparison
#
# Design (72 cells): as Study 2 with impact 0 only (impact 0.5 changed the
# Study 2 error rates by at most .027):
#   n 150, 300, 600 (mean group size); equal or unequal (1.5 n, 0.5 n) sizes;
#   10 items x 4 categories or 20 items x 3 categories; scenarios null,
#   unif1, unif3, sym1, nonsym1, nonsym3 (shift 0.6).
# Data are generated with the same seeds as Study 2, so the datasets are
# the impact-0 datasets of that run.
#
# Usage (from the easyRasch2 package root, on the host):
#   nohup Rscript dev/dif_anchor_consensus/dif_consensus_sim.R 1000 10 > dev/dif_anchor_consensus/run.log 2>&1 &
# Arguments: R (replications per cell, default 1000), cores (default 10),
#   out_dir (default dev/dif_anchor_consensus), B (Monte-Carlo draws for the
#   single-step adjustment, default 5000), S_max (seeds per set, default 200).
# Needs psychotools, multcomp and mvtnorm. Runs with parallel::mclapply.

args <- commandArgs(trailingOnly = TRUE)
R <- if (length(args) >= 1L) as.integer(args[1L]) else 1000L
cores <- if (length(args) >= 2L) as.integer(args[2L]) else 10L
out_dir <- if (length(args) >= 3L) args[3L] else "dev/dif_anchor_consensus"
B <- if (length(args) >= 4L) as.integer(args[4L]) else 5000L
S_max <- if (length(args) >= 5L) as.integer(args[5L]) else 200L
S_grid <- c(1L, 50L, 100L, 200L)
S_grid <- S_grid[S_grid <= S_max]

suppressMessages({
  library(psychotools)
  library(multcomp)
})
source(file.path("dev", "dif_anchor_threshold", "dif_thresh.R"))
source(file.path("dev", "dif_anchor_consensus", "dif_consensus.R"))

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
  impact = 0,
  scenario = c("null", "unif1", "unif3", "sym1", "nonsym1", "nonsym3"),
  stringsAsFactors = FALSE
)
design$cell <- seq_len(nrow(design))

# Study 2's cell numbers, so that the data seeds match that run.
design_s2 <- expand.grid(
  n = c(150L, 300L, 600L), ratio = c("equal", "unequal"),
  structure = names(structures), impact = c(0, 0.5),
  scenario = c("null", "unif1", "unif3", "sym1", "nonsym1", "nonsym3"),
  stringsAsFactors = FALSE
)
design_s2$cell <- seq_len(nrow(design_s2))
design <- merge(design, design_s2, by = c("n", "ratio", "structure", "impact", "scenario"),
                suffixes = c("", "_s2"))
design <- design[order(design$cell), ]
rownames(design) <- NULL

one_rep <- function(cell, rep) {
  d <- design[cell, ]
  st <- structures[[d$structure]]
  set.seed(d$cell_s2 * 100000L + rep)
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
  mc_seed <- 1000000000L + d$cell_s2 * 100000L + rep

  res <- tryCatch({
    fit <- dif_thr_fit(rbind(x1, x2), rep(1:2, c(n1, n2)))
    s_thr <- unlist(shift)
    loc_null <- vapply(shift, function(s) mean(s) == 0, logical(1))
    thr_null <- s_thr == 0
    # Two disjoint seed sets of S_max seeds each.
    seedsA <- mc_seed + seq_len(S_max)
    seedsB <- mc_seed + S_max + seq_len(S_max)
    As_all <- dif_thr_anchor_ties(fit, "forward", c(seedsA, seedsB))
    AsA <- As_all[seq_len(S_max)]
    AsB <- As_all[S_max + seq_len(S_max)]
    cons <- function(As, S) {
      if (S == 1L) As[[1L]] else dif_consensus_from(As[seq_len(S)], k)$anchor
    }
    cache <- list()
    flags <- function(A, drop_first = TRUE) {
      key <- paste(drop_first, A[1L], paste(sort(A), collapse = ","))
      if (is.null(cache[[key]])) {
        set.seed(mc_seed)
        r <- dif_thr_anchored(fit, A, drop_first = drop_first, B = B)
        cache[[key]] <<- list(
          loc = !is.na(r$locations$padj_ss) & r$locations$padj_ss < alpha,
          loc_tested = !is.na(r$locations$padj_ss),
          thr = !is.na(r$thresholds$padj_ss) & r$thresholds$padj_ss < alpha,
          thr_tested = !is.na(r$thresholds$padj_ss)
        )
      }
      cache[[key]]
    }
    row <- function(method, A, f) {
      data.frame(
        method = method, family = c("locations", "thresholds"),
        n_false = c(sum(f$loc & loc_null), sum(f$thr & thr_null)),
        n_alt = c(sum(!loc_null), sum(!thr_null)),
        n_alt_tested = c(sum(f$loc_tested & !loc_null), sum(f$thr_tested & !thr_null)),
        n_found = c(sum(f$loc & !loc_null), sum(f$thr & !thr_null)),
        anchor_len = length(A),
        anchor_locdif = sum(!loc_null[A])
      )
    }
    out <- list()
    agree <- list()
    for (S in S_grid) {
      a <- cons(AsA, S)
      b <- cons(AsB, S)
      fa <- flags(a)
      fb <- flags(b)
      out[[length(out) + 1L]] <- row(if (S == 1L) "single" else paste0("cons", S), a, fa)
      fr <- if (S == 1L) rep(NA_real_, k) else dif_consensus_from(AsA[seq_len(S)], k)$freq
      agree[[length(agree) + 1L]] <- data.frame(
        S = S,
        same_set = identical(sort(a), sort(b)),
        same_ref = a[1L] == b[1L],
        loc_same = identical(fa$loc, fb$loc),
        thr_same = identical(fa$thr, fb$thr),
        found_diff_loc = sum(fa$loc & !loc_null) != sum(fb$loc & !loc_null),
        found_diff_thr = sum(fa$thr & !thr_null) != sum(fb$thr & !thr_null),
        n_near = sum(fr > 0.4 & fr < 0.6)
      )
    }
    a200 <- cons(AsA, max(S_grid))
    b200 <- cons(AsB, max(S_grid))
    drop <- length(a200) < 2L
    fa <- flags(a200, drop_first = drop)
    fb <- flags(b200, drop_first = length(b200) < 2L)
    out[[length(out) + 1L]] <- row("cons_all", a200, fa)
    agree[[length(agree) + 1L]] <- data.frame(
      S = -max(S_grid), same_set = identical(sort(a200), sort(b200)),
      same_ref = NA, loc_same = identical(fa$loc, fb$loc),
      thr_same = identical(fa$thr, fb$thr),
      found_diff_loc = sum(fa$loc & !loc_null) != sum(fb$loc & !loc_null),
      found_diff_thr = sum(fa$thr & !thr_null) != sum(fb$thr & !thr_null),
      n_near = NA_integer_
    )
    g <- dif_thr_anchor(fit, "gini")
    out[[length(out) + 1L]] <- row("gini", g, flags(g))
    fr200 <- dif_consensus_from(AsA, k)$freq
    list(rows = do.call(rbind, out), agree = do.call(rbind, agree),
         freq = data.frame(item = seq_len(k), freq = fr200, loc_null = loc_null),
         n_distinct = length(unique(vapply(As_all, paste, character(1), collapse = ","))))
  }, error = function(e) conditionMessage(e))
  if (is.character(res)) return(list(status = paste("error:", res)))
  res$rows$cell <- res$agree$cell <- res$freq$cell <- cell
  res$rows$rep <- res$agree$rep <- res$freq$rep <- rep
  res$agree$n_distinct <- res$n_distinct
  res$status <- "ok"
  res
}

# --- Run ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cat(sprintf("Design: %d cells x %d replications on %d cores, B = %d, %d seeds per set\n",
            nrow(design), R, cores, B, S_max))
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
    agree = do.call(rbind, lapply(reps, `[[`, "agree")),
    freq = do.call(rbind, lapply(reps, `[[`, "freq")),
    status = st_rows
  )
  d <- design[cell, ]
  cat(sprintf(
    "[%s] cell %2d/%d  n=%d %-7s %s %-7s  ok %d/%d  %.1f min\n",
    format(Sys.time(), "%H:%M"), cell, nrow(design), d$n, d$ratio, d$structure,
    d$scenario, sum(st_rows$status == "ok"), R,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
}
per_rep <- do.call(rbind, lapply(all_res, `[[`, "per_rep"))
agree <- do.call(rbind, lapply(all_res, `[[`, "agree"))
freq <- do.call(rbind, lapply(all_res, `[[`, "freq"))
status <- do.call(rbind, lapply(all_res, `[[`, "status"))
status <- merge(status, design, by = "cell")

# --- Summary -----------------------------------------------------------------
# fwer      share of replications with any flag among units without DIF
# power     mean share of DIF units flagged (NA without DIF)
# untested  share of replications where a DIF unit was not tested
grp <- split(per_rep, per_rep[c("cell", "method", "family")], drop = TRUE)
summ <- do.call(rbind, lapply(grp, function(x) data.frame(
  x[1L, c("cell", "method", "family")], reps = nrow(x),
  fwer = mean(x$n_false > 0), mean_false = mean(x$n_false),
  power = if (x$n_alt[1L] > 0) mean(x$n_found / x$n_alt) else NA_real_,
  untested = if (x$n_alt[1L] > 0) mean(x$n_alt_tested < x$n_alt) else NA_real_,
  anchor_len = mean(x$anchor_len), anchor_locdif = mean(x$anchor_locdif > 0)
)))
summ <- merge(design, summ, by = "cell")
summ <- summ[order(summ$family, summ$scenario, summ$ratio, summ$structure,
                   summ$n, summ$method), ]

# Reproducibility: share of replications where seed sets A and B agree.
# S = -200 is cons_all (every item tested).
asum <- do.call(rbind, lapply(split(agree, agree[c("cell", "S")], drop = TRUE),
  function(x) data.frame(
    cell = x$cell[1L], S = x$S[1L], reps = nrow(x),
    same_set = mean(x$same_set), same_ref = mean(x$same_ref),
    loc_same = mean(x$loc_same), thr_same = mean(x$thr_same),
    found_diff_loc = mean(x$found_diff_loc), found_diff_thr = mean(x$found_diff_thr),
    any_near = mean(x$n_near > 0), n_distinct = mean(x$n_distinct)
  )))
asum <- merge(design, asum, by = "cell")
asum <- asum[order(asum$S, asum$scenario, asum$ratio, asum$structure, asum$n), ]
fails <- as.data.frame(table(cell = status$cell, status = status$status))
fails <- fails[fails$Freq > 0 & fails$status != "ok", ]

out <- list(
  design = design, per_rep = per_rep, agree = agree, freq = freq,
  status = status, summary = summ, agree_summary = asum, failures = fails,
  R = R, B = B, S_max = S_max, S_grid = S_grid, alpha = alpha, size = size,
  session = utils::sessionInfo(),
  run_minutes = as.numeric(difftime(Sys.time(), t_start, units = "mins"))
)
stamp <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(out, file.path(out_dir, paste0("dif_consensus_sim_", stamp, ".rds")))
saveRDS(out, file.path(out_dir, "dif_consensus_sim.rds"))
utils::write.csv(summ, file.path(out_dir, "dif_consensus_sim_summary.csv"), row.names = FALSE)
utils::write.csv(asum, file.path(out_dir, "dif_consensus_sim_agree.csv"), row.names = FALSE)

options(width = 200)
cat("\nFailures (cell, status, count):\n")
print(fails, row.names = FALSE)
cat("\nFWER and power by method (single-step), mean over cells per scenario:\n")
print(aggregate(cbind(fwer, power) ~ family + scenario + method, data = summ,
                FUN = function(x) round(mean(x), 3), na.action = na.pass))
cat("\nReproducibility (seed set A vs B), mean over cells:\n")
print(aggregate(cbind(same_set, same_ref, loc_same, thr_same) ~ S, data = asum,
                FUN = function(x) round(mean(x), 3)))
cat(sprintf("\nTotal run time: %.1f min\n", out$run_minutes))
