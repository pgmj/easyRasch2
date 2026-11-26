# Check A2: code checks for the generic band code and the adjacent and
# residual layouts. Run from the package root:
#   LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/check_A2.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/catfit_dif/R/perm_functions.R")

ok <- function(cond, label) {
  cat(sprintf("%-4s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  invisible(cond)
}

set.seed(202)
d <- gen_null(n = 400, k = 6, impact = 1, props = c(0.7, 0.3))
pr <- prep_dif(d$X, d$grp)

# 1. Generic code on the calibration layout reproduces perm_bands().
cs <- calib_setup(pr)
set.seed(5); b_old <- perm_bands(pr, B = 199L)
set.seed(5); b_new <- perm_bands_generic(pr, cs$lay, cs$fun, B = 199L)
ok(isTRUE(all.equal(b_old$curves$p_sim, b_new$curves$p_sim)) &&
     isTRUE(all.equal(b_old$points$p_point, b_new$points$p_point)),
   "generic code reproduces perm_bands() for calibration")

# 2. Adjacent: observed shares per group equal the package's
#    .catplot_adjacent() run on that group's respondents with pooled P.
as <- adjacent_setup(pr)
P <- .cond_cat_probs(pr$thr)
adj_ok <- vapply(seq_len(pr$G), function(g) {
  rows <- pr$grp == g
  ref <- .catplot_adjacent(pr$X[rows, , drop = FALSE], pr$r[rows], P, 0.05)
  mine <- as$lay[as$lay$group == g, ]
  mine$obs <- as$fun(pr$grp)[as$lay$group == g]
  ref$pair_j <- as.integer(sub(" vs.*", "", ref$pair))
  mm <- merge(ref, mine, by.x = c("item", "pair_j", "score"),
              by.y = c("item", "pair", "score"))
  nrow(mm) == nrow(ref) && nrow(mm) == nrow(mine) &&
    isTRUE(all.equal(mm$observed, mm$obs))
}, logical(1))
ok(all(adj_ok), "adjacent shares match .catplot_adjacent() per group")

# 3. Residuals with one group equal the package's .catplot_residuals().
pr1 <- pr
pr1$grp <- rep(1L, length(pr$grp))
pr1$G <- 1L
rs <- residual_setup(pr1, class_intervals = 5L)
bin <- .cicc_class_bins(pr1$r, "quantile", 5L, NULL)
sg <- .score_group_labels(pr1$r, bin)
ref <- .catplot_residuals(pr1$X, pr1$r, P, sg)
ref <- ref[!is.na(ref$z), ]
mine <- rs$lay
mine$z <- rs$fun(pr1$grp)
mine$score_group <- levels(sg)[mine$q]
mm <- merge(ref, mine, by = c("item", "score_group", "category"))
ok(nrow(mm) == nrow(ref) && isTRUE(all.equal(mm$z.x, mm$z.y)),
   "one-group residuals match .catplot_residuals()")

# 4. Residuals per group: E summed over groups equals the pooled E.
rs2 <- residual_setup(pr, 5L)
Eg <- stats::aggregate(E ~ item + q + category, data = rs2$lay, FUN = sum)
E1 <- rs$lay[, c("item", "q", "category", "E")]
mm <- merge(Eg, E1, by = c("item", "q", "category"))
ok(isTRUE(all.equal(mm$E.x, mm$E.y)), "group E values sum to pooled E")

# 5. Constant-within-score label: every band has zero width (p_sim = 1 by
#    convention of ties, all points fixed).
prc <- pr
prc$grp <- as.integer(pr$r %% 2L) + 1L
for (nm in c("adjacent", "residuals")) {
  st <- if (nm == "adjacent") adjacent_setup(prc) else residual_setup(prc, 5L)
  bb <- perm_bands_generic(prc, st$lay, st$fun, B = 49L)
  ok(!any(bb$points$varies) && all(bb$curves$p_sim == 1),
     paste0(nm, ": label constant within scores gives fixed points"))
}

# 6. Timing at the largest smoke cell.
set.seed(3)
d3 <- gen_null(500, 10, 0.5, rep(1 / 3, 3))
pr3 <- prep_dif(d3$X, d3$grp)
for (nm in c("adjacent", "residuals")) {
  st <- if (nm == "adjacent") adjacent_setup(pr3) else residual_setup(pr3, 5L)
  tm <- system.time(bb <- perm_bands_generic(pr3, st$lay, st$fun, B = 999L))
  cat(sprintf("     %s: %.1f s, %d curves, %.1f varying points per curve\n",
              nm, tm[["elapsed"]], nrow(bb$curves),
              mean(bb$curves$n_points_var)))
}
