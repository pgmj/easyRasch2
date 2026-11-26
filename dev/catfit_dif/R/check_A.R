# Check A: deterministic checks of the permutation-band code.
# Run from the package root:
#   LANG=en_US.UTF-8 Rscript dev/catfit_dif/R/check_A.R
suppressMessages(devtools::load_all(".", quiet = TRUE))
source("dev/catfit_dif/R/perm_functions.R")

ok <- function(cond, label) {
  cat(sprintf("%-4s %s\n", if (isTRUE(cond)) "PASS" else "FAIL", label))
  invisible(cond)
}

set.seed(101)
d <- gen_null(n = 400, k = 6, impact = 1, props = c(0.7, 0.3))
pr <- prep_dif(d$X, d$grp)
lay <- curve_layout(pr)
S <- length(pr$sc)

# 1. Predictions never decrease with the score, and ties occur only where
#    the event is certain (P = 0 or 1, at extreme scores). There the
#    outcome is fixed, so keeping tied scores as separate points cannot
#    change the PAV values.
tie_ok <- apply(pr$Pev, 2L, function(p) {
  d <- diff(p)
  all(d >= -1e-12) &&
    all(p[-1][abs(d) <= 1e-12] < 1e-12 | p[-1][abs(d) <= 1e-12] > 1 - 1e-12)
})
ok(all(tie_ok),
   "P(X >= k | r) non-decreasing in r, ties only at P = 0 or 1")
cert <- pr$Pev[cbind(pr$s_idx, rep(seq_len(ncol(pr$Y)), each = nrow(pr$Y)))]
cert <- matrix(cert, ncol = ncol(pr$Y))
ok(all(pr$Y[cert > 1 - 1e-12] == 1) && all(pr$Y[cert < 1e-12] == 0),
   "  ... and responses at certain events agree with them")

# 2. The fast aggregated PAV agrees with the package's .pav() on
#    individual data, for every group x event curve (matched on p).
fast <- curve_values(pr, pr$grp, lay)
ref <- unlist(lapply(unique(lay$curve), function(cv) {
  li <- lay[lay$curve == cv, ]
  rows <- pr$grp == li$group[1L]
  d <- .pav(pr$Pev[pr$s_idx[rows], li$event[1L]], pr$Y[rows, li$event[1L]])
  d$cep[match(li$predicted, d$p)]
}))
ok(isTRUE(all.equal(fast, ref, tolerance = 1e-12)),
   "aggregated PAV equals .pav() on individual data")

# 3. Every permutation keeps group counts per score, hence the layout.
n0 <- agg_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)$n
same <- vapply(1:200, function(b) {
  gp <- permute_within(pr$grp, pr$s_idx)
  identical(agg_counts(pr$Y, pr$s_idx, gp, pr$G, S)$n, n0) &&
    identical(curve_layout(modifyList(pr, list(grp = gp))), lay)
}, logical(1))
ok(all(same), "200 permutations keep counts per score and the curve layout")

# 4. Permutations do move labels (not the identity).
moved <- mean(vapply(1:50, function(b) {
  any(permute_within(pr$grp, pr$s_idx) != pr$grp)
}, logical(1)))
ok(moved == 1, "permutations change the labelling")

# 5. A DIF variable that is constant within every score: zero-width band.
pr_c <- pr
pr_c$grp <- as.integer(pr$r %% 2L) + 1L
bc <- perm_bands(pr_c, B = 99L, keep_perm = TRUE)
ok(all(abs(bc$perm - bc$points$observed) < 1e-12),
   "label constant within scores: every permutation equals the observed")
ok(all(bc$points$p_point == 1), "  ... and every pointwise p is 1")

# 6. Same seed, same result.
set.seed(7); b1 <- perm_bands(pr, B = 99L)
set.seed(7); b2 <- perm_bands(pr, B = 99L)
ok(identical(b1, b2), "same seed gives identical bands and p-values")

# 7. Two groups: raw per-score shares mirror around the pooled share,
#    n1 * (y1 - ybar) = -n2 * (y2 - ybar). (PAV is not linear, so the mirror
#    holds for raw shares, not for PAV curves.)
raw <- curve_raw(pr, pr$grp, lay)
a <- agg_counts(pr$Y, pr$s_idx, pr$grp, pr$G, S)
dev_ok <- vapply(seq_len(ncol(pr$Y)), function(e) {
  both <- which(a$n[1, ] > 0 & a$n[2, ] > 0)
  ybar <- (a$sums[both, e] + a$sums[S + both, e]) / (a$n[1, both] + a$n[2, both])
  y1 <- a$sums[both, e] / a$n[1, both]
  y2 <- a$sums[S + both, e] / a$n[2, both]
  isTRUE(all.equal(a$n[1, both] * (y1 - ybar), -a$n[2, both] * (y2 - ybar)))
}, logical(1))
ok(all(dev_ok), "two groups: weighted raw deviations from pooled mirror")

# 8. Permutation distribution of a group's success count at one score is
#    hypergeometric (mean and variance), checked at the busiest score.
s_star <- which.max(tabulate(pr$s_idx))
e_star <- 3L
N <- sum(pr$s_idx == s_star)
Ks <- sum(pr$Y[pr$s_idx == s_star, e_star])
n1 <- sum(pr$s_idx == s_star & pr$grp == 1L)
set.seed(11)
cnt <- vapply(1:20000, function(b) {
  gp <- permute_within(pr$grp, pr$s_idx)
  sum(pr$Y[pr$s_idx == s_star & gp == 1L, e_star])
}, numeric(1))
hm <- n1 * Ks / N
hv <- n1 * (Ks / N) * (1 - Ks / N) * (N - n1) / (N - 1)
cat(sprintf("     score %d: N = %d, successes = %d, n1 = %d\n",
            pr$sc[s_star], N, Ks, n1))
cat(sprintf("     mean %.3f vs %.3f, var %.3f vs %.3f\n",
            mean(cnt), hm, var(cnt), hv))
ok(abs(mean(cnt) - hm) < 4 * sqrt(hv / 20000) &&
     abs(var(cnt) / hv - 1) < 0.05,
   "within-score permutation gives hypergeometric counts")

# 9. Pointwise p on the grid of (B + 1), and simultaneous p in (0, 1].
b <- perm_bands(pr, B = 199L)
ok(all(b$points$p_point > 0 & b$points$p_point <= 1) &&
     all(b$curves$p_sim >= 1 / 200 & b$curves$p_sim <= 1),
   "p-values in range")
cat(sprintf("     %d curves, %d points; pointwise outside at .05: %.3f; ",
            nrow(b$curves), nrow(b$points), mean(b$points$p_point <= .05)),
    sprintf("curves any-point: %.3f; simultaneous: %.3f\n",
            mean(b$curves$min_p_point <= .05), mean(b$curves$p_sim <= .05)))
