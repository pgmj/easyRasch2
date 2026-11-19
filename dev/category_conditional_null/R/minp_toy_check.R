# Toy check of the single-step min-p combination of four events (2026-11-18).
# Observed and B replicates exchangeable (four correlated normal statistics,
# correlation .6), so a valid item p has P(p <= a) <= a. "old" is the first
# version in analyse_calibration() (replicate p on denominator B, observed on
# B + 1), "pooled" the current .minp_pooled().
#   Rscript R/minp_toy_check.R   (about 10 min)
# Output 2026-11-18: B = 199, old .0168 / .0564, pooled .0013 / .0447 at
# .01 / .05. B = 999 pooled .0088 / .0497.


set.seed(1); Bu <- 199; m <- 4; R <- 20000
S <- matrix(.6, m, m); diag(S) <- 1; L <- chol(S)
old <- new <- numeric(R)
for (r in 1:R) {
  T <- matrix(rnorm((Bu + 1) * m), Bu + 1) %*% L     # row 1 = observed
  obs <- T[1, ]; Tb <- t(T[-1, ])                    # m x Bu
  # old (as in analyse_calibration)
  p_raw <- (1 + rowSums(Tb >= obs)) / (1 + Bu)
  pb <- apply(Tb, 1, function(t) (1 + rowSums(outer(t, t, "<=")) - 1) / Bu)
  minp <- apply(matrix(pb, ncol = m), 1, min)
  old[r] <- min(vapply(p_raw, function(p) (1 + sum(minp <= p)) / (1 + Bu), 0))
  # pooled: observed and replicates ranked in one pool of Bu + 1
  P <- apply(T, 2, function(t) rowSums(outer(t, t, "<=")) / (Bu + 1))  # (Bu+1) x m
  mp <- apply(P, 1, min)
  new[r] <- mean(mp <= mp[1])
}
for (a in c(.01, .025, .05, .1)) cat(a, " old:", round(mean(old <= a + 1e-9), 4), " pooled:", round(mean(new <= a + 1e-9), 4), "\n")

set.seed(2); m <- 4
S <- matrix(.6, m, m); diag(S) <- 1; L <- chol(S)
pooled_p <- function(T) {           # T: (B+1) x m, row 1 observed; upper tail
  N <- nrow(T)
  P <- apply(T, 2, function(t) (N - rank(t, ties.method = "min") + 1) / N)
  mp <- apply(P, 1, min); mean(mp <= mp[1])
}
for (Bu in c(199, 999)) {
  R <- if (Bu == 199) 20000 else 10000
  p <- replicate(R, pooled_p(matrix(rnorm((Bu + 1) * m), Bu + 1) %*% L))
  cat("B =", Bu, " min attainable seen:", round(min(p), 4), " ")
  for (a in c(.01, .025, .05, .1)) cat(a, ":", round(mean(p <= a + 1e-9), 4), " ")
  cat("\n")
}
