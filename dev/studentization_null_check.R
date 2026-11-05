# Synthetic exchangeable null check for joint studentization in
# .bootstrap_pvalues() (R/utils-multiplicity.R).
#
# Each replication draws B + 1 iid rows of k statistics. Row 1 plays the
# observed statistic, rows 2..B+1 the bootstrap null, so observed and simulated
# are exchangeable by construction and any valid Monte Carlo test must reject
# at most floor(alpha * (B + 1) - eps) / (B + 1) of the time.
#
# Two studentizations are compared on the same draws:
#   standard: mean/SD of the B simulated values only (current package code)
#   joint:    mean/SD of all B + 1 values, observed included (proposed)
# for both tails ("two.sided", "upper"). Reported per cell: familywise error
# of the Westfall-Young step-down (any padj < alpha) and the per-comparison
# rejection rate of the marginal p-value (p < alpha).
#
# Under the complete null the step-down rejects anything only if its first
# step rejects, and the first-step adjusted p is the smallest padj, so FWER is
# computed from the first step alone. The assertion block checks that shortcut
# against the package helper before the main run.
#
# Columns get different locations and scales, and the skewed arm uses
# chi-square(4) noise, so that studentization has something to do.

library(parallel)

ALPHA <- 0.05
R <- 10000
CORES <- 10
KS <- c(5, 10, 20)
BS <- c(100, 200, 400, 1000, 2000)
DISTS <- c("normal", "skewed")

here <- "dev"
stopifnot(file.exists(file.path(here, "studentization_null_check.R")))
helpers <- new.env()
sys.source("R/utils-multiplicity.R", envir = helpers)
stopifnot(is.function(helpers$.bootstrap_pvalues), is.function(helpers$.wy_stepdown))

attainable <- function(B) floor(ALPHA * (B + 1) - 1e-9) / (B + 1)

draw <- function(n, k, dist) {
  z <- if (dist == "normal") {
    matrix(stats::rnorm(n * k), n, k)
  } else {
    (matrix(stats::rchisq(n * k, df = 4), n, k) - 4) / sqrt(8)
  }
  mu <- seq(-1, 1, length.out = k)
  sg <- exp(seq(-1, 1, length.out = k))
  sweep(sweep(z, 2L, sg, `*`), 2L, mu, `+`)
}

row_max <- function(x) do.call(pmax, lapply(seq_len(ncol(x)), function(j) x[, j]))

# Candidate package code: the joint variant of .bootstrap_pvalues(). Only the
# moments change, everything downstream is the current helper.
bootstrap_pvalues_joint <- function(observed, sim, correction = "fwer",
                                    tail = c("two.sided", "upper")) {
  tail <- match.arg(tail)
  nm <- names(observed)
  sim <- sim[, nm, drop = FALSE]
  B <- nrow(sim)
  all <- rbind(observed, sim)
  m <- colMeans(all, na.rm = TRUE)
  s <- apply(all, 2L, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- NA_real_
  t_obs <- (observed - m) / s
  t_sim <- sweep(sweep(sim, 2L, m, `-`), 2L, s, `/`)
  stat_obs <- if (tail == "two.sided") abs(t_obs) else t_obs
  stat_sim <- if (tail == "two.sided") abs(t_sim) else t_sim
  p <- vapply(seq_along(nm), function(i) {
    (1 + sum(stat_sim[, i] >= stat_obs[i], na.rm = TRUE)) / (B + 1)
  }, numeric(1))
  data.frame(name = nm, p = p, padj = helpers$.wy_stepdown(stat_obs, stat_sim, B))
}

# Fast per-replication statistics for one variant and one tail.
fast_stats <- function(obs, sim, joint, tail) {
  B <- nrow(sim)
  src <- if (joint) rbind(obs, sim) else sim
  m <- colMeans(src)
  s <- sqrt(colSums(sweep(src, 2L, m)^2) / (nrow(src) - 1))
  t_obs <- (obs - m) / s
  t_sim <- sweep(sweep(sim, 2L, m), 2L, s, `/`)
  if (tail == "two.sided") {
    t_obs <- abs(t_obs)
    t_sim <- abs(t_sim)
  }
  p_marg <- (1 + colSums(sweep(t_sim, 2L, t_obs, `>=`))) / (B + 1)
  p_first <- (1 + sum(row_max(t_sim) >= max(t_obs))) / (B + 1)
  c(fwer = p_first < ALPHA, marg = mean(p_marg < ALPHA), p_first = p_first,
    p_marg1 = unname(p_marg[1]))
}

one_rep <- function(k, B, dist) {
  x <- draw(B + 1, k, dist)
  obs <- x[1, ]
  sim <- x[-1, , drop = FALSE]
  out <- c()
  for (tail in c("two.sided", "upper")) {
    st <- fast_stats(obs, sim, joint = FALSE, tail = tail)
    jt <- fast_stats(obs, sim, joint = TRUE, tail = tail)
    tag <- if (tail == "two.sided") "two" else "up"
    out[paste0(tag, "_fwer_std")] <- st[["fwer"]]
    out[paste0(tag, "_fwer_joint")] <- jt[["fwer"]]
    out[paste0(tag, "_marg_std")] <- st[["marg"]]
    out[paste0(tag, "_marg_joint")] <- jt[["marg"]]
  }
  out
}

# ---- assertions: fast path == package helper / candidate --------------------
set.seed(20261023)
for (a in seq_len(40)) {
  k <- sample(KS, 1)
  B <- sample(c(100, 200), 1)
  dist <- sample(DISTS, 1)
  x <- draw(B + 1, k, dist)
  colnames(x) <- paste0("I", seq_len(k))
  obs <- x[1, ]
  sim <- x[-1, , drop = FALSE]
  for (tail in c("two.sided", "upper")) {
    pk <- helpers$.bootstrap_pvalues(obs, sim, correction = "fwer", tail = tail)
    pj <- bootstrap_pvalues_joint(obs, sim, correction = "fwer", tail = tail)
    fs <- fast_stats(obs, sim, joint = FALSE, tail = tail)
    fj <- fast_stats(obs, sim, joint = TRUE, tail = tail)
    stopifnot(
      isTRUE(all.equal(min(pk$padj), fs[["p_first"]])),
      isTRUE(all.equal(min(pj$padj), fj[["p_first"]])),
      isTRUE(all.equal(pk$p[1], fs[["p_marg1"]])),
      isTRUE(all.equal(pj$p[1], fj[["p_marg1"]]))
    )
    # Upper tail: the marginal p is a rank of the raw values within a column,
    # so the two studentizations must give identical marginal p-values.
    if (tail == "upper") stopifnot(identical(pk$p, pj$p))
  }
}
message("Assertions passed: fast path matches the package helper and the candidate.")

# ---- main run ----------------------------------------------------------------
grid <- expand.grid(k = KS, B = BS, dist = DISTS, stringsAsFactors = FALSE)
RNGkind("L'Ecuyer-CMRG")
set.seed(1023)
t0 <- Sys.time()
res <- do.call(rbind, lapply(seq_len(nrow(grid)), function(g) {
  k <- grid$k[g]
  B <- grid$B[g]
  dist <- grid$dist[g]
  reps <- mclapply(seq_len(R), function(r) one_rep(k, B, dist),
                   mc.cores = CORES, mc.set.seed = TRUE)
  m <- do.call(rbind, reps)
  stopifnot(nrow(m) == R, all(is.finite(m)))
  data.frame(
    k = k, B = B, dist = dist, R = R, attainable = 100 * attainable(B),
    two_fwer_std = 100 * mean(m[, "two_fwer_std"]),
    two_fwer_joint = 100 * mean(m[, "two_fwer_joint"]),
    two_fwer_diff_se = 100 * sd(m[, "two_fwer_std"] - m[, "two_fwer_joint"]) / sqrt(R),
    up_fwer_std = 100 * mean(m[, "up_fwer_std"]),
    up_fwer_joint = 100 * mean(m[, "up_fwer_joint"]),
    up_fwer_diff_se = 100 * sd(m[, "up_fwer_std"] - m[, "up_fwer_joint"]) / sqrt(R),
    two_marg_std = 100 * mean(m[, "two_marg_std"]),
    two_marg_joint = 100 * mean(m[, "two_marg_joint"]),
    up_marg_std = 100 * mean(m[, "up_marg_std"]),
    up_marg_joint = 100 * mean(m[, "up_marg_joint"])
  )
}))
message("Run time: ", format(Sys.time() - t0))

saveRDS(res, file.path(here, "studentization_null_check.rds"))
saveRDS(res, file.path(here, paste0("studentization_null_check_", Sys.Date(), ".rds")))
print(res, digits = 3)
