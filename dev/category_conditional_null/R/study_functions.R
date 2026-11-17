# Study functions, built on R/functions.R.
# Purpose: does the PCM-conditional category-by-score table (Glas's R1c
# table) detect a response scale that does not work as the integer scoring
# says, without flagging threshold disorder that the PCM itself produces?

source("R/functions.R")

# ── Data generation ─────────────────────────────────────────────

# PCM with step-specific discriminations a (Glas & Verhelst 1989, sec. 7):
# P(X = j | theta) proportional to exp(sum_{l <= j} a_l (theta - delta_l)).
# a = 1 for every step is the PCM. a_l = 0 makes categories l-1 and l
# equally likely at every theta (interchangeable).
sim_steps <- function(theta, delta, a = rep(1, length(delta))) {
  m <- length(delta)
  lp <- sapply(0:m, function(j) if (j == 0) 0 * theta else
    rowSums(sweep(outer(theta, delta[1:j], "-"), 2, a[1:j], "*")))
  lp <- matrix(lp, nrow = length(theta))
  p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
  rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
}

# Fixed thresholds for the non-focus items: locations spread over
# [-1, 1], step pattern (-1.5, -0.5, 0.5, 1.5) with a fixed jitter.
# withr::with_seed restores the caller's RNG state afterwards. The first host
# run (2026-11-16) called set.seed() here directly, which reset the global
# stream inside gen_data() and made almost every replicate in a cell the
# same dataset (16 distinct of 500). Archived in results/cells_rngbug/.
other_items <- function(k) {
  withr::with_seed(20261116 + k, .rng_kind = "Mersenne-Twister", {
    loc <- seq(-1, 1, length.out = k - 1)
    lapply(loc, function(l) l + c(-1.5, -0.5, 0.5, 1.5) + stats::rnorm(4, 0, 0.2))
  })
}

WIDE   <- c(-1.5, -0.5, 0.5, 1.5)
NARROW <- c(-1.5, 0.4, -0.4, 1.5)  # middle category rare, thresholds disordered

# Condition table. Focus item is item 1. `type` picks the generator.
conditions <- function() {
  rbind(
    data.frame(cond = "A0_pcm_ordered",    type = "pcm",  null = TRUE),
    data.frame(cond = "A1_pcm_disordered", type = "pcm",  null = TRUE),
    data.frame(cond = "S_outer_a0",        type = "step", null = FALSE),
    data.frame(cond = "S_inner_a0",        type = "step", null = FALSE),
    data.frame(cond = "S_inner_a05",       type = "step", null = FALSE),
    data.frame(cond = "S_inner_a2",        type = "step", null = FALSE),
    data.frame(cond = "SWAP_1_2",          type = "swap", null = FALSE),
    data.frame(cond = "DK_wide_05",        type = "dk",   null = FALSE),
    data.frame(cond = "DK_wide_10",        type = "dk",   null = FALSE),
    data.frame(cond = "DK_wide_25",        type = "dk",   null = FALSE),
    data.frame(cond = "DK_narrow_05",      type = "dk",   null = FALSE),
    data.frame(cond = "DK_narrow_10",      type = "dk",   null = FALSE),
    data.frame(cond = "DK_narrow_25",      type = "dk",   null = FALSE),
    data.frame(cond = "ERS_30",            type = "ers",  null = FALSE)
  )
}

gen_data <- function(cond, n, k) {
  theta <- stats::rnorm(n)
  oth <- other_items(k)
  Xo <- sim_pcm(theta, oth)
  x1 <- switch(cond,
    A0_pcm_ordered    = sim_steps(theta, WIDE),
    A1_pcm_disordered = sim_steps(theta, NARROW),
    S_outer_a0  = sim_steps(theta, WIDE, c(0, 1, 1, 1)),
    S_inner_a0  = sim_steps(theta, WIDE, c(1, 0, 1, 1)),
    S_inner_a05 = sim_steps(theta, WIDE, c(1, 0.5, 1, 1)),
    S_inner_a2  = sim_steps(theta, WIDE, c(1, 2, 1, 1)),
    SWAP_1_2 = {  # categories scored 0-4: swaps the 2nd and 3rd options
      x <- sim_steps(theta, WIDE)
      x[x == 1] <- -1L; x[x == 2] <- 1L; x[x == -1L] <- 2L; x
    },
    {
      if (startsWith(cond, "DK_")) {
        base <- if (grepl("narrow", cond)) NARROW else WIDE
        pi_dk <- as.numeric(sub(".*_", "", cond)) / 100
        x <- sim_steps(theta, base)
        x[stats::runif(n) < pi_dk] <- 2L
        x
      } else if (cond == "ERS_30") {
        NULL
      } else stop("unknown condition ", cond)
    })
  if (cond == "ERS_30") {
    # 30% of respondents use every item's scale with thresholds compressed
    # towards the item location (factor 0.3): more extreme categories.
    ers <- stats::runif(n) < 0.30
    all_d <- c(list(WIDE), oth)
    X <- sapply(all_d, function(d) {
      dc <- mean(d) + 0.3 * (d - mean(d))
      x <- sim_steps(theta, d)
      x[ers] <- sim_steps(theta[ers], dc)
      x
    })
    colnames(X) <- paste0("I", seq_len(ncol(X)))
    return(X)
  }
  X <- cbind(x1, Xo)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  X
}

# ── Statistics ──────────────────────────────────────────────────

fit_cml <- function(X, hessian = FALSE) {
  m <- psychotools::pcmodel(X, hessian = hessian)
  tp <- psychotools::threshpar(m, type = "mode", vcov = hessian)
  list(delta = lapply(tp, unname), vcov = if (hessian) stats::vcov(tp) else NULL)
}

# Glas item index: sum of squared deviates over categories j >= 1 (category
# 0 is redundant given the others, Glas & Verhelst 1989).
.item_index <- function(s) sum(s$Q[s$category >= 1])

all_item_index <- function(X, delta) {
  vapply(seq_along(delta), function(i) .item_index(item_stats(X, delta, i)), 0)
}

# Null 1: AUC of the rest score between adjacent categories of item i, with
# a one-sided Wilcoxon p for AUC > .5.
adjacent_auc <- function(X, i) {
  rest <- rowSums(X[, -i, drop = FALSE])
  x <- X[, i]
  m <- max(x)
  t(vapply(1:m, function(j) {
    a <- rest[x == j - 1]; b <- rest[x == j]
    if (length(a) < 2 || length(b) < 2) return(c(auc = NA, p = NA))
    w <- stats::wilcox.test(b, a, alternative = "greater", exact = FALSE)
    c(auc = unname(w$statistic) / (length(a) * length(b)), p = w$p.value)
  }, c(auc = 0, p = 0)))
}

# Wald check of threshold order for item i: z for each adjacent difference.
threshold_z <- function(fit, i) {
  d <- fit$delta[[i]]
  nm <- paste0("I", i, "-C", seq_along(d))
  V <- fit$vcov[nm, nm]
  vapply(seq_len(length(d) - 1), function(k) {
    g <- numeric(length(d)); g[k] <- -1; g[k + 1] <- 1
    (d[k + 1] - d[k]) / sqrt(drop(t(g) %*% V %*% g))
  }, 0)
}

# One dataset: observed statistics, then a bootstrap that samples patterns
# exactly given the observed totals and refits the CML thresholds.
analyse_dataset <- function(X, B) {
  k <- ncol(X); m <- 4L
  if (any(vapply(seq_len(k), function(i) length(unique(X[, i])) < m + 1, TRUE)))
    return(list(skipped = "empty category in data"))
  fit <- tryCatch(fit_cml(X, hessian = TRUE), error = function(e) NULL)
  if (is.null(fit)) return(list(skipped = "fit failed"))
  K <- k * m
  r <- rowSums(X)
  totals <- r[r > 0 & r < K]
  s1 <- item_stats(X, fit$delta, 1)
  obs_idx <- all_item_index(X, fit$delta)
  obs_cat <- .boot_vec(s1)
  boot <- replicate(B, {
    Xb <- sample_given_totals(fit$delta, totals)
    fb <- tryCatch(fit_cml(Xb), error = function(e) NULL)
    if (is.null(fb) || any(lengths(fb$delta) != m)) {
      rep(NA_real_, k + length(obs_cat))
    } else c(all_item_index(Xb, fb$delta), .boot_vec(item_stats(Xb, fb$delta, 1)))
  })
  boot <- t(boot)
  ok <- stats::complete.cases(boot)
  boot <- boot[ok, , drop = FALSE]
  bi <- boot[, seq_len(k), drop = FALSE]
  bc <- boot[, -seq_len(k), drop = FALSE]
  Bu <- nrow(boot)
  p_idx <- vapply(seq_len(k), function(i) (1 + sum(bi[, i] >= obs_idx[i])) / (1 + Bu), 0)
  is_q <- startsWith(names(obs_cat), "Q_")
  p_cat <- vapply(seq_along(obs_cat), function(j) {
    tb <- bc[, j]
    if (is_q[j]) (1 + sum(tb >= obs_cat[j])) / (1 + Bu)
    else min(1, 2 * min(1 + sum(tb >= obs_cat[j]), 1 + sum(tb <= obs_cat[j])) / (1 + Bu))
  }, 0)
  names(p_cat) <- names(obs_cat)
  list(
    skipped = NA_character_, B_used = Bu,
    thresholds1 = fit$delta[[1]], z1 = threshold_z(fit, 1),
    p_item = p_idx, p_item_holm = stats::p.adjust(p_idx, "holm"),
    stat_cat = obs_cat, p_cat = p_cat,
    auc1 = adjacent_auc(X, 1)
  )
}
