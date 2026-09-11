# ---------------------------------------------------------------------------
# Does the ratio form of marginal reliability really track observed-score
# reliability better than the subtractive (Green/Lord) form?
#
# Part A replicates the Table 1 design of Milanzi, Molenberghs, Alonso,
# Verbeke & De Boeck (2015, BJMSP 68, 43-64) and checks the recomputation that
# was done by inverting their published rho_f values.
#
# Part B repeats the comparison for polytomous Rasch (PCM) data at scale
# lengths easyRasch2 actually sees, using the package's own information
# engine.
#
# No models are fitted anywhere: as in their Section 5.1.1, the generated
# sample is treated as the population and true values are plugged in.
# ---------------------------------------------------------------------------

PKG <- "/Users/magnus.johansson.3/Library/CloudStorage/OneDrive-KarolinskaInstitutet/Claude/easyRasch2"
suppressMessages(devtools::load_all(PKG, quiet = TRUE))

set.seed(20260911)
N_REP_A <- 1000L
N_REP_B <- 500L
N_PERSON <- 400L
N_ITEM <- 24L

# ---------------------------------------------------------------------------
# Part A: binary 1PL / 2PL, Milanzi's design
# ---------------------------------------------------------------------------

#' All six coefficients for one replication of the binary design
#'
#' P is n x I of response probabilities, a the discriminations, theta the
#' person locations, b the difficulties.
binary_coefs <- function(theta, a, b) {
  eta <- outer(theta, b, "-") * rep(a, each = length(theta))
  P <- stats::plogis(eta)
  I_items <- length(a)

  # --- exact expected-sum-score reliability (their Section 3.1) -------------
  mu <- rowSums(P)                       # E[S | theta]
  var_mu <- mean(mu^2) - mean(mu)^2      # Var(mu)
  var_eps <- sum(colMeans(P * (1 - P)))  # sum_i E_theta[P(1-P)]
  rho_S <- var_mu / (var_mu + var_eps)

  # --- Fisher information coefficients --------------------------------------
  info <- as.numeric((P * (1 - P)) %*% (a^2))   # I(theta_j)
  m <- mean(1 / info)                           # E[1 / I(theta)]
  s2 <- mean(theta^2) - mean(theta)^2           # observed variance of theta

  rho_f <- 1 - m / s2                           # their eq. (12), subtractive
  rho_ratio <- s2 / (s2 + m)                    # marginal ratio form
  rho_ratio_pw <- mean(s2 / (s2 + 1 / info))    # average of the ratio CURVE

  # --- latent correlation (their eq. 30 generalised) -------------------------
  sum_a <- sum(a)
  rho_Sl <- s2 * sum_a^2 / (I_items * pi^2 / 3 + s2 * sum_a^2)

  # --- Taylor series approximation (their eq. 30 generalised) ---------------
  # v_i is the Bernoulli variance at the random-effects mean, theta = 0.
  p0 <- stats::plogis(-a * b)
  v <- p0 * (1 - p0)
  num <- s2 * sum(v * a)^2
  rho_SA <- num / (sum(v) + num)

  c(
    exact = rho_S,
    subtractive = rho_f,
    ratio = rho_ratio,
    ratio_pointwise = rho_ratio_pw,
    latent = rho_Sl,
    taylor = rho_SA,
    # identity check on the recomputation: ratio should equal 1 / (2 - rho_f)
    ratio_from_rho_f = 1 / (2 - rho_f)
  )
}

run_binary_cell <- function(model, sigma2, n_rep) {
  out <- vapply(seq_len(n_rep), function(i) {
    b <- stats::runif(N_ITEM, -4, 4)
    a <- if (model == "1PL") rep(1, N_ITEM) else stats::rnorm(N_ITEM, 2, 0.8)
    theta <- stats::rnorm(N_PERSON, 0, sqrt(sigma2))
    binary_coefs(theta, a, b)
  }, numeric(7L))
  t(out)
}

cat("=== Part A: replication of Milanzi et al. (2015) Table 1 ===\n")
cat(sprintf("%d replications per cell, %d persons, %d items\n\n",
            N_REP_A, N_PERSON, N_ITEM))

cells <- expand.grid(
  model = c("1PL", "2PL"),
  sigma2 = c(0.25, 1, 4),
  stringsAsFactors = FALSE
)
cells <- cells[order(cells$model, cells$sigma2), ]

# Published Table 1 values, for side-by-side comparison
published <- data.frame(
  model = c("1PL", "1PL", "1PL", "2PL", "2PL", "2PL"),
  sigma2 = c(0.25, 1, 4, 0.25, 1, 4),
  p_taylor = c(.461, .774, .932, .687, .898, .972),
  p_subtractive = c(-.120, .676, .887, .615, .881, .891),
  p_latent = c(.646, .879, .967, .889, .970, .992),
  p_exact = c(.480, .778, .928, .704, .895, .963),
  stringsAsFactors = FALSE
)

partA <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  res <- run_binary_cell(cells$model[i], cells$sigma2[i], N_REP_A)
  data.frame(
    model = cells$model[i],
    sigma2 = cells$sigma2[i],
    exact = mean(res[, "exact"]),
    subtractive = mean(res[, "subtractive"]),
    ratio = mean(res[, "ratio"]),
    ratio_pointwise = mean(res[, "ratio_pointwise"]),
    latent = mean(res[, "latent"]),
    taylor = mean(res[, "taylor"]),
    ae_subtractive = mean(abs(res[, "subtractive"] - res[, "exact"])),
    ae_ratio = mean(abs(res[, "ratio"] - res[, "exact"])),
    ae_ratio_pw = mean(abs(res[, "ratio_pointwise"] - res[, "exact"])),
    ae_taylor = mean(abs(res[, "taylor"] - res[, "exact"])),
    pct_negative = 100 * mean(res[, "subtractive"] < 0),
    sd_exact = stats::sd(res[, "exact"]),
    sd_subtractive = stats::sd(res[, "subtractive"]),
    max_identity_gap = max(abs(res[, "ratio"] - res[, "ratio_from_rho_f"])),
    exact_lo = stats::quantile(res[, "exact"], 0.025, names = FALSE),
    exact_hi = stats::quantile(res[, "exact"], 0.975, names = FALSE),
    subtr_lo = stats::quantile(res[, "subtractive"], 0.025, names = FALSE),
    subtr_hi = stats::quantile(res[, "subtractive"], 0.975, names = FALSE),
    stringsAsFactors = FALSE
  )
}))

cat("-- Simulated means against the published Table 1 --\n")
cmp <- merge(partA, published, by = c("model", "sigma2"))
cmp <- cmp[order(cmp$model, cmp$sigma2), ]
print(format(
  data.frame(
    model = cmp$model, sigma2 = cmp$sigma2,
    exact = round(cmp$exact, 3), pub_exact = cmp$p_exact,
    subtr = round(cmp$subtractive, 3), pub_subtr = cmp$p_subtractive,
    latent = round(cmp$latent, 3), pub_latent = cmp$p_latent,
    taylor = round(cmp$taylor, 3), pub_taylor = cmp$p_taylor
  ),
  nsmall = 3
), row.names = FALSE)

cat("\n-- The comparison under test --\n")
print(format(
  data.frame(
    model = partA$model, sigma2 = partA$sigma2,
    exact = round(partA$exact, 3),
    subtractive = round(partA$subtractive, 3),
    ratio = round(partA$ratio, 3),
    ratio_pw = round(partA$ratio_pointwise, 3),
    taylor = round(partA$taylor, 3)
  ),
  nsmall = 3
), row.names = FALSE)

cat("\n-- Mean absolute error against the exact sum-score reliability --\n")
print(format(
  data.frame(
    model = partA$model, sigma2 = partA$sigma2,
    subtractive = round(partA$ae_subtractive, 3),
    ratio = round(partA$ae_ratio, 3),
    ratio_pw = round(partA$ae_ratio_pw, 3),
    taylor = round(partA$ae_taylor, 3),
    pct_neg = round(partA$pct_negative, 1)
  ),
  nsmall = 3
), row.names = FALSE)

cat("\n-- Is the published single draw inside the replication spread? --\n")
print(format(
  data.frame(
    model = cmp$model, sigma2 = cmp$sigma2,
    pub_exact = cmp$p_exact,
    exact_95 = paste0("[", round(cmp$exact_lo, 3), ", ",
                      round(cmp$exact_hi, 3), "]"),
    exact_in = cmp$p_exact >= cmp$exact_lo & cmp$p_exact <= cmp$exact_hi,
    pub_subtr = cmp$p_subtractive,
    subtr_95 = paste0("[", round(cmp$subtr_lo, 3), ", ",
                      round(cmp$subtr_hi, 3), "]"),
    subtr_in = cmp$p_subtractive >= cmp$subtr_lo &
      cmp$p_subtractive <= cmp$subtr_hi
  )
), row.names = FALSE)

cat("\nOverall MAE  subtractive:", round(mean(partA$ae_subtractive), 4),
    " ratio:", round(mean(partA$ae_ratio), 4),
    " ratio_pointwise:", round(mean(partA$ae_ratio_pw), 4),
    " taylor:", round(mean(partA$ae_taylor), 4), "\n")
cat("Largest gap between the directly computed ratio form and 1/(2 - rho_f):",
    format(max(partA$max_identity_gap), scientific = TRUE), "\n")

# ---------------------------------------------------------------------------
# Part A2: one fixed item set per model, reused across the three sigma^2 cells
#
# Their Table 1 is a single realization (Tables 2 and 4 show the same 1PL
# difficulties at sigma^2 = 1 and 0.25, so one draw served all three cells).
# Redrawing items every replication, as Part A does, therefore compares
# against an average rather than against their draw. This block matches their
# design: draw the items once, replicate over persons only. The point is to
# check that the ORDERING of the methods does not depend on the item draw.
# ---------------------------------------------------------------------------

cat("\n=== Part A2: item parameters held fixed, replication over persons ===\n")

n_draw <- 200L
cell_lab <- paste0(rep(c("1PL", "2PL"), each = 3L), " sigma2=",
                   rep(c(0.25, 1, 4), times = 2L))

wins <- vapply(seq_len(n_draw), function(d) {
  b1 <- stats::runif(N_ITEM, -4, 4)
  b2 <- stats::runif(N_ITEM, -4, 4)
  a2 <- stats::rnorm(N_ITEM, 2, 0.8)
  res <- vapply(c(0.25, 1, 4), function(s2) {
    r1 <- rowMeans(vapply(seq_len(25L), function(r) {
      binary_coefs(stats::rnorm(N_PERSON, 0, sqrt(s2)), rep(1, N_ITEM), b1)
    }, numeric(7L)))
    r2 <- rowMeans(vapply(seq_len(25L), function(r) {
      binary_coefs(stats::rnorm(N_PERSON, 0, sqrt(s2)), a2, b2)
    }, numeric(7L)))
    c(
      abs(r1["ratio"] - r1["exact"]) < abs(r1["subtractive"] - r1["exact"]),
      abs(r2["ratio"] - r2["exact"]) < abs(r2["subtractive"] - r2["exact"]),
      abs(r1["subtractive"] - r1["exact"]), abs(r1["ratio"] - r1["exact"]),
      abs(r2["subtractive"] - r2["exact"]), abs(r2["ratio"] - r2["exact"])
    )
  }, numeric(6L))
  c(res[1L, ], res[2L, ],
    mean(c(res[3L, ], res[5L, ])),
    mean(c(res[4L, ], res[6L, ])))
}, numeric(8L))

cat("-- Per-cell: how often is the ratio form closer to exact? --\n")
print(data.frame(
  cell = cell_lab,
  ratio_closer_pct = round(100 * rowMeans(wins[1:6, , drop = FALSE]), 1),
  row.names = NULL
))
cat(sprintf(
  "\nMean MAE across %d item draws  subtractive: %.4f  ratio: %.4f\n",
  n_draw, mean(wins[7L, ]), mean(wins[8L, ])
))
cat("The ratio form does not win everywhere. Where both are well behaved the\n",
    "two are close and either can edge ahead. The overall gap comes from the\n",
    "low-information cells, where the subtractive form breaks down entirely.\n",
    sep = "")

# --- Sensitivity: sample variance of theta vs the true sigma^2 --------------
cat("\n-- Sensitivity: sigma^2 taken as the true value rather than Var(theta) --\n")
sens <- do.call(rbind, lapply(c(0.25, 1, 4), function(s2) {
  res <- t(vapply(seq_len(300L), function(r) {
    b <- stats::runif(N_ITEM, -4, 4)
    theta <- stats::rnorm(N_PERSON, 0, sqrt(s2))
    P <- stats::plogis(outer(theta, b, "-"))
    mu <- rowSums(P)
    rho_S <- (mean(mu^2) - mean(mu)^2) /
      ((mean(mu^2) - mean(mu)^2) + sum(colMeans(P * (1 - P))))
    info <- rowSums(P * (1 - P))
    m <- mean(1 / info)
    c(
      exact = rho_S,
      sub_samp = 1 - m / (mean(theta^2) - mean(theta)^2),
      sub_true = 1 - m / s2,
      ratio_samp = (mean(theta^2) - mean(theta)^2) /
        ((mean(theta^2) - mean(theta)^2) + m),
      ratio_true = s2 / (s2 + m)
    )
  }, numeric(5L)))
  data.frame(
    sigma2 = s2,
    exact = round(mean(res[, "exact"]), 3),
    sub_sample_var = round(mean(res[, "sub_samp"]), 3),
    sub_true_var = round(mean(res[, "sub_true"]), 3),
    ratio_sample_var = round(mean(res[, "ratio_samp"]), 3),
    ratio_true_var = round(mean(res[, "ratio_true"]), 3)
  )
}))
print(sens, row.names = FALSE)

# ---------------------------------------------------------------------------
# Part B: polytomous Rasch (PCM), easyRasch2's actual territory
# ---------------------------------------------------------------------------

#' Vectorised PCM category probabilities and score moments over a theta grid
#'
#' Returns E[score | theta] and Var(score | theta) per item, as n x I matrices.
pcm_moments <- function(theta, thr_list) {
  n <- length(theta)
  E <- V <- matrix(0, n, length(thr_list))
  for (i in seq_along(thr_list)) {
    thr <- thr_list[[i]]
    cats <- 0:length(thr)
    # cumulative sums of (theta - tau) across steps, category 0 = 0
    cum <- cbind(0, t(apply(outer(theta, thr, "-"), 1L, cumsum)))
    if (n == 1L) cum <- matrix(cum, nrow = 1L)
    mx <- apply(cum, 1L, max)
    w <- exp(cum - mx)
    P <- w / rowSums(w)
    e <- as.numeric(P %*% cats)
    E[, i] <- e
    V[, i] <- as.numeric(P %*% (cats^2)) - e^2
  }
  list(E = E, V = V)
}

# Confirm the vectorised version agrees with the package's own engine before
# using it in the loop.
chk_thr <- list(c(-1, 0, 1), c(-0.5, 0.5, 1.5), c(0, 0.3, 0.9))
chk_theta <- seq(-3, 3, length.out = 25L)
chk_fast <- rowSums(pcm_moments(chk_theta, chk_thr)$V)
chk_pkg <- easyRasch2:::.test_information(chk_thr, chk_theta)
cat("\n=== Part B: polytomous Rasch (PCM) ===\n")
cat("Vectorised PCM information matches .test_information(): ",
    isTRUE(all.equal(chk_fast, chk_pkg, tolerance = 1e-12)),
    " (max abs diff ", format(max(abs(chk_fast - chk_pkg)), scientific = TRUE),
    ")\n", sep = "")

pcm_coefs <- function(theta, thr_list) {
  mom <- pcm_moments(theta, thr_list)
  mu <- rowSums(mom$E)
  var_mu <- mean(mu^2) - mean(mu)^2
  var_eps <- sum(colMeans(mom$V))
  rho_S <- var_mu / (var_mu + var_eps)

  info <- rowSums(mom$V)           # for a Rasch family, I(theta) = Var(S|theta)
  m <- mean(1 / info)
  s2 <- mean(theta^2) - mean(theta)^2

  c(
    exact = rho_S,
    subtractive = 1 - m / s2,
    ratio = s2 / (s2 + m),
    ratio_pointwise = mean(s2 / (s2 + 1 / info))
  )
}

pcm_cells <- expand.grid(
  n_items = c(5L, 10L, 20L),
  sigma = c(0.5, 1.0, 1.5),
  stringsAsFactors = FALSE
)

partB <- do.call(rbind, lapply(seq_len(nrow(pcm_cells)), function(i) {
  k <- pcm_cells$n_items[i]
  sg <- pcm_cells$sigma[i]
  res <- t(vapply(seq_len(N_REP_B), function(r) {
    centres <- stats::runif(k, -2, 2)
    thr_list <- lapply(centres, function(cc) sort(cc + c(-0.8, 0, 0.8)))
    theta <- stats::rnorm(N_PERSON, 0, sg)
    pcm_coefs(theta, thr_list)
  }, numeric(4L)))
  data.frame(
    n_items = k, sigma = sg,
    exact = mean(res[, "exact"]),
    subtractive = mean(res[, "subtractive"]),
    ratio = mean(res[, "ratio"]),
    ratio_pw = mean(res[, "ratio_pointwise"]),
    ae_subtractive = mean(abs(res[, "subtractive"] - res[, "exact"])),
    ae_ratio = mean(abs(res[, "ratio"] - res[, "exact"])),
    ae_ratio_pw = mean(abs(res[, "ratio_pointwise"] - res[, "exact"])),
    pct_negative = 100 * mean(res[, "subtractive"] < 0),
    stringsAsFactors = FALSE
  )
}))

cat("\n-- PCM, 4 categories per item,", N_REP_B, "replications per cell --\n")
print(format(
  data.frame(
    items = partB$n_items, sigma = partB$sigma,
    exact = round(partB$exact, 3),
    subtractive = round(partB$subtractive, 3),
    ratio = round(partB$ratio, 3),
    ratio_pw = round(partB$ratio_pw, 3),
    ae_sub = round(partB$ae_subtractive, 3),
    ae_ratio = round(partB$ae_ratio, 3),
    pct_neg = round(partB$pct_negative, 1)
  ),
  nsmall = 3
), row.names = FALSE)

cat("\nOverall MAE  subtractive:", round(mean(partB$ae_subtractive), 4),
    " ratio:", round(mean(partB$ae_ratio), 4),
    " ratio_pointwise:", round(mean(partB$ae_ratio_pw), 4), "\n")

saveRDS(list(partA = partA, partB = partB),
        file.path(PKG, "dev", "milanzi_check_results.rds"))
