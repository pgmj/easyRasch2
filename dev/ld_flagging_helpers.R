# =====================================================================
# Shared setup for the two local dependence flagging studies
# (dev/ld_q3_band_vs_wy.R and dev/ld_gamma_band_vs_wy.R).
#
# Ported from the companion simulation repository rasch_q3fwer so that
# the two scripts here run from this package alone. The generator is the
# Gaussian copula arm of that study: item responses are drawn from a
# partial credit model by inverse transform sampling, and a dependent
# pair is created by correlating the two uniforms it draws from. The
# marginal PCM distribution of each item is untouched by construction,
# which is what keeps the dependent and independent conditions
# comparable.
#
# Sourced, not run directly.
# =====================================================================

# Cumulative PCM category probabilities, n rows by (m + 1) columns.
.pcm_cum <- function(theta, deltas) {
  m <- length(deltas)
  num <- matrix(0, length(theta), m + 1L)
  cd <- cumsum(deltas)
  for (x in seq_len(m)) num[, x + 1L] <- x * theta - cd[x]
  num <- num - apply(num, 1L, max) # stabilise before exponentiating
  p <- exp(num)
  p <- p / rowSums(p)
  for (x in 2L:(m + 1L)) p[, x] <- p[, x - 1L] + p[, x]
  p
}

# Inverse transform sampling. Supplying `u` is what makes the copula
# possible: correlate the uniforms and the marginals are unchanged.
.pcm_draw <- function(theta, deltas, u = stats::runif(length(theta))) {
  cp <- .pcm_cum(theta, deltas)
  pmin(rowSums(cp < u), length(deltas))
}

# phq9 item parameters: CML thresholds and the estimated latent SD, from
# the real responses shipped with the package.
phq9_arm <- function() {
  e <- new.env()
  utils::data("phq9", package = "easyRasch2", envir = e)
  obs <- as.matrix(e$phq9[, paste0("q", 1:9)])
  obs <- obs[stats::complete.cases(obs), ]
  thr <- easyRasch2:::.fit_cml_thresholds(obs)
  g <- seq(-6, 6, length.out = 81)
  sd_theta <- easyRasch2:::.estimate_prior_sd(
    easyRasch2:::.grid_loglik(obs, easyRasch2:::.logp_tables(thr, g), g), g, 0
  )
  list(thr = thr, sd_theta = sd_theta)
}

# One dataset. `pair` is NULL for the complete null, or a length-2 vector
# of item indices made locally dependent at Gaussian copula correlation
# `rho`.
gen_data <- function(n, thr, sd_theta, pair = NULL, rho = 0) {
  k <- length(thr)
  theta <- stats::rnorm(n, 0, sd_theta)
  X <- matrix(0L, n, k, dimnames = list(NULL, names(thr)))
  independent <- if (is.null(pair)) seq_len(k) else setdiff(seq_len(k), pair)
  for (i in independent) X[, i] <- .pcm_draw(theta, thr[[i]])
  if (!is.null(pair)) {
    e1 <- stats::rnorm(n)
    e2 <- stats::rnorm(n)
    X[, pair[1L]] <- .pcm_draw(theta, thr[[pair[1L]]], stats::pnorm(e1))
    X[, pair[2L]] <- .pcm_draw(
      theta, thr[[pair[2L]]], stats::pnorm(rho * e1 + sqrt(1 - rho^2) * e2)
    )
  }
  as.data.frame(X)
}

# Binomial standard error, for reporting a rate with its precision.
rate_se <- function(p, n) sqrt(p * (1 - p) / n)
