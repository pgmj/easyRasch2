# Before/after check for the 1.3.1 latent-mean fix.
#
# Through 1.3.0 the latent distribution's mean was held at 0 and only its SD
# was estimated, so the SD absorbed any distance between the sample and the
# item locations. Both paths are computed here from the current code: the old
# behaviour by calling the SD-only estimator with the mean pinned, the new one
# through RMreliabilityCurve().
#
# Same seed, items and sample sizes as the original diagnosis, so the "before"
# column reproduces the numbers that prompted the fix.

pkgload::load_all(".", quiet = TRUE)
er <- asNamespace("easyRasch2")

thr <- list(
  i1 = c(-1.5, -0.3, 1.0), i2 = c(-1.0,  0.1, 1.4),
  i3 = c(-0.6,  0.5, 1.8), i4 = c(-0.2,  0.9, 2.2),
  i5 = c( 0.3,  1.4, 2.7), i6 = c(-2.0, -0.8, 0.6),
  i7 = c(-1.2,  0.2, 1.6), i8 = c(-0.4,  0.7, 2.0)
)

sim_pcm <- function(theta, thr_list) {
  out <- matrix(NA_integer_, length(theta), length(thr_list))
  for (j in seq_along(thr_list)) {
    d <- thr_list[[j]]
    # matrix() rather than t(vapply(...)) directly: with one threshold the
    # vapply result is a plain vector and t() would transpose it the wrong way.
    cs <- matrix(vapply(theta, function(th) cumsum(th - d),
                        numeric(length(d))), nrow = length(d))
    eta <- cbind(0, t(cs))
    p <- exp(eta - apply(eta, 1, max))
    p <- p / rowSums(p)
    out[, j] <- apply(p, 1, function(pr)
      sample.int(length(pr), 1L, prob = pr)) - 1L
  }
  colnames(out) <- names(thr_list)
  as.data.frame(out)
}

# The pre-1.3.1 path, rebuilt from the estimator that still backs the EAP prior
old_sigma <- function(d, thr_list) {
  ge <- seq(-6, 6, length.out = 81L)
  er$.estimate_prior_sd(
    er$.grid_loglik(as.matrix(d), er$.logp_tables(thr_list, ge), ge), ge, 0
  )
}

n <- 1500
true_sd <- 0.9
set.seed(7)

rows <- lapply(c(0, 1, 2), function(shift) {
  th <- rnorm(n, shift, true_sd)
  d  <- sim_pcm(th, thr)

  fitted <- er$.fit_cml_thresholds(as.matrix(d))
  s_old  <- old_sigma(d, fitted)
  m_old  <- er$.marginal_summaries(fitted, s_old, mu = 0)$ratio

  cd <- suppressWarnings(suppressMessages(
    RMreliabilityCurve(d, output = "dataframe")))
  pp <- suppressWarnings(suppressMessages(
    RMpersonParameters(d, output = "dataframe")))
  rr <- suppressWarnings(suppressMessages(
    RMreliability(d, draws = 200, rmu_iter = 10, boot = FALSE,
                  parallel = FALSE, seed = 1, output = "dataframe")))

  data.frame(
    shift        = shift,
    wle_mean     = mean(pp$theta),
    latent_mean  = attr(cd, "latent_mean", exact = TRUE),
    sigma_old    = s_old,
    sigma_new    = attr(cd, "sigma", exact = TRUE),
    marginal_old = m_old,
    marginal_new = attr(cd, "marginal_ratio", exact = TRUE),
    psi          = rr$estimate[rr$metric == "PSI"],
    alpha        = rr$estimate[grepl("alpha", rr$metric)]
  )
})

out <- do.call(rbind, rows)
cat("\nTrue person SD", true_sd, "at every shift; n =", n, "; 8 items, 4 categories.\n\n")
show <- data.frame(
  `shift`          = sprintf("%+.0f", out$shift),
  `WLE mean`       = sprintf("%.3f", out$wle_mean),
  `latent mean`    = sprintf("%.3f", out$latent_mean),
  `sigma 1.3.0`    = sprintf("%.3f", out$sigma_old),
  `sigma 1.3.1`    = sprintf("%.3f", out$sigma_new),
  `marginal 1.3.0` = sprintf("%.3f", out$marginal_old),
  `marginal 1.3.1` = sprintf("%.3f", out$marginal_new),
  `PSI`            = sprintf("%.3f", out$psi),
  `alpha`          = sprintf("%.3f", out$alpha),
  check.names = FALSE
)
print(show, row.names = FALSE)


# ---------------------------------------------------------------------------
# Where the latent mean and the mean of the WLE estimates part company.
#
# They are different objects. The latent mean is a parameter of the assumed
# normal population, fitted by integrating each response vector over theta, so
# no respondent is ever assigned a location. mean(theta_hat) averages n point
# estimates. A WLE estimate at the minimum or maximum score is a finite, and
# therefore bounded, extrapolation, so respondents whose true location lies
# past that bound are pulled back to it and the mean of the estimates is
# dragged inward. The marginal fit has no such ceiling.
#
# The two agree whenever extreme scores are rare, which is most well-targeted
# data, and separate sharply when they are not.
# ---------------------------------------------------------------------------

extreme_check <- function() {
  run <- function(lbl, thr_list, n, mu, sd) {
    th <- stats::rnorm(n, mu, sd)
    d  <- sim_pcm(th, thr_list)
    if (any(vapply(d, function(x) min(x) != 0, logical(1)))) {
      cat(sprintf("%-28s (skipped: an item lost its lowest category)\n", lbl))
      return(invisible())
    }
    fitted <- er$.fit_cml_thresholds(as.matrix(d))
    lat <- er$.latent_moments(as.matrix(d), fitted)
    w <- suppressWarnings(suppressMessages(
      RMpersonParameters(d, method = "WLE", output = "dataframe")))
    e <- suppressWarnings(suppressMessages(
      RMpersonParameters(d, method = "EAP", output = "dataframe")))
    cat(sprintf("%-28s %7.3f %8.3f %8.3f %8.3f %6.1f%%\n",
                lbl, mu, lat$mean, mean(w$theta), mean(e$theta),
                100 * mean(w$extreme)))
  }

  cat(sprintf("\n%-28s %7s %8s %8s %8s %7s\n",
              "case", "true mu", "latent", "WLE", "EAP", "%extr"))
  set.seed(21)
  long <- setNames(
    lapply(seq(-1.5, 1.5, length.out = 12), function(b) b + c(-0.8, 0, 0.8)),
    paste0("i", 1:12))
  short <- setNames(
    lapply(seq(-0.8, 0.8, length.out = 6), function(b) b), paste0("i", 1:6))
  run("12 poly items, on target",   long,  2000, 0.0, 0.9)
  run("6 dich items, on target",    short, 2000, 0.0, 1.0)
  run("6 dich items, wide spread",  short, 2000, 0.0, 2.2)
  run("6 dich items, off +1.5",     short, 2000, 1.5, 1.5)
  run("6 dich items, off +2, wide", short, 2000, 2.0, 2.0)
}

extreme_check()
