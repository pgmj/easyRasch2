# Does the EAP prior carry the same fixed-mean problem as the latent SD?
#
# RMpersonParameters(method = "EAP") has prior_mean = 0 by default, documented
# as 'the natural choice when item parameters are centred at mean difficulty
# zero'. That is the same reasoning the 1.3.1 fix rejected: centring the items
# fixes the ITEM mean and says nothing about where the respondents sit.
#
# Answer: yes, but prior_sd is estimated with the mean pinned to 0, so it comes
# out inflated, and a very wide prior barely shrinks. The two errors cancel.
# Consequence: fixing ONE of them makes EAP worse than fixing neither. See the
# 'fix SD only' row.
#
# Item parameters are supplied at their true values so the scale origin is
# fixed; re-fitting by CML re-centres on whichever items survive and makes the
# bias column incomparable across rows.

pkgload::load_all(".", quiet = TRUE)
er <- asNamespace("easyRasch2")
sim_pcm <- function(theta, thr_list) {
  out <- matrix(NA_integer_, length(theta), length(thr_list))
  for (j in seq_along(thr_list)) {
    d <- thr_list[[j]]
    cs <- matrix(vapply(theta, function(th) cumsum(th - d), numeric(length(d))),
                 nrow = length(d))
    eta <- cbind(0, t(cs)); p <- exp(eta - apply(eta,1,max)); p <- p/rowSums(p)
    out[,j] <- apply(p,1,function(pr) sample.int(length(pr),1L,prob=pr))-1L
  }
  colnames(out) <- names(thr_list); as.data.frame(out)
}
# Wide item spread so an off-target sample still uses every category and no
# item has to be dropped; dropping items re-centres the scale and makes the
# bias column incomparable across rows.
thr <- setNames(lapply(seq(-2, 2, length.out = 12),
                       function(b) b + c(-0.8, 0, 0.8)), paste0("i", 1:12))

report <- function(shift, sd_true = 0.9, n = 2000, seeds = c(31, 77)) {
  cat(sprintf("\n--- true theta ~ N(%.1f, %.1f^2), n = %d ---\n", shift, sd_true, n))
  acc <- list()
  for (sd_i in seeds) {
    set.seed(sd_i)
    th <- rnorm(n, shift, sd_true); d <- sim_pcm(th, thr)
    # Item parameters are supplied at their true values rather than
    # re-estimated, so the scale origin is fixed and the bias column means
    # what it says. Re-fitting by CML re-centres on whichever items survive.
    fitted <- thr
    lat <- er$.latent_moments(as.matrix(d), fitted)
    ge <- seq(-6, 6, length.out = 81L)
    sd0 <- er$.estimate_prior_sd(
      er$.grid_loglik(as.matrix(d), er$.logp_tables(fitted, ge), ge), ge, 0)
    grab <- function(...) suppressWarnings(suppressMessages(
      RMpersonParameters(d, method = "EAP", item_params = thr,
                         output = "dataframe", ...)))$theta
    w <- suppressWarnings(suppressMessages(
      RMpersonParameters(d, method = "WLE", item_params = thr,
                         output = "dataframe")))$theta
    cat(sprintf("  seed %d: latent mean %.3f, latent sd %.3f | sd with mean pinned to 0: %.3f\n",
                sd_i, lat$mean, lat$sd, sd0))
    cs <- list("WLE" = w,
               "EAP default (mean 0, sd estimated)" = grab(),
               "EAP mean 0, sd = latent sd (fix SD only)" = grab(prior_sd = lat$sd),
               "EAP mean and sd both from latent fit" = grab(prior_mean = lat$mean,
                                                             prior_sd = lat$sd))
    acc[[as.character(sd_i)]] <- vapply(cs, function(e)
      c(bias = mean(e - th), rmse = sqrt(mean((e - th)^2))), numeric(2))
  }
  if (!length(acc)) return(invisible())
  m <- Reduce(`+`, acc) / length(acc)
  cat(sprintf("  %-42s %8s %8s\n", "estimator (mean over seeds)", "bias", "RMSE"))
  for (j in seq_len(ncol(m)))
    cat(sprintf("  %-42s %8.3f %8.3f\n", colnames(m)[j], m[1, j], m[2, j]))
}
report(0); report(1.5); report(3)
