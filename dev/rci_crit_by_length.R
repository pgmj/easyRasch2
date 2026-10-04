## Exact RCI critical value by test length, up to 600 items.
##
## Companion to dev/rci_exact_null.qmd, which tabulates 4 to 40 items. This
## extends the table to find where the critical value reaches 1.96. It does
## not reach it by 600 items: about 1.955 at four categories and 1.953 for
## dichotomous items.
##
## Person locations are integrated over N(0, 1.4) by quadrature rather than
## sampled, so the result is deterministic. Up to 40 items it matches the
## note's sampled table to the third decimal.
##
## Run: source() this file from easyRasch2/dev. Around one minute up to 200
## items, a further three minutes for 300 to 600.

suppressMessages(devtools::load_all("..", quiet = TRUE))

## ---- thresholds -------------------------------------------------------
# The note's constructor: k items over [-1.2, 1.2], ncat response categories
mk <- function(k, ncat) {
  lapply(seq(-1.2, 1.2, length.out = k),
         function(b) b + seq(-0.9, 0.9, length.out = ncat - 1L))
}
# Dichotomous items over the same range. mk(k, 2) would shift every item
# down by 0.9 logits, which confounds test length with targeting.
mk2 <- function(k) as.list(seq(-1.2, 1.2, length.out = k))

## ---- Lord-Wingersky: exact sum-score distribution given theta ----------
score_dist <- function(thr_list, theta) {
  d <- 1
  for (thr in thr_list) {
    p <- easyRasch2:::.pcm_cat_probs(theta, thr)
    nd <- numeric(length(d) + length(p) - 1L)
    for (a in seq_along(d)) {
      nd[a + seq_along(p) - 1L] <- nd[a + seq_along(p) - 1L] + d[a] * p
    }
    d <- nd
  }
  d                                   # index is score + 1
}

## ---- exact upper critical value, latent density by quadrature ----------
exact_crit_quad <- function(thr_list, mu = 0, sigma = 1.4, nq = 141L,
                            alpha = 0.05) {
  steps <- vapply(thr_list, length, integer(1L))
  lut <- vapply(0:sum(steps), function(r) {
    easyRasch2:::.theta_wle(
      easyRasch2:::.score_pattern(r, steps), thr_list, c(-10, 10)
    )
  }, numeric(2L))
  th <- lut[1L, ]
  se <- lut[2L, ]
  rci <- outer(seq_along(th), seq_along(th),
               function(i, j) (th[j] - th[i]) / sqrt(se[i]^2 + se[j]^2))
  q <- seq(mu - 5 * sigma, mu + 5 * sigma, length.out = nq)
  qw <- dnorm(q, mu, sigma)
  w <- Reduce(`+`, Map(function(t, wt) {
    d <- score_dist(thr_list, t)
    wt * outer(d, d)
  }, q, qw))
  w <- w / sum(w)
  o <- order(as.numeric(rci))
  v <- as.numeric(rci)[o]
  cw <- cumsum(as.numeric(w)[o])
  v[which(cw >= 1 - alpha / 2)[1L]]
}

## ---- sweep --------------------------------------------------------------
ks <- c(4, 6, 10, 20, 30, 40, 60, 80, 100, 120, 150, 200, 300, 400, 600)
res <- data.frame(
  items = ks,
  cat4  = vapply(ks, function(k) exact_crit_quad(mk(k, 4)), numeric(1L)),
  cat2  = vapply(ks, function(k) exact_crit_quad(mk2(k)), numeric(1L))
)
print(res, digits = 4, row.names = FALSE)

## Result, 2026-09-23:
##  items  cat4  cat2
##      4 1.644 1.271
##      6 1.725 1.602
##     10 1.815 1.734
##     20 1.868 1.720
##     30 1.896 1.854
##     40 1.915 1.855
##     60 1.933 1.901
##     80 1.939 1.890
##    100 1.938 1.926
##    120 1.941 1.918
##    150 1.948 1.937
##    200 1.949 1.931
##    300 1.953 1.940
##    400 1.957 1.947
##    600 1.956 1.953
