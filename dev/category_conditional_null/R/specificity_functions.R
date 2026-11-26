# Specificity study (2026-11-23): does RMitemCatFit() flag items whose
# categories work as scored but which misfit in other ways, and does
# conditional infit tell the two apart? Builds on the main study
# (study_functions.R, calibration_functions.R), which is not changed.
#
# Four new conditions for the focus item (item 1), one level each, calibrated
# against conditional infit in R/infit_calibration_stylised.R to the size of
# the real phq9 misfit (underfit about 1.28, overfit about 0.83):
#   U_multi_r75:  driven by theta2, correlation .75 with theta (infit 1.22-1.31)
#   U_slope_055:  generalised PCM slope 0.55              (infit 1.23-1.33)
#   O_ld_d1:      response dependence on the other item nearest location 0,
#                 location shift d * z with d = 1 and z = (2 * x - m) / m,
#                 as in the Q3 study                     (infit 0.80-0.83)
#   O_slope_16:   generalised PCM slope 1.6               (infit 0.80-0.84)

source("R/study_functions.R")
source("R/calibration_functions.R")

NUISANCE <- data.frame(
  cond = c("U_multi_r75", "U_slope_055", "O_ld_d1", "O_slope_16"),
  type = "nuisance", null = FALSE)

spec_conditions <- function() rbind(conditions(), NUISANCE)

# Thresholds of the focus item under the model the condition departs from,
# which is the null model for its infit reference.
base_of <- function(cond) {
  if (cond == "A1_pcm_disordered" || grepl("^DK_narrow", cond)) "NARROW" else "WIDE"
}

gen_spec <- function(cond, n, k) {
  if (!cond %in% NUISANCE$cond) return(gen_data(cond, n, k))
  # Same order of draws as gen_data(): theta, other items, focus item.
  theta <- stats::rnorm(n)
  oth <- other_items(k)
  Xo <- sim_pcm(theta, oth)
  m <- length(WIDE)
  partner <- which.min(abs(vapply(oth, mean, 0)))
  x1 <- switch(cond,
    U_multi_r75 = sim_steps(0.75 * theta + sqrt(1 - 0.75^2) * stats::rnorm(n), WIDE),
    U_slope_055 = sim_steps(theta, WIDE, rep(0.55, m)),
    O_ld_d1     = sim_steps(theta + 1 * (2 * Xo[, partner] - m) / m, WIDE),
    O_slope_16  = sim_steps(theta, WIDE, rep(1.6, m)))
  X <- cbind(x1, Xo)
  colnames(X) <- paste0("I", seq_len(ncol(X)))
  X
}

# Conditional infit of every item (iarm, through easyRasch2::RMitemInfit()
# without a cutoff). NULL when the data cannot be fitted (empty category).
infit_vec <- function(X) {
  d <- tryCatch(easyRasch2::RMitemInfit(as.data.frame(X), output = "dataframe"),
                error = function(e) NULL)
  if (is.null(d)) NULL else d$Infit_MSQ
}

# Infit reference set for one cell: infit of every item in `R` datasets
# simulated from the PCM with the focus item at the given base thresholds.
# The true parameters are used (an oracle reference, computed once per cell
# instead of a bootstrap per dataset).
infit_reference <- function(n, k, base, R, cores) {
  cond <- if (base == "NARROW") "A1_pcm_disordered" else "A0_pcm_ordered"
  res <- parallel::mclapply(seq_len(R), function(r) infit_vec(gen_data(cond, n, k)),
                            mc.cores = cores, mc.set.seed = TRUE)
  do.call(rbind, Filter(Negate(is.null), res))
}

# Westfall-Young step-down on statistics where larger = more extreme.
# s_obs: length k, s_sim: B x k.
stepdown <- function(s_obs, s_sim) {
  B <- nrow(s_sim)
  ord <- order(s_obs, decreasing = TRUE)
  run <- rep(-Inf, B)
  p <- numeric(length(ord))
  for (j in rev(seq_along(ord))) {
    run <- pmax(run, s_sim[, ord[j]])
    p[j] <- (1 + sum(run >= s_obs[ord[j]])) / (B + 1)
  }
  out <- numeric(length(ord))
  out[ord] <- pmin(cummax(p), 1)
  out
}

# Infit p-values from a reference set: each item's infit studentised by the
# mean and SD of the observed and reference values together, two-sided
# (|z|), item p = share of the pool at or above the observed |z|, and the
# Westfall-Young step-down across items, floored at the item p.
infit_p <- function(obs, ref) {
  pool <- rbind(obs, ref)
  z <- sweep(sweep(pool, 2, colMeans(pool)), 2, apply(pool, 2, stats::sd), "/")
  a <- abs(z)
  N <- nrow(a)
  p <- vapply(seq_along(obs), function(i) mean(a[, i] >= a[1, i]), 0)
  padj <- pmax(stepdown(a[1, ], a[-1, , drop = FALSE]), p)
  list(p = p, padj = padj)
}

# RMitemCatFit()'s combination of the events (studentised max-T) and the
# Westfall-Young step-down across items, from analyse_calibration()'s stored
# MCB values (as in compare_combination.R).
catfit_p <- function(obs, boot) {
  k <- nrow(obs); m <- ncol(obs); N <- dim(boot)[3] + 1L
  pool <- array(c(obs, boot), c(k, m, N))
  T_max <- matrix(NA_real_, N, k)
  for (i in seq_len(k)) {
    E <- matrix(pool[i, , ], nrow = m)
    Z <- t(apply(E, 1, function(t) {
      s <- stats::sd(t)
      if (!is.finite(s) || s < 1e-12) rep(0, length(t)) else (t - mean(t)) / s
    }))
    T_max[, i] <- apply(matrix(Z, nrow = m), 2, max)
  }
  p <- vapply(seq_len(k), function(i) mean(T_max[, i] >= T_max[1, i]), 0)
  list(p = p, padj = pmax(stepdown(T_max[1, ], T_max[-1, , drop = FALSE]), p))
}
