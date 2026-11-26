# Rough calibration of the misfit conditions planned for the specificity study
# (does RMitemCatFit() flag items whose categories work but which misfit in
# other ways?). Maps each generator parameter to the conditional infit it
# produces, with phq9-style data, n = 600, for comparison with the infit
# values of the real phq9 items (0.78 to 1.31).
#
# Item parameters: CML thresholds from the packaged phq9 responses (9 items,
# 4 categories), latent SD 1.39 estimated from the same responses (as in
# infit_width_iterations_poly_weak.qmd). Focus item q5 (location -0.15, on
# target). Three arms:
#   slope:  q5 follows a generalised PCM with slope a (a = 1 is the PCM).
#   multi:  q5 is driven by theta2, correlated rho with theta (the mechanism
#           of the infit studies, same SD).
#   ld:     response dependence (as in the Q3 study, prestudy_copulas.R):
#           q5's location is shifted by d * z, z = (2 * x_q1 - m) / m, where
#           x_q1 is the response to q1 (location -0.13).
# Conditional infit (iarm, via RMitemInfit without a cutoff) for all items.
#
# Run from easyRasch2/:
#   Rscript dev/category_conditional_null/R/infit_calibration.R [reps] [cores]
# Container: reps = 100, cores = 4 (a few minutes).

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 100L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L
N <- 600L
FOCUS <- 5L
PARTNER <- 1L

suppressMessages(devtools::load_all(".", quiet = TRUE))
out_dir <- "dev/category_conditional_null/results"

X0 <- as.matrix(phq9[, paste0("q", 1:9)])
X0 <- X0[stats::complete.cases(X0), ]
THR <- .fit_cml_thresholds(X0)
g <- seq(-6, 6, length.out = 81)
SD_THETA <- .estimate_prior_sd(.grid_loglik(X0, .logp_tables(THR, g), g), g, 0)
REAL_INFIT <- RMitemInfit(as.data.frame(X0), output = "dataframe")

# Draw one PCM item for locations theta (vector) with slope a.
draw_item <- function(theta, t, a = 1) {
  lp <- a * (outer(theta, 0:length(t)) -
    matrix(c(0, cumsum(t)), length(theta), length(t) + 1L, byrow = TRUE))
  p <- exp(lp - apply(lp, 1, max)); p <- p / rowSums(p)
  rowSums(stats::runif(length(theta)) > t(apply(p, 1, cumsum)))
}

gen <- function(arm, par) {
  theta <- stats::rnorm(N, 0, SD_THETA)
  X <- vapply(THR, function(t) draw_item(theta, t), numeric(N))
  m <- length(THR[[FOCUS]])
  if (arm == "slope") {
    X[, FOCUS] <- draw_item(theta, THR[[FOCUS]], par)
  } else if (arm == "multi") {
    theta2 <- SD_THETA * (par * theta / SD_THETA +
                            sqrt(1 - par^2) * stats::rnorm(N))
    X[, FOCUS] <- draw_item(theta2, THR[[FOCUS]])
  } else if (arm == "ld") {
    z <- (2 * X[, PARTNER] - m) / m
    X[, FOCUS] <- draw_item(theta + par * z, THR[[FOCUS]])
  }
  colnames(X) <- names(THR)
  as.data.frame(X)
}

grid <- rbind(
  data.frame(arm = "slope", par = c(0.5, 0.6, 0.7, 0.85, 1, 1.25, 1.5, 2)),
  data.frame(arm = "multi", par = c(0.7, 0.8, 0.85, 0.9, 0.95)),
  data.frame(arm = "ld",    par = c(0.5, 1, 1.5, 2, 2.5))
)

RNGkind("L'Ecuyer-CMRG")
set.seed(20261123)
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(grid)), function(gi) {
  do.call(rbind, lapply(seq_len(REPS), function(rp) {
    X <- gen(grid$arm[gi], grid$par[gi])
    inf <- RMitemInfit(X, output = "dataframe")
    data.frame(arm = grid$arm[gi], par = grid$par[gi], rep = rp,
               item = inf$Item, infit = inf$Infit_MSQ)
  }))
}, mc.cores = CORES)
bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE)
if (any(bad)) stop(res[bad][[1]])
reps <- do.call(rbind, res)

foc <- names(THR)[FOCUS]
par_ <- names(THR)[PARTNER]
summ <- do.call(rbind, lapply(split(reps, list(reps$arm, reps$par), drop = TRUE),
  function(d) {
    f <- d$infit[d$item == foc]
    p <- d$infit[d$item == par_]
    o <- d$infit[!d$item %in% c(foc, par_)]
    data.frame(arm = d$arm[1], par = d$par[1], reps = length(f),
               focus_mean = mean(f), focus_lo = stats::quantile(f, .025, names = FALSE),
               focus_hi = stats::quantile(f, .975, names = FALSE),
               partner_mean = mean(p), others_mean = mean(o))
  }))
# Rough null range per item from the PCM cell (slope a = 1).
null <- reps[reps$arm == "slope" & reps$par == 1, ]
null_range <- do.call(rbind, lapply(split(null, null$item), function(d)
  data.frame(item = d$item[1],
             lo = stats::quantile(d$infit, .025, names = FALSE),
             hi = stats::quantile(d$infit, .975, names = FALSE))))

out <- list(settings = list(reps = REPS, n = N, focus = foc, partner = par_,
                            sd_theta = SD_THETA, thr = THR),
            real_infit = REAL_INFIT, summary = summ[order(summ$arm, summ$par), ],
            null_range = null_range, reps = reps)
saveRDS(out, file.path(out_dir, "infit_calibration.rds"))
saveRDS(out, file.path(out_dir, paste0("infit_calibration_", Sys.Date(), ".rds")))
message("done, ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
