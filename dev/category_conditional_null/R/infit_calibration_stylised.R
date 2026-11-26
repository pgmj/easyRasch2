# Calibration of the two nuisance conditions for the specificity study on the
# stylised design of the main study (study_functions.R): focus item 1 at
# location 0 with steps (-1.5, -0.5, 0.5, 1.5), k - 1 PCM items with
# locations over [-1, 1], theta ~ N(0, 1), 5 categories.
#
# Targets, from the real phq9 items (infit_calibration.R): underfit about
# 1.28 (q3, q8, q9 at 1.23 to 1.31), overfit about 0.83 (q2, q4 at 0.78 to
# 0.83). Each target gets a mechanism and a pure slope with matched infit:
#   underfit: second dimension (focus item driven by theta2, correlation rho)
#             and generalised PCM slope a < 1;
#   overfit:  response dependence on the other item nearest location 0
#             (location shift d * z, z = (2 * x_partner - m) / m, as in the
#             Q3 study) and slope a > 1.
# Conditional infit (iarm, via RMitemInfit without a cutoff).
#
# Run from easyRasch2/dev/category_conditional_null/:
#   Rscript R/infit_calibration_stylised.R [reps] [cores]

args <- commandArgs(trailingOnly = TRUE)
REPS  <- if (length(args) >= 1) as.integer(args[1]) else 100L
CORES <- if (length(args) >= 2) as.integer(args[2]) else 4L

suppressMessages(devtools::load_all("../..", quiet = TRUE))
source("R/study_functions.R")

gen <- function(n, k, arm, par) {
  theta <- stats::rnorm(n)
  oth <- other_items(k)
  X <- cbind(0, sim_pcm(theta, oth))
  m <- length(WIDE)
  partner <- 1L + which.min(abs(vapply(oth, mean, 0)))
  X[, 1] <- switch(arm,
    pcm   = sim_steps(theta, WIDE),
    slope = sim_steps(theta, WIDE, rep(par, m)),
    multi = sim_steps(par * theta + sqrt(1 - par^2) * stats::rnorm(n), WIDE),
    ld    = sim_steps(theta + par * (2 * X[, partner] - m) / m, WIDE))
  colnames(X) <- paste0("I", seq_len(k))
  list(X = as.data.frame(X), partner = partner)
}

grid <- rbind(
  data.frame(arm = "pcm",   par = 1),
  data.frame(arm = "slope", par = c(0.55, 0.6, 0.65, 0.7, 1.3, 1.4, 1.5, 1.6)),
  data.frame(arm = "multi", par = c(0.75, 0.8, 0.85, 0.9)),
  data.frame(arm = "ld",    par = c(0.75, 1, 1.25, 1.5, 2))
)
grid <- merge(grid, expand.grid(k = c(5L, 10L), n = c(300L, 600L)))

RNGkind("L'Ecuyer-CMRG")
set.seed(20261124)
t0 <- Sys.time()
res <- parallel::mclapply(seq_len(nrow(grid)), function(gi) {
  gr <- grid[gi, ]
  do.call(rbind, lapply(seq_len(REPS), function(rp) {
    g <- gen(gr$n, gr$k, gr$arm, gr$par)
    # A dataset with an empty category is skipped, as in the main study.
    inf <- tryCatch(RMitemInfit(g$X, output = "dataframe")$Infit_MSQ,
                    error = function(e) NULL)
    if (is.null(inf)) return(NULL)
    data.frame(arm = gr$arm, par = gr$par, k = gr$k, n = gr$n, rep = rp,
               focus = inf[1], partner = inf[g$partner],
               others = mean(inf[-c(1, g$partner)]))
  }))
}, mc.cores = CORES)
bad <- vapply(res, function(x) inherits(x, "try-error"), TRUE)
if (any(bad)) stop(res[bad][[1]])
reps <- do.call(rbind, res)
n_skipped <- nrow(grid) * REPS - nrow(reps)
summ <- stats::aggregate(cbind(focus, partner, others) ~ arm + par + k + n,
                         reps, mean)
out <- list(settings = list(reps = REPS, skipped = n_skipped), summary = summ, reps = reps)
saveRDS(out, "results/infit_calibration_stylised.rds")
saveRDS(out, paste0("results/infit_calibration_stylised_", Sys.Date(), ".rds"))
message("done, ", round(as.numeric(Sys.time() - t0, units = "mins"), 1), " min")
