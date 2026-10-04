# Calibration for restscore_infit_power.qmd: which misfit parameters inject a
# realistic infit level into a 7-item, 5-category polytomous scale.
#
# Johansson (2026, FWER paper) anchors misfit severity on the infit scale
# rather than on the generating parameter, because the parameter is not
# portable across item formats or scale lengths (rho = .85 to .90 at nine phq9
# items, .75 at five items). The anchors are the phq9 items flagged at n = 600
# with simulation-based cutoffs:
#   underfit  q3, q8, q9   infit 1.23 to 1.31
#   overfit   q2, q4       infit 0.78 to 0.83
#
# Underfit is multidimensional, as in both earlier papers: the misfitting item
# responds to a second dimension correlated rho with the first. Overfit is a
# discrimination above 1 (generalised partial credit model, slope a). The
# scan uses a single misfitting item and observed conditional infit only, no
# bootstrap, at n = 600 as the anchor was.
#
# Run from easyRasch2/dev:  Rscript restscore_power_calibration.R
# Results are cached to restscore_power_calibration.rds.

suppressMessages({
  devtools::load_all("..", quiet = TRUE)
  library(parallel)
  library(dplyr)
})
options(rgl.useNULL = TRUE)

source("restscore_power_generator.R")

RES    <- "restscore_power_calibration.rds"
RHOS   <- c(0.65, 0.70, 0.75, 0.80, 0.85, 0.90)
SLOPES <- c(1.4, 1.6, 1.8, 2.0, 2.2)
N      <- 600
R_REPS <- 80
NCORES <- 10

grid <- rbind(
  expand.grid(mechanism = "underfit", par = RHOS, target = c(0, -2),
              rep = seq_len(R_REPS), stringsAsFactors = FALSE),
  expand.grid(mechanism = "overfit", par = SLOPES, target = c(0, -2),
              rep = seq_len(R_REPS), stringsAsFactors = FALSE)
)

if (!file.exists(RES)) {
  t0 <- Sys.time()
  vals <- mclapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    d <- sim_power_data(
      n = N, n_misfit = 1, target = g$target, direction = g$mechanism,
      rho = if (g$mechanism == "underfit") g$par else 1,
      slope = if (g$mechanism == "overfit") g$par else 1,
      seed = 7e5 + i
    )
    x <- try(suppressMessages(RMitemInfit(d, output = "dataframe")),
             silent = TRUE)
    if (inherits(x, "try-error")) NA_real_ else x$Infit_MSQ[1]
  }, mc.cores = NCORES)
  grid$infit <- unlist(vals)
  saveRDS(list(grid = grid, n = N, reps = R_REPS,
               elapsed_min = as.numeric(Sys.time() - t0, units = "mins")),
          RES)
}

cal <- readRDS(RES)
cal$grid |>
  group_by(mechanism, par, target) |>
  summarise(sd = round(sd(infit, na.rm = TRUE), 3),
            infit = round(mean(infit, na.rm = TRUE), 3),
            .groups = "drop") |>
  as.data.frame() |>
  print()
