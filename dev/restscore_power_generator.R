# Data generator shared by restscore_power_calibration.R and
# restscore_infit_power.qmd, so the calibrated parameters mean the same thing
# in both.
#
# 7 polytomous items, 5 categories, thresholds at location +/- 0.5 and 1.5.
# theta ~ N(0, 1.5^2). Misfitting items are I1 (and I2), placed at `target`
# (+/- 0.25 when there are two); fitting items are evenly spread over
# [-1.5, 1.5].
#
# underfit: the misfitting items respond to a second dimension correlated
#           `rho` with the first, as in Johansson (2025, 2026). With two
#           misfitting items they share that dimension.
# overfit:  the misfitting items have slope `slope` > 1 on the first
#           dimension (generalised partial credit model).

THETA_SD <- 1.5
K_ITEMS  <- 7

item_setup <- function(n_misfit, target) {
  loc <- numeric(K_ITEMS)
  mis <- seq_len(n_misfit)
  loc[mis] <- if (n_misfit == 1) target else target + c(-0.25, 0.25)
  loc[-mis] <- seq(-1.5, 1.5, length.out = K_ITEMS - n_misfit)
  list(thr = lapply(loc, function(l) l + c(-1.5, -0.5, 0.5, 1.5)),
       misfit = mis)
}

# Generalised partial credit responses for one item; slope 1 is the PCM.
sim_gpcm_item <- function(theta, d, a = 1) {
  n <- length(theta)
  m <- length(d)
  eta <- a * (outer(theta, 0:m) -
                matrix(c(0, cumsum(d)), n, m + 1, byrow = TRUE))
  p <- exp(eta - apply(eta, 1, max))
  p <- p / rowSums(p)
  cp <- p %*% upper.tri(diag(m + 1), diag = TRUE)
  rowSums(runif(n) > cp)
}

# One dataset with every category of every item observed (redrawn otherwise).
sim_power_data <- function(n, n_misfit, target, direction, rho = 1,
                           slope = 1, seed) {
  setup <- item_setup(n_misfit, target)
  set.seed(seed)
  repeat {
    th1 <- rnorm(n, 0, THETA_SD)
    th2 <- rho * th1 + sqrt(1 - rho^2) * rnorm(n, 0, THETA_SD)
    d <- sapply(seq_len(K_ITEMS), function(j) {
      if (j %in% setup$misfit) {
        if (direction == "underfit") {
          sim_gpcm_item(th2, setup$thr[[j]])
        } else {
          sim_gpcm_item(th1, setup$thr[[j]], slope)
        }
      } else {
        sim_gpcm_item(th1, setup$thr[[j]])
      }
    })
    if (all(apply(d, 2, function(x) all(tabulate(x + 1, 5) > 0)))) break
  }
  d <- as.data.frame(d)
  names(d) <- paste0("I", seq_len(K_ITEMS))
  d
}
