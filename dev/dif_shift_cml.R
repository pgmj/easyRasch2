# Prototype: conditional ML for the partial credit model with uniform DIF
# shifts between two groups. Development code for dev/dif-impact-design.md,
# not package code.
#
# Model. Item j has cumulative category parameters psi_j1, ..., psi_jm
# (pcmodel()'s parameterisation, psi_jc = sum of the first c Andrich
# thresholds). For a DIF item i, group 2 has psi_ic + c * d_i, which shifts
# every threshold of the item by d_i. All other items are common to both
# groups and act as the anchor set. With two groups the conditional likelihood
# is the sum of the two groups' PCM conditional likelihoods, each built from
# that group's category counts and score distribution.
#
# Identification follows pcmodel(): the first category parameter of the first
# item is fixed at 0. Thresholds are returned centred to grand mean zero for
# the reference group, as .center_thresholds() does in the package.

# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

# Downcode categories that are unobserved in the pooled data, as
# pcmodel(nullcats = "downcode") does. A category empty in only one group is
# left alone: under the shift model its thresholds are shared, so it needs no
# special handling.
.shift_downcode <- function(X) {
  recode <- vector("list", ncol(X))
  for (j in seq_len(ncol(X))) {
    u <- sort(unique(X[, j]))
    recode[[j]] <- u
    X[, j] <- match(X[, j], u) - 1L
  }
  attr(X, "recode") <- recode
  X
}

.shift_prep <- function(data, group, dif_items) {
  X <- as.matrix(data)
  storage.mode(X) <- "integer"
  if (anyNA(X)) stop("Prototype requires complete data.", call. = FALSE)
  if (is.null(colnames(X))) colnames(X) <- paste0("I", seq_len(ncol(X)))
  g <- factor(group)
  if (nlevels(g) != 2L || anyNA(g)) {
    stop("`group` must have exactly two levels and no NA.", call. = FALSE)
  }
  X <- .shift_downcode(X)
  m <- apply(X, 2, max)
  if (any(m < 1L)) stop("Constant item(s): ", paste(colnames(X)[m < 1L], collapse = ", "), call. = FALSE)
  if (is.character(dif_items)) dif_items <- match(dif_items, colnames(X))
  dif_items <- sort(unique(as.integer(dif_items)))
  P <- sum(m)
  item_of <- rep.int(seq_along(m), m)
  cat_of <- unlist(lapply(m, seq_len))
  # S maps the DIF shifts onto the cumulative parameters: psi_2 = psi + S d.
  S <- matrix(0, P, length(dif_items))
  for (q in seq_along(dif_items)) {
    pos <- which(item_of == dif_items[q])
    S[pos, q] <- cat_of[pos]
  }
  stats <- lapply(levels(g), function(lev) {
    Xg <- X[g == lev, , drop = FALSE]
    ctot <- unlist(lapply(seq_along(m), function(j) {
      tabulate(Xg[, j], nbins = m[j])
    }))
    ptot <- tabulate(rowSums(Xg) + 1L, nbins = P + 1L)
    list(ctot = ctot, ptot = ptot)
  })
  list(X = X, g = g, m = m, P = P, item_of = item_of, cat_of = cat_of,
       dif_items = dif_items, S = S, stats = stats)
}

# ---------------------------------------------------------------------------
# Likelihood and gradient
# ---------------------------------------------------------------------------

.shift_unpack <- function(par, prep) {
  P <- prep$P
  psi <- c(0, par[seq_len(P - 1L)])
  d <- if (length(prep$dif_items)) par[P - 1L + seq_along(prep$dif_items)] else numeric(0)
  psi2 <- if (length(d)) psi + drop(prep$S %*% d) else psi
  list(psi = psi, d = d, psi_g = list(psi, psi2))
}

.shift_negll <- function(par, prep) {
  u <- .shift_unpack(par, prep)
  ll <- 0
  for (k in 1:2) {
    st <- prep$stats[[k]]
    esf <- psychotools::elementary_symmetric_functions(
      split(u$psi_g[[k]], prep$item_of), order = 0L
    )[[1]]
    keep <- st$ptot > 0
    ll <- ll - sum(st$ctot * u$psi_g[[k]]) - sum(st$ptot[keep] * log(esf[keep]))
  }
  if (!is.finite(ll)) return(.Machine$double.xmax)
  -ll
}

# Gradient of the negative log-likelihood. Per group, with gamma0 the ESF and
# gamma1 its derivative from elementary_symmetric_functions(order = 1):
#   d(-ll)/d(psi_g) = ctot_g - sum_r ptot_gr * gamma1[r, ] / gamma0[r]
# which is the aggregated form of pcmodel()'s agrad().
.shift_grad <- function(par, prep) {
  u <- .shift_unpack(par, prep)
  gpsi <- vector("list", 2)
  for (k in 1:2) {
    st <- prep$stats[[k]]
    esf <- psychotools::elementary_symmetric_functions(
      split(u$psi_g[[k]], prep$item_of), order = 1L
    )
    keep <- which(st$ptot > 0)
    g0 <- esf[[1]][keep]
    g1 <- esf[[2]][keep, , drop = FALSE]
    gpsi[[k]] <- st$ctot - colSums(st$ptot[keep] * g1 / g0)
  }
  g_psi <- (gpsi[[1]] + gpsi[[2]])[-1]
  g_d <- if (length(prep$dif_items)) drop(crossprod(prep$S, gpsi[[2]])) else numeric(0)
  c(g_psi, g_d)
}

# ---------------------------------------------------------------------------
# Fit
# ---------------------------------------------------------------------------

#' @param data Complete response matrix scored from 0.
#' @param group Two-level grouping variable; the first level is the reference.
#' @param dif_items Items given a group-2 shift (names or column indices).
#'   Empty gives the ordinary pooled PCM.
#' @param hessian Compute the covariance matrix (numerical Hessian of the
#'   analytic gradient).
shift_cml <- function(data, group, dif_items = integer(0), hessian = TRUE,
                      reltol = 1e-10, maxit = 1000L) {
  prep <- .shift_prep(data, group, dif_items)
  P <- prep$P
  nd <- length(prep$dif_items)
  # Start from the pooled pcmodel() fit on complete data (unaffected by the
  # NA bug), with zero shifts.
  pooled <- suppressWarnings(psychotools::pcmodel(prep$X, hessian = FALSE))
  start <- c(unname(coef(pooled)), rep(0, nd))
  opt <- stats::optim(
    start, .shift_negll, .shift_grad, prep = prep, method = "BFGS",
    control = list(reltol = reltol, maxit = maxit)
  )
  u <- .shift_unpack(opt$par, prep)
  vc <- NULL
  if (hessian) {
    H <- stats::optimHess(opt$par, .shift_negll, .shift_grad, prep = prep)
    vc <- tryCatch(solve(H), error = function(e) NULL)
  }
  se_d <- if (nd && !is.null(vc)) sqrt(diag(vc)[P - 1L + seq_len(nd)]) else rep(NA_real_, nd)
  # Andrich thresholds, centred on the reference group's grand mean.
  thr <- lapply(split(u$psi, prep$item_of), function(p) diff(c(0, p)))
  shift <- mean(unlist(thr))
  thr <- lapply(thr, function(x) x - shift)
  names(thr) <- colnames(prep$X)
  thr2 <- thr
  for (q in seq_len(nd)) thr2[[prep$dif_items[q]]] <- thr[[prep$dif_items[q]]] + u$d[q]
  structure(list(
    d = stats::setNames(u$d, colnames(prep$X)[prep$dif_items]),
    se_d = stats::setNames(se_d, colnames(prep$X)[prep$dif_items]),
    thresholds = list(ref = thr, focal = thr2),
    psi = u$psi,
    loglik = -opt$value,
    convergence = opt$convergence,
    counts = opt$counts,
    vcov = vc,
    groups = levels(prep$g),
    recode = attr(prep$X, "recode"),
    prep = prep
  ), class = "shift_cml")
}

# ---------------------------------------------------------------------------
# Helpers used by the validation script
# ---------------------------------------------------------------------------

# Simulate PCM responses. thr_list: list of Andrich threshold vectors.
sim_pcm <- function(theta, thr_list) {
  n <- length(theta)
  sapply(thr_list, function(thr) {
    eta <- cbind(0, outer(theta, thr, "-"))
    eta <- t(apply(eta, 1, cumsum))
    p <- exp(eta - apply(eta, 1, max))
    p <- p / rowSums(p)
    cp <- p %*% upper.tri(diag(ncol(p)), diag = TRUE)
    rowSums(stats::runif(n) > cp)
  })
}

# Two-group data with uniform DIF: group 2 has every threshold of the items in
# `dif` shifted by the matching value.
sim_dif <- function(n1, n2, thr_list, dif = c(), mu2 = 0, sd = 1) {
  theta <- c(stats::rnorm(n1, 0, sd), stats::rnorm(n2, mu2, sd))
  g <- factor(rep(c("g1", "g2"), c(n1, n2)))
  thr2 <- thr_list
  for (nm in names(dif)) thr2[[nm]] <- thr2[[nm]] + dif[[nm]]
  X <- rbind(sim_pcm(theta[g == "g1"], thr_list), sim_pcm(theta[g == "g2"], thr2))
  colnames(X) <- names(thr_list)
  list(X = X, g = g, theta = theta)
}

# Item-split data: column `item` becomes one column per group.
split_items <- function(X, g, items) {
  out <- as.data.frame(X)
  for (it in items) {
    a <- ifelse(g == levels(g)[1], X[, it], NA)
    b <- ifelse(g == levels(g)[2], X[, it], NA)
    out[[it]] <- NULL
    out[[paste0(it, "_", levels(g)[1])]] <- a
    out[[paste0(it, "_", levels(g)[2])]] <- b
  }
  out
}

# psychotools pcmodel() with the SIMPLIFY = FALSE fix from
# dev/psychotools-na-bug/, applied at run time. Needed for free-split fits
# until a fixed psychotools is released.
pcmodel_fixed <- local({
  src <- deparse(psychotools::pcmodel)
  h <- grep("mapply(function(x, ", src, fixed = TRUE)
  h <- h[grepl("\\)\\)$", src[h + 1])]
  src[h + 1] <- sub("\\)\\)$", ", SIMPLIFY = FALSE), SIMPLIFY = FALSE)", src[h + 1])
  f <- eval(parse(text = src))
  environment(f) <- asNamespace("psychotools")
  f
})
