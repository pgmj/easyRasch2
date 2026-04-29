# Internal patched copies of two iarm functions used by RMscoreSE().
#
# The upstream `iarm::person_estimates()` and `iarm::persons.mle()` (iarm 0.4.x)
# converge poorly for extreme raw scores: the default `sthetarange` and `maxval`
# can leave the boundary scores with NA / non-finite SEs. The functions below
# are direct copies of the iarm code with these defaults widened so that score-
# to-logit tables produced for typical Rasch analyses come out clean at the
# tails.
#
# Original source: https://github.com/cran/iarm/blob/master/R/Person-Fit.R
# License: GPL (>= 2), as iarm.

#' Patched iarm::person_estimates()
#'
#' Score-to-theta lookup with MLE and WLE columns plus optional precision
#' properties (SEM, bias, RMSE, score-SEM). Used internally by [RMscoreSE()].
#'
#' @param object An eRm `Rm` object or a psychotools `raschmodel` / `pcmodel`.
#' @param properties Logical. Return precision properties.
#' @param allperson Logical. If TRUE, return one row per respondent rather than
#'   per raw score.
#' @param sthetarange Numeric length 2. Theta search range for boundary scores.
#'
#' @return A matrix or list, depending on `properties` / `allperson`.
#'
#' @keywords internal
#' @noRd
iarm_person_estimates <- function(object, properties = FALSE, allperson = FALSE,
                                  sthetarange = c(-10, 10)) {

  if (!any("Rm" %in% class(object),
           class(object) %in% c("raschmodel", "pcmodel"))) {
    stop("`object` must be of class Rm, raschmodel, or pcmodel.",
         call. = FALSE)
  }
  if (class(object)[1L] == "pcmodel")    object$model <- "pcmodel"
  if (class(object)[1L] == "raschmodel") object$model <- "raschmodel"

  X <- if (object$model %in% c("raschmodel", "pcmodel")) object$data else object$X

  if (object$model %in% c("RM", "raschmodel")) {
    k <- dim(X)[2L]
    coeff <- if (object$model == "RM") (-1) * stats::coef(object) else psychotools::itempar(object)
    m <- k
    respm <- rbind(rep(0, k), lower.tri(matrix(1, k, k)) + diag(k))
  } else {
    if (object$model == "PCM") {
      thr <- eRm::thresholds(object)[[3L]][[1L]][, -1L]
      coeff <- thr - mean(thr, na.rm = TRUE)
    } else {
      coeff <- stats::coef(psychotools::threshpar(object), type = "matrix")
    }
    k  <- dim(X)[2L]
    mi <- apply(X, 2L, max, na.rm = TRUE)
    m  <- sum(mi)
    respm <- matrix(0, ncol = k, nrow = m + 1L)
    respm[, 1L] <- c(0:mi[1L], rep(mi[1L], nrow(respm) - mi[1L] - 1L))
    for (i in 2:k) {
      respm[, i] <- c(rep(0, cumsum(mi)[i - 1L] + 1L),
                      1:mi[i],
                      rep(mi[i], nrow(respm) - cumsum(mi)[i] - 1L))
    }
  }

  mode <- if (object$model == "pcmodel") {
    "PCM"
  } else if (object$model == "raschmodel") {
    "RM"
  } else {
    object$model
  }

  mm <- cbind(
    0:m,
    iarm_persons_mle(respm, coeff, model = mode, type = "MLE")[, 1L],
    iarm_persons_mle(respm, coeff, model = mode, type = "WLE")[, 1L]
  )
  rownames(mm) <- rep(" ", m + 1L)
  colnames(mm) <- c("Raw Score", "MLE", "WLE")

  if (allperson) {
    rv <- rowSums(X, na.rm = TRUE)
    return(mm[rv + 1L, ])
  }

  if (!properties) return(mm)

  if (object$model %in% c("RM", "raschmodel")) {
    koeff <- as.list(coeff)
  } else {
    koeff <- lapply(as.list(as.data.frame(t(coeff))),
                    function(x) cumsum(stats::na.omit(x)))
  }
  gr <- psychotools::elementary_symmetric_functions(koeff)[[1L]]

  s.theta <- function(r) {
    function(x) {
      ((exp(x * (0:m)) * gr) / as.vector(exp(x * (0:m)) %*% gr)) %*% (0:m) - r
    }
  }

  if (object$model %in% c("pcmodel", "raschmodel")) {
    mm[1L, 2L] <- NA
  } else {
    mm[1L, 2L] <- eRm::person.parameter(object)$pred.list[[1L]]$y[1L]
  }
  try(mm[1L, 2L]      <- stats::uniroot(s.theta(0.25),     sthetarange)$root,
      silent = TRUE)
  try(mm[m + 1L, 2L]  <- stats::uniroot(s.theta(m - 0.25), sthetarange)$root,
      silent = TRUE)

  rvec <- 0:m
  pers_prop <- function(x, persons) {
    pr        <- (exp(x[2L] * rvec) * gr) / as.vector(exp(x[2L] * rvec) %*% gr)
    bias      <- pr %*% persons - x[2L]
    sem       <- sqrt((persons - as.vector(pr %*% persons))^2 %*% pr)
    rsem      <- sqrt((persons - x[2L])^2 %*% pr)
    scoresem  <- sqrt((rvec - x[1L])^2 %*% pr)
    c(SEM = sem, Bias = bias, RMSE = rsem, Score.SEM = scoresem)
  }

  list(
    cbind(mm[, 1:2], t(apply(mm[, c(1L, 2L)], 1L, pers_prop, persons = mm[, 2L]))),
    cbind(mm[, c(1L, 3L)], t(apply(mm[, c(1L, 3L)], 1L, pers_prop, persons = mm[, 3L])))
  )
}

#' Patched iarm::persons.mle()
#'
#' Newton-style WLE / MLE estimation given a response matrix and item-threshold
#' matrix. `maxval = 9` retains finite WLE values across a wide raw-score range.
#'
#' @keywords internal
#' @noRd
iarm_persons_mle <- function(respm, thresh,
                             model = c("RM", "PCM"),
                             theta = rep(0, dim(respm)[1L]),
                             type  = c("MLE", "WLE"),
                             extreme = TRUE,
                             maxit = 20L, maxdelta = 3, tol = 1e-04, maxval = 9) {

  n  <- dim(respm)[1L]
  k  <- dim(respm)[2L]
  rv <- rowSums(respm, na.rm = TRUE)
  mode <- match.arg(model)
  typ  <- match.arg(type)

  cll.rasch <- function(theta) {
    ksi <- exp(theta)
    mm  <- outer(ksi, 1 / exp(thresh))
    mm[is.na(respm)] <- 0
    dll  <- rv - rowSums(mm / (1 + mm))
    d2ll <- -rowSums(mm / (1 + mm)^2)
    d3ll <- 0
    if (typ == "WLE") {
      d3ll <- -rowSums((mm * (1 - mm)) / (1 + mm)^3)
      if (!extreme) {
        d3ll[rv == 0L]    <- 0
        d3ll[rv == maxr]  <- 0
      }
    }
    list(dll = dll, d2ll = d2ll, d3ll = d3ll)
  }

  cll.pcm <- function(theta) {
    dlogki <- function(i) {
      mmn <- exp(outer(theta, 1:mi[i]) +
                   matrix(psi.l[[i]], ncol = mi[i], nrow = n, byrow = TRUE))
      kd  <- 1 + rowSums(mmn)
      kd1 <- rowSums(matrix(1:mi[i],         ncol = mi[i], nrow = n, byrow = TRUE) * mmn)
      kd2 <- rowSums(matrix((1:mi[i])^2,     ncol = mi[i], nrow = n, byrow = TRUE) * mmn)
      kd3 <- rowSums(matrix((1:mi[i])^3,     ncol = mi[i], nrow = n, byrow = TRUE) * mmn)
      cbind(
        dlli  = kd1 / kd,
        d2lli = kd2 / kd - (kd1 / kd)^2,
        d3lli = -kd3 / kd + 3 * kd2 * kd1 / (kd^2) - 2 * (kd1 / kd)^3
      )
    }
    mm <- sapply(1:k, dlogki)
    mm[is.na(rbind(respm, respm, respm))] <- 0
    dll  <- rv -  rowSums(mm[1:n, ])
    d2ll <- -rowSums(mm[(n + 1):(2 * n), ])
    d3ll <- 0
    if (typ == "WLE") {
      d3ll <- rowSums(mm[(2 * n + 1):(3 * n), ])
      if (!extreme) {
        d3ll[rv == 0L]   <- 0
        d3ll[rv == maxr] <- 0
      }
    }
    list(dll = dll, d2ll = d2ll, d3ll = d3ll)
  }

  iter <- 1L
  conv <- 1
  if (mode == "RM") {
    maxr <- apply(respm, 1L, function(x) sum(!is.na(x)))
    clog <- cll.rasch
  } else {
    thresh.l <- apply(thresh, 1L, function(x) as.vector(stats::na.omit(x)),
                      simplify = FALSE)
    psi.l    <- lapply(thresh.l, function(x) (-1) * cumsum(x))
    mi       <- sapply(psi.l, length)
    maxr     <- apply(respm, 1L, function(x) sum(mi[!is.na(x)]))
    clog     <- cll.pcm
  }

  while ((conv > tol) & (iter <= maxit)) {
    theta0 <- theta
    fn     <- clog(theta)
    delta  <- -fn[[1L]] / fn[[2L]]
    if (typ == "WLE") {
      delta <- -fn[[1L]] / fn[[2L]] - fn[[3L]] / (2 * fn[[2L]]^2)
    }
    maxdelta <- maxdelta / 1.05
    delta    <- ifelse(abs(delta) > maxdelta, sign(delta) * maxdelta, delta)
    theta    <- theta + delta
    theta    <- ifelse(abs(theta) > maxval,   sign(theta) * maxval,    theta)
    conv     <- max(abs(theta - theta0))
    iter     <- iter + 1L
  }

  se    <- sqrt(abs(-1 / fn[[2L]]))
  se    <- ifelse(abs(theta) == maxval, NA, se)
  theta <- ifelse(theta ==  maxval,  Inf,
           ifelse(theta == -maxval, -Inf, theta))

  structure(data.frame(est = theta, se = se),
            model = mode, type = typ)
}
