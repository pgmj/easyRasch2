# ---------------------------------------------------------------------------
# PROTOTYPE: partial gamma profile across the score range, with a uniform-LD
# reference ribbon.
#
# Not part of the package yet. Source this file after `library(easyRasch2)`
# and call RMlocdepGammaProfile(). Package internals are reached with `:::`;
# those qualifiers come off if and when this moves into R/.
#
# The local dependence analogue of RMitemICCPlot(): where the CICC plot shows
# item fit across the score range, this shows local dependence across the
# score range.
# ---------------------------------------------------------------------------

# The package imports rlang's `.data` pronoun via NAMESPACE. A sourced file has
# to bind it itself; this line comes out when the function moves into R/.
.data <- rlang::.data

#' Partial gamma profile across the score range
#'
#' Davis's partial gamma pools concordant and discordant pair counts over
#' strata of the rest score, and the pooling is what discards any information
#' about *where* on the scale a pair's dependence lives. This plots the
#' stratum-wise coefficients against the latent continuum, with a reference
#' ribbon simulated under **uniform** local dependence at the pair's observed
#' strength.
#'
#' @param data A data.frame or matrix of item responses, items scored from 0.
#' @param pairs A length-2 vector of item names or column indices, or a list of
#'   such vectors. Intended as a follow-up on pairs already flagged by
#'   [RMlocdepQ3()] or [RMlocdepGamma()], not as a screen over all pairs.
#' @param class_intervals Integer >= 2. Number of class intervals along the
#'   rest score. Default `5`. Quantile grouping can form fewer bins than
#'   requested when rest scores tie at the quantile boundaries.
#' @param direction One of `"auto"` (default), `"1"` or `"2"`. Partial gamma
#'   has two conditioning directions, one per item removed from the rest
#'   score. `"auto"` uses whichever gives the larger absolute pooled
#'   coefficient in the observed data, and holds that direction fixed for the
#'   simulated datasets.
#' @param band One of `"uniform"` (default) or `"none"`. The `"uniform"` band
#'   is simulated under constant dependence calibrated to the pair's observed
#'   pooled gamma, so departures read as *more or less dependence than
#'   uniform*. There is deliberately no independence band here: see Details.
#' @param band_type One of `"simultaneous"` (default) or `"pointwise"`.
#' @param conf_level Numeric in (0, 1). Coverage of the ribbon. Default `0.95`.
#' @param iterations Number of simulated datasets for the ribbon. Default
#'   `400`.
#' @param cal_iterations Number of simulated datasets per calibration grid
#'   point. Default `40`.
#' @param min_n Integer. A class interval needs at least this many respondents
#'   to be drawn. Default `10`.
#' @param x_axis One of `"theta"` (default) or `"restscore"`.
#' @param parallel,n_cores,seed,verbose Parallel and reproducibility control.
#' @param output One of `"ggplot"` (default) or `"dataframe"`.
#'
#' @details
#' \strong{Bin the results, not the data.} For a CICC plot a class interval is
#' display aggregation and is harmless. For partial gamma the stratification
#' *is* the method: conditioning on the rest score is what removes the
#' part-whole confounding. Collapse five rest-score levels into one bin and the
#' rest score still varies inside it, both items correlate with it, and they
#' reacquire a spurious positive association. Every coefficient here is
#' therefore computed within a single rest-score level and aggregated
#' afterwards,
#' \deqn{\gamma_B = \frac{\sum_{k \in B} (C_k - D_k)}{\sum_{k \in B} (C_k + D_k)},}
#' a weighted mean of the stratum coefficients with weights \eqn{C_k + D_k}.
#' Pooling over all bins returns the coefficient [RMlocdepGamma()] reports, so
#' the plot cannot contradict the printed number.
#'
#' \strong{Why the ribbon is not flat.} Under the Rasch model the stratum
#' coefficients are centred on zero, so a no-dependence profile is a flat line
#' at zero. A *uniformly dependent* profile is not flat at its pooled value.
#' It arches, because comparable pairs run out at the score extremes, and the
#' ends are attenuated for that reason rather than because the dependence
#' weakens there. A horizontal reference at the pooled gamma would therefore
#' make every dependent pair look non-uniform at both ends. The ribbon is
#' simulated instead, and it is the direct counterpart of averaging the
#' model expectation over the people in a CICC bin rather than evaluating it
#' at the bin's mean theta.
#'
#' \strong{Why there is no independence band.} Heterogeneity of the stratum
#' coefficients *decreases* as dependence strengthens, so a band drawn under
#' independence runs the wrong way for this question: a pair with strong
#' uniform dependence would look more homogeneous than that null. An
#' independence band answers "is there dependence", which is what
#' [RMlocdepGammaCutoff()] and [RMlocdepGammaPlot()] already do.
#'
#' \strong{Simultaneous versus pointwise.} With five bins, the chance that at
#' least one observed point falls outside a 95% *pointwise* ribbon under
#' uniformity is well above 5%. The default ribbon is simultaneous: a common
#' multiplier \eqn{c} is chosen so that `conf_level` of the simulated profiles
#' lie entirely inside \eqn{\mathrm{median}_b \pm c\,s_b}. Read the profile as
#' a whole rather than pointing at one bin.
#'
#' \strong{Sample size.} This is a large-sample diagnostic. The null SD of a
#' pooled partial gamma at n = 600 is around .067, a single rest-score stratum
#' holds roughly 1/25 of the information, and five bins leave a per-bin SD near
#' .13. Below n around 600 the ribbon will be wide enough that only gross
#' departures are readable, which is the honest answer rather than a defect.
#'
#' \strong{Dichotomous items.} A dichotomous pair yields far fewer comparable
#' pairs per stratum than a polytomous one, so the ribbon is frequently too
#' wide to read, and the end bins often carry no comparable pair at all. Those
#' bins are dropped rather than drawn. Read the `n` and `drawn` columns of the
#' `"dataframe"` output before concluding anything from a sparse figure.
#'
#' \strong{Calibration.} The dependence strength is found by inverting a
#' simulated grid, so the ribbon is centred on the observed coefficient only to
#' within Monte Carlo error, on the order of .01 at the default
#' `cal_iterations`. Raise `cal_iterations` if a pair's `expected` at the
#' middle bins is visibly off its pooled coefficient.
#'
#' \strong{What it is good at.} Not the uniform-or-not verdict, which is
#' poorly powered at realistic n, but showing that a pair's dependence is
#' concentrated at one end of the scale. That is actionable: it points at a
#' floor or ceiling artefact, or at a response process shared only among low
#' scorers. A profile crossing zero is the mixed-sign case where the pooled
#' coefficient is not a meaningful measure of partial association at all.
#'
#' @return A ggplot (faceted when several pairs are supplied), or a data.frame
#'   with one row per pair and class interval. Both carry a `"pairs"`
#'   attribute holding the pooled coefficient, the chosen direction, the
#'   calibrated dependence strength and the weight-weighted slope of stratum
#'   gamma on rest score.
RMlocdepGammaProfile <- function(
  data,
  pairs,
  class_intervals = 5,
  direction = c("auto", "1", "2"),
  band = c("uniform", "none"),
  band_type = c("simultaneous", "pointwise"),
  conf_level = 0.95,
  iterations = 400,
  cal_iterations = 40,
  min_n = 10,
  x_axis = c("theta", "restscore"),
  parallel = TRUE,
  n_cores = NULL,
  seed = NULL,
  verbose = FALSE,
  output = c("ggplot", "dataframe")
) {
  direction <- match.arg(direction)
  band      <- match.arg(band)
  band_type <- match.arg(band_type)
  x_axis    <- match.arg(x_axis)
  output    <- match.arg(output)

  if (output == "ggplot" && !requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\".", call. = FALSE)
  }
  if (class_intervals < 2L) {
    stop("`class_intervals` must be at least 2.", call. = FALSE)
  }

  easyRasch2:::validate_response_data(data)

  n_total      <- nrow(data)
  has_na       <- anyNA(data)
  complete_idx <- stats::complete.cases(data)
  data         <- data[complete_idx, , drop = FALSE]
  if (nrow(data) == 0L) {
    stop("No complete cases in data after removing rows with NA.", call. = FALSE)
  }

  X <- as.matrix(data)
  storage.mode(X) <- "integer"
  if (is.null(colnames(X))) colnames(X) <- paste0("V", seq_len(ncol(X)))
  item_names <- colnames(X)
  n <- nrow(X)

  pair_list <- .lp_normalise_pairs(pairs, item_names)

  # `m` is held at the observed number of categories for every dataset, so the
  # observed and simulated coefficients come out of one code path.
  m  <- max(X) + 1L
  G  <- outer(seq_len(m), seq_len(m), function(a, b) as.numeric(b > a))
  tG <- t(G)

  if (!is.null(seed)) set.seed(seed)

  if (is.null(n_cores)) n_cores <- getOption("mc.cores", 1L)
  cores <- if (isTRUE(parallel)) max(1L, as.integer(n_cores)) else 1L

  # Generating model for the ribbon: CML thresholds plus WLE person locations,
  # the same pair the rest of the package bootstraps from. Row-aligned thetas
  # are kept separately for the x axis, since the pool drops non-finite ones.
  need_sim <- band == "uniform"
  thr_list <- easyRasch2:::.fit_cml_thresholds(X)
  theta_row <- easyRasch2:::.estimate_thetas(X, thr_list, method = "WLE")$theta
  theta_row[!is.finite(theta_row)] <- NA_real_
  theta_pool <- theta_row[is.finite(theta_row)]

  rows <- list()
  meta <- list()

  for (pn in names(pair_list)) {
    pr <- pair_list[[pn]]
    a  <- pr[1L]
    b  <- pr[2L]

    # --- direction -----------------------------------------------------------
    st1 <- .lp_strata(X, a, b, m, G, tG)   # rest score excludes item b
    st2 <- .lp_strata(X, b, a, m, G, tG)   # rest score excludes item a
    use_1 <- switch(
      direction,
      auto = isTRUE(abs(attr(st1, "gamma")) >= abs(attr(st2, "gamma"))),
      `1`  = TRUE,
      `2`  = FALSE
    )
    st <- if (use_1) st1 else st2
    i  <- if (use_1) a else b
    j  <- if (use_1) b else a
    g_obs <- attr(st, "gamma")
    if (!is.finite(g_obs)) {
      warning("Pair ", pn, " has no comparable pairs; skipped.", call. = FALSE)
      next
    }

    z_row <- rowSums(X) - X[, j]

    # --- class intervals -----------------------------------------------------
    # Breaks are formed once, on the observed rest score, and reused for every
    # simulated dataset so that observed and ribbon are binned identically.
    # The outer breaks are open so a simulated rest score outside the observed
    # range is kept rather than dropped.
    qs <- unique(stats::quantile(
      z_row,
      probs = seq(0, 1, length.out = class_intervals + 1L),
      type = 1L, names = FALSE
    ))
    if (length(qs) < 3L) {
      warning("Pair ", pn, " has too few distinct rest scores to bin; skipped.",
              call. = FALSE)
      next
    }
    brk <- c(-Inf, qs[-c(1L, length(qs))], Inf)
    labs <- levels(cut(z_row, breaks = brk, include.lowest = TRUE))
    nb   <- length(labs)

    obs_gamma <- .lp_bin(st, brk, labs)
    bin_row   <- cut(z_row, breaks = brk, include.lowest = TRUE)
    n_bin     <- as.integer(table(bin_row)[labs])
    x_bin <- if (x_axis == "theta") {
      as.numeric(tapply(theta_row, bin_row, mean, na.rm = TRUE)[labs])
    } else {
      as.numeric(tapply(z_row, bin_row, mean, na.rm = TRUE)[labs])
    }

    # Weight-weighted slope of stratum gamma on rest score. Not drawn, but the
    # inferential version of this diagnostic is a test on exactly this number,
    # and it is free to return.
    ok_k  <- is.finite(st$gamma_k) & st$weight_k > 0
    slope <- if (sum(ok_k) >= 3L) {
      stats::coef(stats::lm(gamma_k ~ rest, data = st[ok_k, ],
                            weights = st$weight_k[ok_k]))[["rest"]]
    } else NA_real_

    d_cal <- NA_real_
    clamped <- FALSE
    lo <- hi <- med <- rep(NA_real_, nb)

    if (need_sim) {
      if (verbose) message("Calibrating uniform-LD strength for ", pn, " ...")
      cal <- .lp_calibrate(thr_list, theta_pool, n, c(a, b), i, j, m, G, tG,
                           target = g_obs, R = cal_iterations, cores = cores)
      d_cal   <- cal$d
      clamped <- cal$clamped
      if (clamped) {
        warning("Pair ", pn, ": observed gamma (", round(g_obs, 3),
                ") lies outside the range the calibration grid reached; ",
                "the ribbon uses the nearest attainable strength.",
                call. = FALSE)
      }

      if (verbose) message("Simulating ", iterations, " uniform-LD datasets ...")
      prof <- .lp_par(seq_len(iterations), function(r) {
        Xs <- .lp_sim_uniform(thr_list, theta_pool, n, c(a, b), d_cal)
        ss <- .lp_strata(Xs, i, j, m, G, tG)
        .lp_bin(ss, brk, labs)
      }, cores)
      M <- do.call(rbind, prof)

      env <- .lp_envelope(M, band_type, conf_level)
      med <- env$med
      lo  <- env$lo
      hi  <- env$hi
    }

    # A class interval is drawable only when it holds enough respondents, the
    # observed coefficient is estimable (a stratum can hold people yet no
    # comparable pair, which is common at the score extremes) and the
    # reference is not degenerate.
    keep <- n_bin >= min_n & is.finite(obs_gamma) &
      (band == "none" | (is.finite(lo) & is.finite(hi) & hi - lo < 1.98))
    rows[[pn]] <- data.frame(
      pair      = pn,
      interval  = labs,
      x         = x_bin,
      n         = n_bin,
      gamma     = obs_gamma,
      expected  = med,
      band_lo   = lo,
      band_hi   = hi,
      outside   = is.finite(obs_gamma) & is.finite(lo) &
                  (obs_gamma < lo | obs_gamma > hi),
      drawn     = keep,
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    meta[[pn]] <- data.frame(
      pair = pn,
      item1 = item_names[a], item2 = item_names[b],
      # The rest score conditioned on excludes one item of the pair, and which
      # one is the conditioning direction.
      rest_excludes = item_names[j],
      direction = if (use_1) 1L else 2L,
      gamma = g_obs, slope = slope, d = d_cal, clamped = clamped,
      n_drawn = sum(keep),
      stringsAsFactors = FALSE, row.names = NULL
    )
  }

  if (length(rows) == 0L) stop("No pair could be profiled.", call. = FALSE)

  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  info <- do.call(rbind, meta)
  attr(out, "pairs") <- info

  sample_clause <- easyRasch2:::.n_caption(
    n, n_total, if (has_na) "complete cases" else character()
  )

  if (output == "dataframe") {
    attr(out, "sample_clause") <- sample_clause
    return(out)
  }

  # An omitted interval is a finding, not a tidying step: it says the scale
  # ran out of comparable pairs there. Name the count rather than leaving a
  # gap the reader has to notice.
  n_drop <- sum(!out$drawn)
  drop_clause <- if (n_drop > 0L) {
    sprintf(paste0(
      " %d of %d class intervals are omitted, holding too few respondents ",
      "or no comparable response pairs."), n_drop, nrow(out))
  } else ""

  .lp_plot(out, info, band, band_type, conf_level, iterations, x_axis,
           sample_clause, drop_clause)
}


# ---------------------------------------------------------------------------
# Internals
# ---------------------------------------------------------------------------

#' Normalise the `pairs` argument to a named list of column indices
#' @noRd
.lp_normalise_pairs <- function(pairs, item_names) {
  if (!is.list(pairs)) pairs <- list(pairs)
  out <- lapply(pairs, function(p) {
    if (length(p) != 2L) {
      stop("Each pair must have exactly two items.", call. = FALSE)
    }
    idx <- if (is.character(p)) match(p, item_names) else as.integer(p)
    if (anyNA(idx) || any(idx < 1L) || any(idx > length(item_names))) {
      stop("Unknown item in pair: ", paste(p, collapse = ", "), call. = FALSE)
    }
    if (idx[1L] == idx[2L]) {
      stop("A pair must name two different items.", call. = FALSE)
    }
    idx
  })
  names(out) <- vapply(out, function(p) {
    paste(item_names[p], collapse = " - ")
  }, character(1L))
  out
}

#' Stratum-wise partial gamma for one pair and one conditioning direction
#'
#' Wraps the package internal so the strata come back labelled by rest score.
#' `.partgam_one()` indexes strata from `min(z)`, so the labels are recovered
#' from that origin.
#' @noRd
.lp_strata <- function(X, i, j, m, G, tG) {
  z  <- rowSums(X) - X[, j]
  st <- easyRasch2:::.partgam_one(X[, i], X[, j], z, m, G, tG, strata = TRUE)
  out <- data.frame(
    rest     = seq.int(min(z), max(z)),
    gamma_k  = st$gamma_k,
    weight_k = st$weight_k,
    n_k      = st$n_k
  )
  attr(out, "gamma") <- st$gamma
  out
}

#' Aggregate stratum coefficients into class intervals
#'
#' The weighted form, never a coefficient recomputed on collapsed data.
#' @noRd
.lp_bin <- function(st, brk, labs) {
  bin <- cut(st$rest, breaks = brk, include.lowest = TRUE)
  num <- tapply(st$gamma_k * st$weight_k, bin, sum, na.rm = TRUE)
  den <- tapply(st$weight_k, bin, sum, na.rm = TRUE)
  g <- as.numeric(num[labs]) / as.numeric(den[labs])
  g[!is.finite(g)] <- NA_real_
  g
}

#' PCM draw by inverse transform, vectorised over persons
#'
#' Accepts a person-specific location vector, which the dependence mechanism
#' needs, and consumes one uniform per person.
#' @noRd
.lp_pcm_draw <- function(theta, deltas, u = stats::runif(length(theta))) {
  m  <- length(deltas)
  cd <- cumsum(deltas)
  nn <- length(theta)
  num <- matrix(0, nn, m + 1L)
  for (x in seq_len(m)) num[, x + 1L] <- x * theta - cd[x]
  num <- num - apply(num, 1L, max)          # stabilise before exponentiating
  pr <- exp(num)
  pr <- pr / rowSums(pr)
  for (x in 2L:(m + 1L)) pr[, x] <- pr[, x - 1L] + pr[, x]
  as.integer(pmin(rowSums(pr < u), m))
}

#' Simulate one dataset with uniform local dependence in one pair
#'
#' Response dependence in the shift form: the second item is answered at a
#' location displaced by `d * z`, where `z` maps the first item's response onto
#' [-1, 1]. The shift is linear in the observed category and does not depend on
#' theta, so the association is constant both across the score range and across
#' the category grid. That is what "uniform" means here, and `d = 0` reduces to
#' the plain Rasch/PCM null.
#' @noRd
.lp_sim_uniform <- function(thr_list, theta_pool, n, pair, d) {
  th <- sample(theta_pool, n, replace = TRUE)
  k  <- length(thr_list)
  X  <- matrix(0L, n, k)
  for (a in setdiff(seq_len(k), pair)) {
    X[, a] <- .lp_pcm_draw(th, thr_list[[a]])
  }
  i <- pair[1L]
  j <- pair[2L]
  X[, i] <- .lp_pcm_draw(th, thr_list[[i]])
  mi <- length(thr_list[[i]])
  z  <- (2 * X[, i] - mi) / mi
  X[, j] <- .lp_pcm_draw(th + d * z, thr_list[[j]])
  X
}

#' Find the dependence strength reproducing the observed pooled gamma
#'
#' A grid is evaluated, made monotone with `isoreg()` so Monte Carlo wobble
#' cannot produce a non-invertible relation, then inverted at the target. The
#' grid spans both signs, since a negatively associated pair is a real case and
#' the mechanism is not quite symmetric in `d`.
#' @noRd
.lp_calibrate <- function(thr_list, theta_pool, n, pair, i, j, m, G, tG,
                          target, R, cores, max_d = 2.5, n_grid = 9L) {
  grid <- seq(-max_d, max_d, length.out = n_grid)
  means <- vapply(grid, function(d) {
    g <- unlist(.lp_par(seq_len(R), function(r) {
      Xs <- .lp_sim_uniform(thr_list, theta_pool, n, pair, d)
      easyRasch2:::.partgam_one(Xs[, i], Xs[, j],
                                rowSums(Xs) - Xs[, j], m, G, tG)
    }, cores))
    mean(g, na.rm = TRUE)
  }, numeric(1L))

  yf   <- stats::isoreg(grid, means)$yf
  keep <- !duplicated(yf)
  if (target <= min(yf)) return(list(d = grid[1L], means = means, clamped = TRUE))
  if (target >= max(yf)) {
    return(list(d = grid[length(grid)], means = means, clamped = TRUE))
  }
  list(d = stats::approx(yf[keep], grid[keep], xout = target)$y,
       means = means, clamped = FALSE)
}

#' Pointwise or simultaneous envelope over simulated profiles
#'
#' The simultaneous band widens the pointwise quantiles to a common level until
#' `conf_level` of the simulated profiles lie entirely inside it. Working on the
#' empirical CDF scale rather than as median plus or minus a scale keeps the
#' band inside [-1, 1] and lets it be asymmetric, which matters in the sparse
#' end bins where the simulated coefficients pile up against a bound. Same idea
#' as a Westfall-Young max statistic, taken across class intervals rather than
#' across item pairs.
#' @noRd
.lp_envelope <- function(M, band_type, conf_level) {
  med <- apply(M, 2L, stats::median, na.rm = TRUE)
  a   <- (1 - conf_level) / 2
  if (band_type == "pointwise") {
    plo <- a
    phi <- 1 - a
  } else {
    U <- apply(M, 2L, function(v) {
      ok <- is.finite(v)
      u  <- rep(NA_real_, length(v))
      if (any(ok)) {
        u[ok] <- (rank(v[ok], ties.method = "average") - 0.5) / sum(ok)
      }
      u
    })
    dev <- suppressWarnings(apply(abs(U - 0.5), 1L, max, na.rm = TRUE))
    dev <- dev[is.finite(dev)]
    q <- if (length(dev)) {
      stats::quantile(dev, conf_level, names = FALSE)
    } else 0.5
    plo <- max(0, 0.5 - q)
    phi <- min(1, 0.5 + q)
  }
  lo <- apply(M, 2L, stats::quantile, probs = plo, na.rm = TRUE, names = FALSE)
  hi <- apply(M, 2L, stats::quantile, probs = phi, na.rm = TRUE, names = FALSE)
  list(med = med, lo = pmax(lo, -1), hi = pmin(hi, 1))
}

#' mclapply where available, lapply otherwise, with a parallel-safe RNG
#' @noRd
.lp_par <- function(X, FUN, cores) {
  if (cores > 1L && .Platform$OS.type != "windows" &&
      requireNamespace("parallel", quietly = TRUE)) {
    old <- RNGkind("L'Ecuyer-CMRG")
    on.exit(RNGkind(old[1L]), add = TRUE)
    parallel::mclapply(X, FUN, mc.cores = cores, mc.set.seed = TRUE)
  } else {
    lapply(X, FUN)
  }
}

#' Assemble the figure
#' @noRd
.lp_plot <- function(out, info, band, band_type, conf_level, iterations,
                     x_axis, sample_clause, drop_clause = "") {
  df <- out[out$drawn, , drop = FALSE]
  lab <- vapply(seq_len(nrow(info)), function(r) {
    sprintf("%s  (gamma = %.2f, rest score excludes %s)",
            info$pair[r], info$gamma[r], info$rest_excludes[r])
  }, character(1L))
  names(lab) <- info$pair
  df$panel <- factor(lab[df$pair], levels = lab)

  p <- ggplot2::ggplot()
  if (band == "uniform") {
    p <- p + ggplot2::geom_ribbon(
      data = df,
      ggplot2::aes(x = .data$x, ymin = .data$band_lo, ymax = .data$band_hi),
      fill = "grey85", alpha = 0.6
    ) +
      ggplot2::geom_line(
        data = df, ggplot2::aes(x = .data$x, y = .data$expected),
        colour = "black", linewidth = 0.8
      )
  }
  p <- p +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_line(
      data = df, ggplot2::aes(x = .data$x, y = .data$gamma),
      colour = "sienna2", linewidth = 0.5
    ) +
    ggplot2::geom_point(
      data = df, ggplot2::aes(x = .data$x, y = .data$gamma),
      shape = 18, size = 2.9, colour = "sienna2"
    ) +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::labs(
      x = if (x_axis == "theta") "Person location (logits)" else "Rest score",
      y = "Partial gamma",
      caption = easyRasch2:::er2_caption(paste0(
        if (band == "uniform") {
          paste0(
            "Ribbon: ", round(100 * conf_level), "% ",
            if (band_type == "simultaneous") "simultaneous" else "pointwise",
            " reference from ", iterations,
            " datasets simulated under uniform local dependence, ",
            "calibrated to the observed partial gamma of the pair. ",
            "The line is the simulated median. "
          )
        } else "",
        "Coefficients are computed within single rest-score levels and ",
        "aggregated into class intervals. ", sample_clause, ".",
        drop_clause
      ))
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey95",
                                               colour = "grey70"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.spacing = ggplot2::unit(0.7, "cm")
    ) +
    easyRasch2:::er2_axis_margins() +
    easyRasch2:::er2_plot_caption()
  p
}
