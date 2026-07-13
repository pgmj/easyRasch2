#' Person-Fit Statistics for a Rasch / Partial Credit Model
#'
#' Computes per-respondent person-fit statistics and resampling-based
#' p-values. Three statistics are reported: conditional **infit** and
#' **outfit** mean-squares (MSQ) -- the magnitude (effect-size) measures
#' familiar from Rasch analysis -- and the standardized log-likelihood
#' statistic **lz**. The conditional MSQ statistics use response
#' probabilities conditional on the total score, which require no person
#' estimate and are therefore unbiased (Kreiner & Christensen, 2011); they are
#' computed on each person's observed response pattern, so partial
#' missingness is handled directly.
#'
#' Statistical significance is assessed by Monte-Carlo resampling under
#' the fitted model rather than by an assumed asymptotic distribution.
#' The asymptotic null distributions of MSQ and lz are known to be
#' unreliable -- the Wilson-Hilferty (ZSTD) transformation of MSQ in
#' particular adds nothing once conditional estimation is used (Müller,
#' 2020) -- so resampling provides the valid reference (Sinharay, 2016).
#'
#' @param data A data.frame or matrix of item responses. Items must be
#'   scored starting at 0 (non-negative integers). Missing values (`NA`)
#'   are allowed.
#' @param statistics Character vector. Which statistics to report; any of
#'   `"infit"`, `"outfit"`, `"lz"`. Defaults to all three.
#' @param estimator Character. How item parameters are estimated:
#'   `"CML"` (default) or `"MML"` (via \pkg{mirt}).
#' @param theta_method Character. Person-location estimator used for the
#'   lz statistic: `"WLE"` (default) or `"EAP"`. Ignored if `"lz"` is not
#'   requested.
#' @param iterations Integer. Number of Monte-Carlo replications per
#'   person. `0` skips resampling and returns statistics only. Default
#'   `500`.
#' @param flag_alpha Numeric in (0, 1). Significance level for the `flagged`
#'   column and the plot highlighting. Default `0.05`.
#' @param flag Character. Which misfit direction drives flagging for the MSQ
#'   statistics: `"both"` (default) flags either direction with a two-sided
#'   p-value; `"underfit"` flags only underfit (MSQ > 1, noisy responding)
#'   with a one-sided upper p-value, which is more powerful for detecting it
#'   and ignores (benign) overfit. The reported `p_infit` / `p_outfit` follow
#'   this choice. lz is one-sided regardless. Overfit can occasionally be
#'   suspicious (e.g. fabricated or copied response patterns), so `"both"` is
#'   the default.
#' @param zstd Logical. If `TRUE`, also report the Wilson-Hilferty
#'   standardized (ZSTD) versions of the MSQ statistics. These are
#'   provided only as a familiar cross-person comparability metric and
#'   are **not** used for inference (see Müller, 2020). Default `FALSE`.
#' @param parallel,n_cores Logical / integer. Parallelise the resampling
#'   across persons via \pkg{mirai} when available. Default sequential.
#' @param seed Optional integer for reproducible resampling.
#' @param output Character. `"kable"` (default), `"dataframe"`,
#'   `"ggplot"`, or `"list"`. For `"ggplot"`, a **named list** of
#'   person-fit maps is returned -- one per requested statistic (e.g.
#'   `$infit`, `$outfit`, `$lz`) -- each plotting the statistic against
#'   person location, with
#'   respondents flagged by *that* statistic highlighted and a caption
#'   reporting the assessed sample size and the proportion flagged. For the
#'   MSQ maps the flagged count is split into underfit (MSQ > 1, noisy/erratic
#'   responding) and overfit (MSQ < 1, overly deterministic responding).
#'
#' @return
#' For `output = "dataframe"`, a data.frame with one row per respondent
#' (input order): `id`, `n_answered`, `sum_score`, the requested
#' statistics (`infit_msq`, `outfit_msq`, `lz`), their resampled
#' p-values (`p_infit`, `p_outfit`, `p_lz`) when `iterations > 0`,
#' `flagged`, and -- if `zstd = TRUE` -- `infit_zstd`, `outfit_zstd`.
#' The `flagged` column is `TRUE` when the marginal (uncorrected) p-value of
#' **any** requested statistic is below `flag_alpha` for that respondent; it
#' is therefore a per-person screening flag, not corrected for the number of
#' respondents tested (so under fit expect about `flag_alpha` of respondents
#' to be flagged by chance). Each `output = "ggplot"` map instead colours and
#' counts respondents flagged by its **own** statistic alone. Extreme scorers
#' (minimum or maximum possible given the items they answered) receive `NA`
#' statistics and are not assessed. For `output = "kable"` the same content as
#' a `knitr_kable`; for `output = "ggplot"` the named list of maps described
#' under that argument; for `output = "list"`, a list with both views from a
#' single computation: `fit` (the data.frame) and `plots` (the named list of
#' maps).
#'
#' @details
#' \strong{Conditional MSQ.} For person \eqn{v} with answered items
#' \eqn{A_v} and total score \eqn{r_v}, the conditional residual is
#' \eqn{z_{vi} = (x_{vi} - E(X_{vi} \mid R_v = r_v)) /
#' \sqrt{\mathrm{Var}(X_{vi} \mid R_v = r_v)}}, with conditional moments
#' obtained from elementary symmetric functions of the answered items'
#' parameters. Outfit is the unweighted mean of \eqn{z_{vi}^2}; infit is
#' the information-weighted mean. Because the conditional moments use only
#' item parameters, no (biased) person estimate enters the residual.
#'
#' \strong{Interpreting MSQ direction.} Values near 1 indicate fit. MSQ > 1
#' (\emph{underfit}) means the responses are noisier / more erratic than the
#' model expects -- the validity-relevant direction, indicating careless or
#' aberrant responding. MSQ < 1 (\emph{overfit}) means responses are overly
#' deterministic (Guttman-like); usually benign for score validity, though a
#' suspiciously perfect pattern can occasionally signal copied or fabricated
#' data. Use `flag = "underfit"` to flag only the former. lz is one-sided:
#' low (negative) values flag the aberrant (underfit-like) direction.
#'
#' \strong{lz.} The standardized log-likelihood of the response pattern
#' evaluated at the estimated person location (Drasgow, Levine &
#' Williams, 1985); small (negative) values indicate misfit. Because the
#' location is estimated rather than known, the variance of lz is below 1
#' and its asymptotic standard-normal null is invalid, which makes a
#' naive test conservative (Snijders, 2001; Sinharay, 2016). Instead of
#' applying the analytic lz\* standardization -- derived by Snijders (2001)
#' for dichotomous items and extended to polytomous / mixed-format items by
#' Sinharay (2016, 2026) -- easyRasch2 obtains the reference distribution by
#' resampling with the person location \strong{re-estimated for every
#' simulated pattern} (see \strong{Resampling}). This reproduces the
#' ability-estimation effect and is the resampling analogue of lz\*, so the
#' resampled p-value is well calibrated.
#'
#' \strong{Resampling.} Each statistic is referenced against the null that
#' matches it, removing the need to choose a scheme. The conditional MSQ
#' statistics use patterns sampled **conditional on the person's total
#' score** -- the exact Rasch-native null, requiring no person estimate and
#' fully consistent with the (conditional) statistic. lz, which is defined at
#' the estimated location, uses patterns **simulated at that location and then
#' scored with the location re-estimated from each simulated pattern**, so the
#' null carries the same ability-estimation effect as the observed lz (for a
#' Rasch / PCM model the location is a function of the total score, so this
#' re-estimation is a cheap score-based lookup). Both schemes are computed per
#' person over the items actually answered, so partial missingness is handled
#' by either. The p-value is the proportion of the `iterations` replicates at
#' least as extreme as the observed value. This follows the resampling-based
#' person-fit approach (Sinharay, 2016) and the bootstrap recommendation of
#' Müller (2020).
#'
#' @references
#' de la Torre, J., & Deng, W. (2008). Improving person-fit assessment by
#' correcting the ability estimate and its reference distribution.
#' *Journal of Educational Measurement, 45*(2), 159-177.
#' \doi{10.1111/j.1745-3984.2008.00058.x}
#'
#' Drasgow, F., Levine, M. V., & Williams, E. A. (1985). Appropriateness
#' measurement with polychotomous item response models and standardized
#' indices. *British Journal of Mathematical and Statistical Psychology,
#' 38*(1), 67-86. \doi{10.1111/j.2044-8317.1985.tb00817.x}
#'
#' Kreiner, S., & Christensen, K. B. (2011). Exact evaluation of bias in
#' Rasch model residuals. *Advances in Mathematics Research, 12*,
#' 19-40.
#'
#' Müller, M. (2020). Item fit statistics for Rasch analysis: can we
#' trust them? *Journal of Statistical Distributions and Applications,
#' 7*(5). \doi{10.1186/s40488-020-00108-7}
#'
#' Sinharay, S. (2016). Assessment of person fit using resampling-based
#' approaches. *Journal of Educational Measurement, 53*(1), 63-85.
#' \doi{10.1111/jedm.12101}
#'
#' Sinharay, S. (2016). Asymptotically correct standardization of
#' person-fit statistics beyond dichotomous items. *Psychometrika,
#' 81*(4), 992-1013. \doi{10.1007/s11336-015-9465-x}
#'
#' Sinharay, S. (2026). Refining the asymptotically correct
#' standardization of person-fit statistics for mixed-format tests.
#' *British Journal of Mathematical and Statistical Psychology*.
#' \doi{10.1111/bmsp.70049}
#'
#' Snijders, T. A. B. (2001). Asymptotic null distribution of person fit
#' statistics with estimated person parameter. *Psychometrika, 66*(3),
#' 331-342. \doi{10.1007/BF02294440}
#'
#' @seealso [RMpersonParameters()], [RMitemInfit()], [RMitemParameters()]
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' dat <- as.data.frame(
#'   matrix(sample(0:2, 200 * 8, replace = TRUE), nrow = 200, ncol = 8)
#' )
#' colnames(dat) <- paste0("Item", 1:8)
#'
#' # Conditional infit/outfit MSQ + lz with resampled p-values
#' RMpersonFit(dat, iterations = 200, output = "dataframe") |> head()
#'
#' # Person-fit maps: a named list with one plot per statistic
#' if (requireNamespace("ggplot2", quietly = TRUE)) {
#'   plots <- RMpersonFit(dat, iterations = 200, output = "ggplot")
#'   plots$infit    # infit map (each plot's caption reports the % flagged)
#' }
#'
#' # Flag only underfit (noisy responding), the validity-relevant direction
#' RMpersonFit(dat, iterations = 200, flag = "underfit",
#'             output = "dataframe") |> head()
#' }
RMpersonFit <- function(
  data,
  statistics = c("infit", "outfit", "lz"),
  estimator = c("CML", "MML"),
  theta_method = c("WLE", "EAP"),
  iterations = 500L,
  flag_alpha = 0.05,
  flag = c("both", "underfit"),
  zstd = FALSE,
  parallel = FALSE,
  n_cores = NULL,
  seed = NULL,
  output = c("kable", "dataframe", "ggplot", "list")
) {
  statistics <- match.arg(
    statistics,
    c("infit", "outfit", "lz"),
    several.ok = TRUE
  )
  estimator <- match.arg(estimator)
  theta_method <- match.arg(theta_method)
  flag <- match.arg(flag)
  output <- match.arg(output)
  # Underfit-only flagging uses a one-sided (upper) MSQ test; lz is always
  # one-sided (low = aberrant) and is unaffected.
  msq_tail <- if (flag == "underfit") "upper" else "two.sided"

  validate_response_data(data)
  if (
    !is.numeric(flag_alpha) ||
      length(flag_alpha) != 1L ||
      flag_alpha <= 0 ||
      flag_alpha >= 1
  ) {
    stop("`flag_alpha` must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }

  data <- as.data.frame(data)
  n_total <- nrow(data)
  has_na <- anyNA(data)
  data <- .drop_empty_respondents(data)

  # --- Item parameters --------------------------------------------------------
  fit <- if (estimator == "CML") {
    .rasch_fit_cml(data, se = FALSE)
  } else {
    .rasch_fit_mml(data, se = FALSE)
  }
  # Grand-mean-zero identification so the lz location matches
  # RMpersonParameters(); conditional MSQ is invariant to this shift.
  thr_list <- .center_thresholds(fit$thr_list)
  data <- .align_items(data, thr_list)
  data_mat <- as.matrix(data)
  n <- nrow(data_mat)

  # --- Person bookkeeping -----------------------------------------------------
  n_answered <- rowSums(!is.na(data_mat))
  sum_score <- rowSums(data_mat, na.rm = TRUE)
  max_score <- vapply(
    seq_len(n),
    function(p) {
      ans <- !is.na(data_mat[p, ])
      sum(vapply(thr_list[ans], length, integer(1)))
    },
    numeric(1)
  )
  extreme <- n_answered > 0L & (sum_score == 0 | sum_score == max_score)

  # --- Person location -------------------------------------------------------
  # Needed for lz (both the statistic and its theta-based resampling) and for
  # the ggplot x-axis. The MSQ statistics and their conditional resampling
  # need no person estimate.
  need_theta <- ("lz" %in% statistics) || (output %in% c("ggplot", "list"))
  theta <- rep(NA_real_, n)
  if (need_theta) {
    # Shared engine (utils-theta.R): WLE caches by (answered-set, score); EAP
    # uses its quadrature grid and estimates the prior SD from the data. The
    # wider c(-10, 10) range locates extreme scorers (person-reporting scale).
    theta <- .estimate_thetas(
      data_mat,
      thr_list,
      method = theta_method,
      theta_range = c(-10, 10)
    )$theta
  }

  # --- Observed statistics ----------------------------------------------------
  obs <- .person_fit_stats(data_mat, thr_list, theta, statistics)

  out <- data.frame(
    id = if (!is.null(rownames(data))) rownames(data) else seq_len(n),
    n_answered = n_answered,
    sum_score = sum_score,
    stringsAsFactors = FALSE
  )
  if ("infit" %in% statistics) {
    out$infit_msq <- obs$infit
  }
  if ("outfit" %in% statistics) {
    out$outfit_msq <- obs$outfit
  }
  if ("lz" %in% statistics) {
    out$lz <- obs$lz
  }

  # --- Resampling p-values ----------------------------------------------------
  if (iterations > 0L) {
    pv <- .person_fit_resample(
      data_mat,
      thr_list,
      theta,
      statistics,
      obs,
      iterations = iterations,
      theta_method = theta_method,
      parallel = parallel,
      n_cores = n_cores,
      msq_tail = msq_tail
    )
    if ("infit" %in% statistics) {
      out$p_infit <- pv$infit
    }
    if ("outfit" %in% statistics) {
      out$p_outfit <- pv$outfit
    }
    if ("lz" %in% statistics) {
      out$p_lz <- pv$lz
    }

    p_cols <- intersect(c("p_infit", "p_outfit", "p_lz"), names(out))
    out$flagged <- apply(out[p_cols], 1L, function(p) {
      any(p < flag_alpha, na.rm = TRUE)
    })
  }

  # --- Optional ZSTD (non-inferential) ----------------------------------------
  if (zstd) {
    if ("infit" %in% statistics) {
      out$infit_zstd <- .msq_to_zstd(obs$infit, n_answered)
    }
    if ("outfit" %in% statistics) {
      out$outfit_zstd <- .msq_to_zstd(obs$outfit, n_answered)
    }
  }

  # Extreme scorers carry no fit information
  stat_cols <- intersect(
    c(
      "infit_msq",
      "outfit_msq",
      "lz",
      "p_infit",
      "p_outfit",
      "p_lz",
      "infit_zstd",
      "outfit_zstd"
    ),
    names(out)
  )
  out[extreme, stat_cols] <- NA
  if ("flagged" %in% names(out)) {
    out$flagged[extreme] <- NA
  }

  .person_fit_output(
    out,
    theta,
    extreme,
    statistics,
    output,
    flag_alpha,
    iterations,
    flag,
    n_total = n_total,
    has_na = has_na
  )
}

# ---------------------------------------------------------------------
# Internal: conditional moments by total score (missingness-aware core)
# ---------------------------------------------------------------------

#' Conditional expected score and variance by total score
#'
#' For a set of items (a single missingness pattern), returns lookup
#' matrices of \eqn{E(X_i \mid R = r)} and \eqn{\mathrm{Var}(X_i \mid
#' R = r)} for every attainable total score \eqn{r}, computed from
#' elementary symmetric functions of the item parameters.
#'
#' @param thr_sub List of threshold vectors for the answered items.
#' @return List with `E` and `V` matrices (rows = total score 0..max,
#'   cols = items) and `max_score`.
#' @keywords internal
#' @noRd
.cond_moments <- function(thr_sub) {
  steps <- vapply(thr_sub, length, integer(1))
  is_poly <- any(steps > 1L)
  k <- length(thr_sub)
  max_sc <- sum(steps)

  pars <- if (is_poly) lapply(thr_sub, cumsum) else unlist(thr_sub)
  esf <- psychotools::elementary_symmetric_functions(pars, order = 1L)
  g <- esf[[1L]]
  dg <- esf[[2L]]

  E <- matrix(NA_real_, max_sc + 1L, k)
  V <- matrix(NA_real_, max_sc + 1L, k)

  if (!is_poly) {
    for (r in seq_len(max_sc - 1L)) {
      # 1..max-1 (extremes degenerate)
      e <- dg[r + 1L, ] / g[r + 1L]
      E[r + 1L, ] <- e
      V[r + 1L, ] <- e * (1 - e)
    }
  } else {
    col0 <- c(0L, cumsum(steps))
    for (r in seq_len(max_sc - 1L)) {
      for (i in seq_len(k)) {
        cols <- (col0[i] + 1L):col0[i + 1L]
        p_ge <- dg[r + 1L, cols] / g[r + 1L] # P(X_i = c | R = r), c = 1..m_i
        p <- c(1 - sum(p_ge), p_ge)
        cats <- 0:steps[i]
        ev <- sum(cats * p)
        E[r + 1L, i] <- ev
        V[r + 1L, i] <- sum(cats^2 * p) - ev^2
      }
    }
  }
  list(E = E, V = V, max_score = max_sc)
}

#' Conditional residuals for one response matrix (grouped by NA pattern)
#'
#' @param data_mat Numeric response matrix.
#' @param thr_list Threshold list aligned with the columns.
#' @return List of length nrow with, per person, the conditional
#'   residuals `z` and variance weights `w` for their answered items;
#'   `NULL` for extreme/empty rows.
#' @keywords internal
#' @noRd
.conditional_residuals <- function(data_mat, thr_list) {
  n <- nrow(data_mat)
  res <- vector("list", n)

  pattern <- apply(!is.na(data_mat), 1L, function(a) {
    paste0(which(a), collapse = ",")
  })
  for (pat in unique(pattern)) {
    rows <- which(pattern == pat)
    items <- as.integer(strsplit(pat, ",", fixed = TRUE)[[1L]])
    if (length(items) == 0L) {
      next
    }
    mom <- .cond_moments(thr_list[items])

    for (p in rows) {
      x <- data_mat[p, items]
      r <- sum(x)
      if (r == 0L || r == mom$max_score) {
        next
      } # extreme: no information
      e <- mom$E[r + 1L, ]
      v <- mom$V[r + 1L, ]
      z <- (x - e) / sqrt(v)
      res[[p]] <- list(z = z, w = v)
    }
  }
  res
}

# ---------------------------------------------------------------------
# Internal: observed statistics
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.person_fit_stats <- function(data_mat, thr_list, theta, statistics) {
  n <- nrow(data_mat)
  infit <- outfit <- lz <- rep(NA_real_, n)

  if (any(c("infit", "outfit") %in% statistics)) {
    resid <- .conditional_residuals(data_mat, thr_list)
    for (p in seq_len(n)) {
      r <- resid[[p]]
      if (is.null(r)) {
        next
      }
      outfit[p] <- mean(r$z^2)
      infit[p] <- sum(r$w * r$z^2) / sum(r$w)
    }
  }
  if ("lz" %in% statistics) {
    for (p in seq_len(n)) {
      lz[p] <- .person_lz(data_mat[p, ], thr_list, theta[p])
    }
  }
  list(infit = infit, outfit = outfit, lz = lz)
}

#' Standardized log-likelihood (lz) for one response pattern at theta
#' @keywords internal
#' @noRd
.person_lz <- function(resp, thr_list, theta) {
  if (is.na(theta) || !is.finite(theta)) {
    return(NA_real_)
  }
  ok <- !is.na(resp)
  if (!any(ok)) {
    return(NA_real_)
  }
  resp <- resp[ok]
  thr_list <- thr_list[ok]

  l <- 0
  El <- 0
  Vl <- 0
  for (i in seq_along(thr_list)) {
    p <- .pcm_cat_probs(theta, thr_list[[i]])
    lp <- log(p)
    l <- l + lp[resp[i] + 1L]
    m1 <- sum(p * lp)
    El <- El + m1
    Vl <- Vl + sum(p * lp^2) - m1^2
  }
  if (Vl <= 0) {
    return(NA_real_)
  }
  (l - El) / sqrt(Vl)
}

#' Wilson-Hilferty (ZSTD) transform of an MSQ value (non-inferential)
#' @keywords internal
#' @noRd
.msq_to_zstd <- function(msq, df) {
  q <- pmax(df, 1L)
  (msq^(1 / 3) - 1) * (3 / sqrt(2 / q)) + (sqrt(2 / q) / 3)
}

# ---------------------------------------------------------------------
# Internal: resampling
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.person_fit_resample <- function(
  data_mat,
  thr_list,
  theta,
  statistics,
  obs,
  iterations,
  theta_method = "WLE",
  parallel,
  n_cores,
  msq_tail = "two.sided"
) {
  n <- nrow(data_mat)
  want_msq <- any(c("infit", "outfit") %in% statistics)
  want_lz <- "lz" %in% statistics

  # .cond_moments() depends only on the answered-item set (not the person or
  # the observed score), so compute it once per distinct answered-set and reuse
  # it -- on complete data a single build instead of one per respondent. It is
  # deterministic, so caching does not touch the resampling RNG stream.
  ans_sets <- lapply(seq_len(n), function(p) which(!is.na(data_mat[p, ])))
  ans_key <- vapply(
    ans_sets,
    function(a) paste0(a, collapse = ","),
    character(1L)
  )
  mom_cache <- NULL
  if (want_msq) {
    mom_cache <- list()
    for (k in unique(ans_key)) {
      a <- ans_sets[[match(k, ans_key)]]
      if (length(a) > 0L) mom_cache[[k]] <- .cond_moments(thr_list[a])
    }
  }

  # Each statistic uses its matching null: the conditional MSQ statistics are
  # referenced against patterns sampled conditional on the total score (no
  # person estimate, consistent with the statistic); lz, which is defined at
  # the estimated location, is referenced against patterns simulated there and
  # scored with a location *re-estimated* from each simulated pattern (see the
  # lz block below).
  one_person <- function(p) {
    ans <- ans_sets[[p]]
    res <- c(infit = NA_real_, outfit = NA_real_, lz = NA_real_)
    if (length(ans) == 0L) {
      return(res)
    }
    thr_sub <- thr_list[ans]
    r_obs <- sum(data_mat[p, ans])

    # MSQ: conditional-on-total-score resampling (moments cached per answered-set)
    if (want_msq) {
      mom <- mom_cache[[ans_key[p]]]
      if (r_obs > 0L && r_obs < mom$max_score) {
        sims <- .sim_conditional(thr_sub, r_obs, iterations)
        if (!is.null(sims)) {
          stat <- t(apply(sims, 1L, function(xs) {
            rr <- sum(xs)
            if (rr == 0L || rr == mom$max_score) {
              return(c(NA, NA))
            }
            e <- mom$E[rr + 1L, ]
            v <- mom$V[rr + 1L, ]
            z <- (xs - e) / sqrt(v)
            c(mean(z^2), sum(v * z^2) / sum(v)) # outfit, infit
          }))
          if ("outfit" %in% statistics) {
            res["outfit"] <- .mc_pvalue(stat[, 1L], obs$outfit[p], msq_tail)
          }
          if ("infit" %in% statistics) {
            res["infit"] <- .mc_pvalue(stat[, 2L], obs$infit[p], msq_tail)
          }
        }
      }
    }

    # lz: simulate at the estimated location, then *re-estimate* the location
    # from each simulated pattern and score lz there (lower tail = misfit).
    # Re-estimation is essential: the observed lz carries the ability-estimation
    # effect (its variance is < 1; Snijders 2001), so the null must reproduce it
    # or the test is conservative (Sinharay 2016). For a Rasch/PCM model the WLE
    # / EAP location is a function of the total score, so this re-estimation is a
    # cheap score-based lookup (cached inside .estimate_thetas()), and it is the
    # resampling analogue of the analytic lz* standardization.
    if (want_lz && is.finite(theta[p])) {
      sims <- .sim_theta(thr_sub, theta[p], iterations)
      if (!is.null(sims)) {
        full <- matrix(NA_real_, iterations, ncol(data_mat))
        full[, ans] <- sims
        th_sim <- .estimate_thetas(
          full,
          thr_list,
          method = theta_method,
          theta_range = c(-10, 10)
        )$theta
        lz_sim <- vapply(
          seq_len(iterations),
          function(b) {
            .person_lz(full[b, ], thr_list, th_sim[b])
          },
          numeric(1L)
        )
        res["lz"] <- .mc_pvalue(lz_sim, obs$lz[p], tail = "lower")
      }
    }
    res
  }

  rows <- .maybe_parallel_lapply(seq_len(n), one_person, parallel, n_cores)
  m <- do.call(rbind, rows)
  list(infit = m[, "infit"], outfit = m[, "outfit"], lz = m[, "lz"])
}

#' Empirical p-value for an MSQ-type statistic against its simulated null
#'
#' @param sim Numeric vector of simulated null values.
#' @param obs Observed value.
#' @param tail `"two.sided"` (deviation in either direction), `"upper"`
#'   (underfit only -- the validity-relevant direction), or `"lower"` (e.g.
#'   the lz statistic, where small values indicate misfit).
#' @keywords internal
#' @noRd
.mc_pvalue <- function(sim, obs, tail = "two.sided") {
  sim <- sim[!is.na(sim)]
  B <- length(sim)
  if (B == 0L || is.na(obs)) {
    return(NA_real_)
  }
  # (1 + count) / (B + 1) Monte-Carlo p-value, consistent with
  # .bootstrap_pvalues() in utils-multiplicity.R: never exactly 0, and the
  # denominator counts only the valid (non-NA) simulations.
  if (tail == "upper") {
    return((1 + sum(sim >= obs)) / (B + 1))
  }
  if (tail == "lower") {
    return((1 + sum(sim <= obs)) / (B + 1))
  }
  p_upper <- (1 + sum(sim >= obs)) / (B + 1)
  p_lower <- (1 + sum(sim <= obs)) / (B + 1)
  min(1, 2 * min(p_upper, p_lower))
}

#' Simulate response patterns at a fixed theta (per-person MC scheme)
#' @keywords internal
#' @noRd
.sim_theta <- function(thr_sub, theta, B) {
  if (is.na(theta) || !is.finite(theta)) {
    return(NULL)
  }
  k <- length(thr_sub)
  out <- matrix(0L, B, k)
  for (i in seq_len(k)) {
    p <- .pcm_cat_probs(theta, thr_sub[[i]])
    out[, i] <- sample.int(length(p), B, replace = TRUE, prob = p) - 1L
  }
  out
}

#' Simulate response patterns conditional on a fixed total score
#'
#' Sequential exact sampling via suffix elementary symmetric functions:
#' draws each item in turn from its conditional distribution given the
#' remaining target score and the remaining items.
#' @keywords internal
#' @noRd
.sim_conditional <- function(thr_sub, total, B) {
  k <- length(thr_sub)
  steps <- vapply(thr_sub, length, integer(1))
  if (total <= 0L || total >= sum(steps)) {
    return(NULL)
  }

  # per-item category weights w_{i,c} = exp(-eta_{i,c}), c = 0..m_i
  wlist <- lapply(thr_sub, function(thr) exp(-c(0, cumsum(thr))))
  # suffix ESF: gamma of items i..k, for i = 1..k+1 (k+1 = empty -> gamma_0 = 1)
  suffix_g <- vector("list", k + 1L)
  suffix_g[[k + 1L]] <- 1 # empty product
  for (i in k:1) {
    suffix_g[[i]] <- .esf_convolve(suffix_g[[i + 1L]], wlist[[i]])
  }

  out <- matrix(0L, B, k)
  for (b in seq_len(B)) {
    s <- total
    for (i in seq_len(k)) {
      w_i <- wlist[[i]]
      g_rest <- suffix_g[[i + 1L]] # gamma of items i+1..k
      cats <- 0:steps[i]
      # P(x_i = c) propto w_{i,c} * gamma_{s-c}(rest)
      probs <- vapply(
        cats,
        function(c) {
          idx <- s - c
          if (idx < 0L || idx >= length(g_rest)) {
            return(0)
          }
          w_i[c + 1L] * g_rest[idx + 1L]
        },
        numeric(1)
      )
      if (sum(probs) <= 0) {
        out[b, i] <- 0L
        next
      }
      xc <- sample.int(length(probs), 1L, prob = probs) - 1L
      out[b, i] <- xc
      s <- s - xc
    }
  }
  out
}

#' Convolve an ESF vector with one item's category weights
#'
#' Given gamma (ESF of a set of items, indexed by score 0..S) and an
#' item's weights w (score 0..m), returns the ESF of the combined set.
#' @keywords internal
#' @noRd
.esf_convolve <- function(gamma, w) {
  S <- length(gamma) - 1L
  m <- length(w) - 1L
  out <- numeric(S + m + 1L)
  for (s in 0:S) {
    gs <- gamma[s + 1L]
    if (gs == 0) {
      next
    }
    for (c in 0:m) {
      out[s + c + 1L] <- out[s + c + 1L] + gs * w[c + 1L]
    }
  }
  out
}

#' Optionally parallel lapply via mirai
#' @keywords internal
#' @noRd
.maybe_parallel_lapply <- function(x, fun, parallel, n_cores) {
  if (parallel && requireNamespace("mirai", quietly = TRUE)) {
    workers <- if (is.null(n_cores)) min(2L, length(x)) else n_cores
    mirai::daemons(workers)
    on.exit(mirai::daemons(0L), add = TRUE)
    res <- mirai::mirai_map(x, fun)[]
    return(res)
  }
  lapply(x, fun)
}

# ---------------------------------------------------------------------
# Internal: output formatting
# ---------------------------------------------------------------------

#' @keywords internal
#' @noRd
.person_fit_output <- function(
  out,
  theta,
  extreme,
  statistics,
  output,
  flag_alpha,
  iterations,
  flag = "both",
  n_total = nrow(out),
  has_na = FALSE
) {
  # Assessed = respondents that received a fit statistic (extreme scorers get
  # NA). n_total is the raw supplied row count (before all-NA rows were
  # dropped), so the caption anchors to the full sample. Each parenthetical
  # qualifier appears only when it applies (missingness / extreme scorers).
  n_assessed <- sum(!extreme)
  pf_quals <- c(
    if (has_na) "incomplete responses retained",
    if (any(extreme)) "extreme scores excluded"
  )
  n_clause <- paste0(
    " ",
    .n_caption(n_assessed, n_total, pf_quals, noun = "respondents assessed"),
    "."
  )
  if (output == "dataframe") {
    # Unrounded: rounding is a display concern (the kable path rounds a
    # display copy below).
    return(out)
  }

  if (output %in% c("ggplot", "list")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop(
        "Package 'ggplot2' is required for output = \"", output, "\". ",
        "Install it with: install.packages(\"ggplot2\")",
        call. = FALSE
      )
    }
    # One plot per requested statistic, returned as a named list.
    stat_map <- list(
      infit = list(
        col = "infit_msq",
        pcol = "p_infit",
        label = "Infit MSQ",
        ref = 1,
        split = TRUE
      ),
      outfit = list(
        col = "outfit_msq",
        pcol = "p_outfit",
        label = "Outfit MSQ",
        ref = 1,
        split = TRUE
      ),
      lz = list(col = "lz", pcol = "p_lz", label = "lz", ref = 0, split = FALSE)
    )
    plots <- lapply(statistics, function(s) {
      m <- stat_map[[s]]
      .person_fit_plot(
        out,
        theta,
        m$col,
        m$pcol,
        m$label,
        m$ref,
        m$split,
        flag_alpha,
        iterations,
        flag,
        n_total = n_total,
        qualifiers = pf_quals
      )
    })
    names(plots) <- statistics
    if (output == "list") {
      # One computation, both views: the per-person results and the maps.
      return(list(fit = out, plots = plots))
    }
    return(plots)
  }

  display <- out
  num <- vapply(display, is.numeric, logical(1))
  display[num] <- lapply(display[num], function(x) round(x, 3))
  caption <- paste0(
    "Person fit (conditional MSQ + lz). ",
    if (iterations > 0L) {
      msq_desc <- if (flag == "underfit") {
        "one-sided (underfit) MSQ p-values"
      } else {
        "two-sided MSQ p-values"
      }
      paste0(
        msq_desc,
        " from ",
        iterations,
        " Monte-Carlo replications (MSQ conditional on the total score, ",
        "lz simulated at the estimated location); ",
        "flagged at alpha = ",
        flag_alpha,
        "."
      )
    } else {
      "No resampling (iterations = 0)."
    },
    n_clause
  )
  knitr::kable(display, format = "pipe", caption = caption, row.names = FALSE)
}

#' Build one person-fit map (statistic vs person location)
#'
#' @param out The assembled person-fit data.frame.
#' @param theta Person locations (x-axis).
#' @param ycol,pcol Column names for the statistic and its marginal p-value.
#' @param ylab Axis label for the statistic.
#' @param ref Horizontal reference line (expected value), or `NA` for none.
#' @param split Logical. If `TRUE` (the MSQ statistics), flagged respondents
#'   are split into underfit (above `ref`) and overfit (below `ref`) for both
#'   the point colour and the caption. `FALSE` for one-sided statistics (lz).
#' @param flag_alpha Flagging threshold on the marginal p-value.
#' @param iterations Number of resampling iterations (0 = no p-values).
#' @param flag `"both"` (two-sided MSQ test) or `"underfit"` (one-sided);
#'   controls the caption wording for the MSQ maps.
#' @return A `ggplot`. Caption reports the assessed sample size and the
#'   proportion flagged by this statistic.
#' @keywords internal
#' @noRd
.person_fit_plot <- function(
  out,
  theta,
  ycol,
  pcol,
  ylab,
  ref,
  split,
  flag_alpha,
  iterations,
  flag = "both",
  n_total = NULL,
  qualifiers = character()
) {
  has_p <- iterations > 0L && pcol %in% names(out)
  df <- data.frame(
    theta = theta,
    y = out[[ycol]],
    flagged = if (has_p) out[[pcol]] < flag_alpha else FALSE
  )
  df <- df[is.finite(df$theta) & is.finite(df$y), , drop = FALSE]
  n_assessed <- nrow(df)
  n_lead <- .n_caption(
    n_assessed,
    if (is.null(n_total)) n_assessed else n_total,
    qualifiers,
    noun = "respondents assessed"
  )

  # Classify each respondent for colouring.
  df$status <- "Not flagged"
  if (split && is.finite(ref)) {
    df$status[df$flagged & df$y > ref] <- "Underfit"
    df$status[df$flagged & df$y < ref] <- "Overfit"
  } else {
    df$status[df$flagged] <- "Flagged"
  }
  pal <- c(
    "Not flagged" = "grey55",
    "Underfit" = "#e41a1c",
    "Overfit" = "#377eb8",
    "Flagged" = "#e41a1c"
  )
  present <- intersect(names(pal), unique(df$status))
  df$status <- factor(df$status, levels = present)

  if (has_p) {
    n_flag <- sum(df$flagged, na.rm = TRUE)
    pct <- if (n_assessed > 0L) {
      round(100 * n_flag / n_assessed, 1)
    } else {
      NA_real_
    }
    if (split && is.finite(ref)) {
      n_under <- sum(df$flagged & df$y > ref, na.rm = TRUE) # MSQ > 1: noisy
      n_over <- sum(df$flagged & df$y < ref, na.rm = TRUE) # MSQ < 1: muted
      caption <- if (flag == "underfit") {
        paste0(
          n_lead,
          "; ",
          n_flag,
          " (",
          pct,
          "%) flagged as underfit (MSQ > ",
          ref,
          "; one-sided p < ",
          flag_alpha,
          ")."
        )
      } else {
        paste0(
          n_lead,
          "; ",
          n_flag,
          " (",
          pct,
          "%) flagged: ",
          n_under,
          " underfit (MSQ > ",
          ref,
          "), ",
          n_over,
          " overfit (MSQ < ",
          ref,
          ") (p < ",
          flag_alpha,
          ")."
        )
      }
    } else {
      caption <- paste0(
        n_lead,
        "; ",
        n_flag,
        " (",
        pct,
        "%) flagged as misfitting (p < ",
        flag_alpha,
        ")."
      )
    }
  } else {
    caption <- paste0(
      n_lead,
      " (no resampling; iterations = 0)."
    )
  }

  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$theta, y = .data$y))
  if (is.finite(ref)) {
    p <- p +
      ggplot2::geom_hline(
        yintercept = ref,
        linetype = "dashed",
        colour = "grey60"
      )
  }
  if (has_p && length(present) > 1L) {
    p <- p +
      ggplot2::geom_point(ggplot2::aes(colour = .data$status), alpha = 0.7) +
      ggplot2::scale_colour_manual(
        values = pal,
        name = "Person fit",
        breaks = present
      )
  } else {
    p <- p + ggplot2::geom_point(colour = "grey45", alpha = 0.7)
  }
  p +
    ggplot2::labs(
      x = "Person location (logits)",
      y = ylab,
      caption = er2_caption(caption)
    ) +
    ggplot2::theme_bw() +
    er2_axis_margins() +
    er2_plot_caption()
}
