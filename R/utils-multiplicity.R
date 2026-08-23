# Internal helpers for bootstrap p-values and multiple-comparison correction.
#
# Shared by the simulation-cutoff functions that compare an observed per-item
# (or per-pair) statistic against a parametric-bootstrap null distribution.
# The per-comparison statistic is the residual studentised by the *bootstrap*
# mean and SD, putting all comparisons on equal footing for the family-wise
# maximum (Ferreira, 2024) without invoking the Wilson-Hilferty / ZSTD
# transform criticised by Müller (2020).

#' Bootstrap p-values with optional multiplicity correction
#'
#' @param observed Named numeric vector of observed statistics (one per item
#'   or pair).
#' @param sim Numeric matrix of simulated null values, `iterations` rows by
#'   `length(observed)` columns, with column names matching `names(observed)`.
#' @param correction One of `"fwer"` (Westfall-Young studentised-max
#'   step-down), `"fdr_bh"`, `"fdr_by"`, or `"none"`.
#' @param tail `"two.sided"` (default; deviations in either direction are
#'   extreme, e.g. infit over/underfit) or `"upper"` (only large positive
#'   values are extreme, e.g. excess Q3 local dependence).
#' @return A data.frame with columns `name`, `p` (marginal Monte-Carlo
#'   p-value), and `padj` (corrected p-value). Comparisons whose simulated
#'   values have no variance return `NA`.
#' @keywords internal
#' @noRd
.bootstrap_pvalues <- function(
  observed,
  sim,
  correction = c("fwer", "fdr_bh", "fdr_by", "none"),
  tail = c("two.sided", "upper")
) {
  correction <- match.arg(correction)
  tail <- match.arg(tail)
  nm <- names(observed)
  sim <- sim[, nm, drop = FALSE]
  B <- nrow(sim)

  m <- colMeans(sim, na.rm = TRUE)
  s <- apply(sim, 2L, stats::sd, na.rm = TRUE)
  s[!is.finite(s) | s < 1e-8] <- NA_real_ # guard zero-variance comparisons

  # Studentise observed and simulated by the bootstrap mean / SD; the test
  # statistic is |t| (two-sided) or t (upper tail).
  t_obs <- (observed - m) / s
  t_sim <- sweep(sweep(sim, 2L, m, `-`), 2L, s, `/`)
  if (tail == "two.sided") {
    stat_obs <- abs(t_obs)
    stat_sim <- abs(t_sim)
  } else {
    stat_obs <- t_obs
    stat_sim <- t_sim
  }

  # Marginal Monte-Carlo p-value, (1 + count) / (B + 1).
  p_marg <- vapply(
    seq_along(nm),
    function(i) {
      if (is.na(stat_obs[i])) {
        return(NA_real_)
      }
      (1 + sum(stat_sim[, i] >= stat_obs[i], na.rm = TRUE)) / (B + 1)
    },
    numeric(1)
  )

  padj <- switch(
    correction,
    none = p_marg,
    fdr_bh = stats::p.adjust(p_marg, method = "BH"),
    fdr_by = stats::p.adjust(p_marg, method = "BY"),
    fwer = .wy_stepdown(stat_obs, stat_sim, B)
  )

  data.frame(
    name = nm,
    p = p_marg,
    padj = padj,
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Format a p-value for prose captions
#'
#' Plain-number formatting: exact to three decimals down to 0.001, then the
#' conventional "p < 0.001". Avoids the scientific notation ("<1e-04") that
#' `format.pval()` emits, which most readers are not used to.
#'
#' @param p Numeric p-value (length 1).
#' @return A string starting with "p", e.g. "p = 0.043" or "p < 0.001".
#' @keywords internal
#' @noRd
.format_p <- function(p) {
  if (is.na(p)) {
    return("p = NA")
  }
  if (p < 0.001) "p < 0.001" else paste0("p = ", sprintf("%.3f", p))
}

#' Human-readable label for a multiplicity-correction method
#'
#' Shared by the kable captions of the p-value-reporting functions.
#'
#' @param correction One of `"fwer"`, `"fdr_bh"`, `"fdr_by"`, `"none"`.
#' @return A character label.
#' @keywords internal
#' @noRd
.correction_label <- function(correction) {
  switch(
    correction,
    fwer = "Westfall-Young step-down (FWER)",
    fdr_bh = "Benjamini-Hochberg (FDR)",
    fdr_by = "Benjamini-Yekutieli (FDR)",
    none = "uncorrected"
  )
}

# ---------------------------------------------------------------------------
# The interval width as an error-rate parameter
#
# Flagging every item whose statistic falls outside a width-w interval tests
# k hypotheses at once, so the family-wise error rate is 1 - w^k (Sidak). The
# helpers below turn that relation around, and translate a width into the
# number of iterations it needs before the interval means what it says.
# Johansson (2026) measures all three.
# ---------------------------------------------------------------------------

#' Family-wise error rate implied by an interval width
#'
#' @param w Interval width in (0, 1).
#' @param k Number of comparisons.
#' @return Numeric family-wise error rate.
#' @keywords internal
#' @noRd
.sidak_fwe <- function(w, k) 1 - w^k

#' Interval width implying a target family-wise error rate
#'
#' @param k Number of comparisons.
#' @param alpha Target family-wise error rate (default `0.05`).
#' @return Numeric width.
#' @keywords internal
#' @noRd
.sidak_width <- function(k, alpha = 0.05) (1 - alpha)^(1 / k)

#' Iterations an interval of a given width needs to converge
#'
#' Ten simulated values per tail, `B = 20 / (1 - w)`. Verified for
#' `w <= .999` (Johansson, 2026).
#'
#' @param w Interval width in (0, 1).
#' @return Numeric iteration count.
#' @keywords internal
#' @noRd
.width_iterations <- function(w) 20 / (1 - w)

#' Capitalise the first letter of a string
#'
#' @param x Character string.
#' @return `x` with its first character upper-cased.
#' @keywords internal
#' @noRd
.capitalise <- function(x) {
  if (!nzchar(x)) {
    return(x)
  }
  paste0(toupper(substr(x, 1L, 1L)), substr(x, 2L, nchar(x)))
}

#' Format a width for display, dropping the leading zero
#'
#' @param w Numeric width.
#' @return A string such as `".95"` or `".99432"`.
#' @keywords internal
#' @noRd
.fmt_width <- function(w) {
  s <- sub("0+$", "", sprintf("%.5f", w))
  sub("^0", "", sub("[.]$", "", s))
}

#' Clause describing the error rate of interval-based flagging
#'
#' Used in the kable caption whenever items are flagged against the interval
#' rather than against a corrected p-value. Returns `NULL` when the width is
#' unknown, which happens when only the bare `$item_cutoffs` data.frame was
#' passed and the metadata went with it.
#'
#' @param width Interval width, or `NULL`.
#' @param k Number of comparisons.
#' @param unit Plural noun for what is being compared. `"items"` for item fit,
#'   `"item pairs"` for local dependence.
#' @return A single string, or `NULL`.
#' @keywords internal
#' @noRd
.band_error_clause <- function(width, k, unit = "items") {
  if (is.null(width) || !is.finite(width) || width <= 0 || width >= 1) {
    return(NULL)
  }
  paste0(
    " Flagging against the interval tests all ",
    k,
    " ",
    unit,
    " at once, which implies a family-wise error rate of about ",
    sprintf("%.0f", 100 * .sidak_fwe(width, k)),
    "% (Johansson, 2026)."
  )
}

#' One-time console notice that flagging is interval-based
#'
#' Emitted once per session when items are flagged against the interval. The
#' interval is a description of where a fitting item's statistic is expected
#' to fall, and using it as a decision rule sets an error rate implicitly and
#' usually far above .05. A `message()` rather than a `warning()`, because the
#' call is legitimate and a hard warning on a supported path trains users to
#' ignore warnings.
#'
#' @param width Interval width, or `NULL` when unknown.
#' @param k Number of comparisons.
#' @param unit Plural noun for what is being compared, e.g. `"items"` or
#'   `"item pairs"`.
#' @param fn Name of the cutoff function to point the reader at.
#' @param id Frequency id, so that each family of functions gets to say this
#'   once per session rather than the first one silencing the rest.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.notify_band_flagging <- function(
  width,
  k,
  unit = "items",
  fn = "RMitemInfitCutoff()",
  id = "easyRasch2_band_flagging"
) {
  head <- paste0(
    .capitalise(unit),
    " are flagged against the interval, not against a corrected p-value."
  )
  body <- if (is.null(width) || !is.finite(width) || width <= 0 || width >= 1) {
    c(i = paste(
      "The interval tests all", k, unit, "at once, so its width sets a",
      "family-wise error rate of 1 - width^k."
    ))
  } else {
    need_w <- .sidak_width(k)
    c(
      i = sprintf(
        "The %s interval over %d %s implies a family-wise error rate of about %.0f%%.",
        .fmt_width(width),
        k,
        unit,
        100 * .sidak_fwe(width, k)
      ),
      i = sprintf(
        "A rate of .05 would need width %s, which needs roughly %s iterations to converge.",
        .fmt_width(need_w),
        format(
          signif(.width_iterations(need_w), 2),
          big.mark = " ",
          scientific = FALSE,
          trim = TRUE
        )
      )
    )
  }
  rlang::inform(
    c(
      head,
      body,
      i = paste0(
        "Pass the full ",
        fn,
        " object and leave `p_value = NULL` to flag on the Westfall-Young ",
        "corrected p-value, which targets .05 directly at 400 iterations."
      ),
      i = "See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2."
    ),
    .frequency = "once",
    .frequency_id = id
  )
  invisible(NULL)
}

#' One-time console notice that the bootstrap ran too few iterations
#'
#' Below 400 iterations the Westfall-Young rule is mildly liberal under the
#' null (Johansson, 2026, measured 6.8% at 100 iterations against a nominal
#' 5%, and 4.3% at 400). Between 400 and 1000 the error rate is trustworthy
#' and only reproducibility improves further, which the table caption covers
#' instead.
#'
#' @param n_iter Number of completed iterations.
#' @param requested Number of iterations asked for, or `NULL`. When the two
#'   differ the notice says so, since someone who asked for 400 and landed
#'   below it needs to know that iterations were discarded rather than that
#'   they chose too few.
#' @param fn Name of the cutoff function whose `iterations` should be raised.
#' @param id Frequency id, one per family of functions.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.notify_low_iterations <- function(
  n_iter,
  requested = NULL,
  fn = "RMitemInfitCutoff()",
  id = "easyRasch2_low_iterations"
) {
  lost <- !is.null(requested) && is.finite(requested) && n_iter < requested
  rlang::inform(
    c(
      sprintf(
        "Bootstrap p-values are based on %d iterations, below the calibrated floor of 400.",
        n_iter
      ),
      i = paste(
        "Below 400 the Westfall-Young correction is mildly liberal under the",
        "null, so the family-wise error rate is above the nominal level."
      ),
      if (lost) {
        c(i = sprintf(
          paste(
            "%d of the %d simulated datasets could not be refitted, usually",
            "because an item ended up with an unused response category or",
            "almost no variation."
          ),
          requested - n_iter,
          requested
        ))
      },
      i = "See Johansson (2026), doi:10.31234/osf.io/7pqz4_v2.",
      i = paste0("Raise `iterations` in ", fn, ".")
    ),
    .frequency = "once",
    .frequency_id = id
  )
  invisible(NULL)
}

#' Caption sentence about the iteration count
#'
#' The second tier of the two-tier treatment. Below 400 the correction itself
#' is mildly liberal, which the console also reports through
#' `.notify_low_iterations()`. Between 400 and 1000 the error rate is
#' trustworthy and only reproducibility keeps improving, which is a caption
#' matter rather than something to interrupt over. At 1000 and above the
#' caption says nothing.
#'
#' The measured seed-disagreement percentages in `RMitemInfit()`'s caption are
#' for items and are not repeated here, since the equivalent quantity for item
#' pairs has not been measured.
#'
#' @param n_iter Completed iterations, or `NULL`.
#' @return A single string, empty when there is nothing to say.
#' @keywords internal
#' @noRd
.iteration_note <- function(n_iter) {
  if (is.null(n_iter) || !is.finite(n_iter)) {
    return("")
  }
  if (n_iter < 400L) {
    return(paste0(
      " This is below the calibrated floor of 400, where the correction is",
      " mildly liberal and the family-wise error rate sits above the nominal",
      " level (Johansson, 2026)."
    ))
  }
  if (n_iter < 1000L) {
    return(paste0(
      " Error rates are calibrated at this many iterations, but decisions are",
      " still somewhat seed-dependent, so use 1000 to 2000 for a final",
      " analysis (Johansson, 2026)."
    ))
  }
  ""
}

#' Iterations a false discovery rate procedure needs to reject anything
#'
#' A Monte Carlo p-value cannot fall below `1 / (B + 1)`. When `s` of the `m`
#' p-values sit at that floor, the smallest attainable adjusted value is
#' `m / (s * (B + 1))` for Benjamini-Hochberg and `c(m) * m / (s * (B + 1))`
#' for Benjamini-Yekutieli, with `c(m) = sum(1 / seq_len(m))`. Both are exact
#' against [stats::p.adjust()]. The requirement is therefore
#' `B >= m / (s * alpha) - 1`, which relaxes as more comparisons reach the
#' floor.
#'
#' `s = 1` is the binding case, since a single genuinely dependent comparison
#' is the situation an analyst cannot rule out in advance, and it is what this
#' function returns. Note that the familiar `B >= m / alpha - 1` is this case
#' and not a general bound. It is also one iteration short, because flagging
#' tests `padj < alpha` strictly and `B = m / alpha - 1` lands the adjusted
#' value exactly on `alpha`. The count is therefore found by search rather
#' than from the closed form, which also keeps it exact when `m / alpha` is
#' not representable in binary.
#'
#' The Westfall-Young alternative needs `1 / alpha - 1` iterations whatever
#' `m` is, which is why it is the default correction.
#'
#' @param m Number of comparisons.
#' @param correction One of `"fdr_bh"`, `"fdr_by"`, `"fwer"` or `"none"`.
#' @param alpha Significance level.
#' @return Integer iterations needed for a single standout comparison to be
#'   rejectable, or `NULL` when the correction has no such floor beyond the
#'   trivial one.
#' @keywords internal
#' @noRd
.fdr_min_iterations <- function(m, correction, alpha = 0.05) {
  if (!correction %in% c("fdr_bh", "fdr_by")) {
    return(NULL)
  }
  if (!is.finite(m) || m < 1 || !is.finite(alpha) || alpha <= 0) {
    return(NULL)
  }
  penalty <- if (identical(correction, "fdr_by")) sum(1 / seq_len(m)) else 1
  # Mirrors how p.adjust() forms the value, `q * n/i * p` with `p = 1/(b+1)`,
  # so that the two agree in the last bit when `m / alpha` is a whole number.
  reachable <- function(b) (penalty * m) * (1 / (b + 1)) < alpha
  b <- max(1, floor(penalty * m / alpha) - 1)
  while (!reachable(b)) b <- b + 1
  while (b > 1 && reachable(b - 1)) b <- b - 1
  as.integer(b)
}

#' Warn when the bootstrap cannot reach the false discovery rate threshold
#'
#' Fires only for `"fdr_bh"` and `"fdr_by"`, which the user has to select
#' deliberately, so it never sounds on a default call. The wording is careful
#' not to claim that nothing can be flagged, since several comparisons at the
#' Monte Carlo floor lower the requirement.
#'
#' @param n_iter Completed iterations.
#' @param m Number of comparisons.
#' @param correction The correction in force.
#' @param alpha Significance level.
#' @param unit Plural noun for what is being compared.
#' @param fn Name of the cutoff function whose `iterations` should be raised.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.warn_fdr_floor <- function(
  n_iter,
  m,
  correction,
  alpha = 0.05,
  unit = "comparisons",
  fn = "RMitemInfitCutoff()"
) {
  need <- .fdr_min_iterations(m, correction, alpha)
  if (is.null(need) || is.null(n_iter) || !is.finite(n_iter) || n_iter >= need) {
    return(invisible(NULL))
  }
  label <- if (identical(correction, "fdr_by")) {
    "Benjamini-Yekutieli"
  } else {
    "Benjamini-Hochberg"
  }
  warning(
    "The ",
    label,
    " threshold is out of reach of ",
    n_iter,
    " iterations. A bootstrap p-value cannot be smaller than 1/(",
    n_iter,
    "+1) = ",
    signif(1 / (n_iter + 1), 3),
    ", and over ",
    m,
    " ",
    unit,
    " a single standout ",
    sub("s$", "", unit),
    " needs a raw p-value below ",
    signif(alpha / m, 3),
    " to clear the threshold, so it cannot be flagged below ",
    need,
    " iterations. Several extreme ",
    unit,
    " lower that requirement, since the smallest attainable adjusted p-value ",
    "is m/(s(B+1)) for s of them at the floor. Raise `iterations` in ",
    fn,
    " to at least ",
    need,
    ", or use the default correction = \"fwer\", which needs ",
    ceiling(1 / alpha) - 1,
    " iterations whatever the number of ",
    unit,
    ".",
    call. = FALSE
  )
  invisible(NULL)
}

#' Westfall-Young studentised-max step-down adjusted p-values
#'
#' Rejects the most extreme comparison against the maximum over all
#' comparisons, then recomputes the maximum over the remaining comparisons,
#' enforcing monotonicity. Controls the family-wise error rate while using the
#' bootstrap dependence among the statistics (Westfall & Young, 1993;
#' Ferreira, 2024). Validity of the step-down rests on subset pivotality.
#'
#' @param stat_obs Numeric vector of observed test statistics (already
#'   |studentised| or signed-studentised depending on the tail).
#' @param stat_sim Numeric matrix (`iterations` x comparisons) of the
#'   simulated test statistics.
#' @param B Number of iterations.
#' @return Numeric vector of adjusted p-values (NA for invalid comparisons).
#' @keywords internal
#' @noRd
.wy_stepdown <- function(stat_obs, stat_sim, B) {
  k <- length(stat_obs)
  padj <- rep(NA_real_, k)

  valid <- which(is.finite(stat_obs))
  if (length(valid) == 0L) {
    return(padj)
  }

  # Process valid comparisons in order of decreasing observed evidence.
  ord <- valid[order(stat_obs[valid], decreasing = TRUE)]
  prev <- 0
  for (j in seq_along(ord)) {
    remaining <- ord[j:length(ord)]
    max_star <- apply(
      stat_sim[, remaining, drop = FALSE],
      1L,
      max,
      na.rm = TRUE
    )
    p_j <- (1 + sum(max_star >= stat_obs[ord[j]], na.rm = TRUE)) / (B + 1)
    prev <- max(prev, p_j) # enforce monotonicity (step-down)
    padj[ord[j]] <- min(prev, 1)
  }
  padj
}
