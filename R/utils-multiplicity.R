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
#' @param k Number of items.
#' @return A single string, or `NULL`.
#' @keywords internal
#' @noRd
.band_error_clause <- function(width, k) {
  if (is.null(width) || !is.finite(width) || width <= 0 || width >= 1) {
    return(NULL)
  }
  paste0(
    " Flagging against the interval tests all ",
    k,
    " items at once, which implies a family-wise error rate of about ",
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
#' @param k Number of items.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.notify_band_flagging <- function(width, k) {
  head <- "Items are flagged against the interval, not against a corrected p-value."
  body <- if (is.null(width) || !is.finite(width) || width <= 0 || width >= 1) {
    c(i = paste(
      "The interval tests all", k, "items at once, so its width sets a",
      "family-wise error rate of 1 - width^k."
    ))
  } else {
    need_w <- .sidak_width(k)
    c(
      i = sprintf(
        "The %s interval over %d items implies a family-wise error rate of about %.0f%%.",
        .fmt_width(width),
        k,
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
      i = paste(
        "Pass the full RMitemInfitCutoff() object and leave `p_value = NULL`",
        "to flag on the Westfall-Young corrected p-value, which targets .05",
        "directly at 400 iterations."
      ),
      i = "See Johansson (2026), doi:10.31234/osf.io/7pqz4_v1."
    ),
    .frequency = "once",
    .frequency_id = "easyRasch2_band_flagging"
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
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.notify_low_iterations <- function(n_iter, requested = NULL) {
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
      i = "See Johansson (2026), doi:10.31234/osf.io/7pqz4_v1.",
      i = "Raise `iterations` in RMitemInfitCutoff()."
    ),
    .frequency = "once",
    .frequency_id = "easyRasch2_low_iterations"
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
