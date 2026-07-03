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
