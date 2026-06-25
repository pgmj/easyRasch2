#' Raw-Score to Logit Score Transformation Table
#'
#' For a given set of items, returns the score-to-theta lookup that maps each
#' possible raw sum score to a person-location estimate (in logits) and its
#' standard error. Useful when reporting a scale's measurement properties or
#' converting raw totals to interval-scaled scores for downstream analysis.
#'
#' @param data A data.frame or matrix of item responses. Items must be scored
#'   starting at 0 (non-negative integers). Missing values (`NA`) are allowed;
#'   the underlying model fit handles them.
#' @param method Character string. Either `"WLE"` (default) for Warm's Weighted
#'   Likelihood Estimator computed from a CML-fitted Rasch / Partial Credit
#'   Model via `psychotools`, or `"EAP"` for Expected A Posteriori sum-score
#'   estimates from an MML-fitted model via `mirt`.
#' @param output Character string controlling the return value: `"kable"`
#'   (default) for a formatted `knitr::kable()` table, `"dataframe"` for the
#'   underlying data.frame, or `"ggplot"` for a `ggplot2` figure showing each
#'   raw score's logit estimate with `ci_multiplier`-scaled error bars.
#' @param ci_multiplier Numeric. Multiplier applied to the standard error to
#'   draw error bars on the figure. Default `1.96` (\eqn{\approx}95% CI under a Gaussian
#'   approximation). Ignored when `output != "ggplot"`.
#' @param point_size Numeric. Point size for the figure. Default `3`.
#' @param error_width Numeric. Cap width for error bars on the figure.
#'   Default `0.5`.
#' @param theta_range Numeric length 2. Theta search range used for boundary
#'   raw scores under WLE estimation. Default `c(-10, 10)`. Ignored when
#'   `method = "EAP"`.
#'
#' @return
#' * If `output = "kable"`: a `knitr_kable` object with columns "Ordinal sum
#'   score", "Logit score", and "Logit std.error", and a caption noting the
#'   estimation method.
#' * If `output = "dataframe"`: a data.frame with columns `raw_score`,
#'   `logit_score`, and `logit_se` (one row per possible raw sum score from 0
#'   to the theoretical maximum).
#' * If `output = "ggplot"`: a `ggplot` object — points at each
#'   (`logit_score`, `raw_score`) with horizontal error bars at
#'   ± `ci_multiplier × logit_se`.
#'
#' @details
#' The function automatically detects whether the data is dichotomous
#' (max score 1) or polytomous (max score > 1) and selects the appropriate
#' Rasch / Partial Credit model.
#'
#' **`method = "WLE"`** fits the model by CML with `psychotools::pcmodel()`,
#' centres the item thresholds to grand-mean-zero, and solves Warm's
#' weighted-likelihood equation for each raw score with the same engine used
#' by [RMpersonParameters()]; the two functions therefore report identical
#' locations and standard errors. Warm's bias correction yields finite
#' locations even at the minimum and maximum scores (only a root outside
#' `theta_range` is clamped to the boundary with `NA` SE). The reported
#' `logit_se` is the information-based standard error `1 / sqrt(I(theta))`,
#' matching `catR`, `TAM` and most Rasch software.
#'
#' **`method = "EAP"`** fits the model with `mirt::mirt(..., itemtype =
#' "Rasch")` (MML) and obtains sum-score-based EAP estimates and posterior
#' SDs via `mirt::fscores(method = "EAPsum", full.scores = FALSE,
#' full.scores.SE = TRUE)`. EAP estimates are finite at all score boundaries
#' (the prior shrinks them inward), but they depend on the assumed normal
#' prior on theta. Item parameters from MML differ slightly from the CML
#' values used by the WLE path; for well-behaved data the difference is small.
#'
#' @references
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item
#' response theory. *Psychometrika, 54*(3), 427-450.
#' \doi{10.1007/BF02294627}
#'
#' Bock, R. D., & Mislevy, R. J. (1982). Adaptive EAP estimation of ability
#' in a microcomputer environment. *Applied Psychological Measurement, 6*(4),
#' 431-444. \doi{10.1177/014662168200600405}
#'
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(42)
#' sim_data <- as.data.frame(
#'   matrix(sample(0:3, 200 * 6, replace = TRUE), nrow = 200, ncol = 6)
#' )
#' colnames(sim_data) <- paste0("Item", 1:6)
#'
#' # Default kable output, WLE
#' RMscoreSE(sim_data)
#'
#' # Underlying data.frame
#' RMscoreSE(sim_data, output = "dataframe")
#'
#' # ggplot figure
#' RMscoreSE(sim_data, output = "ggplot")
#'
#' # EAP via mirt
#' RMscoreSE(sim_data, method = "EAP")
#' }
RMscoreSE <- function(data,
                      method        = "WLE",
                      output        = "kable",
                      ci_multiplier = 1.96,
                      point_size    = 3,
                      error_width   = 0.5,
                      theta_range   = c(-10, 10)) {

  method <- match.arg(method, c("WLE", "EAP"))
  output <- match.arg(output, c("kable", "dataframe", "ggplot"))

  validate_response_data(data)

  if (!is.numeric(theta_range) || length(theta_range) != 2L ||
      theta_range[1L] >= theta_range[2L]) {
    stop("`theta_range` must be a numeric vector of length 2 with ",
         "theta_range[1] < theta_range[2].", call. = FALSE)
  }

  # --- Compute the score-to-theta table ---------------------------------------
  score_table <- if (method == "WLE") {
    .scoreSE_wle(data, theta_range)
  } else {
    .scoreSE_eap(data)
  }

  # --- Return: dataframe ------------------------------------------------------
  if (output == "dataframe") {
    return(score_table)
  }

  # --- Return: kable ----------------------------------------------------------
  if (output == "kable") {
    display <- data.frame(
      raw_score    = score_table$raw_score,
      logit_score  = round(score_table$logit_score, 3),
      logit_se     = round(score_table$logit_se,    3),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
    caption <- paste0(
      if (method == "WLE") {
        "Person locations via Warm's WLE (CML item parameters)."
      } else {
        "Person locations via EAPsum (MML item parameters from mirt; depends on a normal theta prior)."
      }
    )
    return(knitr::kable(
      display,
      format    = "pipe",
      col.names = c("Ordinal sum score", "Logit score", "Logit std.error"),
      caption   = caption
    ))
  }

  # --- Return: ggplot ---------------------------------------------------------
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for output = \"ggplot\". ",
         "Install it with: install.packages(\"ggplot2\")", call. = FALSE)
  }

  fig_df <- score_table
  fig_df$lower <- fig_df$logit_score - ci_multiplier * fig_df$logit_se
  fig_df$upper <- fig_df$logit_score + ci_multiplier * fig_df$logit_se
  # Drop boundary rows where the SE is NA so the figure stays clean
  fig_df <- fig_df[is.finite(fig_df$logit_score) & !is.na(fig_df$logit_se), ,
                   drop = FALSE]

  ggplot2::ggplot(fig_df,
                  ggplot2::aes(x = .data$logit_score, y = .data$raw_score)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$lower, xmax = .data$upper),
      width = error_width, colour = "darkgrey",
      orientation = "y"
    ) +
    ggplot2::geom_point(size = point_size, shape = 18) +
    ggplot2::labs(
      x = "Logit interval score",
      y = "Ordinal sum score"
    ) +
    ggplot2::theme_bw() +
    er2_axis_margins()
}

# ---------------------------------------------------------------------------
# Internal: WLE path (shared Warm-WLE engine)
# ---------------------------------------------------------------------------

#' WLE score-to-theta lookup via the shared Warm-WLE solver
#'
#' Uses the same CML item parameters (centred to grand-mean-zero) and the
#' same `.theta_wle()` engine as [RMpersonParameters()], so the two
#' functions report identical locations and standard errors. The SE is
#' the information-based value `1 / sqrt(I(theta))`.
#'
#' @param data Validated response data.
#' @param theta_range Boundary-score search interval.
#' @return data.frame with columns `raw_score`, `logit_score`, `logit_se`.
#' @keywords internal
#' @noRd
.scoreSE_wle <- function(data, theta_range) {
  fit       <- .rasch_fit_cml(data, se = FALSE)
  thr_list  <- .center_thresholds(fit$thr_list)
  steps     <- vapply(thr_list, length, integer(1))
  max_score <- sum(steps)

  est <- vapply(0:max_score, function(r) {
    .theta_wle(.score_pattern(r, steps), thr_list, theta_range)
  }, numeric(2))

  data.frame(
    raw_score   = 0:max_score,
    logit_score = est[1L, ],
    logit_se    = est[2L, ],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

#' Build a representative response pattern summing to a given raw score
#'
#' For complete data the raw total is the sufficient statistic, so any
#' pattern with that sum yields the same WLE; this greedily fills items
#' up to their category maximum.
#'
#' @param r Target raw score.
#' @param steps Integer vector of per-item threshold counts (category max).
#' @return Integer response vector summing to `r`.
#' @keywords internal
#' @noRd
.score_pattern <- function(r, steps) {
  resp <- integer(length(steps))
  rem  <- r
  for (i in seq_along(steps)) {
    take    <- min(steps[i], rem)
    resp[i] <- take
    rem     <- rem - take
  }
  resp
}

# ---------------------------------------------------------------------------
# Internal: EAP path (mirt EAPsum)
# ---------------------------------------------------------------------------

#' EAPsum score-to-theta lookup via mirt
#'
#' @param data Validated response data.
#' @return data.frame with columns `raw_score`, `logit_score`, `logit_se`.
#' @keywords internal
#' @noRd
.scoreSE_eap <- function(data) {
  if (!requireNamespace("mirt", quietly = TRUE)) {
    stop("Package 'mirt' is required for method = \"EAP\". ",
         "Install it with: install.packages(\"mirt\")", call. = FALSE)
  }

  mirt_fit <- suppressMessages(
    mirt::mirt(
      data       = as.data.frame(data),
      model      = 1,
      itemtype   = "Rasch",
      verbose    = FALSE,
      accelerate = "squarem"
    )
  )

  sscores <- mirt::fscores(
    mirt_fit,
    method         = "EAPsum",
    full.scores    = FALSE,
    full.scores.SE = TRUE,
    verbose        = FALSE
  )
  sscores <- as.data.frame(sscores)

  raw_col   <- intersect(c("Sum.Scores", "Sum.Score"), names(sscores))[1L]
  theta_col <- intersect(c("F1", "Theta", "EAP"),       names(sscores))[1L]
  se_col    <- intersect(c("SE_F1", "SE", "SE_Theta"),  names(sscores))[1L]

  if (is.na(raw_col) || is.na(theta_col) || is.na(se_col)) {
    stop("Unexpected mirt::fscores() output structure; cannot locate ",
         "score / theta / SE columns.", call. = FALSE)
  }

  data.frame(
    raw_score   = as.integer(sscores[[raw_col]]),
    logit_score = as.numeric(sscores[[theta_col]]),
    logit_se    = as.numeric(sscores[[se_col]]),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}
