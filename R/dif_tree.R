# Tree-based DIF analysis with effect-size classification and optional
# stability assessment.
#
# Acknowledgement
# ----------------
# Several internal helpers below (Mantel-Haenszel effect-size with
# iterative purification, partial-gamma effect-size with iterative
# purification, ETS A/B/C classification, partial-gamma A/B/C
# classification, walking the partykit tree to obtain per-split group
# assignments) follow the algorithms implemented by
# Mirka Henninger and Jan Radek in the GitHub-only packages
#   * raschtreeMH (https://github.com/mirka-henninger/raschtreeMH)
#   * effecttree  (https://github.com/mirka-henninger/effecttree)
# both released under the MIT licence:
#
#   Copyright (c) Mirka Henninger and Jan Radek
#
#   Permission is hereby granted, free of charge, to any person
#   obtaining a copy of this software and associated documentation
#   files (the "Software"), to deal in the Software without
#   restriction, including without limitation the rights to use, copy,
#   modify, merge, publish, distribute, sublicense, and/or sell copies
#   of the Software, and to permit persons to whom the Software is
#   furnished to do so, subject to the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND.
#
# The MIT-licensed code from those packages has been adapted here for
# easyRasch2's API and combined with the rest of the package, which is
# distributed under GPL (>= 3) (the combined work is therefore GPL).
#
# References
#  - Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping
#    criterion for Rasch trees based on the Mantel-Haenszel effect size
#    measure for DIF. Educational Psychological Measurement, 83, 181-212.
#    doi:10.1177/00131644221077135
#  - Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial
#    credit trees meet the partial gamma coefficient for quantifying DIF
#    and DSF in polytomous items. Behaviormetrika, 52, 221-257.
#  - Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new
#    method for detecting differential item functioning in the Rasch
#    model. Psychometrika, 80, 289-316.
#  - Komboz, B., Strobl, C., & Zeileis, A. (2018). Tree-based global
#    model tests for polytomous Rasch models. Educational and
#    Psychological Measurement, 78, 128-166.
#  - Philipp, M., Rusch, T., Hornik, K., & Strobl, C. (2018). Measuring
#    the stability of results from supervised statistical learning.
#    Journal of Computational and Graphical Statistics, 27, 685-700.

#' Tree-based DIF analysis with effect-size classification
#'
#' Detects Differential Item Functioning (DIF) by recursively splitting
#' the sample on one or more covariates using \code{psychotree::raschtree()}
#' (dichotomous data) or \code{psychotree::pctree()} (polytomous data),
#' then computes per-split effect-size measures for every item -- the
#' Mantel-Haenszel odds-ratio (in the ETS Delta scale) for dichotomous
#' data, or the partial gamma coefficient for polytomous data -- and
#' classifies them into ETS A/B/C categories. Optionally, an iterative
#' purification step (Holland & Thayer, 1988; Bjorner et al., 1998) is
#' applied, and tree stability across resamples can be assessed via
#' \code{stablelearner::stabletree()}.
#'
#' Continuous covariates (e.g., age in years) and interactions among
#' multiple covariates are handled natively by the model-based
#' recursive-partitioning machinery of \code{partykit} / \code{psychotree}.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item. Items only -- pass
#'   covariates separately via \code{covariates}.
#' @param covariates A data.frame with \code{nrow(data)} rows giving the
#'   DIF covariates. Columns may be numeric (continuous), factor, ordered
#'   factor, or logical. At least one column is required.
#' @param model One of \code{"auto"} (default), \code{"PCM"}, or
#'   \code{"RM"}. \code{"auto"} picks RM when the data are dichotomous
#'   (max response = 1) and PCM otherwise.
#' @param effect_size One of \code{"auto"} (default), \code{"MH"}
#'   (Mantel-Haenszel), or \code{"pgamma"} (partial gamma). \code{"auto"}
#'   uses MH for RM and partial gamma for PCM.
#' @param purification \code{"none"} (default) or \code{"iterative"}.
#'   Iterative purification recomputes the effect size while excluding
#'   items already classified as DIF from the matching score (Bjorner et
#'   al., 1998).
#' @param p_adj Multiple-testing adjustment for the partial-gamma
#'   classification: \code{"none"} (default), \code{"fdr"} (Benjamini &
#'   Hochberg), or \code{"bonferroni"}. Ignored for MH (which uses the
#'   ETS test directly).
#' @param thresholds Numeric length-2 vector with the partial-gamma
#'   B/C boundaries (default \code{c(0.21, 0.31)}). Ignored for MH (the
#'   ETS Delta-scale boundaries 1.0 / 1.5 are used).
#' @param alpha Significance level for the A/C tests. Default
#'   \code{0.05}.
#' @param prune_negligible Logical. If \code{TRUE}, splits whose every
#'   item is classified as A are pruned from the tree. Default
#'   \code{FALSE}.
#' @param stability Logical. If \code{TRUE}, runs
#'   \code{stablelearner::stabletree()} on the fitted tree to assess
#'   variable-selection and cutpoint stability across resamples.
#'   Default \code{FALSE}.
#' @param stability_B Integer. Number of resamples for stability
#'   assessment. Default \code{100}.
#' @param stability_sampler One of \code{"subsampling"} (default) or
#'   \code{"bootstrap"}. Subsampling is recommended for tree stability
#'   because bootstrap resamples can drop levels of categorical
#'   covariates, breaking the model fit on small subsamples.
#' @param min_n_per_level Integer. Minimum count required for any factor
#'   level in \code{covariates}. Default \code{20}. Set to \code{0} to
#'   disable.
#' @param on_rescale One of \code{"message"} (default), \code{"warning"},
#'   or \code{"stop"}. Controls how the function reports the situation
#'   in which one or more items have no responses in their lowest
#'   category within a terminal node of the tree (psychotree silently
#'   rescales those items per node, which leaves the per-node item
#'   parameters on a different metric). \emph{The MH and partial-gamma
#'   effect sizes reported here are computed on the raw responses and
#'   are not affected.} The diagnostic concerns only the per-node item
#'   parameters that are visible via \code{plot(tree)}.
#' @param output One of \code{"kable"} (default), \code{"dataframe"},
#'   \code{"tree"}, or \code{"plot"}. \code{"kable"} renders the
#'   per-split per-item effect-size table as a \code{knitr::kable()};
#'   \code{"dataframe"} returns the underlying tidy data.frame;
#'   \code{"tree"} returns the augmented partykit tree object;
#'   \code{"plot"} returns the partykit tree plot with item names on
#'   the terminal-node x-axis (\code{tp_args = list(names = TRUE)}).
#'   For full control over the plot, use \code{output = "tree"} and
#'   call \code{plot()} on the result with your own \code{tp_args}.
#'   When
#'   \code{stability = TRUE}, the stability summary is attached as
#'   \code{attr(result, "stability")} (a small data.frame) and
#'   \code{attr(result, "stability_kable")} (pre-rendered kable).
#' @param ... Additional arguments forwarded to
#'   \code{psychotree::raschtree()} or \code{psychotree::pctree()}, e.g.
#'   \code{minsize}, \code{alpha} for the parameter-instability test.
#'   Note: this \code{alpha} (the tree-fitting argument) is independent
#'   of the \code{alpha} for ETS classification above; pass it via
#'   \code{...}.
#'
#' @return Depending on \code{output}:
#' \describe{
#'   \item{\code{"kable"}}{A \code{knitr_kable} object with one row per
#'     (split node x item), grouped by node via \code{pack_rows}-style
#'     section headers. The caption summarises model, effect-size
#'     measure, purification, and (if requested) stability.}
#'   \item{\code{"dataframe"}}{A data.frame with one row per (split node
#'     x item): columns \code{NodeID}, \code{Split} (human-readable
#'     description of the split), \code{Variable}, \code{Direction},
#'     \code{Item}, \code{EffectSize}, \code{SE}, \code{Class} (A/B/C),
#'     \code{Flagged} (TRUE for B or C), \code{n_left}, \code{n_right}.}
#'   \item{\code{"tree"}}{The fitted partykit tree object with class
#'     \code{c("RMdifTree", ...)}, with effect-size results stored at
#'     \code{tree$info$effectsize}.}
#'   \item{\code{"plot"}}{A plotted partykit tree.}
#' }
#'
#' Stability results (when \code{stability = TRUE}) are attached to the
#' return value as \code{attr(result, "stability")} (data.frame) and
#' \code{attr(result, "stability_kable")} (pre-rendered kable).
#'
#' @details
#' \strong{ETS Delta-scale Mantel-Haenszel.} The MH common odds ratio
#' \eqn{\hat{\alpha}_{MH}} is mapped to the ETS Delta scale via
#' \eqn{\Delta_{MH} = -2.35 \log \hat{\alpha}_{MH}}. The classification
#' rules (Holland & Thayer, 1988; Zwick, 2012) are:
#' \itemize{
#'   \item Class \strong{A}: \eqn{|\Delta_{MH}| < 1} or test of
#'     \eqn{\Delta_{MH} = 0} not rejected at \code{alpha}.
#'   \item Class \strong{C}: \eqn{|\Delta_{MH}| \ge 1.5} and test of
#'     \eqn{|\Delta_{MH}| \le 1} rejected at \code{alpha}.
#'   \item Class \strong{B}: otherwise.
#' }
#' easyRasch2 uses the sign convention of \code{effecttree}: positive
#' \eqn{\Delta} indicates that the item is more difficult for the second
#' (reference) group.
#'
#' \strong{Partial gamma.} \code{iarm::partgam_DIF()} provides the gamma
#' estimate and SE per item. ETS-style classification follows Bjorner et
#' al. (1998): A if \eqn{|\gamma| < 0.21} or the test of \eqn{\gamma = 0}
#' is not rejected at \code{alpha}; C if \eqn{|\gamma| > 0.31} and the
#' test of \eqn{|\gamma| \le 0.21} is rejected at \code{alpha}; B
#' otherwise.
#'
#' \strong{Continuous and interaction effects.} The recursive
#' partitioning step automatically handles continuous covariates (the
#' parameter-instability test searches for the optimal cutpoint) and
#' multivariate interactions (later splits are conditional on earlier
#' ones). To assess sensitivity of variable selection and cutpoints to
#' the particular sample, set \code{stability = TRUE}.
#'
#' \strong{Stability assessment.} When \code{stability = TRUE}, the same
#' tree-fitting call is replayed on \code{stability_B} resamples of the
#' data. The returned \code{stabletree} object reports, per covariate,
#' the proportion of resamples in which it was selected at any split,
#' and (for continuous covariates) the empirical distribution of
#' cutpoints. Stability is independent of effect-size classification --
#' a stable split with a negligible effect is still negligible, and a
#' large effect on an unstable split should be interpreted with care.
#' Cost is roughly \code{stability_B} times the original fitting time.
#'
#' \strong{Acknowledgement.} The effect-size machinery -- MH and partial
#' gamma per node, iterative purification, ETS A/B/C classification --
#' follows the implementations by Henninger & Radek in the GitHub
#' packages \emph{raschtreeMH} and \emph{effecttree}. The relevant code
#' has been adapted here under the MIT licence (see file header).
#'
#' @references
#' Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech,
#' P. (1998). Differential item functioning in the Danish translation
#' of the SF-36. \emph{Journal of Clinical Epidemiology, 51}(11),
#' 1189-1202.
#'
#' Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping
#' criterion for Rasch trees based on the Mantel-Haenszel effect size
#' measure for DIF. \emph{Educational and Psychological Measurement,
#' 83}, 181-212. \doi{10.1177/00131644221077135}
#'
#' Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial
#' credit trees meet the partial gamma coefficient for quantifying DIF
#' and DSF in polytomous items. \emph{Behaviormetrika, 52}, 221-257.
#'
#' Philipp, M., Rusch, T., Hornik, K., & Strobl, C. (2018). Measuring
#' the stability of results from supervised statistical learning.
#' \emph{Journal of Computational and Graphical Statistics, 27},
#' 685-700.
#'
#' @seealso \code{\link{RMdifLR}}, \code{\link{RMdifGamma}}
#'
#' @examples
#' \donttest{
#' if (requireNamespace("psychotree", quietly = TRUE) &&
#'     requireNamespace("partykit", quietly = TRUE) &&
#'     requireNamespace("iarm", quietly = TRUE)) {
#'   data("DIFSimPC", package = "psychotree")
#'
#'   items <- as.data.frame(as.matrix(DIFSimPC$resp))
#'   covs  <- DIFSimPC[, c("age", "gender", "motivation")]
#'
#'   # Default: kable of per-split effect sizes
#'   RMdifTree(items, covariates = covs)
#'
#'   # Tidy data.frame -- one row per (split node x item)
#'   df <- RMdifTree(items, covariates = covs, output = "dataframe")
#'   df[df$Flagged, ]
#'
#'   # Tree object (for plotting via partykit::plot)
#'   tree <- RMdifTree(items, covariates = covs, output = "tree")
#'   plot(tree)
#' }
#' }
#' \dontrun{
#' # Stability assessment refits the tree on B resamples -- slow.
#' kbl <- RMdifTree(items, covariates = covs,
#'                  purification = "iterative", p_adj = "fdr",
#'                  stability = TRUE, stability_B = 100)
#' kbl                          # main effect-size kable
#' attr(kbl, "stability_kable") # pre-rendered stability kable
#' attr(kbl, "stability")       # raw stability data.frame
#' }
#'
#' @export
RMdifTree <- function(data,
                      covariates,
                      model             = c("auto", "PCM", "RM"),
                      effect_size       = c("auto", "MH", "pgamma"),
                      purification      = c("none", "iterative"),
                      p_adj             = c("none", "fdr", "bonferroni"),
                      thresholds        = c(0.21, 0.31),
                      alpha             = 0.05,
                      prune_negligible  = FALSE,
                      stability         = FALSE,
                      stability_B       = 100L,
                      stability_sampler = c("subsampling", "bootstrap"),
                      min_n_per_level   = 20L,
                      on_rescale        = c("message", "warning", "stop"),
                      output            = c("kable", "dataframe", "tree", "plot"),
                      ...) {

  # Capture the unevaluated `covariates` expression so we can recover a
  # readable column name when the user passes a single vector like
  # RMdifTree(d, dif$age) instead of a data.frame.
  covariates_expr <- substitute(covariates)

  model             <- match.arg(model)
  effect_size       <- match.arg(effect_size)
  purification      <- match.arg(purification)
  p_adj             <- match.arg(p_adj)
  stability_sampler <- match.arg(stability_sampler)
  on_rescale        <- match.arg(on_rescale)
  output            <- match.arg(output)

  # ---------------------------------------------------------------------
  # Package availability
  # ---------------------------------------------------------------------
  if (!requireNamespace("psychotree", quietly = TRUE)) {
    stop("Package 'psychotree' is required for RMdifTree(). ",
         "Install it with: install.packages(\"psychotree\")",
         call. = FALSE)
  }
  if (!requireNamespace("partykit", quietly = TRUE)) {
    stop("Package 'partykit' is required for RMdifTree().", call. = FALSE)
  }
  if (!requireNamespace("difR", quietly = TRUE)) {
    stop("Package 'difR' is required for RMdifTree() (Mantel-Haenszel ",
         "effect size). Install it with: install.packages(\"difR\")",
         call. = FALSE)
  }
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop("Package 'iarm' is required for RMdifTree() (partial-gamma ",
         "effect size). Install it with: install.packages(\"iarm\")",
         call. = FALSE)
  }
  if (stability && !requireNamespace("stablelearner", quietly = TRUE)) {
    stop("Package 'stablelearner' is required for stability = TRUE. ",
         "Install it with: install.packages(\"stablelearner\")",
         call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Validate items
  # ---------------------------------------------------------------------
  validate_response_data(data)
  data <- as.data.frame(data)

  # ---------------------------------------------------------------------
  # Validate covariates
  # ---------------------------------------------------------------------
  if (missing(covariates) || is.null(covariates)) {
    stop("`covariates` is required: a data.frame of one or more DIF ",
         "covariates with nrow(data) rows.", call. = FALSE)
  }
  if (!is.data.frame(covariates)) {
    if (is.atomic(covariates) || is.list(covariates)) {
      derived_name <- derive_covariate_name(covariates_expr)
      covariates <- data.frame(x = covariates, stringsAsFactors = FALSE)
      names(covariates) <- derived_name
    } else {
      stop("`covariates` must be a data.frame.", call. = FALSE)
    }
  }
  if (nrow(covariates) != nrow(data)) {
    stop("`covariates` must have the same number of rows as `data` (",
         nrow(data), "). Got ", nrow(covariates), ".", call. = FALSE)
  }
  if (ncol(covariates) < 1L) {
    stop("`covariates` must contain at least one column.", call. = FALSE)
  }
  validate_covariates(covariates, min_n_per_level = min_n_per_level)

  # ---------------------------------------------------------------------
  # Resolve model and effect_size
  # ---------------------------------------------------------------------
  is_polytomous <- max(as.matrix(data), na.rm = TRUE) > 1L

  resolved_model <- if (model == "auto") {
    if (is_polytomous) "PCM" else "RM"
  } else {
    if (model == "RM" && is_polytomous) {
      stop("`model = \"RM\"` requires dichotomous data; use \"PCM\" or ",
           "\"auto\".", call. = FALSE)
    }
    model
  }

  resolved_es <- if (effect_size == "auto") {
    if (resolved_model == "RM") "MH" else "pgamma"
  } else {
    if (effect_size == "MH" && resolved_model == "PCM") {
      stop("`effect_size = \"MH\"` is only defined for dichotomous data ",
           "(model = \"RM\"). Use \"pgamma\" or \"auto\".", call. = FALSE)
    }
    effect_size
  }

  # ---------------------------------------------------------------------
  # Validate effect-size args
  # ---------------------------------------------------------------------
  if (!is.numeric(thresholds) || length(thresholds) != 2L ||
      any(!is.finite(thresholds)) || thresholds[1L] >= thresholds[2L] ||
      thresholds[1L] <= 0) {
    stop("`thresholds` must be a length-2 numeric vector with ",
         "0 < thresholds[1] < thresholds[2].", call. = FALSE)
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 ||
      alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(stability_B) || length(stability_B) != 1L ||
      stability_B < 1) {
    stop("`stability_B` must be a positive integer.", call. = FALSE)
  }
  stability_B <- as.integer(stability_B)

  # ---------------------------------------------------------------------
  # Build the combined data.frame in the form psychotree expects:
  # one column called `.items` of class matrix on the LHS, covariates on
  # the RHS. Drop rows where any covariate is NA (psychotree handles
  # item NAs itself).
  # ---------------------------------------------------------------------
  na_cov <- !stats::complete.cases(covariates)
  if (any(na_cov)) {
    message(sum(na_cov), " row(s) with NA in `covariates` dropped.")
    data       <- data[!na_cov, , drop = FALSE]
    covariates <- covariates[!na_cov, , drop = FALSE]
  }
  if (nrow(data) < 30L) {
    stop("Need at least 30 complete cases after NA handling. Got ",
         nrow(data), ".", call. = FALSE)
  }

  items_mat <- as.matrix(data)
  storage.mode(items_mat) <- "integer"
  combined  <- covariates
  combined$.items <- items_mat

  fmla <- stats::reformulate(
    termlabels = names(covariates),
    response   = ".items"
  )

  # ---------------------------------------------------------------------
  # Fit the tree. Use do.call so stablelearner::stabletree() can later
  # update() the call cleanly.
  # ---------------------------------------------------------------------
  fit_fun <- if (resolved_model == "RM") {
    psychotree::raschtree
  } else {
    psychotree::pctree
  }

  tree_call <- as.call(c(
    list(if (resolved_model == "RM")
      quote(psychotree::raschtree)
    else
      quote(psychotree::pctree),
      fmla, quote(combined)),
    list(...)
  ))
  names(tree_call)[2:3] <- c("formula", "data")

  # The "Minimum score is not zero" warning fires repeatedly during the
  # tree's parameter-instability search (psychotree refits pcmodel many
  # times). Capturing it here just to muffle the noise -- we inspect the
  # *final* tree's terminal nodes below to detect the actual condition.
  tree <- withCallingHandlers(
    eval(tree_call, envir = list(combined = combined)),
    warning = function(w) {
      if (grepl("Minimum score is not zero",
                conditionMessage(w), fixed = TRUE)) {
        invokeRestart("muffleWarning")
      }
      # Non-matching warnings: let them propagate normally.
    }
  )

  # ---------------------------------------------------------------------
  # Compute per-split effect sizes
  # ---------------------------------------------------------------------
  effect_info <- compute_tree_effectsize(
    tree         = tree,
    items        = items_mat,
    es_method    = resolved_es,
    purification = purification,
    p_adj        = p_adj,
    thresholds   = thresholds,
    alpha        = alpha
  )

  tree$info$effectsize <- effect_info

  # ---------------------------------------------------------------------
  # Optional pruning of negligible (all-A) splits
  # ---------------------------------------------------------------------
  if (prune_negligible && length(effect_info$nodes) > 0L) {
    to_prune <- vapply(effect_info$nodes, function(nd) {
      all(nd$Class == "A", na.rm = TRUE)
    }, logical(1L))
    prune_ids <- as.integer(names(effect_info$nodes)[to_prune])
    if (length(prune_ids) > 0L) {
      tree <- partykit::nodeprune(tree, ids = prune_ids)
      # Recompute on the pruned tree -- splits and ids may have shifted
      effect_info <- compute_tree_effectsize(
        tree         = tree,
        items        = items_mat,
        es_method    = resolved_es,
        purification = purification,
        p_adj        = p_adj,
        thresholds   = thresholds,
        alpha        = alpha
      )
      tree$info$effectsize <- effect_info
    }
  }

  # ---------------------------------------------------------------------
  # Diagnostic: detect items rescaled within terminal nodes
  # ---------------------------------------------------------------------
  rescale_info <- detect_rescaled_items(tree, items_mat)
  if (rescale_info$n_terminal_cells > 0L) {
    msg <- format_rescale_message(rescale_info)
    switch(on_rescale,
           message = message(msg),
           warning = warning(msg, call. = FALSE),
           stop    = stop(msg, call. = FALSE))
  }

  # ---------------------------------------------------------------------
  # Stability assessment
  # ---------------------------------------------------------------------
  # Resample-level fits inside stabletree() can produce three flavours of
  # noise that we don't want to leak to the user (each fires repeatedly
  # during the parameter-instability search but is not separately
  # actionable -- the user has already been told about rescaling in the
  # *original* tree above):
  #
  #   1. "Minimum score is not zero ..."   -- pcmodel rescaling warning
  #   2. "items with null categories ..."  -- interior categories empty
  #   3. "Error in chol.default(meat) ..." -- M-fluctuation test cannot
  #      invert a rank-deficient score covariance (printed via try()
  #      inside partykit::mob; non-fatal, the affected split is skipped).
  #
  # We muffle the warnings via withCallingHandlers and capture stderr
  # via sink() to swallow the printed try() error text. After the call
  # returns we emit a single summary message so the user knows what
  # happened during the resample fits.
  stability_result <- NULL
  if (stability) {
    sampler <- if (stability_sampler == "subsampling") {
      stablelearner::subsampling
    } else {
      stablelearner::bootstrap
    }

    n_min_zero  <- 0L
    n_null_cat  <- 0L

    err_file <- tempfile(fileext = ".log")
    err_con  <- file(err_file, open = "wt")
    sink(err_con, type = "message")

    stability_result <- tryCatch(
      withCallingHandlers(
        stablelearner::stabletree(tree, B = stability_B, sampler = sampler),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (grepl("Minimum score is not zero", msg, fixed = TRUE)) {
            n_min_zero <<- n_min_zero + 1L
            invokeRestart("muffleWarning")
          } else if (grepl("items with null categories", msg, fixed = TRUE)) {
            n_null_cat <<- n_null_cat + 1L
            invokeRestart("muffleWarning")
          }
        }
      ),
      finally = {
        sink(NULL, type = "message")
        close(err_con)
      }
    )

    err_lines <- tryCatch(readLines(err_file, warn = FALSE),
                          error = function(e) character())
    unlink(err_file)
    n_resample_errors <- sum(grepl("^Error", err_lines))

    # Single user-facing summary, only if anything noteworthy happened.
    summary_parts <- character()
    if (n_min_zero > 0L) {
      summary_parts <- c(summary_parts,
        paste0(n_min_zero, " 'minimum score not zero' rescaling warning(s)"))
    }
    if (n_null_cat > 0L) {
      summary_parts <- c(summary_parts,
        paste0(n_null_cat, " 'items with null categories' warning(s)"))
    }
    if (n_resample_errors > 0L) {
      summary_parts <- c(summary_parts,
        paste0(n_resample_errors,
               " resample fit(s) errored out (rank-deficient score ",
               "covariance in the M-fluctuation test); affected splits ",
               "skipped"))
    }
    if (length(summary_parts) > 0L) {
      message(
        "Stability assessment notes (suppressed during fitting): ",
        paste(summary_parts, collapse = "; "),
        ". These reflect sparsity in some resamples and are normal ",
        "sampling variability; selection-frequency and cutpoint summaries ",
        "are based on resamples that completed successfully."
      )
    }

    attr(tree, "stability") <- stability_result
  }

  class(tree) <- unique(c("RMdifTree", class(tree)))

  # ---------------------------------------------------------------------
  # Build the per-split per-item dataframe (used by all non-tree outputs)
  # ---------------------------------------------------------------------
  df <- effectsize_to_dataframe(
    effect_info,
    item_names         = colnames(items_mat),
    tree               = tree,
    rescaled_terminals = rescale_info$per_terminal
  )

  stab_df <- NULL
  stab_kbl <- NULL
  if (!is.null(stability_result)) {
    stab_df  <- stabletree_summary_df(stability_result, covariates)
    stab_kbl <- stabletree_summary_kable(stab_df)
  }

  attr(df, "model")       <- resolved_model
  attr(df, "effect_size") <- resolved_es
  if (!is.null(stab_df))   attr(df, "stability")              <- stab_df
  if (!is.null(stab_kbl))  attr(df, "stability_kable")        <- stab_kbl
  if (rescale_info$n_terminal_cells > 0L) {
    attr(df, "rescaled_terminal_items") <- rescale_info$per_terminal
  }

  # ---------------------------------------------------------------------
  # Dispatch on output
  # ---------------------------------------------------------------------
  if (output == "tree") {
    if (!is.null(stab_df))  attr(tree, "stability")       <- stab_df
    if (!is.null(stab_kbl)) attr(tree, "stability_kable") <- stab_kbl
    if (rescale_info$n_terminal_cells > 0L) {
      attr(tree, "rescaled_terminal_items") <- rescale_info$per_terminal
    }
    return(tree)
  }
  if (output == "plot") {
    # Show item names (rather than 1..I numeric indices) on the
    # terminal-node x-axis. Users wanting full control can use
    # `output = "tree"` and call `plot()` themselves with custom
    # `tp_args`.
    return(graphics::plot(tree, tp_args = list(names = TRUE)))
  }
  if (output == "dataframe") {
    return(df)
  }

  # output == "kable"
  kbl <- render_effectsize_kable(
    df,
    n_persons        = nrow(combined),
    n_items          = ncol(items_mat),
    es_method        = resolved_es,
    purification     = purification,
    p_adj            = p_adj,
    thresholds       = thresholds,
    alpha            = alpha,
    is_polytomous    = is_polytomous,
    has_stability    = !is.null(stab_df),
    n_rescaled_items = length(unique(rescale_info$per_terminal$Item))
  )
  if (!is.null(stab_df))  attr(kbl, "stability")       <- stab_df
  if (!is.null(stab_kbl)) attr(kbl, "stability_kable") <- stab_kbl
  if (rescale_info$n_terminal_cells > 0L) {
    attr(kbl, "rescaled_terminal_items") <- rescale_info$per_terminal
  }
  attr(kbl, "model")       <- resolved_model
  attr(kbl, "effect_size") <- resolved_es
  kbl
}

# =====================================================================
# Internal helpers
# =====================================================================

#' Derive a readable column name from an unevaluated `covariates` expression
#'
#' When the user passes a single vector -- e.g.
#' \code{RMdifTree(d, dif$age)} -- \code{as.data.frame(covariates)} would
#' name the column "covariates" (the parameter name). To keep the plot,
#' kable, and dataframe outputs informative we recover the user's
#' expression here and produce a clean column name:
#' * a bare symbol \code{age} -> "age"
#' * \code{dif$age} or \code{dif[["age"]]} -> "age"
#' * \code{dif[, "age"]} -> "age"
#' * anything else -> the deparsed expression as-is.
#'
#' If the resulting name would collide with R reserved syntax, a
#' fallback "covariate" is used.
#'
#' @keywords internal
#' @noRd
derive_covariate_name <- function(expr) {
  bad_fallback <- "covariate"
  if (is.symbol(expr)) {
    nm <- as.character(expr)
    return(if (nzchar(nm)) nm else bad_fallback)
  }
  if (is.call(expr) && length(expr) >= 3L) {
    op <- expr[[1L]]
    # x$name
    if (identical(op, as.name("$"))) {
      return(as.character(expr[[3L]]))
    }
    # x[["name"]] or x[[i]]
    if (identical(op, as.name("[["))) {
      key <- expr[[3L]]
      if (is.character(key)) return(key)
      if (is.symbol(key))    return(as.character(key))
    }
    # x[, "name"] or x[i]
    if (identical(op, as.name("["))) {
      # Last index argument is the column selector
      n <- length(expr)
      key <- expr[[n]]
      if (is.character(key)) return(key)
      if (is.symbol(key))    return(as.character(key))
    }
  }
  out <- tryCatch(deparse(expr, width.cutoff = 60L)[1L],
                  error = function(e) bad_fallback)
  if (is.null(out) || !nzchar(out)) bad_fallback else out
}

#' Validate the covariates data.frame
#'
#' Catches all-NA columns, zero-variance columns, and (for factors) too
#' few cases per level.
#'
#' @keywords internal
#' @noRd
validate_covariates <- function(covariates, min_n_per_level) {
  for (nm in names(covariates)) {
    col <- covariates[[nm]]
    if (all(is.na(col))) {
      stop("Covariate `", nm, "` is all NA.", call. = FALSE)
    }
    if (is.character(col)) {
      covariates[[nm]] <- as.factor(col)
      col <- covariates[[nm]]
    }
    if (is.logical(col)) {
      covariates[[nm]] <- as.factor(col)
      col <- covariates[[nm]]
    }
    if (is.factor(col)) {
      col <- droplevels(col)
      tab <- table(col, useNA = "no")
      if (length(tab) < 2L) {
        stop("Covariate `", nm, "` has fewer than 2 non-NA levels.",
             call. = FALSE)
      }
      if (min_n_per_level > 0L && any(tab < min_n_per_level)) {
        small <- names(tab)[tab < min_n_per_level]
        stop("Covariate `", nm, "` has level(s) with fewer than ",
             min_n_per_level, " cases: ",
             paste(sprintf("%s (n = %d)", small, tab[small]),
                   collapse = ", "),
             ". Combine sparse levels or set `min_n_per_level = 0`.",
             call. = FALSE)
      }
    } else if (is.numeric(col)) {
      if (length(unique(stats::na.omit(col))) < 2L) {
        stop("Covariate `", nm, "` has fewer than 2 unique values.",
             call. = FALSE)
      }
    } else {
      stop("Covariate `", nm, "` has unsupported type ",
           paste(class(col), collapse = "/"),
           ". Use numeric, factor, ordered factor, or logical.",
           call. = FALSE)
    }
  }
  invisible(TRUE)
}

#' Compute effect sizes for every inner (split) node of a tree
#'
#' Walks the partykit tree, gets the binary group assignment defined by
#' the descendants of each inner node, and computes per-item effect
#' sizes on that split.
#'
#' Adapted from \code{effecttree::get_effectsize()} (MIT, Henninger &
#' Radek), with simplifications: we only support the binary-split case
#' that \code{psychotree::raschtree}/\code{pctree} produce, and we
#' return a flat list keyed by node id.
#'
#' @keywords internal
#' @noRd
compute_tree_effectsize <- function(tree, items, es_method, purification,
                                    p_adj, thresholds, alpha) {

  inner_ids <- setdiff(
    partykit::nodeids(tree, terminal = FALSE),
    partykit::nodeids(tree, terminal = TRUE)
  )

  if (length(inner_ids) == 0L) {
    return(list(
      es_method    = es_method,
      purification = purification,
      p_adj        = p_adj,
      thresholds   = thresholds,
      alpha        = alpha,
      nodes        = list(),
      n_splits     = 0L
    ))
  }

  dat_fitted <- partykit::data_party(tree)
  fitted_id  <- dat_fitted[["(fitted)"]]
  root_node  <- partykit::node_party(tree)

  split_groups <- lapply(inner_ids, function(id) {
    nd <- find_node_by_id(root_node, id)
    if (is.null(nd)) {
      stop("Internal: could not locate node ", id, " in tree.",
           call. = FALSE)
    }
    children <- partykit::kids_node(nd)
    left_terminals  <- collect_terminal_ids(children[[1L]])
    right_terminals <- collect_terminal_ids(children[[2L]])
    grp <- rep(NA_integer_, length(fitted_id))
    grp[fitted_id %in% left_terminals]  <- 0L
    grp[fitted_id %in% right_terminals] <- 1L
    grp
  })
  names(split_groups) <- as.character(inner_ids)

  sums <- rowSums(items, na.rm = TRUE)

  per_node <- lapply(seq_along(inner_ids), function(i) {
    nid <- inner_ids[i]
    grp <- split_groups[[i]]
    split_meta <- describe_split(tree, nid)

    result <- if (es_method == "MH") {
      one_split_mh(items = items, group = grp, sums = sums,
                   purification = purification, alpha = alpha)
    } else {
      one_split_pgamma(items = items, group = grp,
                       purification = purification,
                       p_adj = p_adj, thresholds = thresholds,
                       alpha = alpha)
    }

    list(
      NodeID     = nid,
      Variable   = split_meta$variable,
      Direction  = split_meta$direction,
      Split      = split_meta$split,
      EffectSize = result$es,
      SE         = result$se,
      Class      = result$class,
      n_left     = sum(grp == 0L, na.rm = TRUE),
      n_right    = sum(grp == 1L, na.rm = TRUE),
      purification_iterations = result$purification_iterations
    )
  })
  names(per_node) <- as.character(inner_ids)

  list(
    es_method    = es_method,
    purification = purification,
    p_adj        = p_adj,
    thresholds   = thresholds,
    alpha        = alpha,
    nodes        = per_node,
    n_splits     = length(per_node)
  )
}

#' Recursively collect terminal node ids descended from a node
#'
#' Adapted from \code{effecttree::get_terminal_nodes} (MIT).
#'
#' @keywords internal
#' @noRd
collect_terminal_ids <- function(node) {
  kids <- partykit::kids_node(node)
  if (is.null(kids)) return(partykit::id_node(node))
  unlist(lapply(kids, collect_terminal_ids), use.names = FALSE)
}

#' Detect which (terminal node x item) cells had no responses in
#' category 0 in the final fitted tree
#'
#' \code{psychotree::pcmodel} silently rescales such items to start at
#' the new lowest observed category, producing per-node parameters on a
#' shifted metric. We detect this post-hoc by scanning each terminal
#' node's data.
#'
#' @return A list with three components:
#'   * \code{per_terminal}: data.frame with columns \code{NodeID}
#'     (terminal-node id) and \code{Item}, one row per affected cell.
#'   * \code{per_item}: named integer vector -- for each affected item,
#'     how many terminal nodes had it rescaled.
#'   * \code{n_terminal_cells}: total count of (terminal x item)
#'     affected cells.
#'
#' @keywords internal
#' @noRd
detect_rescaled_items <- function(tree, items_mat) {
  dat_fitted <- partykit::data_party(tree)
  fitted_id  <- dat_fitted[["(fitted)"]]
  terminal_ids <- partykit::nodeids(tree, terminal = TRUE)
  item_names <- colnames(items_mat)
  if (is.null(item_names)) item_names <- paste0("V", seq_len(ncol(items_mat)))

  affected_rows <- list()
  for (tid in terminal_ids) {
    rows <- which(fitted_id == tid)
    if (length(rows) == 0L) next
    sub <- items_mat[rows, , drop = FALSE]
    item_min <- apply(sub, 2L, function(x) {
      x <- x[!is.na(x)]
      if (length(x) == 0L) NA_real_ else min(x)
    })
    bad_items <- item_names[!is.na(item_min) & item_min > 0]
    if (length(bad_items) > 0L) {
      affected_rows[[length(affected_rows) + 1L]] <- data.frame(
        NodeID = rep(tid, length(bad_items)),
        Item   = bad_items,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(affected_rows) == 0L) {
    return(list(
      per_terminal     = data.frame(NodeID = integer(), Item = character(),
                                    stringsAsFactors = FALSE),
      per_item         = integer(0),
      n_terminal_cells = 0L
    ))
  }

  per_terminal <- do.call(rbind, affected_rows)
  rownames(per_terminal) <- NULL
  per_item <- table(per_terminal$Item)

  list(
    per_terminal     = per_terminal,
    per_item         = stats::setNames(as.integer(per_item), names(per_item)),
    n_terminal_cells = nrow(per_terminal)
  )
}

#' Build a single user-facing message about rescaled items
#'
#' Item-level: "Item X was rescaled in K of T terminal nodes".
#'
#' @keywords internal
#' @noRd
format_rescale_message <- function(rescale_info) {
  pi  <- rescale_info$per_item
  T_n <- length(unique(rescale_info$per_terminal$NodeID))
  n_aff_terms <- length(unique(rescale_info$per_terminal$NodeID))
  per_item_lines <- vapply(seq_along(pi), function(i) {
    paste0("  - ", names(pi)[i],
           ": rescaled in ", pi[i], " terminal node",
           if (pi[i] != 1L) "s" else "")
  }, character(1L))

  paste0(
    "Some items have no responses in their lowest category within ",
    n_aff_terms,
    " terminal node",
    if (n_aff_terms != 1L) "s" else "",
    " of the fitted tree. ",
    "psychotree silently rescales those items per node, so the ",
    "per-node item parameters (visible via plot(tree)) are on a ",
    "shifted metric for the affected cells. The MH / partial-gamma ",
    "effect sizes reported here are computed on raw responses and ",
    "are NOT affected.\n",
    "Affected items:\n",
    paste(per_item_lines, collapse = "\n"),
    "\nDiagnose with: RMplotTile(data, group = <covariate>) for each ",
    "covariate. To suppress, increase `minsize`. To downgrade or ",
    "escalate this notice see `?RMdifTree` argument `on_rescale`."
  )
}

#' Walk a partykit node tree to find the node with the given id
#'
#' \code{tree[[id]]} on a fitted partykit tree returns a *renumbered*
#' subtree (terminals start at 2, 3, ...), so we cannot use it to map
#' fitted ids back to subtree membership. Walking the root node
#' directly preserves the original ids.
#'
#' @keywords internal
#' @noRd
find_node_by_id <- function(node, id) {
  if (partykit::id_node(node) == id) return(node)
  kids <- partykit::kids_node(node)
  if (is.null(kids)) return(NULL)
  for (kid in kids) {
    found <- find_node_by_id(kid, id)
    if (!is.null(found)) return(found)
  }
  NULL
}

#' Describe a split (variable name + cutpoint or level set)
#'
#' @keywords internal
#' @noRd
describe_split <- function(tree, node_id) {
  na_out <- list(variable = NA_character_, direction = NA_character_,
                 split = NA_character_)
  nd <- find_node_by_id(partykit::node_party(tree), node_id)
  if (is.null(nd)) return(na_out)
  sp <- partykit::split_node(nd)
  if (is.null(sp)) return(na_out)

  varid   <- sp$varid
  varname <- names(partykit::data_party(tree))[varid]
  bp      <- sp$breaks
  idx     <- sp$index

  if (!is.null(bp)) {
    cut_str <- signif(bp, 4)
    return(list(
      variable  = varname,
      direction = paste0("<= ", cut_str),
      split     = paste0(varname, " <= ", cut_str)
    ))
  }
  if (!is.null(idx)) {
    var_col <- partykit::data_party(tree)[[varid]]
    levs    <- levels(var_col)
    if (is.null(levs)) {
      return(list(variable = varname, direction = "split",
                  split = paste0(varname, ": split")))
    }
    left  <- levs[idx == 1L]
    right <- levs[idx == 2L]
    left_str  <- if (length(left)  == 1L) left
                 else paste0("{", paste(left,  collapse = ", "), "}")
    right_str <- if (length(right) == 1L) right
                 else paste0("{", paste(right, collapse = ", "), "}")
    return(list(
      variable  = varname,
      direction = paste0("= ", left_str),
      split     = paste0(varname, ": ", left_str, " vs ", right_str)
    ))
  }
  na_out
}

#' Mantel-Haenszel effect-size for one binary split
#'
#' Adapted from \code{effecttree::calculate_mantelhaenszel} (MIT) with
#' the easyRasch2-style alpha argument.
#'
#' @keywords internal
#' @noRd
one_split_mh <- function(items, group, sums, purification, alpha) {
  keep  <- !is.na(group)
  dat   <- as.data.frame(items[keep, , drop = FALSE])
  grp   <- group[keep]
  smm   <- sums[keep]

  res     <- difR::mantelHaenszel(data = dat, member = grp, match = smm)
  delta   <- -2.35 * log(res$resAlpha)
  delta_sd<- 2.35 * sqrt(res$varLambda)
  iters   <- 0L

  if (purification == "iterative") {
    which_dif <- which(abs(delta) >= 1)
    if (length(which_dif) > 0L) {
      repeat {
        iters <- iters + 1L
        pur <- purify_mh_step(dat, grp, delta)
        delta    <- pur$delta
        delta_sd <- pur$delta_sd
        new_dif  <- which(abs(delta) >= 1)
        if (setequal(which_dif, new_dif) || iters >= ncol(dat)) break
        which_dif <- new_dif
      }
    }
  }

  cls <- ets_classify(delta, delta_sd, alpha = alpha)
  list(es                       = unname(delta),
       se                       = unname(delta_sd),
       class                    = cls,
       purification_iterations  = iters)
}

#' One step of Holland-Thayer iterative purification for MH
#'
#' Items currently classified as DIF are excluded from the matching
#' total when re-computing MH for the non-DIF items. For DIF items,
#' the purified MH includes the item itself in the matching total
#' (Bjorner et al., 1998, p. 1191; Henninger et al., 2023).
#'
#' Adapted from \code{effecttree::purify_MH} (MIT).
#'
#' @keywords internal
#' @noRd
purify_mh_step <- function(dat, group, delta) {
  no_dif  <- which(abs(delta) < 1)
  dif_set <- which(abs(delta) >= 1)

  if (length(no_dif) < 2L) {
    return(list(delta = delta,
                delta_sd = rep(NA_real_, length(delta))))
  }

  base_match <- rowSums(dat[, no_dif, drop = FALSE], na.rm = TRUE)
  res        <- difR::mantelHaenszel(data = dat, member = group,
                                     match = base_match)
  new_delta  <- -2.35 * log(res$resAlpha)
  new_sd     <- 2.35 * sqrt(res$varLambda)

  # Recompute DIF items adding their own score back in
  for (i in dif_set) {
    smm  <- dat[, i] + base_match
    res_i <- difR::mantelHaenszel(data = dat[, i, drop = FALSE],
                                  member = group, match = smm)
    new_delta[i] <- -2.35 * log(res_i$resAlpha)
    new_sd[i]    <- 2.35 * sqrt(res_i$varLambda)
  }

  list(delta = new_delta, delta_sd = new_sd)
}

#' ETS A/B/C classification for MH Delta-scale effect sizes
#'
#' Two non-central chi-square tests on \eqn{(\Delta/SD)^2}:
#' tau = 0 for the A-test, tau = 1 for the C-test (Paek & Holland,
#' 2015).
#'
#' Adapted from \code{effecttree::get_ETS_classification} (MIT).
#'
#' @keywords internal
#' @noRd
ets_classify <- function(delta, delta_sd, alpha = 0.05) {
  pval_at <- function(d, sd, tau) {
    lambda <- (tau / sd)^2
    psi    <- (d / sd)^2
    1 - stats::pchisq(psi, df = 1, ncp = lambda)
  }
  p_A <- pval_at(delta, delta_sd, tau = 0)
  p_C <- pval_at(delta, delta_sd, tau = 1)
  ifelse(abs(delta) < 1 | p_A >= alpha, "A",
         ifelse(abs(delta) >= 1.5 & p_C < alpha, "C", "B"))
}

#' Partial-gamma effect size for one binary split
#'
#' Adapted from \code{effecttree::calculate_pgamma} (MIT).
#'
#' @keywords internal
#' @noRd
one_split_pgamma <- function(items, group, purification, p_adj,
                             thresholds, alpha) {
  keep <- !is.na(group)
  dat  <- as.data.frame(items[keep, , drop = FALSE])
  grp  <- group[keep]

  res   <- quiet_call(iarm::partgam_DIF(dat.items = dat, dat.exo = grp))
  gam   <- res$gamma
  sse   <- res$se
  iters <- 0L

  if (purification == "iterative") {
    which_dif <- which(abs(gam) >= thresholds[1L])
    if (length(which_dif) > 0L) {
      repeat {
        iters <- iters + 1L
        pur <- purify_pgamma_step(dat, grp, gam, sse, thresholds)
        if (is.null(pur)) break
        gam <- pur$gamma
        sse <- pur$se
        new_dif <- which(abs(gam) >= thresholds[1L])
        if (setequal(which_dif, new_dif) || iters >= ncol(dat)) break
        which_dif <- new_dif
      }
    }
  }

  cls <- pgamma_classify(gam, sse, p_adj = p_adj,
                         thresholds = thresholds, alpha = alpha)
  list(es                      = unname(gam),
       se                      = unname(sse),
       class                   = cls,
       purification_iterations = iters)
}

#' One step of iterative purification for partial gamma
#'
#' Adapted from \code{effecttree::purify_pgamma} (MIT). Returns NULL if
#' fewer than 2 non-DIF items remain (purification ill-defined).
#'
#' @keywords internal
#' @noRd
purify_pgamma_step <- function(dat, group, gam, sse, thresholds) {
  no_dif  <- which(abs(gam) < thresholds[1L])
  dif_set <- which(abs(gam) >= thresholds[1L])
  if (length(no_dif) < 2L) return(NULL)

  new_gam <- gam
  new_sse <- sse

  base <- quiet_call(
    iarm::partgam_DIF(dat[, no_dif, drop = FALSE], group)
  )
  new_gam[no_dif] <- base$gamma
  new_sse[no_dif] <- base$se

  for (i in dif_set) {
    cols <- sort(unique(c(i, no_dif)))
    pos  <- match(i, cols)
    res_i <- quiet_call(iarm::partgam_DIF(dat[, cols, drop = FALSE], group))
    new_gam[i] <- res_i$gamma[pos]
    new_sse[i] <- res_i$se[pos]
  }

  list(gamma = new_gam, se = new_sse)
}

#' A/B/C classification for partial gamma
#'
#' Adapted from \code{effecttree::get_ABC_classification} (MIT).
#'
#' @keywords internal
#' @noRd
pgamma_classify <- function(gam, sse, p_adj, thresholds, alpha = 0.05) {
  z_A <- gam / sse
  p_A <- ifelse(gam > 0,
                2 * (1 - stats::pnorm(z_A)),
                2 * stats::pnorm(z_A))
  shifted <- ifelse(gam >= 0,
                    gam - thresholds[1L],
                    gam + thresholds[1L])
  z_C <- shifted / sse
  p_C <- ifelse(shifted > 0,
                2 * (1 - stats::pnorm(z_C)),
                2 * stats::pnorm(z_C))
  if (p_adj != "none") {
    p_A <- stats::p.adjust(p_A, method = p_adj)
    p_C <- stats::p.adjust(p_C, method = p_adj)
  }
  ifelse(abs(gam) < thresholds[1L] | p_A >= alpha, "A",
         ifelse(abs(gam) > thresholds[2L] & p_C < alpha, "C", "B"))
}

#' Convert per-node effect-size results to a tidy data.frame
#'
#' One row per (node x item).
#'
#' @keywords internal
#' @noRd
effectsize_to_dataframe <- function(effect_info, item_names, tree = NULL,
                                    rescaled_terminals = NULL) {
  empty_cols <- c("NodeID","Split","Variable","Direction","Item",
                  "EffectSize","SE","Class","Flagged","Rescaled",
                  "n_left","n_right")
  if (effect_info$n_splits == 0L) {
    out <- data.frame(matrix(ncol = length(empty_cols), nrow = 0L))
    names(out) <- empty_cols
    return(out)
  }

  rows <- lapply(effect_info$nodes, function(nd) {
    n_items <- length(nd$EffectSize)

    # Rescaled at this split for this item: TRUE if item was rescaled
    # in any terminal node descended from EITHER child of this inner
    # node. Computed only when we have the tree + per-terminal table.
    rescaled <- rep(FALSE, n_items)
    if (!is.null(tree) && !is.null(rescaled_terminals) &&
        nrow(rescaled_terminals) > 0L) {
      nd_obj <- find_node_by_id(partykit::node_party(tree), nd$NodeID)
      if (!is.null(nd_obj)) {
        descend_terms <- collect_terminal_ids(nd_obj)
        bad <- rescaled_terminals$Item[
          rescaled_terminals$NodeID %in% descend_terms
        ]
        rescaled <- item_names %in% bad
      }
    }

    data.frame(
      NodeID     = rep(nd$NodeID, n_items),
      Split      = rep(nd$Split, n_items),
      Variable   = rep(nd$Variable, n_items),
      Direction  = rep(nd$Direction, n_items),
      Item       = item_names,
      EffectSize = round(nd$EffectSize, 4),
      SE         = round(nd$SE, 4),
      Class      = nd$Class,
      Flagged    = nd$Class %in% c("B", "C"),
      Rescaled   = rescaled,
      n_left     = rep(nd$n_left, n_items),
      n_right    = rep(nd$n_right, n_items),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' Render the per-split effect-size dataframe as a kable
#'
#' Flagged rows (Class B/C) render in bold. One section per inner node,
#' grouped via \code{kableExtra::pack_rows()} when the package is
#' available, and via plain section headers otherwise.
#'
#' @keywords internal
#' @noRd
render_effectsize_kable <- function(df, n_persons, n_items, es_method,
                                    purification, p_adj, thresholds,
                                    alpha, is_polytomous, has_stability,
                                    n_rescaled_items = 0L) {

  caption <- build_es_caption(
    df               = df,
    n_persons        = n_persons,
    n_items          = n_items,
    es_method        = es_method,
    purification     = purification,
    p_adj            = p_adj,
    thresholds       = thresholds,
    alpha            = alpha,
    is_polytomous    = is_polytomous,
    has_stability    = has_stability,
    n_rescaled_items = n_rescaled_items
  )

  if (nrow(df) == 0L) {
    empty <- data.frame(Note = "No splits -- model is invariant across covariates.",
                        stringsAsFactors = FALSE)
    return(knitr::kable(empty, format = "pipe", row.names = FALSE,
                        caption = caption))
  }

  display <- df
  flag    <- isTRUE_logical(display$Flagged)
  display$Flagged <- ifelse(flag, "yes", "no")
  body_cols <- c("Item", "EffectSize", "SE", "Class", "Flagged")
  if (any(flag)) {
    for (cc in body_cols) {
      col <- display[[cc]]
      if (is.numeric(col)) col <- format(col, nsmall = 3, trim = TRUE)
      display[[cc]] <- ifelse(flag, paste0("**", col, "**"),
                              as.character(col))
    }
  }

  # Build one section per node with a header row showing the split.
  # Pull n_left / n_right / Split from the original `df` to avoid the
  # bold-decoration applied to `display`.
  ids <- unique(df$NodeID)
  parts <- lapply(ids, function(nid) {
    raw <- df[df$NodeID == nid, , drop = FALSE]
    split_lbl <- raw$Split[1L]
    n_l <- raw$n_left[1L]; n_r <- raw$n_right[1L]
    header <- paste0("Node ", nid, " -- ", split_lbl,
                     "  (n left = ", n_l, ", n right = ", n_r, ")")
    body <- display[display$NodeID == nid, body_cols, drop = FALSE]
    list(header = header, body = body)
  })

  # Concatenate as a single kable by interleaving header lines using
  # pipe-format manually.
  header_line <- "| Item | EffectSize | SE | Class | Flagged |"
  divider     <- "|:-----|-----------:|---:|:------|:--------|"

  body_chunks <- lapply(parts, function(p) {
    bl <- knitr::kable(p$body, format = "pipe", row.names = FALSE)
    # Drop the kable's own header line (first 2 lines) and prepend our
    # section header instead.
    bl_lines <- as.character(bl)
    if (length(bl_lines) >= 3L) {
      data_lines <- bl_lines[3:length(bl_lines)]
    } else {
      data_lines <- bl_lines
    }
    c(paste0("**", p$header, "**"), "", header_line, divider, data_lines, "")
  })

  out_lines <- c(unlist(body_chunks, use.names = FALSE),
                 "", paste0(": ", caption))

  structure(out_lines,
            format = "pipe",
            class  = "knitr_kable")
}

#' Build the caption string for the effect-size kable
#'
#' @keywords internal
#' @noRd
build_es_caption <- function(df, n_persons, n_items, es_method,
                             purification, p_adj, thresholds, alpha,
                             is_polytomous, has_stability,
                             n_rescaled_items = 0L) {
  model_name <- if (is_polytomous) "Partial Credit Tree" else "Rasch Tree"
  es_name    <- if (es_method == "MH")
                  "Mantel-Haenszel effect size (ETS Delta scale)"
                else
                  paste0("partial gamma (B/C thresholds = ",
                         thresholds[1L], " / ", thresholds[2L], ")")
  pur_part   <- if (purification == "iterative")
                  " Iterative purification applied (Bjorner et al. 1998)."
                else ""
  padj_part  <- if (es_method == "pgamma" && p_adj != "none")
                  paste0(" p-values adjusted by ", p_adj, ".")
                else ""
  flag_part  <- if (nrow(df) > 0L)
                  paste0(" Flagged (Class B or C): ",
                         sum(df$Flagged), " / ", nrow(df),
                         " (item x split combinations).",
                         " **Bold** rows are flagged.")
                else ""
  stab_part  <- if (has_stability)
                  " Stability summary attached as `attr(result, \"stability\")`."
                else ""
  rescale_part <- if (n_rescaled_items > 0L)
                    paste0(" ", n_rescaled_items,
                           " item(s) were rescaled in at least one terminal ",
                           "node -- see `attr(result, \"rescaled_terminal_items\")`. ",
                           "Effect sizes are unaffected.")
                  else ""
  paste0(
    model_name, " (n = ", n_persons, " persons, ", n_items, " items). ",
    "Effect size: ", es_name, ". alpha = ", alpha, ".",
    pur_part, padj_part, flag_part, stab_part, rescale_part
  )
}

#' Robust isTRUE for character/logical Flagged columns
#'
#' @keywords internal
#' @noRd
isTRUE_logical <- function(x) {
  if (is.logical(x)) return(!is.na(x) & x)
  vapply(x, isTRUE, logical(1L))
}

#' Summarise a stabletree object as a per-variable data.frame
#'
#' Selection frequency = proportion of resamples in which the variable
#' was selected at any split. For continuous variables, also report
#' mean and SD of cutpoint when used.
#'
#' @keywords internal
#' @noRd
stabletree_summary_df <- function(stb, covariates) {
  vsel <- stb$vs
  if (is.null(vsel)) {
    out <- data.frame(
      Variable           = character(),
      Type               = character(),
      Selected_orig      = logical(),
      Selection_freq_pct = numeric(),
      Cutpoint_mean      = numeric(),
      Cutpoint_sd        = numeric(),
      stringsAsFactors = FALSE
    )
    attr(out, "B")       <- if (!is.null(stb$B)) stb$B else NA_integer_
    attr(out, "sampler") <- "unknown"
    attr(out, "caption") <- "No variables available for stability assessment."
    return(out)
  }

  vars <- colnames(vsel)
  freq <- colMeans(vsel > 0)
  vs0  <- if (!is.null(stb$vs0)) stb$vs0[vars] > 0 else rep(NA, length(vars))

  # Per-variable cutpoint summary across resamples (numeric covariates only;
  # for factors the "cutpoint" is a level-index pair that does not have a
  # meaningful mean/SD interpretation, so we report NA).
  cp_mean <- rep(NA_real_, length(vars))
  cp_sd   <- rep(NA_real_, length(vars))
  br_list <- stb$br
  if (!is.null(br_list)) {
    for (i in seq_along(vars)) {
      v <- vars[i]
      if (!is.numeric(covariates[[v]]) || is.factor(covariates[[v]])) next
      brks <- br_list[[v]]
      if (is.null(brks) || length(brks) == 0L) next
      if (is.list(brks)) {
        all_breaks <- unlist(lapply(brks, function(rec) {
          if (is.list(rec) && !is.null(rec$breaks) && is.numeric(rec$breaks)) {
            as.numeric(rec$breaks)
          } else if (is.numeric(rec)) {
            as.numeric(rec)
          } else {
            numeric(0)
          }
        }), use.names = FALSE)
      } else if (is.numeric(brks)) {
        all_breaks <- as.numeric(brks)
      } else {
        all_breaks <- numeric(0)
      }
      all_breaks <- all_breaks[is.finite(all_breaks)]
      if (length(all_breaks) > 0L) cp_mean[i] <- mean(all_breaks)
      if (length(all_breaks) >= 2L) cp_sd[i]   <- stats::sd(all_breaks)
    }
  }

  type_of <- function(v) {
    col <- covariates[[v]]
    if (is.ordered(col)) "ordered"
    else if (is.factor(col)) "factor"
    else if (is.logical(col)) "logical"
    else if (is.numeric(col)) "numeric"
    else class(col)[1L]
  }

  out <- data.frame(
    Variable           = vars,
    Type               = vapply(vars, type_of, character(1L)),
    Selected_orig      = unname(vs0),
    Selection_freq_pct = round(unname(freq) * 100, 1),
    Cutpoint_mean      = round(cp_mean, 4),
    Cutpoint_sd        = round(cp_sd, 4),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  # Sort: selected-in-original first, then by frequency descending
  out <- out[order(!out$Selected_orig, -out$Selection_freq_pct), , drop = FALSE]
  rownames(out) <- NULL

  attr(out, "B")       <- nrow(vsel)
  attr(out, "sampler") <- if (!is.null(stb$sampler))
                            attr(stb$sampler, "method") %||% "subsampling"
                          else "subsampling"
  attr(out, "caption") <- paste0(
    "Tree-stability assessment across ", nrow(vsel),
    " resamples (", attr(out, "sampler"), "). ",
    "Selection_freq_pct: percent of resamples in which the covariate was ",
    "selected at any split. Selected_orig: TRUE if selected in the original tree. ",
    "Cutpoint_mean / Cutpoint_sd: mean and SD of the chosen cutpoint across ",
    "resamples (numeric covariates only)."
  )

  out
}

#' Render the stability summary as a kable
#'
#' @keywords internal
#' @noRd
stabletree_summary_kable <- function(stab_df) {
  if (nrow(stab_df) == 0L) return(NULL)
  cap <- attr(stab_df, "caption")
  knitr::kable(stab_df, format = "pipe", row.names = FALSE,
               caption = cap)
}

# null-coalesce helper
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Suppress messages and printed output from a single call
#'
#' Equivalent to SimDesign::quiet() but without the dependency.
#'
#' @keywords internal
#' @noRd
quiet_call <- function(expr) {
  utils::capture.output({
    val <- suppressMessages(suppressWarnings(force(expr)))
  }, file = nullfile())
  val
}
