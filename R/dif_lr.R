#' DIF analysis via Andersen's likelihood-ratio test
#'
#' Tests whether the item parameters differ between the groups of an
#' external variable, with Andersen's likelihood-ratio (LR) test. The model
#' is fitted by conditional maximum likelihood (CML) in the full sample and
#' in each group. The test answers whether there is any DIF, not which items
#' or thresholds have it. For that, use \code{\link{RMdifAnchor}}.
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person. Person
#'   IDs and grouping variables must not be included -- pass the grouping
#'   variable separately via \code{dif_var}.
#' @param dif_var Vector of length \code{nrow(data)} (factor, character,
#'   or numeric) defining the DIF grouping variable. Coerced to factor;
#'   unused levels are dropped. Rows where \code{dif_var} is \code{NA} are
#'   dropped with a message. Must result in at least 2 groups after
#'   cleaning.
#' @param model One of \code{"auto"} (default), \code{"PCM"}, or
#'   \code{"RM"}. \code{"auto"} reports the Rasch model when the data are
#'   dichotomous (max response = 1) and the partial credit model otherwise
#'   (the two coincide for dichotomous data). \code{"RM"} errors on
#'   polytomous data.
#' @param output \code{"kable"} (default) or \code{"dataframe"}.
#' @param ... Not used. Arguments that earlier versions accepted for the
#'   per-group location table and plots (\code{level}, \code{cutoff},
#'   \code{conf}, \code{sort}, \code{p_value}, \code{correction},
#'   \code{alpha}, \code{plot_type}, \code{reference}) give an error that
#'   points to \code{\link{RMdifAnchor}}.
#'
#' @return A one-row data.frame, or a \code{knitr_kable} of it, with
#'   \code{Chisq} (the LR statistic), \code{df}, \code{p}, \code{Model}
#'   (\code{"PCM"} or \code{"RM"}), \code{Items}, \code{Groups}, \code{N}
#'   (respondents used), \code{N_excluded} (rows dropped for a missing
#'   \code{dif_var} or no responses), \code{N_incomplete} (respondents used
#'   who skipped some items) and one \code{n_<group>} column per group.
#'
#' @details
#' The LR statistic is twice the difference between the sum of the group
#' conditional log-likelihoods and the full-sample one, with df equal to the
#' extra free parameters of the group fits (Andersen, 1973). It equals
#' \code{eRm::LRtest()} on the same data. Each fit is a CML partial credit
#' model from \code{psychotools::pcmodel()} (a dichotomous item is a
#' 2-category PCM). Respondents with missing responses are kept (CML uses
#' the items they answered).
#'
#' The test does not depend on how each group's scale is centred, but any
#' comparison of single items or thresholds between groups does. Centring
#' each group on its own mean makes real DIF in some items appear as
#' artificial DIF in the others (Andrich & Hagquist, 2012), so this function
#' no longer reports per-group locations. \code{\link{RMdifAnchor}} links the
#' groups through anchor items and tests each item, threshold or threshold
#' spacing, with a correction for multiple testing.
#'
#' Like other global tests, the LR test becomes significant for small
#' amounts of DIF in large samples, and has little power in small ones.
#'
#' \strong{Response categories.} Every category an item uses must be used in
#' every group, or the group fits have different thresholds and the test is
#' not defined. The same holds when a category is used in a group only by
#' respondents with a zero or perfect score there, since CML leaves them
#' out. The function then stops and names the items, categories and groups
#' (\code{eRm::LRtest()} instead left items with an unused category out of
#' the test without saying so). With small groups this is common for sparse
#' categories: merge them with an adjacent one, combine groups, or leave the
#' item out.
#'
#' @references
#' Andersen, E. B. (1973). A goodness of fit test for the Rasch model.
#' *Psychometrika, 38*(1), 123-140. \doi{10.1007/BF02291180}
#'
#' Andrich, D., & Hagquist, C. (2012). Real and artificial differential item
#' functioning. *Journal of Educational and Behavioral Statistics, 37*(3),
#' 387-416. \doi{10.3102/1076998611411913}
#'
#' @seealso \code{\link{RMdifAnchor}} for anchored tests of single items,
#'   thresholds and threshold spacing, \code{\link{RMdifGamma}} and
#'   \code{\link{RMdifTree}}.
#'
#' @examples
#' \donttest{
#' if (requireNamespace("eRm", quietly = TRUE)) { # for the example data
#'   set.seed(1)
#'   data("pcmdat2", package = "eRm")
#'   grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
#'
#'   RMdifLR(pcmdat2, dif_var = grp)
#'   RMdifLR(pcmdat2, dif_var = grp, output = "dataframe")
#' }
#' }
#'
#' @export
RMdifLR <- function(
  data,
  dif_var,
  model = c("auto", "PCM", "RM"),
  output = c("kable", "dataframe"),
  ...
) {
  # In 1.4.0 the fourth argument was `level`, so a positional call written
  # for it puts "item" or "threshold" into `output`.
  if (is.character(output) && length(output) == 1L &&
        output %in% c("item", "threshold")) {
    stop(
      "RMdifLR() now reports Andersen's LR test only. Its fourth argument ",
      "is `output` (\"kable\" or \"dataframe\"), and \"", output,
      "\" is a value of the removed `level`. Use RMdifAnchor() for anchored ",
      "tests of items, thresholds and threshold spacing.",
      call. = FALSE
    )
  }
  dots <- list(...)
  if (length(dots) > 0L) {
    removed <- c("level", "cutoff", "conf", "sort", "p_value", "correction",
                 "alpha", "plot_type", "reference")
    nm <- names(dots)
    if (is.null(nm)) nm <- rep("", length(dots))
    if (any(nm %in% removed)) {
      stop(
        "RMdifLR() now reports Andersen's LR test only. Argument(s) ",
        paste0("`", intersect(nm, removed), "`", collapse = ", "),
        " were removed with the per-group location table and plots. Use ",
        "RMdifAnchor() for anchored tests of items, thresholds and ",
        "threshold spacing.",
        call. = FALSE
      )
    }
    stop("Unused argument(s) in RMdifLR(): ",
         paste(ifelse(nzchar(nm), nm, "<unnamed>"), collapse = ", "), ".",
         call. = FALSE)
  }
  if (identical(output, "ggplot")) {
    stop(
      "RMdifLR() no longer draws plots. Use RMdifAnchor(..., output = ",
      "\"ggplot\") for anchored differences between two groups.",
      call. = FALSE
    )
  }
  model <- match.arg(model)
  output <- match.arg(output)
  if (output == "kable" && !requireNamespace("knitr", quietly = TRUE)) {
    stop("Package 'knitr' is required for output = \"kable\".", call. = FALSE)
  }
  rlang::inform(
    paste0(
      "RMdifLR() now reports Andersen's LR test only, without the ",
      "per-group location table and plots. For DIF in single items, ",
      "thresholds and threshold spacing, use RMdifAnchor()."
    ),
    .frequency = "once",
    .frequency_id = "easyRasch2_difLR_anchor"
  )

  validate_response_data(data)
  prep <- .prepare_dif_data(data, dif_var)
  data <- prep$data
  dif_factor <- prep$dif_factor
  n_total <- prep$n_total
  has_na <- prep$has_na
  if (nlevels(dif_factor) < 2L) {
    stop(
      "`dif_var` must have at least 2 distinct non-missing levels.",
      call. = FALSE
    )
  }
  groups <- levels(dif_factor)

  data_mat <- as.matrix(data)
  storage.mode(data_mat) <- "integer"
  is_polytomous <- max(data_mat, na.rm = TRUE) > 1L
  resolved_model <- if (model == "auto") {
    if (is_polytomous) "PCM" else "RM"
  } else {
    if (model == "RM" && is_polytomous) {
      stop(
        "`model = \"RM\"` requires dichotomous (0/1) data; use ",
        "\"PCM\" or \"auto\" for polytomous responses.",
        call. = FALSE
      )
    }
    model
  }

  # Every category an item uses in the full data must occur in every group,
  # or the group fits would have different thresholds.
  .check_dif_categories(data, dif_factor, "RMdifLR")
  fit_full <- .cml_fit_or_stop(data, "the full sample")
  fits <- lapply(groups, function(g) {
    .cml_fit_or_stop(data[dif_factor == g, , drop = FALSE],
                     paste0("group '", g, "'"))
  })

  # LR = 2 * (sum of group conditional log-likelihoods - full one), with df
  # the extra free parameters of the group fits (Andersen, 1973).
  ll <- function(f) as.numeric(stats::logLik(f))
  npar <- function(f) length(stats::coef(f))
  lr_stat <- 2 * (sum(vapply(fits, ll, numeric(1L))) - ll(fit_full))
  lr_df <- sum(vapply(fits, npar, integer(1L))) - npar(fit_full)
  lr_p <- stats::pchisq(lr_stat, lr_df, lower.tail = FALSE)
  n_per_group <- stats::setNames(as.integer(table(dif_factor)), groups)

  # Everything is in the columns (no attributes): the test, the model, and
  # the respondents behind it, overall and per group.
  out <- data.frame(
    Chisq = lr_stat,
    df = lr_df,
    p = lr_p,
    Model = resolved_model,
    Items = ncol(data),
    Groups = length(groups),
    N = nrow(data),
    N_excluded = n_total - nrow(data),
    N_incomplete = sum(!stats::complete.cases(data)),
    stringsAsFactors = FALSE
  )
  for (g in groups) out[[paste0("n_", g)]] <- n_per_group[[g]]
  if (output == "dataframe") {
    return(out)
  }

  count_cols <- setdiff(names(out), c("Chisq", "p", "Model"))
  display <- out
  display$p <- .format_p_table(out$p)
  display <- .round_display(display, c(
    list(Chisq = 3L),
    stats::setNames(as.list(rep(0L, length(count_cols))), count_cols)
  ))
  n_clause <- .n_caption(
    nrow(data), n_total,
    c(if (has_na) "incomplete responses retained",
      paste(paste0(groups, ": ", n_per_group), collapse = ", "))
  )
  # N_excluded: rows dropped for a missing DIF group or no responses.
  caption <- paste0(
    "Andersen's likelihood-ratio test of DIF, ",
    if (resolved_model == "PCM") "partial credit" else "Rasch",
    " model fitted by CML in each of ", length(groups), " groups. ",
    n_clause, ", ", ncol(data), " items. A significant test means that ",
    "some item parameters differ between groups. RMdifAnchor() tests ",
    "single items and thresholds."
  )
  knitr::kable(display, format = "pipe", row.names = FALSE, caption = caption)
}
