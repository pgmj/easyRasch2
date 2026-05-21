#' Conditional Item Characteristic Curves
#'
#' Faceted panel of Conditional Item Characteristic Curves (cICCs) for all
#' items in `data`. Each panel shows the model-expected ICC together with
#' average observed item scores within class intervals (or for every
#' possible total score, with `method = "score"`). When a `dif_var` is
#' supplied, observed averages are computed separately per group, so each
#' panel becomes a visual DIF check.
#'
#' Wraps [iarm::ICCplot()] (Santiago, Mueller, & Müller, 2022) for the
#' per-item computation and composes the panels with [patchwork::wrap_plots()].
#'
#' @param data A data.frame or matrix of item responses (non-negative
#'   integers, 0-based). One column per item, one row per person.
#' @param class_intervals Integer. Number of class intervals used to
#'   group respondents when `method = "cut"`. Default `4`. Each class
#'   interval needs enough respondents to compute a stable average
#'   observed item score; with small samples or many possible total
#'   scores, lower this. Ignored when `method = "score"`.
#' @param method One of `"cut"` (default) or `"score"`. `"cut"` groups
#'   respondents into `class_intervals` adjacent bins; `"score"` uses
#'   every possible total raw score.
#' @param dif_var Optional vector of length `nrow(data)` (factor or
#'   coercible to factor) defining a DIF grouping variable. When
#'   supplied, observed averages are shown separately per group.
#'   Default `NULL` (no DIF).
#' @param output One of `"patchwork"` (default) -- composite patchwork
#'   figure -- or `"list"` -- a named list of per-item `ggplot` objects.
#'
#' @return Either a `ggplot` / `patchwork` composite (default), or a
#'   list of `ggplot` objects (`output = "list"`), one per item.
#'
#' @details
#' \strong{Class intervals vs every total score.} With `method = "cut"`
#' (the default) respondents are divided into `class_intervals` bins of
#' roughly equal size based on their total raw score; each bin's average
#' observed item score is plotted as a point. With `method = "score"`,
#' the average observed item score is computed for every observed total
#' raw score. The cut method is the conventional Rasch-analysis default
#' and works well at moderate sample sizes; the score method is more
#' granular but can be unstable when raw-score totals are sparsely
#' populated.
#'
#' \strong{DIF mode.} With `dif_var` supplied, [iarm::ICCplot()] computes
#' average observed item scores separately for each level of the
#' grouping variable. Each item's panel then shows the model-expected
#' ICC (one black line) plus one set of observed points per DIF group.
#' Stark visual divergence between groups suggests DIF.
#'
#' \strong{Dependencies.} Requires `iarm` for the per-item curve, and
#' `patchwork` + `ggplot2` for the figure composition. All three are in
#' Suggests with runtime guards.
#'
#' @references
#' Mueller, M., & Santiago, P. H. R. (2022). \emph{iarm: Item Analysis
#' in Rasch Models}. R package version 0.4.3.
#' \url{https://CRAN.R-project.org/package=iarm}
#'
#' @seealso [iarm::ICCplot()], [RMtargeting()], [RMtileplot()]
#'
#' @examples
#' \donttest{
#' data("pcmdat2", package = "eRm")
#' RMciccPlot(pcmdat2, class_intervals = 4)
#'
#' # With DIF grouping
#' set.seed(1)
#' grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
#' RMciccPlot(pcmdat2, dif_var = grp)
#' }
#'
#' @export
RMciccPlot <- function(data,
                       class_intervals = 4L,
                       method          = c("cut", "score"),
                       dif_var         = NULL,
                       output          = c("patchwork", "list")) {

  method <- match.arg(method)
  output <- match.arg(output)

  # --- Package availability -------------------------------------------------
  if (!requireNamespace("iarm", quietly = TRUE)) {
    stop("Package 'iarm' is required for RMciccPlot(). ",
         "Install with: install.packages(\"iarm\")",
         call. = FALSE)
  }
  if (output == "patchwork" &&
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for output = \"patchwork\". ",
         "Install with: install.packages(\"patchwork\")",
         call. = FALSE)
  }
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for RMciccPlot(). ",
         "Install with: install.packages(\"ggplot2\")",
         call. = FALSE)
  }

  # --- Data validation ------------------------------------------------------
  validate_response_data(data)
  data <- as.data.frame(data)
  if (ncol(data) < 2L) {
    stop("RMciccPlot() requires at least 2 items.", call. = FALSE)
  }

  # --- class_intervals validation ------------------------------------------
  if (!is.numeric(class_intervals) || length(class_intervals) != 1L ||
      !is.finite(class_intervals) || class_intervals < 2L) {
    stop("`class_intervals` must be a single integer >= 2.", call. = FALSE)
  }
  class_intervals <- as.integer(class_intervals)

  # --- DIF variable + joint NA handling ------------------------------------
  dif_mode   <- !is.null(dif_var)
  dif_factor <- NULL
  if (dif_mode) {
    if (length(dif_var) != nrow(data)) {
      stop("`dif_var` must have the same length as `nrow(data)` (",
           nrow(data), "). Got ", length(dif_var), ".", call. = FALSE)
    }
    dif_factor <- droplevels(as.factor(dif_var))
  }

  # Drop incomplete rows (iarm::ICCplot() calls complete.cases() internally)
  # jointly with dif_factor.
  keep <- stats::complete.cases(data)
  if (dif_mode) keep <- keep & !is.na(dif_factor)
  data <- data[keep, , drop = FALSE]
  if (dif_mode) dif_factor <- droplevels(dif_factor[keep])
  if (nrow(data) == 0L) {
    stop("No complete cases in `data`.", call. = FALSE)
  }
  if (dif_mode && nlevels(dif_factor) < 2L) {
    stop("`dif_var` must have at least 2 distinct non-missing levels ",
         "after NA filtering.", call. = FALSE)
  }

  # --- Build per-item plots -------------------------------------------------
  # iarm::ICCplot() (a) emits ggplot2 `size=` deprecation warnings and (b)
  # prints DIF stats to stdout when difstats = "yes". Wrap with
  # suppressWarnings + capture.output to keep the user-facing output clean.
  iccplots <- vector("list", ncol(data))
  for (i in seq_len(ncol(data))) {
    p <- tryCatch({
      tmp_cap <- character()
      tmp_cap <- utils::capture.output(
        plt <- suppressWarnings(suppressMessages(
          if (dif_mode) {
            iarm::ICCplot(
              data        = data,
              itemnumber  = i,
              method      = method,
              cinumber    = class_intervals,
              axis.rumm   = "yes",
              title       = NULL,
              icclabel    = "yes",
              itemdescrip = names(data)[i],
              dif         = "yes",
              difvar      = dif_factor,
              diflabels   = levels(dif_factor),
              difstats    = "yes"
            )
          } else {
            iarm::ICCplot(
              data        = data,
              itemnumber  = i,
              method      = method,
              cinumber    = class_intervals,
              axis.rumm   = "yes",
              title       = NULL,
              icclabel    = "yes",
              itemdescrip = names(data)[i],
              dif         = "no"
            )
          }
        ))
      )
      plt + ggplot2::theme(legend.direction = "vertical")
    }, error = function(e) {
      stop(paste0(
        "iarm::ICCplot() failed on item ", names(data)[i],
        ": ", conditionMessage(e),
        ". This often means `class_intervals` (", class_intervals,
        ") is too large for the available total scores. ",
        "Try a smaller value or use `method = \"score\"`."),
        call. = FALSE)
    })
    iccplots[[i]] <- p
  }
  names(iccplots) <- names(data)

  if (output == "list") return(iccplots)

  # --- Compose with patchwork ----------------------------------------------
  patchwork::wrap_plots(
    iccplots,
    axes        = "collect",
    guides      = "collect",
    axis_titles = "collect"
  ) +
    patchwork::plot_annotation(
      title = "Conditional Item Characteristic Curves"
    )
}
