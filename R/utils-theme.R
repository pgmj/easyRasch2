# Internal theme helpers shared across easyRasch2 plotting functions.

#' Standard easyRasch2 axis-title margins
#'
#' Returns the ggplot2 theme tweak used throughout easyRasch2 to add a
#' little extra breathing room around the x and y axis titles. Add to
#' any ggplot via `p + er2_axis_margins()`.
#'
#' Implemented as a `ggplot2::theme()` object that sets only the
#' `margin` slot of `axis.title.x` and `axis.title.y` --- other axis-
#' title attributes (size, face, family, colour) inherit from the
#' parent theme and are not disturbed.
#'
#' @return A `ggplot2::theme()` object.
#' @keywords internal
#' @noRd
er2_axis_margins <- function() {
  ggplot2::theme(
    axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 12)),
    axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 12))
  )
}
