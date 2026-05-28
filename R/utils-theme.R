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

#' APA-style plot.caption theme element
#'
#' Returns a `ggplot2::theme()` object that sets `plot.caption` so it
#' renders left-aligned ("Note." prefix in italic, body text in roman
#' face). When `ggtext` is installed (recommended; in Suggests), the
#' caption is rendered via `ggtext::element_markdown()` which lets the
#' italic "*Note.*" prefix coexist with roman body text in a single
#' caption. When `ggtext` is not installed, falls back to a plain
#' `element_text(face = "italic", ...)` --- the whole caption goes
#' italic, but no markdown asterisks leak through to the rendered text.
#'
#' Pair with \code{\link{er2_caption}} when building caption strings, so the
#' "*Note.* " prefix is formatted to match whichever renderer is active.
#'
#' @return A `ggplot2::theme()` object.
#' @keywords internal
#' @noRd
er2_plot_caption <- function() {
  ggplot2::theme(
    plot.caption = if (requireNamespace("ggtext", quietly = TRUE)) {
      ggtext::element_markdown(hjust = 0, size = 9)
    } else {
      ggplot2::element_text(hjust = 0, face = "italic", size = 9)
    }
  )
}

#' APA-style "Note." caption-text prefix
#'
#' Wrap a plot caption with the standard APA prefix. When `ggtext` is
#' installed, returns `"*Note.* <text>"` (markdown italic on "Note.");
#' otherwise returns `"Note. <text>"` (plain text, with the whole
#' caption set italic by \code{\link{er2_plot_caption}}). Centralises the prefix
#' so a single edit here propagates to every plot in the package.
#'
#' @param text Character. The caption body text (everything after
#'   "Note. ").
#' @return A character string ready to pass to
#'   `ggplot2::labs(caption = ...)`.
#' @keywords internal
#' @noRd
er2_caption <- function(text) {
  if (requireNamespace("ggtext", quietly = TRUE)) {
    paste0("*Note.* ", text)
  } else {
    paste0("Note. ", text)
  }
}
