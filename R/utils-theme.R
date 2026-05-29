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
#' renders left-aligned with the APA-conventional italic "Note." prefix
#' followed by roman body text. When `ggtext` is installed
#' (recommended; in Suggests), the caption is rendered via
#' `ggtext::element_markdown()` which lets the italic "*Note.*" prefix
#' coexist with roman body text in a single caption. When `ggtext` is
#' not installed, falls back to a plain `element_text(face = "italic")`
#' --- the whole caption goes italic, and the prefix collapses to a
#' plain "Note. " (no asterisks leaked).
#'
#' Pair with \code{\link{er2_caption}} when building caption strings,
#' so the "Note." prefix matches whichever renderer is active.
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

#' APA-style "Note." caption-text prefix with line wrapping
#'
#' Build a plot caption with the standard "Note." prefix and wrap the
#' body text at `width` characters so long captions don't run off the
#' right edge of the plot. The body text is wrapped first (via
#' [strwrap()]) then prefixed --- so the prefix never gets broken by a
#' line break, and the markdown asterisks (when `ggtext` is installed)
#' stay together.
#'
#' When `ggtext` is installed, returns `"*Note.* <wrapped body>"`
#' (markdown italic on "Note."); otherwise `"Note. <wrapped body>"`
#' (plain text, with the whole caption set italic by
#' \code{\link{er2_plot_caption}}).
#'
#' @param text Character. The caption body text (everything after the
#'   "Note." prefix). Pre-existing `\\n` in `text` is treated as
#'   ordinary whitespace by `strwrap()`; the wrapped output reflows
#'   to fit `width`.
#' @param width Integer. Maximum characters per line, passed to
#'   [strwrap()] after subtracting the prefix length. Default `90L`
#'   (suits the package's default 8-inch plot at 144 dpi).
#' @return A character string ready to pass to
#'   `ggplot2::labs(caption = ...)`.
#' @keywords internal
#' @noRd
er2_caption <- function(text, width = 90L) {
  prefix <- if (requireNamespace("ggtext", quietly = TRUE)) {
    "*Note.* "
  } else {
    "Note. "
  }
  body_width <- max(20L, width - nchar(prefix))
  body <- paste(strwrap(text, width = body_width), collapse = "\n")
  paste0(prefix, body)
}
