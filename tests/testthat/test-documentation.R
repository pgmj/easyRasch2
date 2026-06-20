# Guard against the roxygen "Could not resolve link to topic" warnings that
# devtools::document() emits when a cross-reference link points at an internal
# helper. By package convention every dot-prefixed function (e.g.
# `.estimate_thetas`, `.q3_residual_matrix`) is internal and tagged `@noRd`, so
# it generates no Rd topic and a link to it can never resolve. Reference such
# helpers as plain code -- `.estimate_thetas()` -- never as a link
# (`[.estimate_thetas()]`, `[`.estimate_thetas()`]`, or `\link{.estimate_thetas}`).
#
# This is the blind spot R CMD check does not cover: it validates \link{} targets
# only for documented (exported) topics, so broken links to @noRd helpers slip
# through as document()-time warnings.

test_that("roxygen comments never link to internal (dot-prefixed) helpers", {
  r_dir <- test_path("..", "..", "R")
  skip_if_not(dir.exists(r_dir))
  r_files <- list.files(r_dir, pattern = "\\.R$", full.names = TRUE)

  # Matches a roxygen markdown auto-link or an Rd \link{} whose target begins
  # with a dot: "[.name", "[`.name", or "\\link{.name".
  link_to_dotted <- "\\[`?\\.[A-Za-z]|\\\\link\\{\\.[A-Za-z]"

  offenders <- character()
  for (f in r_files) {
    lines    <- readLines(f, warn = FALSE)
    is_roxy  <- grepl("^\\s*#'", lines)
    hit_rows <- which(is_roxy & grepl(link_to_dotted, lines))
    if (length(hit_rows)) {
      offenders <- c(offenders, sprintf("%s:%d: %s", basename(f), hit_rows,
                                        trimws(lines[hit_rows])))
    }
  }

  if (length(offenders)) {
    fail(paste0(
      "Roxygen links to internal (@noRd) helpers warn in devtools::document() ",
      "(\"Could not resolve link to topic\"). Use plain `code()` instead of a ",
      "[link]:\n", paste(offenders, collapse = "\n")
    ))
  }
  succeed()
})
