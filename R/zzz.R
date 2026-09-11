#' Package startup message
#'
#' Prints the installed version, points at the changelog, and gives the
#' canonical CRAN page so the reader can see which release is current.
#'
#' No network access: the URL is a literal string and nothing is fetched at
#' attach time. Checking CRAN here would consume a shared external resource on
#' every `library()` call, would block on machines without outbound internet
#' (HPC nodes, CI runners, the CRAN check farm), and would be wrong for anyone
#' on a frozen or institutional mirror.
#'
#' `packageStartupMessage()` rather than `message()`, so that
#' `suppressPackageStartupMessages()` silences it.
#'
#' @param libname Library path, supplied by R.
#' @param pkgname Package name, supplied by R.
#' @return Invisibly `NULL`, called for its side effect.
#' @keywords internal
#' @noRd
.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "This is easyRasch2 ", utils::packageVersion("easyRasch2"), ".\n",
    "For recent changes type news(package = 'easyRasch2').\n",
    "Current release: https://CRAN.R-project.org/package=easyRasch2"
  )
  invisible(NULL)
}
