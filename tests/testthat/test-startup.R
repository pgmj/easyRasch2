# The attach-time message. Two things matter beyond its content: that it is a
# startup message rather than a plain one, so `suppressPackageStartupMessages()`
# silences it, and that it reaches for nothing over the network.

test_that(".onAttach reports the installed version and where to read more", {
  msg <- testthat::capture_messages(
    easyRasch2:::.onAttach(libname = NULL, pkgname = "easyRasch2")
  )
  txt <- paste(msg, collapse = "")
  expect_match(txt, "This is easyRasch2", fixed = TRUE)
  # the version is read at run time, never a literal that could go stale
  expect_match(txt, as.character(utils::packageVersion("easyRasch2")),
               fixed = TRUE)
  expect_match(txt, "news(package = 'easyRasch2')", fixed = TRUE)
  # the canonical CRAN form, not the /web/packages/ one
  expect_match(txt, "https://CRAN.R-project.org/package=easyRasch2", fixed = TRUE)
})

test_that("the startup message is suppressible", {
  # packageStartupMessage(), not message(): a plain message would ignore
  # suppressPackageStartupMessages() and there would be no way to quieten
  # library(easyRasch2) in a script.
  expect_silent(
    suppressPackageStartupMessages(
      easyRasch2:::.onAttach(libname = NULL, pkgname = "easyRasch2")
    )
  )
})

test_that("attaching the package touches no network resource", {
  # A version check at attach time was considered and rejected. This pins that
  # decision: the message must stay a literal string.
  #
  # Inspect the installed function rather than R/zzz.R. R CMD check runs the
  # tests against an installed copy of the package, where the source directory
  # is absent and reading it fails with "cannot open the connection". Deparsing
  # without "useSource" also drops comments, so no comment can trip the check.
  code <- deparse(body(easyRasch2:::.onAttach),
                  control = c("keepInteger", "keepNA"))
  expect_false(any(grepl(
    "available\\.packages|old\\.packages|download\\.file|url\\(|curl|httr|readLines\\(\\s*[\"']http",
    code
  )))
})
