# Pre-compile the vignette locally.
# Run from the package root directory:
#   Rscript vignettes/precompile.R
#
# This knits the heavy `.Rmd.orig` source into a lightweight `.Rmd`
# whose code chunks are already evaluated and whose figures are saved
# under vignettes/figures/. The shipped `.Rmd` is what `R CMD build`
# sees, so the vignette installs in seconds even though the analysis
# code that produced it may have taken minutes.

library(knitr)

# Re-install from local source so vignette code calls the current
# package state. Skip if you know the installed version is already
# current and you only want a quick re-render.
if (interactive() ||
    !nzchar(Sys.getenv("EASYRASCH2_SKIP_INSTALL"))) {
  devtools::install(".", upgrade = FALSE, quiet = TRUE)
}

# Knit from within vignettes/ so fig.path = "figures/rasch-" resolves
# relative to the directory the .Rmd lives in.
withr::with_dir("vignettes", {
  knit(
    input  = "easyRasch2.Rmd.orig",
    output = "easyRasch2.Rmd"
  )
})
