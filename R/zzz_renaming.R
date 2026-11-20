#' Function renaming in easyRasch2 0.8.0
#'
#' In version 0.8.0, 25 exported functions were renamed under one scheme:
#' domain prefix, then method, then variant suffix. Only the names
#' changed. Use the table below to update existing scripts.
#'
#' @section Domain prefixes:
#' \itemize{
#'   \item `RMdif`: differential item functioning
#'   \item `RMlocdep`: local dependence
#'   \item `RMitem`: item statistics (fit, restscore, ICC, hierarchy)
#'   \item `RMdim`: dimensionality
#'   \item `RMplot`: descriptive plots of the response distributions
#' }
#'
#' Variant suffixes attached at the end:
#' \itemize{
#'   \item `Cutoff`: simulation-based critical values
#'   \item `Plot`: plots
#'   \item `MI`: version for multiply imputed data (replaces `_mi`)
#' }
#'
#' @section Old to new name map:
#' \tabular{ll}{
#'   \strong{Old name} \tab \strong{New name} \cr
#'   `RMpartgamDIF()`        \tab `RMdifGamma()`             \cr
#'   `RMpgDIFcutoff()`       \tab `RMdifGammaCutoff()`       \cr
#'   `RMpgDIFplot()`         \tab `RMdifGammaPlot()`         \cr
#'   `RMpartgamLD()`         \tab `RMlocdepGamma()`          \cr
#'   `RMpgLDcutoff()`        \tab `RMlocdepGammaCutoff()`    \cr
#'   `RMpgLDplot()`          \tab `RMlocdepGammaPlot()`      \cr
#'   `RMlocdepQ3cutoff()`    \tab `RMlocdepQ3Cutoff()`       \cr
#'   `RMlocdepQ3plot()`      \tab `RMlocdepQ3Plot()`         \cr
#'   `RMiteminfit()`         \tab `RMitemInfit()`            \cr
#'   `RMiteminfit_mi()`      \tab `RMitemInfitMI()`          \cr
#'   `RMinfitcutoff()`       \tab `RMitemInfitCutoff()`      \cr
#'   `RMinfitcutoff_mi()`    \tab `RMitemInfitCutoffMI()`    \cr
#'   `RMinfitcutoffPlot()`   \tab `RMitemInfitCutoffPlot()`  \cr
#'   `RMitemrestscore()`     \tab `RMitemRestscore()`        \cr
#'   `RMbootRestscore()`     \tab `RMitemRestscoreBoot()`    \cr
#'   `RMciccPlot()`          \tab `RMitemICCPlot()`          \cr
#'   `RMresidualPCA()`       \tab `RMdimResidualPCA()`       \cr
#'   `RMpcaCutoff()`         \tab `RMdimResidualPCACutoff()` \cr
#'   `RMcfaCutoff()`         \tab `RMdimCFACutoff()`         \cr
#'   `RMcfaPlot()`           \tab `RMdimCFAPlot()`           \cr
#'   `RMmartinLof()`         \tab `RMdimMartinLof()`         \cr
#'   `RMmartinLofResiduals()`\tab `RMdimMartinLofResiduals()`\cr
#'   `RMtileplot()`          \tab `RMplotTile()`             \cr
#'   `RMbarplot()`           \tab `RMplotBar()`              \cr
#'   `RMstackedbarplot()`    \tab `RMplotStackedbar()`       \cr
#' }
#'
#' Unchanged: `RMreliability()`, `RMUreliability()`,
#' `RMitemHierarchy()`, `RMdifLR()`, `RMdifTree()`, `RMlocdepQ3()`,
#' `RMtargeting()`, `RMscoreSE()`.
#'
#' @section No aliases for the old names:
#' The old names were removed without aliases, to keep the export list
#' clean before the first CRAN release. Search and replace with the table
#' above to update a script. `RMitemInfitCutoffPlot()` has since been
#' replaced by [RMitemInfitPlot()] and remains as a deprecated alias.
#'
#' @name easyRasch2-renaming
#' @aliases easyRasch2-renaming
#' @keywords internal
NULL
