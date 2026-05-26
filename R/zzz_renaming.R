#' Function renaming in easyRasch2 0.8.0
#'
#' In version 0.8.0, 22 exported functions were renamed under a
#' consistent *domain-prefix \eqn{\to} method \eqn{\to} variant-suffix*
#' scheme. The renames have no semantic effect --- only names changed.
#' Use the table below to migrate existing scripts.
#'
#' @section Domain prefixes:
#' \itemize{
#'   \item `RMdif` --- differential item functioning
#'   \item `RMlocdep` --- local dependence
#'   \item `RMitem` --- per-item statistics (fit, restscore, ICC, hierarchy)
#'   \item `RMdim` --- dimensionality assessment
#'   \item `RMplot` --- descriptive response-distribution plots
#' }
#'
#' Variant suffixes attached at the end:
#' \itemize{
#'   \item `Cutoff` --- simulation-based critical values
#'   \item `Plot` --- visualisation
#'   \item `MI` --- multiple-imputation variant (replaces former `_mi`)
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
#' @section No deprecation aliases:
#' No deprecation shims are shipped: each old name is simply gone. This
#' keeps the export list clean ahead of the first CRAN submission. Run a
#' search-and-replace against the table above to migrate.
#'
#' @name easyRasch2-renaming
#' @aliases easyRasch2-renaming
#' @keywords internal
NULL
