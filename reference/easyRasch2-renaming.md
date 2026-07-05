# Function renaming in easyRasch2 0.8.0

In version 0.8.0, 22 exported functions were renamed under a consistent
*domain-prefix \\\to\\ method \\\to\\ variant-suffix* scheme. The
renames have no semantic effect — only names changed. Use the table
below to migrate existing scripts.

## Domain prefixes

- `RMdif` — differential item functioning

- `RMlocdep` — local dependence

- `RMitem` — per-item statistics (fit, restscore, ICC, hierarchy)

- `RMdim` — dimensionality assessment

- `RMplot` — descriptive response-distribution plots

Variant suffixes attached at the end:

- `Cutoff` — simulation-based critical values

- `Plot` — visualisation

- `MI` — multiple-imputation variant (replaces former `_mi`)

## Old to new name map

|  |  |
|----|----|
| **Old name** | **New name** |
| `RMpartgamDIF()` | [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/reference/RMdifGamma.md) |
| `RMpgDIFcutoff()` | [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdifGammaCutoff.md) |
| `RMpgDIFplot()` | [`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/reference/RMdifGammaPlot.md) |
| `RMpartgamLD()` | [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGamma.md) |
| `RMpgLDcutoff()` | [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGammaCutoff.md) |
| `RMpgLDplot()` | [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepGammaPlot.md) |
| `RMlocdepQ3cutoff()` | [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3cutoff.md) |
| `RMlocdepQ3plot()` | [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3plot.md) |
| `RMiteminfit()` | [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/reference/RMiteminfit.md) |
| `RMiteminfit_mi()` | [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitMI.md) |
| `RMinfitcutoff()` | [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoff.md) |
| `RMinfitcutoff_mi()` | [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitCutoffMI.md) |
| `RMinfitcutoffPlot()` | [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/reference/RMitemInfitPlot.md) |
| `RMitemrestscore()` | [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/reference/RMitemrestscore.md) |
| `RMbootRestscore()` | [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/reference/RMitemRestscoreBoot.md) |
| `RMciccPlot()` | [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/reference/RMitemICCPlot.md) |
| `RMresidualPCA()` | [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCA.md) |
| `RMpcaCutoff()` | [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdimResidualPCACutoff.md) |
| `RMcfaCutoff()` | [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/reference/RMdimCFACutoff.md) |
| `RMcfaPlot()` | [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/reference/RMdimCFAPlot.md) |
| `RMmartinLof()` | [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/reference/RMdimMartinLof.md) |
| `RMmartinLofResiduals()` | [`RMdimMartinLofResiduals()`](https://pgmj.github.io/easyRasch2/reference/RMdimMartinLofResiduals.md) |
| `RMtileplot()` | [`RMplotTile()`](https://pgmj.github.io/easyRasch2/reference/RMplotTile.md) |
| `RMbarplot()` | [`RMplotBar()`](https://pgmj.github.io/easyRasch2/reference/RMplotBar.md) |
| `RMstackedbarplot()` | [`RMplotStackedbar()`](https://pgmj.github.io/easyRasch2/reference/RMplotStackedbar.md) |

Unchanged:
[`RMreliability()`](https://pgmj.github.io/easyRasch2/reference/RMreliability.md),
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/reference/RMUreliability.md),
[`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/reference/RMitemHierarchy.md),
[`RMdifLR()`](https://pgmj.github.io/easyRasch2/reference/RMdifLR.md),
[`RMdifTree()`](https://pgmj.github.io/easyRasch2/reference/RMdifTree.md),
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/reference/RMlocdepQ3.md),
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/reference/RMtargeting.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/reference/RMscoreSE.md).

## No deprecation aliases

No deprecation shims are shipped: each old name is simply gone. This
keeps the export list clean ahead of the first CRAN submission. Run a
search-and-replace against the table above to migrate.
