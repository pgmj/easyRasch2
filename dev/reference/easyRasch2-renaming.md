# Function renaming in easyRasch2 0.8.0

In version 0.8.0, 25 exported functions were renamed under one scheme:
domain prefix, then method, then variant suffix. Only the names changed.
Use the table below to update existing scripts.

## Domain prefixes

- `RMdif`: differential item functioning

- `RMlocdep`: local dependence

- `RMitem`: item statistics (fit, restscore, ICC, hierarchy)

- `RMdim`: dimensionality

- `RMplot`: descriptive plots of the response distributions

Variant suffixes attached at the end:

- `Cutoff`: simulation-based critical values

- `Plot`: plots

- `MI`: version for multiply imputed data (replaces `_mi`)

## Old to new name map

|  |  |
|----|----|
| **Old name** | **New name** |
| `RMpartgamDIF()` | [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md) |
| `RMpgDIFcutoff()` | [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md) |
| `RMpgDIFplot()` | [`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md) |
| `RMpartgamLD()` | [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md) |
| `RMpgLDcutoff()` | [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md) |
| `RMpgLDplot()` | [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md) |
| `RMlocdepQ3cutoff()` | [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md) |
| `RMlocdepQ3plot()` | [`RMlocdepQ3Plot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3plot.md) |
| `RMiteminfit()` | [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md) |
| `RMiteminfit_mi()` | [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md) |
| `RMinfitcutoff()` | [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md) |
| `RMinfitcutoff_mi()` | [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md) |
| `RMinfitcutoffPlot()` | [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md) |
| `RMitemrestscore()` | [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md) |
| `RMbootRestscore()` | [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md) |
| `RMciccPlot()` | [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md) |
| `RMresidualPCA()` | [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md) |
| `RMpcaCutoff()` | [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md) |
| `RMcfaCutoff()` | [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md) |
| `RMcfaPlot()` | [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md) |
| `RMmartinLof()` | [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md) |
| `RMmartinLofResiduals()` | [`RMdimMartinLofResiduals()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLofResiduals.md) |
| `RMtileplot()` | [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md) |
| `RMbarplot()` | [`RMplotBar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotBar.md) |
| `RMstackedbarplot()` | [`RMplotStackedbar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotStackedbar.md) |

Unchanged:
[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md),
[`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md),
[`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md),
[`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md),
[`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md),
[`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md),
[`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md),
[`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md).

## No aliases for the old names

The old names were removed without aliases, to keep the export list
clean before the first CRAN release. Search and replace with the table
above to update a script.
[`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
has since been replaced by
[`RMitemInfitPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitPlot.md)
and remains as a deprecated alias.
