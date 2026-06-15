# easyRasch2

`easyRasch2` is an R package for Rasch measurement theory analysis
workflows. It is the successor to
[`easyRasch`](https://pgmj.github.io/easyRasch/), offering a lightweight
and consistent structure with proper namespacing and minimal
dependencies.

A central design choice is **simulation-based critical values** for
various fit statistics. Rather than relying on rule-of-thumb cutoffs,
most diagnostics are paired with a parametric-bootstrap function that
generates an empirical null distribution from the fitted Rasch / PCM
model and the observed sample.

The [Get
Started](https://pgmj.github.io/easyRasch2/articles/easyRasch2.html)
link above contains a short introduction. For broader Rasch-analysis
tutorials, see the
[vignette](https://pgmj.github.io/raschrvignette/RaschRvign.html) for
the archived sibling package
[`easyRasch`](https://pgmj.github.io/easyRasch/).

## Installation

Install from CRAN:

``` r

install.packages("easyRasch2")
```

Install the development version from GitHub:

``` r

# install.packages("remotes") # if needed
remotes::install_github("pgmj/easyRasch2")
```

## Key design principles

- **Estimation**: Conditional Maximum Likelihood (CML) via `eRm`,
  `psychotools`, and `iarm` for item parameters; Weighted Likelihood
  Estimation (WLE) for person parameters; `mirt` (MML) only where CML
  alternatives do not exist (e.g., Q3 residual correlations).
- **Output**:
  [`knitr::kable()`](https://rdrr.io/pkg/knitr/man/kable.html) for
  tables (Quarto-friendly), `ggplot2` for figures, and `"dataframe"`
  output options for downstream use.
- **Naming**: Functions use the `RM` prefix (e.g.,
  [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md)).

## Functions by domain

### Item fit

- [`RMitemInfit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMiteminfit.md)
  — conditional infit MSQ
- [`RMitemInfitCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoff.md) +
  [`RMitemInfitCutoffPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffPlot.md)
  — simulation-based cutoffs and plot
- [`RMitemInfitMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitMI.md) +
  [`RMitemInfitCutoffMI()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemInfitCutoffMI.md)
  — multiple-imputation variants
- [`RMitemRestscore()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemrestscore.md)
  — item-restscore with Goodman-Kruskal’s gamma
- [`RMitemRestscoreBoot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemRestscoreBoot.md)
  — non-parametric bootstrap of item-restscore fit
- [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md) -
  conditional item characteristic curves

### Local dependence

- [`RMlocdepQ3()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3.md) +
  [`RMlocdepQ3Cutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepQ3cutoff.md)
  — Yen’s Q3 residual correlations
- [`RMlocdepGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGamma.md) +
  [`RMlocdepGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaCutoff.md) +
  [`RMlocdepGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMlocdepGammaPlot.md)
  — partial-gamma local dependence

### Dimensionality / unidimensionality

- [`RMdimResidualPCA()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCA.md) +
  [`RMdimResidualPCACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimResidualPCACutoff.md)
  — PCA of standardized residuals, with simulation-based first-contrast
  cutoff (Chou & Wang, 2010)
- [`RMdimMartinLof()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLof.md) +
  [`RMdimMartinLofResiduals()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimMartinLofResiduals.md)
  — Martin-Löf LR test (Christensen & Kreiner, 2007), supports
  polytomous data
- [`RMdimCFACutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFACutoff.md) +
  [`RMdimCFAPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdimCFAPlot.md)
  — posterior-predictive CFA fit-index cutoffs under PCM
  unidimensionality (via `lavaan` WLSMV)

### Differential item functioning

- [`RMdifLR()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifLR.md)
  — Andersen’s likelihood-ratio test
  ([`eRm::LRtest`](https://rdrr.io/pkg/eRm/man/LRtest.html))
- [`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md)
  — Rasch / partial-credit trees (`psychotree`) with Mantel-Haenszel or
  partial-gamma effect sizes per split, optional iterative purification,
  and `stablelearner`-based stability assessment
- [`RMdifGamma()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGamma.md) +
  [`RMdifGammaCutoff()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaCutoff.md) +
  [`RMdifGammaPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifGammaPlot.md)
  — partial-gamma DIF
- [`RMitemICCPlot()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemICCPlot.md) -
  evaluates DIF across class intervals

### Item category threshold ordering

- [`RMitemCatProb()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemCatProb.md)
  for classic item probability function trace plots.
- [`RMitemHierarchy()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemHierarchy.md)
  renders a plot sorting items based on difficulty, showing item
  threshold locations and confidence intervals.

### Reliability, targeting, score conversion

- [`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md) +
  [`RMUreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMUreliability.md)
  — Cronbach’s α, PSI, empirical reliability, and Relative Measurement
  Uncertainty from plausible values
- [`RMtargeting()`](https://pgmj.github.io/easyRasch2/dev/reference/RMtargeting.md)
  — Wright-map style person-item targeting plot
- [`RMscoreSE()`](https://pgmj.github.io/easyRasch2/dev/reference/RMscoreSE.md)
  — raw-score → logit transformation table (WLE / EAP)

### Item & person parameters

- [`RMitemParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMitemParameters.md)
  — item difficulty / threshold locations in long or wide format, with
  optional standard errors and confidence intervals (CML via `eRm`, or
  MML via `mirt` for sparse data)
- [`RMpersonParameters()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonParameters.md)
  — per-respondent person locations (WLE or EAP), estimated on each
  response pattern so partial missingness is handled directly

### Person fit

- [`RMpersonFit()`](https://pgmj.github.io/easyRasch2/dev/reference/RMpersonFit.md)
  — per-respondent conditional infit / outfit MSQ and the standardized
  log-likelihood `lz`, with resampling-based p-values rather than
  unreliable asymptotic nulls (Sinharay, 2016; Müller, 2020)

### Data visualization

- [`RMplotTile()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotTile.md)
  — response-distribution heatmap with optional group faceting
- [`RMplotBar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotBar.md)
  &
  [`RMplotStackedbar()`](https://pgmj.github.io/easyRasch2/dev/reference/RMplotStackedbar.md)

## Example

``` r

library(easyRasch2)
data("pcmdat2", package = "eRm")
options(mc.cores = 4)
set.seed(42)

# Conditional item infit with simulation-based cutoffs
simfit <- RMitemInfitCutoff(pcmdat2, iterations = 250)
RMitemInfit(pcmdat2, cutoff = simfit)

# Test of unidimensionality via posterior-predictive ordinal CFA
cfa_res <- RMdimCFACutoff(pcmdat2, iterations = 250)
cfa_res                # kable: observed vs simulated cutoffs
RMdimCFAPlot(cfa_res)     # histogram + observed diamond

# DIF analysis via Andersen's LR test
grp <- factor(sample(c("A", "B"), nrow(pcmdat2), replace = TRUE))
RMdifLR(pcmdat2, dif_var = grp)

# Rasch-tree DIF with effect-size classification on continuous +
# categorical covariates simultaneously
covs <- data.frame(
  group = grp,
  band  = sample(c("low", "high"), nrow(pcmdat2), replace = TRUE)
)
RMdifTree(pcmdat2, covariates = covs)
```

## References

- Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P.
  (1998). Differential item functioning in the Danish translation of the
  SF-36. *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
  <https://doi.org/10.1016/S0895-4356(98)00111-5>
- Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in
  item-response models with principal component analysis on standardized
  residuals. *Educational and Psychological Measurement, 70*(5),
  717–731. <https://doi.org/10.1177/0013164410379322>
- Christensen, K. B., & Kreiner, S. (2007). A Monte Carlo approach to
  unidimensionality testing in polytomous Rasch models. *Applied
  Psychological Measurement, 31*(1), 20–30.
  <https://doi.org/10.1177/0146621605286204>
- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical
  values for Yen’s Q3: Identification of local dependence in the Rasch
  model using residual correlations. *Applied Psychological Measurement,
  41*(3), 178–194. <https://doi.org/10.1177/0146621616677520>
- Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping
  criterion for Rasch trees based on the Mantel-Haenszel effect size
  measure for DIF. *Educational and Psychological Measurement, 83*,
  181–212. <https://doi.org/10.1177/00131644221077135>
- Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial
  credit trees meet the partial gamma coefficient for quantifying DIF
  and DSF in polytomous items. *Behaviormetrika, 52*, 221–257.
  <https://doi.org/10.1007/s41237-024-00252-3>
- Johansson, M. (2025). Detecting item misfit in Rasch models.
  *Educational Methods & Psychometrics, 3*(18).
  <https://doi.org/10.61186/emp.2025.5>
- Kreiner, S. (2011). A note on item-restscore association in Rasch
  models. *Applied Psychological Measurement, 35*(7), 557–561.
  <https://doi.org/10.1177/0146621611410227>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we
  trust them? *Journal of Statistical Distributions and Applications,
  7*(1), 5. <https://doi.org/10.1186/s40488-020-00108-7>
- Philipp, M., Rusch, T., Hornik, K., & Strobl, C. (2018). Measuring the
  stability of results from supervised statistical learning. *Journal of
  Computational and Graphical Statistics, 27*, 685–700.
  <https://doi.org/10.1080/10618600.2018.1473779>
- Rosseel, Y. (2012). lavaan: An R package for structural equation
  modeling. *Journal of Statistical Software, 48*(2), 1–36.
  <https://doi.org/10.18637/jss.v048.i02>
- Sinharay, S. (2016). Assessment of person fit using resampling-based
  approaches. *Journal of Educational Measurement, 53*(1), 63–85.
  <https://doi.org/10.1111/jedm.12101>
- Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method
  for detecting DIF in the Rasch model. *Psychometrika, 80*, 289–316.
  <https://doi.org/10.1007/s11336-013-9388-3>

## Credits

As mentioned earlier, this is based on my `easyRasch` package, and I am
using Claude to “transfer” functions to this more properly formatted
package. While it uses my earlier code, most of the code in this package
is produced by the LLM and bug fixed by me.

[`RMdifTree()`](https://pgmj.github.io/easyRasch2/dev/reference/RMdifTree.md)
adapts MIT-licensed code from Mirka Henninger and Jan Radek’s
[`raschtreeMH`](https://github.com/mirka-henninger/raschtreeMH) and
[`effecttree`](https://github.com/mirka-henninger/effecttree) packages
for the effect-size and ETS-classification algorithms.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a
licensed psychologist with a PhD in behavior analysis. He works as a
research specialist at [Karolinska
Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research),
Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky:
  [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social)

## License

GPL (\>= 3)
