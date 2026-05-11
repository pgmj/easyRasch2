# easyRasch2

<!-- badges: start -->
<a href="https://buymeacoffee.com/pgmj" target="_blank"><img src="https://cdn.buymeacoffee.com/buttons/default-orange.png" alt="Buy Me A Coffee" height="41" width="174"></a>
<!-- badges: end -->

`easyRasch2` is a CRAN-targeted R package for Rasch measurement theory
analysis workflows. It is the successor to
[easyRasch](https://github.com/pgmj/easyRasch), offering a lightweight,
CRAN-ready structure with proper namespacing and minimal dependencies.

A central design choice is **simulation-based critical values** for various fit
statistics. Rather than relying on rule-of-thumb cutoffs, most diagnostics
are paired with a parametric-bootstrap function that generates an empirical
null distribution from the fitted Rasch / PCM model and the observed
sample.

For broader Rasch-analysis tutorials, see the
[vignette](https://pgmj.github.io/raschrvignette/RaschRvign.html) for the
sibling package [`easyRasch`](https://pgmj.github.io/easyRasch/).

## Key design principles

- **Estimation**: Conditional Maximum Likelihood (CML) via `eRm`,
  `psychotools`, and `iarm` for item parameters; Weighted Likelihood
  Estimation (WLE) for person parameters; `mirt` (MML) only where CML
  alternatives do not exist (e.g., Q3 residual correlations).
- **Output**: `knitr::kable()` for tables (Quarto-friendly), `ggplot2`
  for figures, and `"dataframe"` output options for downstream use.
- **Naming**: Functions use the `RM` prefix (e.g., `RMlocdepQ3()`).

## Functions by domain

### Item fit

- `RMiteminfit()` — conditional infit MSQ
- `RMinfitcutoff()` + `RMinfitcutoffPlot()` — simulation-based cutoffs and plot
- `RMiteminfit_mi()` + `RMinfitcutoff_mi()` — multiple-imputation variants
- `RMitemrestscore()` — item-restscore with Goodman-Kruskal's gamma
- `RMbootRestscore()` — non-parametric bootstrap of item-restscore fit

### Local dependence

- `RMlocdepQ3()` + `RMlocdepQ3cutoff()` — Yen's Q3 residual correlations
- `RMpartgamLD()` + `RMpgLDcutoff()` + `RMpgLDplot()` — partial-gamma local dependence

### Dimensionality / unidimensionality

- `RMresidualPCA()` + `RMpcaCutoff()` — PCA of standardized residuals,
  with simulation-based first-contrast cutoff (Chou & Wang, 2010)
- `RMmartinLof()` + `RMmartinLofResiduals()` — Martin-Löf LR test
  (Christensen & Kreiner, 2007), supports polytomous data
- `RMcfaCutoff()` + `RMcfaPlot()` — posterior-predictive CFA fit-index
  cutoffs under PCM unidimensionality (via `lavaan` WLSMV)

### Differential item functioning

- `RMdifLR()` — Andersen's likelihood-ratio test (`eRm::LRtest`)
- `RMdifTree()` — Rasch / partial-credit trees (`psychotree`) with
  Mantel-Haenszel or partial-gamma effect sizes per split, optional
  iterative purification, and `stablelearner`-based stability assessment
- `RMpartgamDIF()` + `RMpgDIFcutoff()` + `RMpgDIFplot()` — partial-gamma DIF

### Reliability, targeting, score conversion

- `RMreliability()` + `RMUreliability()` — Cronbach's α, PSI, empirical
  reliability, and Bignardi-Kievit-Bürkner RMU from plausible values
- `RMtargeting()` — Wright-map style person-item targeting plot
- `RMscoreSE()` — raw-score → logit transformation table (WLE / EAP)
- `RMtileplot()` — response-distribution heatmap with optional group faceting

## Installation

Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("pgmj/easyRasch2")
```

## Example

```r
library(easyRasch2)
data("pcmdat2", package = "eRm")
options(mc.cores = 4)
set.seed(42)

# Conditional item infit with simulation-based cutoffs
simfit <- RMinfitcutoff(pcmdat2, iterations = 250)
RMiteminfit(pcmdat2, cutoff = simfit)

# Test of unidimensionality via posterior-predictive ordinal CFA
cfa_res <- RMcfaCutoff(pcmdat2, iterations = 250)
cfa_res                # kable: observed vs simulated cutoffs
RMcfaPlot(cfa_res)     # histogram + observed diamond

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

- Bjorner, J. B., Kreiner, S., Ware, J. E., Damsgaard, M. T., & Bech, P. (1998).
  Differential item functioning in the Danish translation of the SF-36.
  *Journal of Clinical Epidemiology, 51*(11), 1189–1202.
- Chou, Y.-T., & Wang, W.-C. (2010). Checking dimensionality in item-response
  models with principal component analysis on standardized residuals.
  *Educational and Psychological Measurement, 70*(5), 717–731.
  <https://doi.org/10.1177/0013164410379322>
- Christensen, K. B., & Kreiner, S. (2007). A Monte Carlo approach to
  unidimensionality testing in polytomous Rasch models.
  *Applied Psychological Measurement, 31*(1), 20–30.
  <https://doi.org/10.1177/0146621605286204>
- Christensen, K. B., Makransky, G., & Horton, M. (2017). Critical values for
  Yen's Q3: Identification of local dependence in the Rasch model using
  residual correlations. *Applied Psychological Measurement, 41*(3), 178–194.
  <https://doi.org/10.1177/0146621616677520>
- Henninger, M., Debelak, R., & Strobl, C. (2023). A new stopping criterion
  for Rasch trees based on the Mantel-Haenszel effect size measure for DIF.
  *Educational and Psychological Measurement, 83*, 181–212.
  <https://doi.org/10.1177/00131644221077135>
- Henninger, M., Radek, J., Debelak, R., & Strobl, C. (2025). Partial credit
  trees meet the partial gamma coefficient for quantifying DIF and DSF in
  polytomous items. *Behaviormetrika, 52*, 221–257.
- Johansson, M. (2025). Detecting item misfit in Rasch models.
  *Educational Methods & Psychometrics, 3*(18).
  <https://doi.org/10.61186/emp.2025.5>
- Kreiner, S. (2011). A note on item-restscore association in Rasch models.
  *Applied Psychological Measurement, 35*(7), 557–561.
  <https://doi.org/10.1177/0146621611410227>
- Müller, M. (2020). Item fit statistics for Rasch analysis: Can we trust them?
  *Journal of Statistical Distributions and Applications, 7*(1), 5.
  <https://doi.org/10.1186/s40488-020-00108-7>
- Philipp, M., Rusch, T., Hornik, K., & Strobl, C. (2018). Measuring the
  stability of results from supervised statistical learning. *Journal of
  Computational and Graphical Statistics, 27*, 685–700.
- Rosseel, Y. (2012). lavaan: An R package for structural equation modeling.
  *Journal of Statistical Software, 48*(2), 1–36.
  <https://doi.org/10.18637/jss.v048.i02>
- Strobl, C., Kopf, J., & Zeileis, A. (2015). Rasch trees: A new method for
  detecting DIF in the Rasch model. *Psychometrika, 80*, 289–316.


## Credits

As stated, this is based on my `easyRasch` package, and I am using Claude
to "transfer" functions to this more properly formatted package. While it uses
my earlier code, most of the code in this package is produced by the LLM and
bug fixed by me.

`RMdifTree()` adapts MIT-licensed code from Mirka Henninger and Jan Radek's
[`raschtreeMH`](https://github.com/mirka-henninger/raschtreeMH) and
[`effecttree`](https://github.com/mirka-henninger/effecttree) packages for
the effect-size and ETS-classification algorithms.

[Magnus Johansson](https://ki.se/en/people/magnus-johansson-3) is a licensed psychologist with a PhD in behavior analysis. He works as a research specialist at [Karolinska Institutet](https://ki.se/en/cns/research/centre-for-psychiatry-research), Department of Clinical Neuroscience, Center for Psychiatry Research.

- ORCID: [0000-0003-1669-592X](https://orcid.org/0000-0003-1669-592X)
- Bluesky: [@pgmj.bsky.social](https://bsky.app/profile/pgmj.bsky.social) 

## License

GPL (>= 3)
