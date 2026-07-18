# Relative Measurement Uncertainty (RMU)

Bayesian-style reliability estimate (Bignardi, Kievit & Bürkner, 2025)
computed from a matrix of posterior or plausible-value draws. The
columns of `input_draws` are split at random into two halves;
reliability is the Pearson correlation across persons of paired columns
from the two halves, summarised across pairs as a posterior mean with
HDCI.

## Usage

``` r
RMUreliability(input_draws, level = 0.95, verbose = FALSE)
```

## Arguments

- input_draws:

  Numeric matrix or data.frame of draws. Rows are subjects; columns are
  draws. Must have at least two columns; ideally many.

- level:

  Numeric in (0, 1). Width of the HDCI returned. Default `0.95`.

- verbose:

  Logical. Print summary information about the input. Default `FALSE`.

## Value

A 1-row data.frame with columns `rmu_estimate`, `hdci_lowerbound`,
`hdci_upperbound`, plus the `.width`/`.point`/`.interval` metadata
columns added by
[`ggdist::mean_hdci()`](https://rdrr.io/pkg/ggdist/man/point_interval.html).

## Details

Adapted (with permission, GPL-2/3) from
<https://github.com/giac01/gbtoolbox/blob/main/R/reliability.R>.

The function silently returns 0 for any column pair where either side
has zero variance (the correlation is undefined there).

Requires the `ggdist` package (Suggests).

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2025). A general method for
estimating reliability using Bayesian Measurement Uncertainty.
*PsyArXiv*.
[doi:10.31234/osf.io/h54k8_v1](https://doi.org/10.31234/osf.io/h54k8_v1)

## See also

[`RMreliability()`](https://pgmj.github.io/easyRasch2/reference/RMreliability.md)
