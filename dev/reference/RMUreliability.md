# Relative measurement uncertainty (RMU)

Reliability estimated from a matrix of posterior or plausible-value
draws (Bignardi, Kievit & Bürkner, 2026). The columns of `input_draws`
are split at random into two halves, and the reliability is the Pearson
correlation across persons between paired columns of the two halves,
summarised over the pairs as a mean with an HDCI.
[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md)
uses it with plausible values from mirt.

## Usage

``` r
RMUreliability(input_draws, level = 0.95, verbose = FALSE)
```

## Arguments

- input_draws:

  Numeric matrix or data.frame of draws, one row per person and one
  column per draw. At least two columns, ideally many.

- level:

  Numeric in (0, 1). Width of the HDCI returned. Default `0.95`.

- verbose:

  Logical. Print summary information about the input. Default `FALSE`.

## Value

A 1-row data.frame with columns `rmu_estimate`, `hdci_lowerbound`,
`hdci_upperbound`, plus the `.width`/`.point`/`.interval` metadata
columns added by
[`ggdist::mean_hdci()`](https://mjskay.github.io/ggdist/reference/point_interval.html).

## Details

Adapted (with permission, GPL-2/3) from `rmu()`, formerly
`reliability()`, in the 'gbtoolbox' package,
<https://github.com/giac01/gbtoolbox>.

A column pair in which either column has zero variance gets a
correlation of 0, without a message, since the correlation is undefined
there.

The function takes no `seed`. The split of the columns is random and
follows the session's random number stream (see
[easyRasch2-reproducibility](https://pgmj.github.io/easyRasch2/dev/reference/easyRasch2-reproducibility.md)).
Needs the ggdist package (in Suggests).

## References

Bignardi, G., Kievit, R., & Bürkner, P. C. (2026). A general method for
estimating reliability using Bayesian Measurement Uncertainty.
*PsyArXiv*. <https://osf.io/preprints/psyarxiv/h54k8>

## See also

[`RMreliability()`](https://pgmj.github.io/easyRasch2/dev/reference/RMreliability.md)
