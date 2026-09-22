# Next dev version: `RMreliabilityCurve()` annotations

Status: **noted 2026-09-14, not implemented.** 1.3.1 is on CRAN and the
development version is not open yet. Three items, all raised from using the
curve through easyRasch2jmv 3.2.0. None of them changes a computed value.

## 1. Name sigma when the benchmark is drawn on the SEM or information axis

`benchmark` is stated on the reliability scale and converted to whichever axis
is plotted, in `.curve_plot()`:

```r
sem_crit <- sigma * sqrt((1 - benchmark) / benchmark)
# information axis:
1 / sem_crit^2
```

which is conditional reliability solved for information,

```
I_crit = benchmark / ((1 - benchmark) * sigma^2)
```

and collapses to `4 / sigma^2` at a benchmark of 0.80.

`sigma` is the fitted latent SD from `.latent_moments()`, so the converted
benchmark moves with the sample. Two real runs at the same benchmark of 0.80
produced `information >= 2.34` and `information >= 1.41`, which back out to
latent SDs of 1.31 and 1.68. Nothing is wrong with either number. A wider
latent distribution needs less information per person to reach the same
reliability, because reliability is a variance ratio.

The caption does not say this. It states sigma only on the reliability axis
(`sprintf("sigma = %.2f", sigma)`), which is the one axis where the benchmark
needs no conversion. On the SEM and information axes, where it does, sigma is
absent and the converted number looks like a property of the instrument.

It also sits next to a line that is true of the curve but not of the dotted
benchmark:

> The curve is a property of the items, not of the sample.

The information and SEM curves are properties of the items. The benchmark
drawn across them is not, because it carries sigma. Two instruments compared
on the information axis at the same reliability benchmark are being held to
two different thresholds.

**Proposed:** name sigma in the benchmark clause whenever the axis is not
reliability, for example

> Dotted line and shaded band mark information >= 2.34, the information
> needed for conditional reliability 0.80 at the estimated latent SD of 1.31,
> reached by 62% of respondents.

and consider qualifying the "property of the items" sentence so it is not
read as covering the benchmark.

## 2. Bring `.curve_kable()` into line with the module's summary table

The package already has a summary table, `.curve_kable()`, reachable through
`output = "kable"`, and it already carries the theta range where the benchmark
is reached. easyRasch2jmv 3.2.0 built its own Conditional Precision Summary
against the attributes rather than that table, and the two have drifted apart.
Neither is a superset of the other.

| Row | package kable | jmv table |
|---|---|---|
| Latent SD (sigma) | yes | yes |
| Latent mean | **no** | yes |
| Mean theta (person estimates) | **no** | yes |
| SD theta (person estimates) | **no** | yes |
| Marginal reliability (curve mean) | yes | yes |
| Marginal reliability (Green, superseded) | yes | yes |
| Average SEM | yes | yes |
| Average test information | **no** | yes |
| Minimum SEM | yes | **no** |
| Theta at minimum SEM | yes | **no** |
| Theta range reaching the benchmark | yes | footnote only |
| Respondents in that range | yes | yes |

**Proposed:** add the four missing rows to `.curve_kable()`.

- `latent_mean` is already an attribute, so that row is a one-liner. Its
  distance from 0 is what says how far off target the sample is, and it is
  the row that makes the 1.3.1 fix visible.
- Mean and SD of the person estimates need `theta_hat`, which
  `.curve_kable()` does not currently receive. They are worth the argument:
  they are the distribution the density overlay draws, and putting them
  beside the latent moments is what shows that the SD of the estimates is
  the wider of the two because it carries measurement error.
- Average test information as `1 / sem_average^2`, which is already the flat
  reference line on the information axis. Deriving it from the average SEM
  rather than averaging over the grid keeps the two rows describing one
  summary.

Going the other way, the jamovi module should pick up minimum SEM and theta
at minimum SEM, and promote the benchmark theta range from a cell note to a
row. Recorded on the module side separately.

## 3. "reliability" in the benchmark parenthesis

Requested as "marginal reliability". **I think that would be wrong and want a
decision before it goes in.**

The benchmark is a *conditional* reliability threshold. It is compared against
the curve pointwise:

```r
rxx_person <- .curve_stats(thr_list, theta_hat, sigma)$reliability
benchmark_percent <- 100 * mean(rxx_person >= benchmark)
benchmark_range <- .curve_runs(grid, curve$reliability >= benchmark)
```

Marginal reliability is the single latent-density-weighted average of that
curve. It is already in the caption, correctly named, as the dashed reference
line: `"the marginal reliability (%.3f)"`. Calling the dotted benchmark
marginal too would give one word to two different quantities in one caption,
and would suggest the shaded band marks where marginal reliability reaches
0.80, which is not a statement that has a location on the scale.

The underlying complaint is real. The caption uses a bare "reliability" for
three things: the y-axis (conditional), the dashed line (marginal), and the
dotted benchmark (conditional). Only the middle one is qualified.

**Proposed:** qualify the benchmark as **conditional** reliability, matching
the y-axis label the reliability axis already uses:

```
sem         = "SEM <= %.2f logits (conditional reliability %.2f)"
information = "information >= %.2f (conditional reliability %.2f)"
```

If the intent was instead that the benchmark should be a marginal-reliability
criterion, that is a different feature, not a wording change, and it would
need a definition of what region of the scale such a criterion picks out.
