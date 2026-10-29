# Predictive information criteria for a Bayesian model

Efficient approximate leave-one-out cross validation (`LOO`), deviance
information criterion (`DIC`) and widely applicable information
criterion (`WAIC`) for Bayesian models, calculated on the observed data.

## Usage

``` r
pic(x, criterion = "dic", cases = "cc", ...)
```

## Arguments

- x:

  A `missingHE` object containing the results of a Bayesian model fitted
  in cost-effectiveness analysis using the function
  [`selection()`](https://angabrio.github.io/missingHE/reference/selection.md),
  [`pattern()`](https://angabrio.github.io/missingHE/reference/pattern.md),
  [`hurdle()`](https://angabrio.github.io/missingHE/reference/hurdle.md)
  or [`lmdm()`](https://angabrio.github.io/missingHE/reference/lmdm.md).

- criterion:

  type of information criteria to be produced. Available choices are
  `dic` for the Deviance Information Criterion, `waic` for the Widely
  Applicable Information Criterion, and `looic` for the Leave-One-Out
  Information Criterion.

- cases:

  group of cases for which information criteria should be computed for:
  either `all` for all cases, `cc` for complete cases, `ac_e` and `ac_c`
  for only the observed effect and cost cases, respectively.

- ...:

  Additional parameters that can be provided to manage the output of
  `pic` when `looic` is selected. For more details see
  [bayesplot](https://mc-stan.org/bayesplot/reference/bayesplot-package.html).

## Value

A named list containing different predictive information criteria
results and quantities according to the value of `criterion`. In all
cases, the measures are computed on the observed data for the specific
value selected in `cases`. The elements contained include:

- **d_bar**: Posterior mean deviance (only if `criterion = "dic"`).

- **pD**: Effective number of parameters calculated with the formula
  used by `JAGS` (only if `criterion = "dic"`).

- **dic**: Deviance Information Criterion calculated with the formula
  used by `JAGS` (only if `criterion = "dic"`).

- **d_hat**: Deviance evaluated at the posterior mean of the parameters
  and calculated with the formula used by `JAGS` (only if
  `criterion = "dic"`).

- **elpd**: Expected log pointwise predictive density calculated on the
  selected cases (only if `criterion = "waic"` or
  `criterion = "looic"`).

- **elpd_se**: `elpd` standard error calculated on the selected cases
  (only if `criterion = "waic"` or `criterion = "looic"`).

- **p**: Effective number of parameters calculated on the selected cases
  (only if `criterion = "waic"` or `criterion = "looic"`).

- **p_se**: `p` standard error calculated on the selected cases (only if
  `criterion = "waic"` or `criterion = "looic"`).

- **looic**: The leave-one-out information criterion calculated on the
  selected cases (only if `criterion = "looic"`).

- **looic_se**: The `looic` standard error calculated on the selected
  cases (only if `criterion = "looic"`).

- **waic**: The widely applicable information criterion calculated on
  the selected cases (only if `criterion = "waic"`).

- **waic_se**: The `waic` standard error calculated on the selected
  cases (only if `criterion = "waic"`).

- **pointwise**: A matrix containing the pointwise contributions of each
  of the above measures calculated on the selected cases (only if
  `criterion = "waic"` or `criterion = "looic"`).

- **diagnostics**: A named list containing additional diagnostic
  measures (only if `criterion = "looic"`). See
  [`loo::loo()`](https://mc-stan.org/loo/reference/loo.html) for details
  about interpreting the list elements.

- **psis_object**: A named list containing the matrix of (smoothed) log
  weights (only if `criterion = "looic"` and `save_psis = TRUE`).

## Details

The Deviance Information Criterion (`DIC`), Leave-One-Out Information
Criterion (`LOOIC`) and the Widely Applicable Information Criterion
(`WAIC`) are methods for estimating out-of-sample predictive accuracy
from a Bayesian model using the log-likelihood evaluated at the
posterior simulations of the parameters. `DIC` is computationally simple
to calculate but it is known to have some problems, arising in part from
it not being fully Bayesian in that it is based on a point estimate.
`LOOIC` can be computationally expensive but can be easily approximated
using importance weights that are smoothed by fitting a generalised
Pareto distribution to the upper tail of the distribution of the
importance weights. `WAIC` is fully Bayesian and closely approximates
Bayesian cross-validation. Unlike `DIC`, `WAIC` is invariant to
parameterisation and also works for singular models. In finite cases,
`WAIC` and `LOOIC` give similar estimates, but for influential
observations `WAIC` underestimates the effect of leaving out one
observation.

## References

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling*. (2003).

Vehtari, A. Gelman, A. Gabry, J. (2016a) *Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC*. Statistics
and Computing. Advance online publication.

Vehtari, A. Gelman, A. Gabry, J. (2016b) *Pareto smoothed importance
sampling*. ArXiv preprint.

Gelman, A. Hwang, J. Vehtari, A. (2014) *Understanding predictive
information criteria for Bayesian models*. Statistics and Computing 24,
997-1016.

Watanable, S. (2010). *Asymptotic equivalence of Bayes cross validation
and widely application information criterion in singular learning
theory*. Journal of Machine Learning Research 11, 3571-3594.

## See also

[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`loo::loo()`](https://mc-stan.org/loo/reference/loo.html),
[`loo::waic()`](https://mc-stan.org/loo/reference/waic.html)

## Author

Andrea Gabrio

## Examples

``` r
# For examples see the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#
#
```
