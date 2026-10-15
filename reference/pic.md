# Predictive information criteria for Bayesian models fitted in `JAGS` using the function [`selection`](https://angabrio.github.io/missingHE/reference/selection.md), [`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md), [`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md) or [`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md)

Efficient approximate leave-one-out cross validation (LOO), deviance
information criterion (DIC) and widely applicable information criterion
(WAIC) for Bayesian models, calculated on the observed data.

## Usage

``` r
pic(x, criterion = "dic", cases = "cc", ...)
```

## Arguments

- x:

  A `missingHE` object containing the results of a Bayesian model fitted
  in cost-effectiveness analysis using the function
  [`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
  [`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
  [`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
  or [`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

- criterion:

  type of information criteria to be produced. Available choices are
  `'dic'` for the Deviance Information Criterion, `'waic'` for the
  Widely Applicable Information Criterion, and `'looic'` for the
  Leave-One-Out Information Criterion.

- cases:

  group of cases for which information criteria should be computed for:
  either `'all'` for all cases, `'cc'` for complete cases, `'ac_e'` and
  `'ac_c'` for only the observed effect and cost cases, respectively.

- ...:

  Additional parameters that can be provided to manage the output of
  `pic` when `'looic'` is selected. For more details see **bayesplot**.

## Value

A named list containing different predictive information criteria
results and quantities according to the value of `criterion`. In all
cases, the measures are computed on the observed data for the specific
modules of the model selected in `module`.

- d_bar:

  Posterior mean deviance (only if `criterion` is `'dic'`).

- pD:

  Effective number of parameters calculated with the formula used by
  `JAGS` (only if `criterion` is `'dic'`)

- dic:

  Deviance Information Criterion calculated with the formula used by
  `JAGS` (only if `criterion` is `'dic'`)

- d_hat:

  Deviance evaluated at the posterior mean of the parameters and
  calculated with the formula used by `JAGS` (only if `criterion` is
  `'dic'`)

- elpd, elpd_se:

  Expected log pointwise predictive density and standard error
  calculated on the observed data for the model nodes indicated in
  `module` (only if `criterion` is `'waic'` or `'looic'`).

- p, p_se:

  Effective number of parameters and standard error calculated on the
  cases indicated in `cases`. (only if `criterion` is `'waic'` or
  `'looic'`).

- looic, looic_se:

  The leave-one-out information criterion and standard error calculated
  on the cases indicated in `cases`. (only if `criterion` is `'looic'`).

- waic, waic_se:

  The widely applicable information criterion and standard error
  calculated on the cases indicated in `cases`. (only if `criterion` is
  `'waic'`).

- pointwise:

  A matrix containing the pointwise contributions of each of the above
  measures calculated on the cases indicated in `cases`. (only if
  `criterion` is `'waic'` or `'looic'`).

- diagnostics:

  A named list containing additional diagnostic measures (only if
  `criterion` is `'looic'`). See
  [`loo`](https://mc-stan.org/loo/reference/loo.html) for details about
  interpreting the list elements.

- psis_object:

  A named list containing the matrix of (smoothed) log weights (only if
  `criterion` is `'looic'` with the optional argument `'save_psis'` is
  set to `TRUE`). See
  [`loo`](https://mc-stan.org/loo/reference/loo.html) for details about
  interpreting the list elements.

## Details

The Deviance Information Criterion (DIC), Leave-One-Out Information
Criterion (LOOIC) and the Widely Applicable Information Criterion (WAIC)
are methods for estimating out-of-sample predictive accuracy from a
Bayesian model using the log-likelihood evaluated at the posterior
simulations of the parameters. DIC is computationally simple to
calculate but it is known to have some problems, arising in part from it
not being fully Bayesian in that it is based on a point estimate. LOOIC
can be computationally expensive but can be easily approximated using
importance weights that are smoothed by fitting a generalised Pareto
distribution to the upper tail of the distribution of the importance
weights. WAIC is fully Bayesian and closely approximates Bayesian
cross-validation. Unlike DIC, WAIC is invariant to parameterisation and
also works for singular models. In finite cases, WAIC and LOOIC give
similar estimates, but for influential observations WAIC underestimates
the effect of leaving out one observation.

## References

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling.* (2003).

Vehtari, A. Gelman, A. Gabry, J. (2016a) Practical Bayesian model
evaluation using leave-one-out cross-validation and WAIC. *Statistics
and Computing*. Advance online publication.

Vehtari, A. Gelman, A. Gabry, J. (2016b) Pareto smoothed importance
sampling. *ArXiv* preprint.

Gelman, A. Hwang, J. Vehtari, A. (2014) Understanding predictive
information criteria for Bayesian models. *Statistics and Computing* 24,
997-1016.

Watanable, S. (2010). Asymptotic equivalence of Bayes cross validation
and widely application information criterion in singular learning
theory. *Journal of Machine Learning Research* 11, 3571-3594.

## See also

[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`loo`](https://mc-stan.org/loo/reference/loo.html),
[`waic`](https://mc-stan.org/loo/reference/waic.html)

## Author

Andrea Gabrio

## Examples

``` r
 
# For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
# \code{\link{hurdle}} or \code{\link{lmdm}}
# 
# 
```
