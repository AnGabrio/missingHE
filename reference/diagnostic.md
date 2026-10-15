# Diagnostic checks for assessing MCMC convergence of Bayesian models fitted in `JAGS` using the function [`selection`](https://angabrio.github.io/missingHE/reference/selection.md), [`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md), [`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md) or [`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

The focus is restricted to full Bayesian models in cost-effectiveness
analyses based on the function
[`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md) and
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md), with
convergence of the MCMC chains that is assessed through graphical checks
of the posterior distribution of the parameters of interest, Examples
are density plots, trace plots, autocorrelation plots, etc. Other types
of posterior checks are related to some summary MCMC statistics that are
able to detect possible issues in the convergence of the algorithm, such
as the potential scale reduction factor or the effective sample size.
Different types of diagnostic tools and statistics are used to assess
model convergence using functions contained in the package **ggmcmc**.
Graphics and plots are managed using functions contained in the package
**ggplot2** and **ggthemes**.

## Usage

``` r
diagnostic(x, type = "denplot", param = "all", theme = NULL, ...)
```

## Arguments

- x:

  An object of class "missingHE" containing the posterior results of a
  full Bayesian model implemented using the function
  [`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
  [`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
  [`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
  or [`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

- type:

  Type of diagnostic check to be plotted for the model parameter
  selected. Available choices include: 'histogram' for histogram plots,
  'denplot' for density plots, 'traceplot' for trace plots, 'acf' for
  autocorrelation plots, 'running' for running mean plots, 'compare' for
  comparing the distribution of the whole chain with only its last part,
  'cross' for cross correlation plots, 'Rhat' for the potential scale
  reduction factor, 'geweke' for the geweke diagnostic, 'pairs' for
  posterior correlation among the parameters,'caterpillar' for
  caterpillar plots.

- param:

  Name of the family of parameters to process, as given by a regular
  expression. For example the mean parameters for the effect and cost
  variables can be specified using 'mu.e' and 'mu.c', respectively.
  Different types of models may have different parameters depending on
  the assumed distributions and missing data assumptions. To see a
  complete list of all possible parameters by types of models assumed
  see details.

- theme:

  Type of ggplot theme among some pre-defined themes, mostly taken from
  the package **ggthemes**. For a full list of available themes see
  details.

- ...:

  Additional parameters that can be provided to manage the graphical
  output of `diagnostic`.

## Value

A **ggplot** object containing the plots specified in the argument
`type`

## Details

Depending on the types of plots specified in the argument `type`, the
output of `diagnostic` can produce different combinations of MCMC visual
posterior checks for the family of parameters indicated in the argument
`param`. For a full list of the available plots see the description of
the argument `type` or see the corresponding plots in the package
**ggmcmc**.

The parameters that can be assessed through `diagnostic` are only those
included in the object `x` (see Arguments). Specific character names
must be specified in the argument `param` according to the specific
model implemented. The available names and the parameters associated
with them are:

- "mu.e" the mean of the effects across treatment arms.

- "mu.c" the mean of the costs across treatment arms.

- "sd.e" the standard deviation of the effects.

- "sd.c" the standard deviation of the costs.

- "alpha" the regression coefficients for the effects.

- "beta" the regression coefficients for the costs.

- "beta.f" the regression coefficients for the costs related to the
  effects predictor.

- "alpha.time" the autoregressive coefficients for the effects (only
  with the function `lmdm`).

- "beta.time" the autoregressive coefficients for the costs (only with
  the function `lmdm`).

- "random.alpha" the regression random effects coefficients for the
  effects.

- "random.beta" the regression random effects coefficients for the
  costs.

- "random.alpha.time" the autoregressive random effects coefficients for
  the effects (only with the function `lmdm`).

- "random.beta.time" the autoregressive random effects coefficients for
  the costs (only with the function `lmdm`).

- "p.e" the probability of missingness or structural values for the
  effects (only with the function `selection`, `hurdle` or `lmdm`).

- "p.c" the probability of missingness or structural values for the
  costs (only with the function `selection`, `hurdle` or `lmdm`).

- "gamma.e" the regression coefficients of missingness or structural
  values for the effects (only with the function `selection`, `hurdle`).

- "gamma.c" the regression coefficientd of missingness or structural
  values for the costs (only with the function `selection`, `hurdle`).

- "random.gamma.e" the random effects regression coefficients of
  missingness or structural values for the effects (only with the
  function `selection`, `hurdle` or `lmdm`).

- "random.gamma.c" the random effects regression coefficients of
  missingness or structural values for the costs (only with the function
  `selection`, `hurdle` or `lmdm`).

- "pattern" the probabilities of the missingness patterns (only with the
  function `pattern`).

- "delta.e" the mnar parameter for the effects (only with the function
  `selection`, `pattern` or `lmdm`).

- "delta.c" the mnar parameters for the costs (only with the function
  `selection`, `pattern` or `lmdm`).

- "random.delta.e" the random effects mnar parameters for the effects
  (only with the function `selection` or `lmdm`).

- "random.delta.c" the random effects mnar parameters for the costs
  (only with the function `selection` or `lmdm`).

- "all" all available parameters stored in `x`.

When the object `x` is created using the function `pattern`,
pattern-specific standard deviation ("sd.e", "sd.c") and regression
coefficient parameters ("alpha", "beta") for both outcomes can be
visualised. The parameters associated with a missingness mechanism can
be accessed only when `x` is created using the function `selection`,
`pattern` or `lmdm`, while the parameters associated with the model for
the structural values mechanism can be accessed only when `x` is created
using the function `hurdle`.

The argument `theme` allows to customise the graphical output of the
plots generated by `diagnostic` and allows to choose among a set of
possible pre-defined themes taken form the package **ggtheme**. For a
complete list of the available character names for each theme, see
**ggthemes**.

## References

Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). *Bayesian Data
Analysis, 2nd edition*, CRC Press.

Brooks, S. Gelman, A. Jones, JL. Meng, XL. (2011). *Handbook of Markov
Chain Monte Carlo*, CRC/Chapman and Hall.

## See also

[`ggs`](https://rdrr.io/pkg/ggmcmc/man/ggs.html)
[`selection`](https://angabrio.github.io/missingHE/reference/selection.md)
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md)
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

## Author

Andrea Gabrio

## Examples

``` r
# For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
# \code{\link{hurdle}} or \code{\link{lmdm}}
#
#
```
