# Extract regression coefficient estimates from objects in the class `missingHE`

Produces a table printout with summary statistics for the regression
coefficients of the health economic evaluation probabilistic model run
using the function
[`selection`](https://angabrio.github.io/missingHE/reference/selection.md),
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md),
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md) or
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

## Usage

``` r
# S3 method for class 'missingHE'
coef(object, prob = c(0.025, 0.975), random = FALSE, digits = 3, ...)
```

## Arguments

- object:

  A `missingHE` object containing the results of the Bayesian modelling
  and the economic evaluation

- prob:

  A numeric vector of probabilities within the range (0,1), representing
  the upper and lower CI sample quantiles to be calculated and returned
  for the estimates.

- random:

  Logical. If `random` is `TRUE`, the estimates of the random effects
  parameters are printed, when available.

- digits:

  Integer indicating the number of decimal places to be used for
  rounding (default = 3).

- ...:

  Additional arguments affecting the summary produced.

## Value

Prints a table with some summary statistics, including posterior mean,
standard deviation and lower and upper quantiles based on the values
specified in `prob`, for the posterior distributions of the regression
coefficients of the effects and costs models run using the function
`selection`, `pattern`, `hurdle` or
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md).

## See also

[`selection`](https://angabrio.github.io/missingHE/reference/selection.md)
[`pattern`](https://angabrio.github.io/missingHE/reference/pattern.md)
[`hurdle`](https://angabrio.github.io/missingHE/reference/hurdle.md)
[`lmdm`](https://angabrio.github.io/missingHE/reference/lmdm.md)
[`diagnostic`](https://angabrio.github.io/missingHE/reference/diagnostic.md)
[`plot.missingHE`](https://angabrio.github.io/missingHE/reference/plot.missingHE.md)

## Author

Andrea Gabrio

## Examples

``` r
 
# For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
# \code{\link{hurdle}} or \code{\link{lmdm}}
# 
# 
```
