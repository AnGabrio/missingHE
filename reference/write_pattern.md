# An internal function to select which type of pattern mixture model to execute. Alternatives vary depending on the type of distribution assumed for the effect and cost variables, type of missingness mechanism assumed and independence or joint modelling This function selects which type of model to execute.

An internal function to select which type of pattern mixture model to
execute. Alternatives vary depending on the type of distribution assumed
for the effect and cost variables, type of missingness mechanism assumed
and independence or joint modelling This function selects which type of
model to execute.

## Usage

``` r
write_pattern(dist_e, dist_c, type, model_txt_info)
```

## Arguments

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
  Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative
  Binomial ('negbin') or Bernoulli ('bern')

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (MAR), Missing Not At Random for the effects (MNAR_eff), Missing Not
  At Random for the costs (MNAR_cost), and Missing Not At Random for
  both (MNAR)

- model_txt_info:

  list containing model specification information used to write the txt
  file of the JAGS model.

## Examples

``` r
#Internal function only
#No examples
#
#
```
