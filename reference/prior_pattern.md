# An internal function to change the hyperprior parameters in the pattern mixture model provided by the user depending on the type of missingness mechanism and outcome distributions assumed

This function modifies default hyper prior parameter values in the type
of pattern mixture model selected according to the type of missingness
mechanism and distributions for the outcomes assumed.

## Usage

``` r
prior_pattern(type, dist_e, dist_c, model_txt_info, model_string_jags)
```

## Arguments

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (MAR), Missing Not At Random for the effects (MNAR_eff), Missing Not
  At Random for the costs (MNAR_cost), and Missing Not At Random for
  both (MNAR). For a complete list of all available hyper parameters and
  types of models see the manual.

- dist_e:

  distribution assumed for the effects. Current available choices are:
  Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
  Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative
  Binomial ('negbin') or Bernoulli ('bern')

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')

- model_txt_info:

  list containing model specification information used to write the txt
  file of the JAGS model.

- model_string_jags:

  text file of the model.

## Examples

``` r
#Internal function only
#no examples
#
#
```
