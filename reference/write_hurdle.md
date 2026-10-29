# A function to write the code for a hurdle model Alternatives vary depending on the type of distribution assumed for the effect and cost variables, type of structural value mechanism assumed and independence or joint modelling This function selects which type of model to execute.

A function to write the code for a hurdle model Alternatives vary
depending on the type of distribution assumed for the effect and cost
variables, type of structural value mechanism assumed and independence
or joint modelling This function selects which type of model to execute.

## Usage

``` r
write_hurdle(dist_e, dist_c, type, model_txt_info)
```

## Arguments

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`).

- type:

  Type of structural value mechanism assumed. Choices are Structural
  Completely At Random (`SCAR`) and Structural At Random (`SAR`).

- model_txt_info:

  list containing model specification information used to write the
  `txt` file of the `JAGS` model.

## Examples

``` r
#Internal function only
#No examples
#
#
```
