# A function to read and re-arrange the data in different ways

This internal function imports the data and outputs only those variables
that are needed to run the model according to the information provided
by the user.

## Usage

``` r
data_read_lmdm(
  data,
  model.eff,
  model.cost,
  model.me,
  model.mc,
  cov_matrix,
  type,
  center,
  fixed_e,
  fixed_c,
  random_e,
  random_c,
  trt_lev,
  trt_pos
)
```

## Arguments

- data:

  A data frame in which to find variables supplied in `model.eff` and
  `model.cost`. Among these, effectiveness, cost and treatment indicator
  variables must always be provided and named 'e', 'c' and 'trt'
  respectively.

- model.eff:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics effectiveness outcome ('e') whose
  name must correspond to that used in `data`, and any covariates are
  given on the right-hand side. If there are no covariates, specify `1`
  on the right hand side. By default, covariates are placed on the
  "location" parameter of the distribution through a linear model.
  Random effects can also be specified for each model parameter.

- model.cost:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economics cost outcome ('c') whose name must
  correspond to that used in `data`, and any covariates are given on the
  right-hand side. If there are no covariates, specify `1` on the right
  hand side. By default, covariates are placed on the "location"
  parameter of the distribution through a linear model. Random effects
  can also be specified for each model parameter.

- model.me:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term 'me'(missing effects) and any
  covariates used to estimate the probability of missing effects are
  given on the right-hand side. If there are no covariates, specify `1`
  on the right hand side. By default, covariates are placed on the
  "probability" parameter for the missing effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter.

- model.mc:

  A formula expression in conventional R linear modelling syntax. The
  response must be indicated with the term 'mc'(missing costs) and any
  covariates used to estimate the probability of missing costs should be
  given on the right-hand side. If there are no covariates, specify `1`
  on the right hand side. By default, covariates are placed on the
  "probability" parameter for the missing costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter.

- cov_matrix:

  Data frame containing the covariate matrix of the model.

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (MAR) and Missing Not At Random (MNAR).

- center:

  Logical. If `center` is `TRUE` all the covariates in the model are
  centered.

- fixed_e:

  Fixed effects variables to be included in the effects model

- fixed_c:

  Fixed effects variables to be included in the costs model

- random_e:

  Random effects variables to be included in the effects model

- random_c:

  Random effects variables to be included in the costs model

- trt_lev:

  Vector of names of each treatment factor level

- trt_pos:

  Vector of name indices of each treatment factor level

## Examples

``` r
#Internal function only
#no examples
#
#
```
