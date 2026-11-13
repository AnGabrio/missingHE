# Bayesian Longitudinal Missing Data Models

Full Bayesian cost-effectiveness models to handle missing data in
longitudinal outcomes under different missing data mechanism
assumptions, using alternative parametric distributions for the effect
and cost variables. The analysis is performed using the `BUGS` language,
which is implemented in the software `JAGS` using the function
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html) The output
is stored in an object of class `missingHE`.

## Usage

``` r
lmdm(
  data,
  model.eff,
  model.cost,
  model.me = me ~ 1,
  model.mc = mc ~ 1,
  dist_e,
  dist_c,
  type,
  time_dep = "AR1",
  prob = c(0.025, 0.975),
  n.chains = 2,
  n.iter = 10000,
  n.burnin = floor(n.iter/2),
  inits = NULL,
  n.thin = 1,
  save_model = FALSE,
  prior = "default",
  center = FALSE,
  ...
)
```

## Arguments

- data:

  A data frame in which to find the longitudinal variables supplied in
  `model.eff`, `model.cost` (model formulas for effects and costs) and
  `model.me`, `model.mc` (model formulas for the missing effect and cost
  models). Among these, effectiveness, cost, time and treatment
  indicator variables must always be provided and named `e`, `c`, `time`
  and `trt`, respectively.

- model.eff:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic effectiveness outcome (`e`) whose
  name must correspond to that used in `data`. Any covariates in the
  model must be provided on the right-hand side of the formula. If there
  are no covariates, `1` should be specified on the right hand side of
  the formula. By default, covariates are placed on the "location"
  parameter of the distribution through a linear model. Random effects
  can also be specified for each model parameter. See details for how
  these can be specified.

- model.cost:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be a health economic cost outcome (`c`) whose name must
  correspond to that used in `data`. Any covariates in the model must be
  provided on the right-hand side of the formula. If there are no
  covariates, `1` should be specified on the right hand side of the
  formula. By default, covariates are placed on the "location" parameter
  of the distribution through a linear model. A joint bivariate
  distribution for effects and costs can be specified by including `e`
  on the right-hand side of the formula for the costs model. Random
  effects can also be specified for each model parameter. See details
  for how these can be specified.

- model.me:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `me` (missing effects) and
  any covariates must be provided on the right-hand side of the formula.
  If there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.mc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `mc` (missing costs) and any
  covariates must be provided on the right-hand side of the formula. If
  there are no covariates, `1` should be specified on the right hand
  side of the formula. By default, covariates are placed on the
  "probability" parameter for the missing costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal ('norm'), Gamma (`gamma`) or LogNormal (`lnorm`).

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (`MAR`) and Missing Not At Random (`MNAR`).

- time_dep:

  Type of dependence structure assumed between effectiveness and cost
  outcomes. Current choices include: autoregressive structure of order
  one (`AR1`) - default, bivariate at each time (`biv`) and independence
  (`none`).

- prob:

  A numeric vector of probabilities within the range `(0, 1)`,
  representing the upper and lower CI sample quantiles to be calculated
  and returned for the imputed values.

- n.chains:

  Number of chains.

- n.iter:

  Number of iterations.

- n.burnin:

  Number of warmup iterations.

- inits:

  A list with elements equal to the number of chains selected; each
  element of the list is itself a list of starting values for the `JAGS`
  model, or a function creating (possibly random) initial values. If
  `inits` is `NULL`, `JAGS` will generate initial values for all the
  model parameters.

- n.thin:

  Thinning interval.

- save_model:

  Logical. If `save_model` is `TRUE`, a `txt` file containing the model
  code is printed in the current working directory.

- prior:

  A list containing the hyperprior values provided by the user. Each
  element of this list must be a vector containing the user-provided
  prior distribution and parameter values and must be named with the
  name of the corresponding parameter. For example, the hyperprior
  values for the standard deviation parameter for the effects can be
  provided using the list
  `prior = list("sigma.prior.e" = c("unif", 0, 10))`. For more
  information about how to provide prior hypervalues for different types
  of parameters and models see details. If `prior` is set to `default`,
  the default values will be used.

- center:

  Logical. If `center` is `TRUE`, all covariates in both the effect and
  cost models will be centred.

- ...:

  Additional arguments that can be provided by the user. Examples are
  the additional arguments that can be provided to the function
  [`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html)
  to summarise the health economic evaluation results.

## Value

An object of the class 'missingHE' containing the following elements:

- **data_set**: A list containing the original data set provided in
  `data` (see Arguments). Additional information is also included about,
  among others, the number of observed and missing individuals, the
  total number of individuals by treatment arm and the indicator vectors
  for the missing values.

- **model_output**: A list containing the output of a `JAGS` model
  generated from the functions
  [`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html), and the
  posterior samples for the main parameters of the model.

- **cea**: A list containing the output of the economic evaluation
  performed using the function
  [`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html).

- **type**: A character variable that indicate which type of missing
  value mechanism used to run the model, either `MAR` or `MNAR` (see
  details).

- **data_format**: A character variable that indicates which type of
  analysis was conducted, either using a `wide` or `long` dataset.

- **time_dep**: A character variable that indicate which type of time
  dependence assumption was made, either `none` or `AR1`.

- **w_e**: A numeric vector containing the weight values used to
  calculate aggregated effectiveness mean parameters used in the `cea`
  output.

- **w_c**: A numeric vector containing the weight values used to
  calculate aggregated cost mean parameters used in the `cea` output.

- **call**: The matched call.

- **prior**: The custom prior choices specified by the user. If no
  custom prior specified, default priors are used (see details).

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of missingness
mechanism specified in the argument `type`, different models are built
and run in the background by the function `lmdm`. These models fit
multinomial logistic regressions to estimate the probability of a given
missing data pattern `k` (1 = completers, 2 = intermittent, 3 = dropout)
in one or both longitudinal outcomes. A simple example can be used to
show how longitudinal missing data models are specified. Consider a
longitudinal data set comprising a response variable \\y\\ measured at
\\S\\ occasions and a set of centered covariate \\X_j\\, for \\i = j,
..., J\\. For each subject in the trial \\i = 1, ..., n\\ and time \\s =
1, ..., S\\ we define an indicator variable \\m_i\\ taking value `k = 1`
if the \\i\\-th individual is associated with no missing value
(completer), a value `k = 2` for intermittent missingness over the study
period, and a value `k = 3` for dropout. This is modelled as:

\$\$m_i ~ Multinomial(\pi^k_i)\$\$ \$\$\pi^k_i = \phi^k_i/\sum\phi_i\$\$
\$\$log(\phi^k_i) = \sum\gamma^k_j X_j + \delta^k y_i\$\$ where

- \\\pi^k_i\\ is the individual probability of a missing value in \\y\\
  for pattern \\k\\ at a given time.

- \\\gamma^k_j\\ represents the impact on the missing data probability
  in \\y\\ of the covariate \\X_j\\ for pattern \\k\\ at a given time.

- \\\delta^k\\ represents the impact on the missing data probability in
  \\y\\ for pattern \\k\\ of missingness itself at a given time.

When \\\delta = 0\\ the model assumes a `MAR` mechanism, while when
\\\delta != 0\\ the mechanism is `MNAR`. For the parameters indexing the
missingness model, the default prior distributions assumed are:

- \\\gamma^k_j ~ Normal(0, 0.01)\\.

- \\\delta^k ~ Normal(0, 1)\\.

When user-defined hyperprior values are supplied via the argument
`prior` in the function `lmdm`, the elements of this list (see
Arguments) must be vectors containing the user-provided distribution
name and hyperprior values and must take specific names according to the
parameters they are associated with. Specifically, the names for the
parameters indexing the model which are accepted by **missingHE** are
the following:

- auxiliary parameters \\\sigma\\: `sigma.prior.e` (effects) and/or
  `sigma.prior.c` (costs).

- covariate parameters \\\alpha_j\\ and \\\beta_j\\: `alpha.prior`
  (effects) and/or `beta.prior` (costs).

- covariate parameters in the missingness model \\\gamma_j\\ (if
  covariate data provided): `gamma.prior.e` (effects) and/or
  `gamma.prior.c` (costs).

- mnar parameter \\\delta\\: `delta.prior.e` (effects) and/or
  `delta.prior.c` (costs).

For simplicity, here we have assumed that the set of covariates \\X_j\\
used in the models for the effects/costs and in the model of the missing
effect/cost values is the same. However, it is possible to specify
different sets of covariates for each model using the arguments in the
function `lmdm` (see Arguments).

For each model, random effects can also be specified for each parameter
by adding the term `+ (x | z)` to each model formula, where `x` is the
fixed regression coefficient for which also the random effects are
desired and `z` is the clustering variable across which the random
effects are specified (must be the name of a factor variable in the
dataset).

For an overview of the general modeling framework and additional package
utilities, see
[missingHE-package](https://angabrio.github.io/missingHE/reference/missingHE-package.md).

## References

Gabrio, A., Daniels, MJ. Baio, G. (2019). *A Bayesian parametric
approach to handle missing longitudinal outcome data in trial-based
health economic evaluations*. Journal of the Royal Statistical Society
Series A: Statistics in Society, 183(2), 607-629 \<10.1111/rssa.12522\>

Mason, AJ. Gomes, M. Carpenter, J. Grieve, R. (2021). *Flexible Bayesian
longitudinal models for cost‐effectiveness analyses with informative
missing data*. Health economics, 30(12), 3138-3158. \<10.1002/hec.4408\>

Daniels, MJ. Hogan, JW. *Missing Data in Longitudinal Studies:
strategies for Bayesian modelling and sensitivity analysis*, CRC/Chapman
Hall. <doi:10.1201/9781420011180>

Baio, G.(2012). *Bayesian Methods in Health Economics*. CRC/Chapman
Hall, London. <doi:10.1201/b13099>

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling*. (2003).
<https://www.r-project.org/conferences/DSC-2003/Proceedings/Plummer.pdf>

## See also

[missingHE-package](https://angabrio.github.io/missingHE/reference/missingHE-package.md)
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html).

## Author

Andrea Gabrio

## Examples

``` r
if (isTRUE(tryCatch(requireNamespace("rjags", quietly = TRUE), error = function(e) FALSE))) {
# Quick example to run using subset of PBS dataset

# Load longitudinal dataset

PBS.long <- PBS

# \donttest{
# Run the model using the long_miss function assuming a MAR mechanism
# Use only 100 iterations to run a quick check
model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 100, time_dep = "none")

# Extract regression coefficient estimates
coef(model.long)
#

# Summarise the CEA information from the model
summary(model.long)

# Further examples which take longer to run
model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt + age,
   model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
   type = "MAR", n.chains = 2, n.iter = 500, time_dep = "none")

# Use looic to assess model fit
pic.looic <- pic(model.long, criterion = "looic")
pic.looic

# Show density plots for all parameters
diag.den <- diagnostic(model.long, type = "denplot", param = "alpha")

# Plots of imputations for all effect data
p1 <- plot(model.long, class = "scatter", outcome = "effects", time.plot = "all")

# Summarise the CEA results
summary(model.long)

# }
}
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 2852
#>    Unobserved stochastic nodes: 112
#>    Total graph size: 16272
#> 
#> Initializing model
#> 
#> $Effects
#>                Mean    SD    QL    QU
#> (Intercept).1 0.475 0.035 0.408 0.533
#> trt2.1        0.090 0.046 0.012 0.190
#> (Intercept).2 0.501 0.031 0.443 0.556
#> trt2.2        0.138 0.047 0.039 0.230
#> (Intercept).3 0.478 0.029 0.419 0.526
#> trt2.3        0.138 0.040 0.064 0.207
#> 
#> $Costs
#>                  Mean     SD      QL      QU
#> (Intercept).1 147.799 48.826  62.582 248.354
#> trt2.1        100.570 39.491  34.235 167.534
#> (Intercept).2  73.607 33.988  17.269 139.259
#> trt2.2         45.972 28.579 -23.516  94.459
#> (Intercept).3  49.996 30.584 -21.641 104.961
#> trt2.3         30.893 31.002 -21.537  95.492
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.491 0.032 0.439 0.554
#> 2 0.615 0.031 0.554 0.673
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD   QL   QU
#> 1  244 435 -555 1108
#> 2  451 450 -467 1340
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 24292 1669 21540 27643
#> 2 30293 1551 27081 33234
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 2852
#>    Unobserved stochastic nodes: 115
#>    Total graph size: 16835
#> 
#> Initializing model
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.489 0.022 0.449 0.531
#> 2 0.612 0.023 0.564 0.656
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD   QL   QU
#> 1 5395 343 4746 6078
#> 2 5834 352 5173 6526
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 19042 1129 16809 21224
#> 2 24779 1207 22298 27224
#
```
