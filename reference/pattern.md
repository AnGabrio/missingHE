# Bayesian Pattern Mixture Models

Full Bayesian cost-effectiveness models to handle missing data in the
outcomes under different missing data mechanism assumptions, using
alternative parametric distributions for the effect and cost variables
and using a pattern mixture model approach to identify the model. The
analysis is performed using the `BUGS` language, which is implemented in
the software `JAGS` using the function
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html) The output
is stored in an object of class `missingHE`.

## Usage

``` r
pattern(
  data,
  model.eff,
  model.cost,
  dist_e,
  dist_c,
  type,
  restriction = "CC",
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

  A data frame in which to find the variables supplied in `model.eff`,
  `model.cost` (model formulas for effects and costs). Among these,
  effectiveness, cost and treatment indicator variables must always be
  provided and named `e`, `c` and `trt`, respectively.

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

- dist_e:

  Distribution assumed for the effects. Current available choices are:
  Normal (`norm`), Beta (`beta`), Gamma (`gamma`), Exponential (`exp`),
  Weibull (`weib`), Logistic (`logis`), Poisson (`pois`), Negative
  Binomial (`negbin`) or Bernoulli (`bern`).

- dist_c:

  Distribution assumed for the costs. Current available choices are:
  Normal (`norm`), Gamma (`gamma`) or LogNormal (`lnorm`).

- type:

  Type of missingness mechanism assumed. Choices are Missing At Random
  (`MAR`) and Missing Not At Random (`MNAR`).

- restriction:

  type of identifying restriction to be imposed to identify the
  distributions of the missing data in each pattern. Available choices
  are: complete case restriction (`CC`) - default - or available case
  restriction (`AC`).

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

- **restriction**: A character name that indicates under which type of
  restrictions the model was identified (see details).

- **call**: The matched call.

- **prior**: The custom prior choices specified by the user. If no
  custom prior specified, default priors are used (see details).

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of missingness
mechanism specified in the argument `type`, different pattern mixture
models are built and run in the background by the function `pattern`.
The model for the outcomes is fitted in each missingness pattern and the
parameters indexing the missing data distributions are identified using:
the corresponding parameters identified from the observed data in other
patterns (under `MAR`); or a combination of the parameters identified by
the observed data and some sensitivity parameters (under `MNAR`). A
simple example can be used to show how pattern mixture models are
specified. Consider a data set comprising a response variable \\y\\ and
a set of centered covariates \\X_j\\. We denote with \\d_i\\ the
patterns' indicator variable for each subject in the trial \\i = 1, ...,
n\\ such that: \\d_i = 1\\ indicates the completers (both e and c
observed), \\d_i = 2\\ and \\d_i = 3\\ indicate that only the costs or
effects are observed, respectively, while \\d_i = 4\\ indicates that
neither of the two outcomes is observed. \\d_i\\ is assigned a
multinomial distribution, which probabilities are modelled using a
Dirichlet prior. Next, the outcomes model is fitted in each pattern and
parameters that cannot be identified in each pattern (d = 2, 3, 4), e.g.
\\mu_e\[d\]\\ and \\mu_c\[d\]\\, are identified using some restrictions
based on the parameters estimated from other patterns. Two choices are
currently available: the complete cases (`CC`) or available cases
(`AC`). For example, using the `CC` restriction, the parameters indexing
the distributions of the missing data are identified as:

\$\$mu_e\[2\] = \mu_e\[4\] = \mu_e\[1\] + \delta_e\$\$ \$\$mu_c\[3\] =
\mu_c\[4\] = \mu_c\[1\] + \delta_c\$\$ where

- \\\mu_e\[1\]\\ is the effects mean for the completers.

- \\\mu_c\[1\]\\ is the costs mean for the completers.

- \\\delta_e\\ is the sensitivity parameters associated with the
  marginal effects mean.

- \\\delta_c\\ is the sensitivity parameters associated with the
  marginal costs mean.

If the `AC` restriction is chosen, only the parameters estimated from
the observed data in pattern 2 (costs) and pattern 3 (effects) are used
to identify those in the other patterns. When \\\delta_e = 0\\ and
\\\delta_c = 0\\ the model assumes a `MAR` mechanism. When \\\delta_e !=
0\\ and/or \\\delta_c != 0\\ `MNAR` departures for the effects and/or
costs are explored with priors on sensitivity parameters that must be
provided by the user (see Arguments).

When user-defined hyperprior values are supplied via the argument
`prior` in the function `pattern`, the elements of this list (see
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

- sensitivity parameters \\\delta\\: `delta.prior.e` (effects) and/or
  `delta.prior.c` (costs).

- missingness patterns' probabilities \\\pi\\: `patterns.prior`.

For each outcome model, random effects can also be specified for each
parameter by adding the term `+ (x | z)` to each model formula, where
`x` is the fixed regression coefficient for which also the random
effects are desired and `z` is the clustering variable across which the
random effects are specified (must be the name of a factor variable in
the dataset).

For an overview of the general modeling framework and additional package
utilities, see
[missingHE-package](https://angabrio.github.io/missingHE/reference/missingHE-package.md).

## References

Molenberghs, G. Verbeke G. *Pattern-Mixture Models*, Linear Mixed Models
for Longitudinal Data, 275-293 \<10.1007/978-1-4419-0300-6_18\>

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
# Quick example to run using subset of MenSS dataset
MenSS.subset <- MenSS[50:100, ]

# Run the model using the pattern function assuming a MAR mechanism
# Use only 100 iterations to run a quick check
model.pattern <- pattern(data = MenSS.subset,model.eff = e ~ trt, model.cost = c ~ trt,
   dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 100)

# Print the results of the JAGS model
print(model.pattern)
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.pattern, criterion = "dic", cases = "cc")
pic.dic
#

# Extract regression coefficient estimates
coef(model.pattern)
#
DONTSHOW({
# Use waic information criterion to assess model fit
pic.waic <- pic(model.pattern, criterion = "waic", cases = "cc")
pic.waic
})

# Assess model convergence using graphical tools
# Produce histograms of the posterior samples for the mean effects
diag.hist <- diagnostic(model.pattern, type = "histogram", param = "mu.e")
#

# Compare observed effect data with imputations from the model
# using plots (posteiror means and credible intervals)
p1 <- plot(model.pattern, class = "scatter", outcome = "effects")
#

# Summarise the CEA information from the model
summary(model.pattern)

# \donttest{
# Further examples which take longer to run
model.pattern <- pattern(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 500)
#
# Print results for all imputed values
print(model.pattern)

# Use looic to assess model fit
pic.looic <- pic(model.pattern, criterion = "looic", cases = "cc")
pic.looic

# Show density plots for mean costs parameters
diag.den <- diagnostic(model.pattern, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.pattern, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.pattern)

# }
}
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 81
#>    Unobserved stochastic nodes: 79
#>    Total graph size: 610
#> 
#> Initializing model
#> 
#> Call:
#> pattern(model.eff = e ~ trt, model.cost = c ~ trt, data = MenSS.subset, 
#>     dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, 
#>     n.iter = 100)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept).1        trt2.1 (Intercept).2        trt2.2 
#>         0.907         0.065         0.907         0.065 
#> 
#> $model.cost
#> (Intercept).1        trt2.1 (Intercept).2        trt2.2 
#>        51.971        21.337        51.971        21.337 
#> 
#> $model.pattern
#>    p1    p2 
#> 0.296 0.704 
#> 
#> $Effects
#>                Mean    SD     QL    QU
#> (Intercept).1 0.907 0.036  0.844 0.983
#> trt2.1        0.907 0.036  0.844 0.983
#> (Intercept).2 0.065 0.072 -0.071 0.195
#> trt2.2        0.065 0.072 -0.071 0.195
#> 
#> $Costs
#>                 Mean     SD      QL      QU
#> (Intercept).1 51.971 26.820   3.901 104.555
#> trt2.1        51.971 26.820   3.901 104.555
#> (Intercept).2 21.337 31.658 -33.686  84.128
#> trt2.2        21.337 31.658 -33.686  84.128
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.907 0.038 0.847 0.991
#> 2 0.971 0.063 0.849 1.103
#> 
#>  
#>  Mean costs by intervention 
#>   Mean   SD     QL  QU
#> 1 51.9 31.3 -0.692 107
#> 2 74.9 43.2  7.060 156
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45306 1898 42299 49474
#> 2 48469 3164 42374 55110
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 251
#>    Unobserved stochastic nodes: 234
#>    Total graph size: 2197
#> 
#> Initializing model
#> 
#> Call:
#> pattern(model.eff = e ~ trt, model.cost = c ~ trt + e, data = MenSS, 
#>     dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, 
#>     n.iter = 500)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept).1        trt2.1 (Intercept).2        trt2.2 
#>         0.904        -0.004         0.904        -0.004 
#> 
#> $model.cost
#> (Intercept).1        trt2.1 (Intercept).2        trt2.2           e.1 
#>        78.678        25.821        78.678        25.821        -7.551 
#>           e.2 
#>        -7.551 
#> 
#> $model.pattern
#>   p1   p2 
#> 0.29 0.71 
#> 
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under MAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.905 0.023 0.859 0.949
#> 2 0.900 0.028 0.847 0.957
#> 
#>  
#>  Mean costs by intervention 
#>    Mean   SD   QL  QU
#> 1  78.5 25.8 25.0 124
#> 2 105.6 34.6 36.5 167
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 45154 1173 42889 47345
#> 2 44917 1418 42233 47780
#
```
