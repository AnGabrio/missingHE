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

## References

Daniels, MJ. Hogan, JW. *Missing Data in Longitudinal Studies:
strategies for Bayesian modelling and sensitivity analysis*, CRC/Chapman
Hall.

Baio, G.(2012). *Bayesian Methods in Health Economics*. CRC/Chapman
Hall, London.

Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). *Bayesian Data
Analysis, 2nd edition*, CRC Press.

Plummer, M. *JAGS: A program for analysis of Bayesian graphical models
using Gibbs sampling*. (2003).

## See also

[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`BCEA::bcea()`](https://n8thangreen.github.io/BCEA/reference/bcea.html).

## Author

Andrea Gabrio

## Examples

``` r
# Quick example to run using subset of MenSS dataset
MenSS.subset <- MenSS[50:100, ]

# Run the model using the pattern function assuming a MAR mechanism
# Use only 100 iterations to run a quick check
model.pattern <- pattern(data = MenSS.subset,model.eff = e ~ trt, model.cost = c ~ trt,
   dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 100)
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

# Print the results of the JAGS model
print(model.pattern)
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
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.pattern, criterion = "dic", cases = "cc")
pic.dic
#> $d_bar
#> [1] 211.8549
#> 
#> $pD
#> [1] 2.242907
#> 
#> $dic
#> [1] 214.0978
#> 
#> $d_hat
#> [1] 209.612
#> 
#

# Extract regression coefficient estimates
coef(model.pattern)
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
#


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

# \donttest{
# Further examples which take longer to run
model.pattern <- pattern(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 500)
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
#
# Print results for all imputed values
print(model.pattern)
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

# Use looic to assess model fit
pic.looic <- pic(model.pattern, criterion = "looic", cases = "cc")
pic.looic
#> $elpd
#> [1] -345.7407
#> 
#> $elpd_se
#> [1] -345.7407
#> 
#> $p
#> [1] 6.954153
#> 
#> $p_se
#> [1] 6.954153
#> 
#> $looic
#> [1] 691.4814
#> 
#> $looic_se
#> [1] 691.4814
#> 
#> $pointwise
#>         elpd_loo mcse_elpd_loo      p_loo    looic influence_pareto_k
#>  [1,]  -6.562130   0.008883579 0.03789003 13.12426        0.057764523
#>  [2,]  -6.604481   0.009111540 0.03976173 13.20896        0.010944416
#>  [3,] -11.160392   0.045313396 0.55037839 22.32078        0.186127403
#>  [4,]  -7.492706   0.010195133 0.05021975 14.98541        0.139698648
#>  [5,]  -6.534033   0.009422550 0.04246772 13.06807        0.016189853
#>  [6,]  -6.658334   0.009944525 0.04701667 13.31667        0.054581421
#>  [7,] -11.036120   0.051216539 0.63006576 22.07224        0.224855575
#>  [8,] -10.850288   0.042966503 0.54664663 21.70058        0.309763882
#>  [9,]  -6.885140   0.011514868 0.06021282 13.77028        0.082990864
#> [10,]  -6.569059   0.008861222 0.03754017 13.13812        0.098319729
#> [11,]  -6.900890   0.011114368 0.05623951 13.80178        0.098868382
#> [12,]  -6.955621   0.011446708 0.05958277 13.91124        0.047021358
#> [13,]  -6.547546   0.008784920 0.03703576 13.09509        0.104572628
#> [14,]  -6.851502   0.011387952 0.05891673 13.70300        0.101056853
#> [15,] -15.619908   0.128293259 2.12375872 31.23982        0.606514298
#> [16,]  -6.898305   0.011121371 0.05631058 13.79661        0.105902246
#> [17,]  -6.900890   0.011114368 0.05623951 13.80178        0.098868382
#> [18,]  -6.733632   0.009029423 0.03914134 13.46726        0.160237198
#> [19,]  -6.525899   0.009260123 0.04095375 13.05180        0.148548505
#> [20,]  -7.698479   0.011509749 0.06175517 15.39696        0.065623687
#> [21,]  -7.767926   0.019345348 0.14565378 15.53585       -0.014891944
#> [22,]  -6.572497   0.008940215 0.03835831 13.14499        0.044254523
#> [23,]  -8.204797   0.016897784 0.11852279 16.40959        0.089581109
#> [24,]  -6.955790   0.009585076 0.04425113 13.91158        0.122910648
#> [25,]  -6.900890   0.011114368 0.05623951 13.80178        0.098868382
#> [26,]  -6.855987   0.011491211 0.05994607 13.71197        0.115741959
#> [27,]  -6.887719   0.011260862 0.05717628 13.77544        0.123840214
#> [28,]  -6.719417   0.010252886 0.05177913 13.43883       -0.146252819
#> [29,]  -6.984449   0.012045713 0.06816386 13.96890       -0.003754172
#> [30,]  -6.891898   0.012038936 0.06918361 13.78380       -0.087976603
#> [31,]  -7.956636   0.014994685 0.09473947 15.91327        0.051681034
#> [32,]  -6.984449   0.012045713 0.06816386 13.96890       -0.003754172
#> [33,]  -8.308023   0.027227884 0.25630601 16.61605        0.170142164
#> [34,]  -6.598976   0.010074687 0.04565546 13.19795        0.163550850
#> [35,]  -7.039610   0.014694096 0.08232715 14.07922        0.285501688
#> [36,]  -6.731107   0.011028501 0.05255016 13.46221        0.131546175
#> [37,]  -8.618053   0.029502308 0.31288641 17.23611        0.125684225
#> [38,]  -6.984449   0.012045713 0.06816386 13.96890       -0.003754172
#> [39,]  -6.924948   0.012066606 0.06959529 13.84990       -0.018803394
#> [40,]  -7.096302   0.011624160 0.05920762 14.19260        0.294243178
#> [41,]  -6.984449   0.012045713 0.06816386 13.96890       -0.003754172
#> [42,]  -7.137420   0.010345735 0.05122120 14.27484        0.030399824
#> [43,]  -7.518433   0.012653948 0.07381119 15.03687        0.095790160
#> [44,]  -7.109135   0.012062805 0.06911781 14.21827        0.115814566
#> [45,]  -6.799422   0.012671442 0.06406920 13.59884        0.228746903
#> [46,]  -8.222555   0.021706640 0.17676663 16.44511        0.202698956
#> 
#> $diagnostics
#> $diagnostics$pareto_k
#>  [1]  0.057764523  0.010944416  0.186127403  0.139698648  0.016189853
#>  [6]  0.054581421  0.224855575  0.309763882  0.082990864  0.098319729
#> [11]  0.098868382  0.047021358  0.104572628  0.101056853  0.606514298
#> [16]  0.105902246  0.098868382  0.160237198  0.148548505  0.065623687
#> [21] -0.014891944  0.044254523  0.089581109  0.122910648  0.098868382
#> [26]  0.115741959  0.123840214 -0.146252819 -0.003754172 -0.087976603
#> [31]  0.051681034 -0.003754172  0.170142164  0.163550850  0.285501688
#> [36]  0.131546175  0.125684225 -0.003754172 -0.018803394  0.294243178
#> [41] -0.003754172  0.030399824  0.095790160  0.115814566  0.228746903
#> [46]  0.202698956
#> 
#> $diagnostics$n_eff
#>  [1] 480.6897 479.7533 245.0729 475.0820 478.2149 476.1087 215.5102 259.7523
#>  [9] 468.3735 480.6936 470.3682 468.6498 481.0953 469.0416  53.5543 470.3540
#> [17] 470.3682 480.0010 479.0488 468.5203 420.6516 480.4606 436.9203 477.7177
#> [25] 470.3682 468.5577 469.5187 474.7300 465.4212 465.5353 448.5352 465.4212
#> [33] 363.8053 475.2287 450.9608 470.4830 347.5013 465.4212 465.4135 468.0859
#> [41] 465.4212 474.1681 462.0778 465.5023 462.2457 404.0810
#> 
#> $diagnostics$r_eff
#>  [1] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
#> [39] 1 1 1 1 1 1 1 1
#> 
#> 
#> $psis_object
#> NULL
#> 

# Show density plots for mean costs parameters
diag.den <- diagnostic(model.pattern, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.pattern, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.pattern)
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

# }
#
#
```
