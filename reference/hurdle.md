# Bayesian Hurdle Models

Full Bayesian cost-effectiveness models to handle missing data in the
outcomes using Hurdle models under a variety of alternative parametric
distributions for the effect and cost variables. Alternative assumptions
about the mechanisms of the structural values are implemented using a
hurdle approach. The analysis is performed using the `BUGS` language,
which is implemented in the software `JAGS` using the function
[`R2jags::jags()`](https://rdrr.io/pkg/R2jags/man/jags.html). The output
is stored in an object of class `missingHE`.

## Usage

``` r
hurdle(
  data,
  model.eff,
  model.cost,
  model.se = se ~ 1,
  model.sc = sc ~ 1,
  se = 1,
  sc = 0,
  dist_e,
  dist_c,
  type,
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
  `model.cost` (model formulas for effects and costs) and `model.se`,
  `model.sc` (model formulas for the structural effect and cost models).
  Among these, effectiveness, cost and treatment indicator variables
  must always be provided and named `e`, `c` and `trt`, respectively.

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

- model.se:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `se` (structural effects).
  Any covariates in the model must be provided on the right-hand side of
  the formula. If there are no covariates, `1` should be specified on
  the right hand side of the formula. By default, covariates are placed
  on the "probability" parameter for the structural effects through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- model.sc:

  A formula expression in conventional `R` linear modelling syntax. The
  response must be indicated with the term `sc` (structural costs). Any
  covariates in the model must be provided on the right-hand side of the
  formula. If there are no covariates, `1` should be specified on the
  right hand side of the formula. By default, covariates are placed on
  the "probability" parameter for the structural costs through a
  logistic-linear model. Random effects can also be specified for each
  model parameter. See details for how these can be specified.

- se:

  Structural value to be found in the effect variables defined in
  `data`. If set to `NULL`, no structural value is chosen and a standard
  model for the effects is run.

- sc:

  Structural value to be found in the cost variables defined in `data`.
  If set to `NULL`, no structural value is chosen and a standard model
  for the costs is run.

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
  Completely At Random (`SCAR`), and Structural At Random (`SAR`).

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

  Logical. If `save_model` is `TRUE` a `txt` file containing the model
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

An object of the class `missingHE` containing the following elements:

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

- **se_prior**: A numeric value that indicates the value assumed for the
  standard deviation associated with the degenerate structural effects
  prior distribution to implement the model in `JAGS` (deafult =
  0.00001).

- **sc_prior**: A numeric value that indicates the value assumed for the
  standard deviation associated with the degenerate structural costs
  prior distribution to implement the model in `JAGS` (deafult =
  0.00001).

- **call**: The matched call.

## Details

Depending on the distributions specified for the outcome variables in
the arguments `dist_e` and `dist_c` and the type of structural value
mechanism specified in the argument `type`, different hurdle models are
built and run in the background by the function `hurdle`. These are
mixture models defined by two components: the first one is a mass
distribution at the spike, while the second is a parametric model
applied to the natural range of the relevant variable. Usually, a
logistic regression is used to estimate the probability of incurring a
"structural" value (e.g. `0` for the costs, or `1` for the effects);
this is then used to weight the mean of the "non-structural" values
estimated in the second component. A simple example can be used to show
how hurdle models are specified. Consider a data set comprising a
response variable \\y\\ and a set of centered covariates \\X_j\\, for
\\i = j, ..., J\\.Specifically, for each subject in the trial \\i = 1,
..., n\\ we define an indicator variable \\s_i\\ taking value `1` if the
\\i\\-th individual is associated with a structural value and `0`
otherwise. This is modelled as:

\$\$s_i ~ Bernoulli(\pi_i)\$\$ \$\$logit(\pi_i) = \sum\gamma_j X_j\$\$
where

- \\\pi_i\\ is the individual probability of a structural value in
  \\y\\.

- \\\gamma_j\\ represents the impact on the probability of a structural
  value in \\y\\ of the covariate \\X_j\\.

When \\\gamma_j = 0\\, the model assumes a `SCAR` mechanism, while when
\\\gamma_j != 0\\ the mechanism is `SAR`. For the parameters indexing
the structural value model, the default prior distributions assumed are:

- \\\gamma_j ~ Normal(0, 0.01)\\.

When user-defined hyperprior values are supplied via the argument
`prior` in the function `hurdle`, the elements of this list (see
Arguments) must be vectors containing the user-provided distribution
name and hyperprior values and must take specific names according to the
parameters they are associated with. Specifically, the names accepted by
**missingHE** are the following:

- auxiliary parameters \\\sigma\\: `sigma.prior.e` (effects) and/or
  `sigma.prior.c` (costs).

- covariate parameters \\\alpha_j\\ and \\\beta_j\\: `alpha.prior`
  (effects) and/or `beta.prior` (costs).

- covariate parameters in the structural values model \\\gamma_j\\ (if
  covariate data provided): `gamma.prior.e` (effects) and/or
  `gamma.prior.c` (costs).

For simplicity, here we assumed that the set of covariates \\X_j\\ used
in the models for the effects/costs and in the model of the structural
effect/cost values is the same. However, it is possible to specify
different sets of covariates for each model using the arguments in the
function `hurdle` (see Arguments).

For each model, random effects can also be specified for each parameter
by adding the term `+ (x | z)` to each model formula, where `x` is the
fixed regression coefficient for which also the random effects are
desired and `z` is the clustering variable across which the random
effects are specified (must be the name of a factor variable in the
dataset).

## References

Ntzoufras I. (2009). *Bayesian Modelling Using WinBUGS*, John Wiley and
Sons.

Daniels, MJ. Hogan, JW. (2008). *Missing Data in Longitudinal Studies:
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
# Quick example to run using the MenSS dataset

# Run the model using the hurdle function assuming a SCAR mechanism
# Use subset of dataset to produce quick results
MenSS_ex <- MenSS[MenSS$age > 25, ]

# Use only 50 iterations to run a quick check
model.hurdle <- hurdle(data = MenSS_ex, model.eff = e ~ trt, model.cost = c ~ trt,
   model.se = se ~ 1, model.sc = sc ~ 1, se = NULL, sc = 0, dist_e = "norm", dist_c = "gamma",
   type = "SCAR", n.chains = 2, n.iter = 50)
#> module glm loaded
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 87
#>    Unobserved stochastic nodes: 220
#>    Total graph size: 1884
#> 
#> Initializing model
#> 
#> Registered S3 method overwritten by 'BCEA':
#>   method     from
#>   plot.evppi voi 

# Print the results of the JAGS model
print(model.hurdle)
#> Call:
#> hurdle(model.eff = e ~ trt, model.cost = c ~ trt, model.se = se ~ 
#>     1, model.sc = sc ~ 1, data = MenSS_ex, se = NULL, sc = 0, 
#>     dist_e = "norm", dist_c = "gamma", type = "SCAR", n.chains = 2, 
#>     n.iter = 50)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       0.885      -0.004 
#> 
#> $model.cost
#> (Intercept)        trt2 
#>       6.835      -0.100 
#> 
#> $model.se
#> [1] NA
#> 
#> $model.sc
#> (Intercept) 
#>      -0.945 
#> 
#

# Use dic information criterion to assess model fit
pic.dic <- pic(model.hurdle, criterion = "dic", cases = "cc")
pic.dic
#> $d_bar
#> [1] 129.5131
#> 
#> $pD
#> [1] 113.3485
#> 
#> $dic
#> [1] 242.8616
#> 
#> $d_hat
#> [1] 16.16456
#> 
#

# Extract regression coefficient estimates
coef(model.hurdle)
#> $Effects
#>               Mean    SD     QL    QU
#> (Intercept)  0.885 0.036  0.822 0.950
#> trt2        -0.004 0.043 -0.064 0.085
#> 
#> $Costs
#>               Mean    SD     QL    QU
#> (Intercept)  6.835 0.980  5.567 8.065
#> trt2        -0.100 0.237 -0.463 0.326
#> 
#


# Assess model convergence using graphical tools
# Produce histograms of the posterior samples for the mean effects
diag.hist <- diagnostic(model.hurdle, type = "histogram", param = "mu.e")
#

# Compare observed effect data with imputations from the model
# using plots (posterior means and credible intervals)
p1 <- plot(model.hurdle, class = "scatter", outcome = "effects")
#

# Summarise the CEA information from the model
summary(model.hurdle)
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under SCAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.885 0.034 0.816 0.966
#> 2 0.883 0.038 0.833 0.969
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD    QL   QU
#> 1  712 582 123.6 1739
#> 2  655 561  92.9 1809
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 43530 1609 40577 46986
#> 2 43510 1972 39934 47496

# \donttest{
# Further examples which take longer to run
model.hurdle <- hurdle(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
   model.se = se ~ age, model.sc = sc ~ age, se = 1, sc = 0, dist_e = "beta", dist_c = "gamma",
   type = "SAR", n.chains = 2, n.iter = 500)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 184
#>    Unobserved stochastic nodes: 463
#>    Total graph size: 6241
#> 
#> Initializing model
#> 
#
# Print results for all imputed values
print(model.hurdle)
#> Call:
#> hurdle(model.eff = e ~ trt, model.cost = c ~ trt + e, model.se = se ~ 
#>     age, model.sc = sc ~ age, data = MenSS, se = 1, sc = 0, dist_e = "beta", 
#>     dist_c = "gamma", type = "SAR", n.chains = 2, n.iter = 500)
#> 
#> Coefficients:
#> $model.eff
#> (Intercept)        trt2 
#>       1.779      -0.166 
#> 
#> $model.cost
#> (Intercept)        trt2           e 
#>       6.788       0.093      -0.583 
#> 
#> $model.se
#> (Intercept)         age 
#>       0.207      -0.026 
#> 
#> $model.sc
#> (Intercept)         age 
#>      -1.200       0.004 
#> 

# Use looic to assess model fit
pic.looic <- pic(model.hurdle, criterion = "looic", cases = "cc")
pic.looic
#> $elpd
#> [1] 25.34684
#> 
#> $elpd_se
#> [1] 25.34684
#> 
#> $p
#> [1] 16.87966
#> 
#> $p_se
#> [1] 16.87966
#> 
#> $looic
#> [1] -50.69369
#> 
#> $looic_se
#> [1] -50.69369
#> 
#> $pointwise
#>        elpd_loo mcse_elpd_loo     p_loo      looic influence_pareto_k
#>  [1,]  9.993008    0.01965004 0.1571902 -19.986016      -0.0094787762
#>  [2,] 10.022891    0.01873503 0.1468215 -20.045782      -0.1314426952
#>  [3,] -9.602025    0.03808458 0.5090225  19.204050       0.3475638514
#>  [4,] -7.245874    0.02265437 0.2303490  14.491748       0.1047152584
#>  [5,] -6.059156    0.02142212 0.2865836  12.118313      -0.2072356735
#>  [6,] -6.325934    0.02634010 0.3777673  12.651868      -0.0589286925
#>  [7,] -9.746028    0.04271015 0.5841980  19.492056       0.3436233911
#>  [8,] -7.995870    0.02344037 0.1681116  15.991739       0.4847931750
#>  [9,]  2.280038    0.02414323 0.3736633  -4.560077      -0.0084560518
#> [10,] -5.394239    0.09960927 1.3188364  10.788478       0.4405804357
#> [11,] 18.594957    0.01861416 0.1251534 -37.189915       0.3371369119
#> [12,]  1.997125    0.02219916 0.3092328  -3.994250       0.1875701629
#> [13,] 10.123251    0.01493537 0.1058958 -20.246501       0.1353121031
#> [14,]  2.829577    0.02233003 0.3020330  -5.659155      -0.0236012606
#> [15,] -9.221282    0.06678752 0.5604582  18.442564       0.5563061611
#> [16,]  4.265830    0.09095245 0.6108660  -8.531659       0.6735394136
#> [17,] 17.823217    0.05496899 0.6685305 -35.646434       0.3226176336
#> [18,] -6.719177    0.02230125 0.3068660  13.438354      -0.0698248243
#> [19,] -6.327481    0.02150575 0.2747684  12.654962       0.0003563834
#> [20,]  1.321836    0.01867530 0.1878799  -2.643672       0.2921721457
#> [21,] -7.455084    0.02955410 0.4006794  14.910169       0.2726255407
#> [22,]  9.967619    0.02112852 0.1772704 -19.935239      -0.0246415407
#> [23,] -8.031869    0.03193744 0.4921050  16.063738       0.2478689311
#> [24,] -6.875947    0.02446344 0.3011887  13.751893       0.0200758658
#> [25,] 18.588772    0.01745666 0.1163794 -37.177545       0.3891124817
#> [26,]  2.462092    0.02249810 0.3222907  -4.924184       0.1325672201
#> [27,] -6.734996    0.02936329 0.4146636  13.469992       0.4005480339
#> [28,] -7.404790    0.07712162 0.8783735  14.809580       0.5371811660
#> [29,] 18.556525    0.01564120 0.1061400 -37.113049       0.1321666454
#> [30,]  2.418572    0.02227680 0.3138075  -4.837145       0.1769719469
#> [31,] -7.459966    0.02684732 0.3097709  14.919931       0.0702051395
#> [32,] 18.305526    0.02715843 0.2472626 -36.611052       0.1488625232
#> [33,]  8.818521    0.02191011 0.1922750 -17.637042       0.2025403884
#> [34,] -6.577861    0.02412779 0.3669596  13.155723      -0.1901684449
#> [35,] -6.707055    0.03962710 0.4067749  13.414109       0.3172510936
#> [36,] -6.755315    0.02521062 0.3950378  13.510631      -0.0071318038
#> [37,] -8.088564    0.03338656 0.6199554  16.177128       0.2087177345
#> [38,] 18.482658    0.01758479 0.1298495 -36.965316       0.0581563649
#> [39,]  2.143892    0.02433593 0.3955351  -4.287784       0.1564820099
#> [40,] -6.988775    0.02620126 0.4393613  13.977550       0.0060505310
#> [41,] 18.602364    0.02133105 0.1486695 -37.204728       0.3688809022
#> [42,] -7.284587    0.03960828 0.5948392  14.569175       0.1935831954
#> [43,]  1.574302    0.02500345 0.3892601  -3.148603       0.1210323198
#> [44,]  1.554940    0.02490583 0.3677816  -3.109881       0.3111400836
#> [45,] -6.562733    0.02185917 0.2520500  13.125466       0.1524205580
#> [46,] -7.816062    0.02902394 0.4971514  15.632124      -0.2100493468
#> 
#> $diagnostics
#> $diagnostics$pareto_k
#>  [1] -0.0094787762 -0.1314426952  0.3475638514  0.1047152584 -0.2072356735
#>  [6] -0.0589286925  0.3436233911  0.4847931750 -0.0084560518  0.4405804357
#> [11]  0.3371369119  0.1875701629  0.1353121031 -0.0236012606  0.5563061611
#> [16]  0.6735394136  0.3226176336 -0.0698248243  0.0003563834  0.2921721457
#> [21]  0.2726255407 -0.0246415407  0.2478689311  0.0200758658  0.3891124817
#> [26]  0.1325672201  0.4005480339  0.5371811660  0.1321666454  0.1769719469
#> [31]  0.0702051395  0.1488625232  0.2025403884 -0.1901684449  0.3172510936
#> [36] -0.0071318038  0.2087177345  0.0581563649  0.1564820099  0.0060505310
#> [41]  0.3688809022  0.1935831954  0.1210323198  0.3111400836  0.1524205580
#> [46] -0.2100493468
#> 
#> $diagnostics$n_eff
#>  [1] 417.7371 424.1623 289.1924 396.8054 406.0460 370.1090 260.5894 392.7801
#>  [9] 386.5842  83.6962 425.7338 400.5298 449.2248 399.8266 155.3483  97.6956
#> [17] 199.0434 399.7262 405.2365 425.3642 346.8815 407.2992 329.6778 384.1689
#> [25] 433.6354 398.5721 348.8596 125.2564 445.0245 400.1980 366.1037 364.1468
#> [33] 402.1030 386.3566 279.2174 378.9655 319.6261 431.8167 385.1329 371.3441
#> [41] 407.5884 278.2101 379.9593 381.8756 403.2223 350.3440
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
diag.den <- diagnostic(model.hurdle, type = "denplot", param = "mu.c")

# Plots of imputations for all data
p1 <- plot(model.hurdle, class = "scatter", outcome = "both")

# Summarise the CEA results
summary(model.hurdle)
#> 
#>  Cost-effectiveness analysis summary 
#>  
#>  CE parameter estimates under SAR assumption
#>  
#>  Mean effects by intervention 
#>    Mean    SD    QL    QU
#> 1 0.940 0.014 0.911 0.965
#> 2 0.934 0.018 0.896 0.964
#> 
#>  
#>  Mean costs by intervention 
#>   Mean  SD    QL   QU
#> 1  731 532  87.3 1675
#> 2  761 546 110.1 1765
#> 
#>  
#>  Mean net monetary benefit by intervention and wtp = 50000 
#>    Mean   SD    QL    QU
#> 1 46254  923 44378 47969
#> 2 45926 1111 43440 47739

# }
#
#
```
