# Fitting MNAR models in missingHE

When analysing partially-observed data, any statistical method makes
either explicit or implicit assumptions about the missing values that
can never be verified from the data at hand. Often, analyses rely on a
*missing at random* (MAR) assumption, that is they assume that the
observed data are able to fully explain missingness. However, there is
always the chance that this assumption is not correct and that
missingness also depends on some values which are not observed, leading
to a *missing not at random* (MNAR) assumption. Thus, it is extremely
important that the robustness of the results to a range of alternative
missingness assumptions is assessed in sensitivity analysis, including
MNAR.

Each missingness model approach in `missingHE`, implemented either
through the function `selection`, `pattern`, or `hurdle`, can be fitted
under MNAR for either or both the effectiveness and cost outcomes. MNAR
assumptions are expressed in terms of some suitably-defined departures
from MAR in terms of mean effectiveness and cost parameters. This
tutorial illustrates how to specify a MNAR model in `missingHE` using
the built-in data set `MenSS`, which is directly available when loading
`missingHE` in the user’s `R` workspace. See [Introduction to
missingHE](https://angabrio.github.io/missingHE/articles/Introduction_to_missingHE.md)
for an introductory tutorial on each type of missingness approach in
`missingHE` and for a general presentation of the `MenSS` dataset.

If you would like to have more information on the package, or would like
to point out potential issues in the current version, feel free to
contact the maintainer at <a.gabrio@maastrichtuniversity.nl>.
Suggestions on how to improve the package are also very welcome. A
development version of the package is also available from the author’s
main GitHub repository at <https://github.com/AnGabrio/missingHE>.

## Handling MNAR using selection models

Selection models encode MNAR assumptions by directly modelling the
missingness mechanism, i.e., the model for the missing data indicator,
as a function of some partially-observed variables. In `missingHE`, MNAR
is specified by including the outcome variables (\\e\\ and/or \\c\\)
inside the model formulae of the corresponding missing indicators
(`model.me` and `model.mc`). For example, a model encoding a MNAR
assumption for the missing effectiveness mechanism and a MAR assumption
for the missing cost mechanism can be specified by setting

``` r

model.me <- me ~  e + ...
model.mc <- mc ~  1 + ...
```

where `me` and `mc` indicate the effectiveness and cost missingness
indicator variables, while `...` is a suitable form for the combination
of covariates to be included in each model. Inclusion of `e` into
`model.me` instructs `missingHE` to include the partially-observed
effectiveness variable into the logistic regression model for its
missingness indicator. Since the regression coefficient capturing the
dependence between missingness and outcome (\\\delta_e\\) cannot be
fully estimated from the observed data, its estimation implies some form
of MNAR assumption for the missing effectiveness mechanism. In selection
models, two sources of external information are used to identify these
parameters: **informative prior distributions** and **distributional
outcome assumptions**.

1.  Informative priors on \\\delta_e\\ (or \\\delta_c\\ for costs) can
    be specified through an `R` list object containing information about
    the desired prior specification, including name of the parametric
    distribution assumed (among a set of pre-defined choices) and its
    corresponding prior values. This list can then be passed to the
    `selection` function using the optional argument `prior` to
    overwrite the default priors (standard normals) assumed by
    `missingHE` for these parameters.

2.  The distribution assumed for each outcome can be selected among a
    pre-defined set of parametric choices (among popular choices in the
    literature) using the arguments `dist_e` and `dist_c` inside the
    `selection` function. Type
    [`help(selection)`](https://angabrio.github.io/missingHE/reference/selection.md)
    to access the full list of available distributions for \\e\\ and
    \\c\\.

As an example, let us fit a MNAR selection model for the effectiveness
variables (QALYs) in the `MenSS` dataset using the following command

``` text
> NN.sel1 <- selection(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, model.me = me ~ e, 
+          model.mc = mc ~ 1, type = "MNAR", 
+          n.iter = 1000, ref = 2, 
+          dist_e = "norm", dist_c = "norm")
```

The MNAR assumption is encoded into the model by: setting the argument
`type = "MNAR"`; adding the term `e` inside `model.me`. As for the
outcome model, both effectiveness and costs are assumed to follow normal
distributions, with the inclusion of the baseline utilities (`u.0`) as
covariates into the model of \\e\\. Since we did not provide custom
priors, default standard normals are specified for \\\delta_e\\.

`missingHE` allows a separate specification of the assumptions on the
effectiveness and cost missingness mechanisms through the arguments
`model.me` and `model.mc`, respectively. When desired, other
fully-observed variables can also be included into the missingness
models (either under MAR or MNAR assumptions), to improve the estimation
of the model parameters. Posterior summaries for the estimates of model
parameters under the selected missingness assumptions can then be
retrieved through the `postsum` function.

``` text
> postsum(NN.sel1)
#>                 mean          sd        2.5%         50%       97.5% Rhat n.eff
#> alpha[1]   0.5343462  0.06618759   0.3947543   0.5358181   0.6615649 1.02   100
#> alpha[2]   0.0370997  0.02562776  -0.0093526   0.0365951   0.0901265 1.00   460
#> alpha[3]   0.4004031  0.07009024   0.2661087   0.3995501   0.5489282 1.02    73
#> beta[1]   79.7270676 26.36415284  29.0872210  79.9310792 131.3966846 1.02  1000
#> beta[2]   24.5820423 28.39262444 -36.5701080  24.7092960  79.4482736 1.00  1000
#> delta_e    0.1331422  0.92386355  -1.7175333   0.1502429   1.9427385 1.00  1000
#> gamma_c    0.9107466  0.17214437   0.5746326   0.8988882   1.2509972 1.00  1000
#> gamma_e    0.7754424  0.84945278  -0.8406131   0.7559512   2.4364242 1.00  1000
#> p_c        0.7118731  0.03501391   0.6398314   0.7107210   0.7774724 1.00  1000
#> p_e        0.7099567  0.03643674   0.6365759   0.7110560   0.7774530 1.00  1000
#> s_c      249.6306598 27.60703667 203.4883185 248.3160497 308.4007179 1.01   200
#> s_e        0.0908020  0.00975253   0.0748604   0.0893822   0.1123818 1.00  1000
#> tau_c      0.0000166  0.00000361   0.0000105   0.0000162   0.0000242 1.01   200
#> tau_e    125.3649638 25.96512298  79.1786937 125.1693715 178.4415759 1.00  1000
#> tmu_c     92.7138070 27.56033024  40.0931795  91.2815899 145.7692647 1.00  1000
#> tmu_e      0.9070680  0.01488000   0.8772706   0.9080101   0.9368123 1.01   110
```

We now consider an alternative MNAR specification where we encode
informative priors on \\\delta_e\\. In general, it is difficult to
attach a specific interpretation to the values for this parameter (and
therefore also to its priors) because its effect may vary depending on
the distributional assumptions made on the outcome. In this applied
example we will therefore only illustrate how some already-defined
suitable prior information can be passed into the model, rather than
focussing on the plausibility of the chosen values within the current
context.

We first define our priors by creating a list object called `my.prior`.
Within this list, we create a vector called `"delta.prior.e"` containing
the following elements: name of the prior distribution (`"norm"`); the
prior values for the selected distribution (mean of \\10\\ and precision
of \\1\\).

``` text
> my.prior <- list(
+   "delta.prior.e" = c("norm", 10, 1)
+ )
```

As a simple exercise, we increase the prior mean of \\\delta_e\\ to
\\10\\ to assess the impact on posterior estimates of a quite
informative prior about this parameter (relatively high positive value
on the logit scale). We then fit the model using the new prior
specification by passing `my.prior` to the argument `prior` inside
`selection`.

``` text
> NN.sel2 <- selection(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, model.me = me ~ e, 
+          model.mc = mc ~ 1, type = "MNAR", 
+          n.iter = 1000, ref = 2, prior = my.prior,
+          dist_e = "norm", dist_c = "norm")
```

We can now check the results

``` text
> postsum(NN.sel2)
#>                 mean          sd         2.5%         50%       97.5% Rhat n.eff
#> alpha[1]   0.6821685  0.08350955   0.52236927   0.6783965   0.8578424 1.11    19
#> alpha[2]   0.0514863  0.02740033  -0.00240917   0.0514828   0.1036527 1.00  1000
#> alpha[3]   0.3073438  0.07972403   0.14671111   0.3073416   0.4649993 1.06    35
#> beta[1]   81.3228550 26.51949175  27.82937557  81.0823258 131.7749880 1.01   310
#> beta[2]   25.8025032 28.30527275 -31.11770668  25.0214884  83.0286349 1.00   850
#> delta_e    9.5793668  1.04353773   7.56348347   9.6113175  11.6736802 1.01   110
#> gamma_c    0.9089760  0.18424669   0.56268410   0.9045462   1.2742280 1.02  1000
#> gamma_e   -8.2711381  1.10000005 -10.44726165  -8.2765462  -6.2225967 1.03    51
#> p_c        0.7113300  0.03758449   0.63707336   0.7118828   0.7814657 1.01  1000
#> p_e        0.7534214  0.03873643   0.67589126   0.7537424   0.8298371 1.01   110
#> s_c      249.6378266 31.98874072 199.33091960 246.6704797 321.6339274 1.01   560
#> s_e        0.1100927  0.01772366   0.08275996   0.1081749   0.1514907 1.24    11
#> tau_c      0.0000168  0.00000398   0.00000967   0.0000164   0.0000252 1.01   560
#> tau_e     88.5235463 26.34317865  43.57436355  85.4569055 146.0022191 1.24    11
#> tmu_c     94.9543661 28.03196161  38.78898758  95.7817159 146.5885065 1.00   780
#> tmu_e      0.9804202  0.02692018   0.93365932   0.9765932   1.0426883 1.21    13
```

and compare them to those obtained under the initial MNAR model. We see
that, with respect to the results from `NN.sel1`, the mean effectiveness
estimates across both groups (`tmu_e`) are on average higher. Average
mean estimates for each treatment group and associated incremental
results for each model can also be inspected by using the `summary`
function. To have a better idea of the impact in terms of
cost-effectiveness conclusions, we can use the function `ceac.plot`,
from the package `BCEA`, to display the cost-effectiveness acceptability
curves based on the results from each model.

![Cost-effectiveness acceptability curve (CEAC) derived from model 1
fitted under default MNAR prior specification for the missing
effectiveness
mechanism.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarsel1-1.png)

Cost-effectiveness acceptability curve (CEAC) derived from model 1
fitted under default MNAR prior specification for the missing
effectiveness mechanism.

![Cost-effectiveness acceptability curve (CEAC) derived from model 2
fitted under custom MNAR prior specification for the missing
effectiveness
mechanism.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarsel2-1.png)

Cost-effectiveness acceptability curve (CEAC) derived from model 2
fitted under custom MNAR prior specification for the missing
effectiveness mechanism.

The comparison between Figure @ref(fig:figplotceamnarsel1) and Figure
@ref(fig:figplotceamnarsel2) shows that CEA conclusions are somewhat
affected by the specific assumptions made about the missingness
mechanism. It is very important that the scenarios explored are informed
based on some external information (e.g., expert opinion) to provide a
range of *plausible* assumptions to assess in sensitivity analysis.

## Handling MNAR using pattern mixture models

Pattern mixture models encode MNAR assumptions through the specification
of two elements: **identifying restrictions** and **sensitivity
parameters**. Identifying restrictions consist in modelling restrictions
imposed to identify those parameters in each missingness pattern that
cannot otherwise be identified from the data. In the context of a
bivariate cross-sectional outcome (\\e,c\\), a total of four missingness
patterns can be observed based on the values of their missingness
indicators \\d^m=(m_i^c,m^e_i)\\. These are: complete cases \\(0,0)\\;
missing costs only \\(1,0)\\; missing effects only \\(0,1)\\; fully
missing cases \\(1,1)\\. A typical example of these restrictions is when
unidentifiable parameters in a given pattern (e.g., fully missing cases)
are identified by matching those from a corresponding identifiable
pattern (e.g., complete cases). Under MAR, these restrictions are the
only element used to achieve the identification of the model. However,
when MNAR assumptions are specified, identifying restrictions can be
combined with additional parameters \\\boldsymbol \delta\\ to identify
the model under a range of departures from MAR. These so-called
“sensitivity parameters” are entirely identified using informative
priors to encode in the model evidence that is external to the observed
data. Unlike the priors for MNAR parameters in selection models,
sensitivity parameters in pattern mixture models can have more intuitive
and easier-to-assess interpretations in terms of their impact on
posterior results.

`missingHE` allows the specification of MNAR assumptions through the
arguments `restriction` and `type` inside the function `pattern`. The
first argument captures the type of identifying restrictions imposed:
either complete case (`"CC"`) or available case (`"AC"`) restrictions.
The second argument specifies the type of mechanism assumed for the
missing values: either `MAR` or `MNAR`. Under MAR, the model is
identified only based on the given identifying restrictions imposed in
each unidentified pattern for both outcomes. For example, when `"CC"`
restrictions are selected, the mean effectiveness and cost parameters in
the fully missing cases \\\mu_e^{(1,1)}\\ and \\\mu_c^{(1,1)}\\ are set
equal to those from the complete cases \\\mu_e^{(0,0)}\\ and
\\\mu_c^{(0,0)}\\, respectively. Under MNAR, location-shift sensitivity
parameters (\\\delta_e\\, \\\delta_c\\) are included within the chosen
identification strategy. Using the same example, under MNAR, the mean
outcomes in the fully missing cases are identified as:

\\\begin{equation} \mu_e^{(1,1)}=\mu_e^{(0,0)}+\delta_e \\\\\\
\text{and} \\\\\\ \mu_c^{(1,1)}=\mu_c^{(0,0)}+\delta_c, (\\eq:pmsp)
\end{equation}\\

where \\\delta_e\\ and \\\delta_c\\ are by default assigned with
minimally informative uniform priors between \\0\\ ans \\1\\ by
`missingHE` (both set to \\0\\ under MAR). Sensitivity analysis to
alternative MNAR scenarios can then be conducted by:

1.  Choose alternative identifying restrictions to be passed to the
    argument `restriction`, i.e., either `"CC"` or `"AC"`.

2.  Specify informative priors on \\\delta_e\\ and/or \\\delta_c\\
    through an `R` list object containing information about the desired
    prior specification, including name of the parametric distribution
    assumed (among a set of pre-defined choices) and its corresponding
    prior values. This list can then be passed to the `pattern` function
    using the optional argument `prior` to overwrite the default priors
    (minimally informative uniforms) assumed by `missingHE` for these
    parameters under MNAR.

Going back to our example, based on external information, let us assume
that the mean effectiveness among the missing cases is on average
between \\0.2\\ and \\0.3\\ lower compared to that among the observed
cases, across both treatment arms. This information can be encoded into
a prior distribution for \\\delta_e\\ within a pattern mixture model in
a similar way to what shown in Section @ref(sec-sm) for selection
models. The following code shows how to specify a uniform distribution
(`"unif"`), with prior boundary values of \\0.2\\ and \\0.3\\, for
\\\delta_e\\.

``` text
> my.prior <- list(
+   "delta.prior.e" = c("unif", -0.3, -0.2)
+ )
```

Next, we fit the desired MNAR model using the function `pattern`,
setting the arguments `restriction = "CC"`, `type = "MNAR"` and
overwriting the default priors for \\\delta_e\\ with those contained in
the object `my.prior`.

``` text
> NN.pat2 <- pattern(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, type = "MNAR", 
+          restriction = "CC", n.iter = 1000, 
+          ref = 2, prior = my.prior,
+          dist_e = "norm", dist_c = "norm")
```

Baseline utilities are included in the model for \\e\\, which is
identified under MNAR using complete case restrictions
(`restriction = "CC"`) and informative priors on the sensitivity
parameters for the mean \\e\\. Economic results in terms of posterior
summaries about mean parameters can be inspected using the `postsum`
function

``` text
> postsum(NN.pat2)
#>                     mean          sd       2.5%         50%       97.5% Rhat n.eff
#> alpha_p[1,1]   0.5516640  0.07156196   0.411827   0.5512801   0.6937966 1.00  1000
#> alpha_p[2,1]   0.0321027  0.02813567  -0.022630   0.0320161   0.0905064 1.00  1000
#> alpha_p[3,1]   0.3823995  0.07645248   0.232283   0.3833580   0.5321849 1.00  1000
#> alpha_p[1,2]   0.5516640  0.07156196   0.411827   0.5512801   0.6937966 1.00  1000
#> alpha_p[2,2]   0.0321027  0.02813567  -0.022630   0.0320161   0.0905064 1.00  1000
#> alpha_p[3,2]   0.3823995  0.07645248   0.232283   0.3833580   0.5321849 1.00  1000
#> delta_e       -0.2502551  0.02852923  -0.296478  -0.2498697  -0.2022185 1.00  1000
#> p_prob[1]      0.2908691  0.03571644   0.223282   0.2906947   0.3649819 1.00  1000
#> p_prob[2]      0.7091309  0.03571644   0.635018   0.7093053   0.7767178 1.00  1000
#> s_c_p[1]     247.0366688 25.94525370 199.343643 245.2197013 301.2716930 1.00  1000
#> s_c_p[2]     247.0366688 25.94525370 199.343643 245.2197013 301.2716930 1.00  1000
#> s_e_p[1]       0.0925849  0.01108866   0.073703   0.0910968   0.1189721 1.01   190
#> s_e_p[2]       0.0925849  0.01108866   0.073703   0.0910968   0.1189721 1.01   190
#> tau_c_p[1]     0.0000169  0.00000356   0.000011   0.0000166   0.0000252 1.00  1000
#> tau_c_p[2]     0.0000169  0.00000356   0.000011   0.0000166   0.0000252 1.00  1000
#> tau_e_p[1]   121.4790256 27.70956494  70.649586 120.5018446 184.0898861 1.01   190
#> tau_e_p[2]   121.4790256 27.70956494  70.649586 120.5018446 184.0898861 1.01   190
#> tmu_c         95.3260010 27.51067788  43.760212  95.1896652 150.9973769 1.00   720
#> tmu_e          0.7284119  0.02539561   0.679252   0.7288158   0.7746406 1.00  1000
```

For comparison, we also fit the same model under an alternative MNAR
specification where we assume a normal distribution for the sensitivity
parameters, roughly defined on the same value range,

``` text
> my.prior2 <- list(
+   "delta.prior.e" = c("norm", -0.25, 1/(0.05^2))
+ )
```

which is then passed to the argument `prior` in the `pattern` function

``` text
> NN.pat3 <- pattern(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, type = "MNAR", 
+          restriction = "CC", n.iter = 1000, 
+          ref = 2, prior = my.prior2,
+          dist_e = "norm", dist_c = "norm")
```

We assess the impact on the cost-effectiveness results by looking at the
acceptability curves derived from each MNAR model using again the
function `ceac.plot`.

![Cost-effectiveness acceptability curves (CEAC) derived from model 1
fitted under uniform MNAR priors for the effectiveness sensitivity
parameter.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarpat1-1.png)

Cost-effectiveness acceptability curves (CEAC) derived from model 1
fitted under uniform MNAR priors for the effectiveness sensitivity
parameter.

![Cost-effectiveness acceptability curves (CEAC) derived from model 2
fitted under normal MNAR priors for the effectiveness sensitivity
parameter.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarpat2-1.png)

Cost-effectiveness acceptability curves (CEAC) derived from model 2
fitted under normal MNAR priors for the effectiveness sensitivity
parameter.

Figure @ref(fig:figplotceamnarpat1) and Figure
@ref(fig:figplotceamnarpat2) indicate a slightly lower chance for the
new intervention to be cost-effective when comparing the results under
an informative normal distribution for \\\delta_e\\ (`NN.pat3`) with
respect to a uniform distribution (`NN.pat3`). Note that, in real
applications, the range of values and choice of distribution for the
sensitivity parameters should be informed based on some external source
of information (e.g., expert opinion) to guide the elicitation process.

## Handling MNAR using hurdle models

Hurdle models consist of two-part mixture models, where the outcome
variable (i.e., \\e\\, \\c\\ or both) are split into two components: one
associated with identical values (called “structural”) which are assumed
to be fixed; and one associated with a set of values (called
non-structural) which are assumed to be generated from a given
distribution. `missingHE` models the non-structural outcome values using
some parametric distributions and uses logistic regressions to estimate
the probability of incurring in a structural value. The function
`hurdle` requires the user to specify separate formulae for: the models
estimating the mean outcome for the non-structural component through the
usual arguments `model.eff` and `model.cost`; the logistic regressions
estimating the structural value probabilities using the arguments
`model.se` and `model.sc`:

``` r

model.se <- se ~  1 + ...
model.sc <- sc ~  1 + ...
```

where `se` and `sc` indicate the effectiveness and cost structural value
indicator variables, while `...` is a suitable form for the combination
of covariates in each logistic regression model. Adopting a nomenclature
similar to the missingness mechanism, it is possible to argue that
exclusion of any covariates from these models corresponds to a
“Structural Completely At Random” (SCAR) assumption, since the
probability of a structural value does not depend on any observed
variable. Similarly, inclusion of some covariates into these models
suggests a “Structural At Random” (SAR) assumption, where the structural
probability can only be correctly estimated conditional on some other
observed variables. In contrast to selection or pattern mixture models,
hurdle models do not provide an explicit framework for encoding
informative missingness assumptions given that the missingness process
is not directly modelled. However, they can be used to assess the impact
of specific MNAR assumptions about the number of missing individuals
assigned with a structural value.

Let us first consider a standard hurdle model specification under MAR.
We specify the model to handle both structural ones and zeros in \\e\\
and \\c\\ from the `MenSS` dataset, where structural values are set by
using the arguments `se = 1` and `sc = 0` inside `hurdle`, respectively.
Since these values are assumed to be “fixed”, their occurrence does not
affect the choice of a distribution for the non-structural component of
the outcomes. As a result, we can assign distributions that are not
defined at \\1\\ and \\0\\ to the non-structural outcome values, such as
Beta and Gamma distributions for \\e\\ and \\c\\, respectively.

``` text
> NN.hur1 <- hurdle(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, model.se = se ~ 1, 
+          model.sc = sc ~ 1, type = "SCAR", 
+          se = 1, sc = 0, n.iter = 1000, ref = 2, 
+          dist_e = "beta", dist_c = "gamma")
```

The model assumes that the structural values do not depend on any
observed covariate, i.e., a *structural completely at random* (SCAR)
mechanism. The function `hurdle` automatically assigns all individuals
with an observed one effect and zero cost to the structural components
of the effectiveness and cost mixture distributions, respectively. In
general, we do not know to which component of the mixture the
individuals with a missing value should be assigned, as this information
cannot be obtained from the data. However, based on some external
information, we can arbitrarily impose this assignment (which
effectively corresponds to a MNAR mechanism).

In `missingHE`, this can be achieved by first creating an indicator
variable (named `s_e`), telling for each individual whether a structural
value is observed (`s_e = 1`), not observed (`s_e = 0`) or unknown
(`s_e = NA`). For example, regarding the effectiveness variables, we can
generate this indicator by typing

``` text
> s_e <- ifelse(MenSS$e == 1, 1, 0)
> #number of ones
> sum(s_e == 1, na.rm = T)
#> [1] 17
```

Next, for all or some of the individuals with a missing effect value, we
set the value of `s_e = 1` to assign them to the structural component of
the model. For example, we may believe that it is likely for all
individuals aged \\\<22\\ to be associated with a perfect health state
(i.e., \\e=1\\). We can encode this information into a new object called
`mys_e` by typing

``` text
> mys_e <- ifelse(is.na(s_e) & MenSS$age < 22, 1, s_e)
> #number of ones
> sum(mys_e == 1, na.rm = T)
#> [1] 41
```

The number of individuals associated with \\e=1\\ has considerably
increased with respect to the observed data alone. We can now proceed to
fit the hurdle model by replacing the original structural indicator
variable of \\e\\ with the newly created one using the argument
`s_e = mys_e`.

``` text
> NN.hur2 <- hurdle(data = MenSS, model.eff = e ~ trt + u.0, 
+          model.cost = c ~ trt, model.se = se ~ 1, 
+          model.sc = sc ~ 1, type = "SCAR", s_e = mys_e,
+          se = 1, sc = 0, n.iter = 1000, ref = 2, 
+          dist_e = "beta", dist_c = "gamma")
```

We can inspect the posterior results by typing

``` text
> postsum(NN.hur2)
#>                 mean      sd      2.5%       50%     97.5% Rhat n.eff
#> alpha[1,1]   0.67818  0.3607  -0.04056   0.67473   1.46014 1.10    22
#> alpha[2,1]   0.10168  0.2393  -0.36037   0.09951   0.57370 1.01   130
#> alpha[3,1]   1.22235  0.4152   0.27378   1.25664   1.93649 1.11    21
#> alpha[1,2]  11.51292  0.0000  11.51292  11.51292  11.51292 1.00     1
#> alpha[2,2]   0.00000  0.0000   0.00000   0.00000   0.00000 1.00     1
#> alpha[3,2]   0.00000  0.0000   0.00000   0.00000   0.00000 1.00     1
#> beta[1,1]    5.45748  0.2090   5.08356   5.44604   5.91384 1.03   140
#> beta[2,1]    0.31837  0.1942  -0.08734   0.32592   0.67780 1.02   120
#> beta[1,2]  -11.51293  0.0000 -11.51293 -11.51293 -11.51293 1.00     1
#> beta[2,2]    0.00000  0.0000   0.00000   0.00000   0.00000 1.00     1
#> gamma_c     -1.05437  0.3311  -1.70637  -1.05165  -0.40194 1.00  1000
#> gamma_e      0.35238  0.2522  -0.13022   0.34226   0.86058 1.00  1000
#> p_c          0.26330  0.0633   0.15363   0.25891   0.40085 1.00  1000
#> p_e          0.58587  0.0601   0.46749   0.58474   0.70278 1.00  1000
#> s_c[1]     244.88136 56.9736 169.61731 234.38722 403.65725 1.02   300
#> s_c[2]       0.00001  0.0000   0.00001   0.00001   0.00001 1.00     1
#> s_e[1]       0.09585  0.0163   0.07149   0.09326   0.13380 1.01   300
#> s_e[2]       0.00001  0.0000   0.00001   0.00001   0.00001 1.00     1
#> tmu_c      207.18242 39.8025 146.22813 202.45663 300.01750 1.02   760
#> tmu_e        0.94127  0.0120   0.91318   0.94250   0.96143 1.01   360
```

and look at how imputations in each treatment group are carried out
based using the `plot` function.

    #> Loading required namespace: ggthemes

![Scatter plots of observed and imputed effectiveness data in the SoC
arm (trt = 1) based on an hurdle model assuming all individuals younger
than 22 years old are associated with a perfect healf
state.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplothm2-1.png)

Scatter plots of observed and imputed effectiveness data in the SoC arm
(trt = 1) based on an hurdle model assuming all individuals younger than
22 years old are associated with a perfect healf state.

As Figure @ref(fig:figplothm2) shows, the model imputes each missing
structural effect in `"SoC"` as one with a very small credible interval.
Finally, we compare the economic results from the two hurdle models
using the `ceac.plot` function from the `BCEA` package.

![Cost-effectiveness acceptability curves (CEAC) derived from the hurdle
model fitted under a MAR assumption about the number of structural ones
in
SoC.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarhm1-1.png)

Cost-effectiveness acceptability curves (CEAC) derived from the hurdle
model fitted under a MAR assumption about the number of structural ones
in SoC.

![Cost-effectiveness acceptability curves (CEAC) derived from the hurdle
models fitted under a MNAR assumption about the number of structural
ones in
SoC.](Fitting_MNAR_models_in_missingHE_files/figure-html/figplotceamnarhm2-1.png)

Cost-effectiveness acceptability curves (CEAC) derived from the hurdle
models fitted under a MNAR assumption about the number of structural
ones in SoC.

According to Figure @ref(fig:figplotceamnarhm1) and Figure
@ref(fig:figplotceamnarhm2), the probability of cost-effectiveness for
the hurdle model under MAR (`NN.hur1`) shows a moderate increase
alongside the willingness to pay threshold up to \\0.7\\, where a
plateau is observed. Conversely, results from the MNAR model (`NN.hur2`)
indicate a systematically lower chance of cost-effectiveness, especially
for lower willingness to pay values.
