# A function to summarise results from a fitted model

This function hides missing data distribution from summary results of
`BUGS` models

## Usage

``` r
jagsresults(
  x,
  params,
  regex = FALSE,
  invert = FALSE,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  signif,
  ...
)
```

## Arguments

- x:

  The `rjags`, `rjags.parallel`, or `mcmc.list` object for which results
  will be printed.

- params:

  Character vector or a regular expression pattern. The parameters for
  which results will be printed (unless `invert` is `FALSE`, in which
  case results for all parameters other than those given in `params`
  will be returned). If `regex` is `FALSE`, only those parameters that
  match `params` exactly will be returned. If `regex` is `TRUE`, `param`
  should be a character string giving the regular expression pattern to
  be matched.

- regex:

  If `regex` is `TRUE`, then `param` is expected to be a single string
  giving a text pattern to be matched. Parameters with names matching
  the pattern will be returned (unless `invert` is `TRUE`, which results
  in all parameters that do not match the pattern being returned). Text
  pattern matching uses regular expressions (`regex`).

- invert:

  Logical. If `invert` is `TRUE`, only those parameters that do not
  match elements of `params` will be returned.

- probs:

  A numeric vector of probabilities, representing the sample quantiles
  to be calculated and returned.

- signif:

  If supplied, all columns other than `n.eff` will have their values
  rounded such that the most extreme value has the specified number of
  significant digits.

- ...:

  Additional arguments accepted by `grep`, e.g. `perl = TRUE`, to allow
  look-around pattern matching.

## Examples

``` r
if (FALSE) { # \dontrun{
## Data
N <- 100
temp <- runif(N)
rain <- runif(N)
wind <- runif(N)
a <- 0.13
beta.temp <- 1.3
beta.rain <- 0.86
beta.wind <- -0.44
sd <- 0.16
y <- rnorm(N, a + beta.temp*temp + beta.rain*rain + beta.wind*wind, sd)
dat <- list(N=N, temp=temp, rain=rain, wind=wind, y=y)

### bugs example
library(R2jags)

## Model
M <- function() {
  for (i in 1:N) {
    y[i] ~ dnorm(y.hat[i], sd^-2)
    y.hat[i] <- a + beta.temp*temp[i] + beta.rain*rain[i] + beta.wind*wind[i]
    resid[i] <- y[i] - y.hat[i]
  }
  sd ~ dunif(0, 100)
  a ~ dnorm(0, 0.0001)
  beta.temp ~ dnorm(0, 0.0001)
  beta.rain ~ dnorm(0, 0.0001)
  beta.wind ~ dnorm(0, 0.0001)
}

## Fit model
jagsfit <- jags(dat, inits=NULL,
                parameters.to.save=c('a', 'beta.temp', 'beta.rain',
                                     'beta.wind', 'sd', 'resid'),
                model.file=M, n.iter=10000)

## Output
# model summary
jagsfit

# Results for beta.rain only
jagsresults(x=jagsfit, param='beta.rain')

# Results for 'a' and 'sd' only
jagsresults(x=jagsfit, param=c('a', 'sd'))
jagsresults(x=jagsfit, param=c('a', 'sd'),
            probs=c(0.01, 0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))

# Results for all parameters including the string 'beta'
jagsresults(x=jagsfit, param='beta', regex=TRUE)

# Results for all parameters not including the string 'beta'
jagsresults(x=jagsfit, param='beta', regex=TRUE, invert=TRUE)

# Note that the above is NOT equivalent to the following, which returns all
# parameters that are not EXACTLY equal to 'beta'.
jagsresults(x=jagsfit, param='beta', invert=TRUE)

# Results for all parameters beginning with 'b' or including 'sd'.
jagsresults(x=jagsfit, param='^b|sd', regex=TRUE)

# Results for all parameters not beginning with 'beta'.
# This is equivalent to using param='^beta' with invert=TRUE and regex=TRUE
jagsresults(x=jagsfit, param='^(?!beta)', regex=TRUE, perl=TRUE)
} # }
#
#
```
