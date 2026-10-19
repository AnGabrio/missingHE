# missingHE: Missing Cost Effectiveness Models

> Fit Bayesian Missing Cost Effectiveness Models in R.

🚀 This is the **development version** of the `R` package `missingHE`
(currently on version 1.6.2.1). The stable version is now release 1.6.2,
on [CRAN](https://cran.r-project.org/web/packages/missingHE/index.html).

## Contents

- [Overview](#introduction)
- [Features](#features)
- [Installation](#installation)
- [Further details](#further-details)

## Overview

Given a individual-level partially-observed cost-effectiveness outcome
data in the form of a dataframe, fits different types of Bayesian models
to handle missing effect and/or cost values under different assumptions
about the missing data mechanisms, assesses model performance and Markov
chain Monte Carlo (MCMC) convergence, and summarises statistical and
economic results based on the posterior distribution of key model
parameters.

## Features

Main features of `missingHE` include:

- Analysis of partially-observed effects/costs data via Bayesian
  missingness models
- Model assessment measures/plots and MCMC convergence checks
- Posterior summaries for key model parameters and cost-effectiveness
  measures
- Standard cost-effectiveness measures/graphs based on model results

## Installation

Install the [released version from
CRAN](https://cran.r-project.org/package=missingHE) with

``` r

install.packages("missingHE")
```

The development version (in this repo, which can be updated more quickly
and more often than the stable one) can be installed using the following
command

``` r

pkgs <- c("R2jags", "ggplot2", "gridExtra", "BCEA", "ggmcmc", "loo", "Rtools", "devtools", "utils")
repos <- c("https://cran.rstudio.com") 
install.packages(pkgs, repos = repos, dependencies = "Depends")
```

before installing the package using `devtools`:

``` r

devtools::install_github("AnGabrio/missingHE", build_vignettes = TRUE)
```

The optional argument `build_vignettes = TRUE` allows to install the
vignettes of the package locally on your computer. These consist in
brief tutorials to guide the user on how to use and customise the models
in `missingHE` using different functions of the package. Once the
package is installed, they can be accessed using the command

``` r

utils::browseVignettes(package = "missingHE")
```

which shows all the vignettes available for the package. Alternatively,
you can install the development version using `remotes`, with the
following command.

``` r

install.packages("remotes")
remotes::install_github("AnGabrio/missingHE")
```

**NB**: On Windows machines, you need to install a few dependencies,
including [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
first, e.g. by running

``` r

pkgs <- c("MASS", "Rtools", "remotes")
repos <- "https://cran.rstudio.com"
install.packages(pkgs, repos = repos, dependencies = "Depends")
```

before installing the package using `remotes`.

All models implemented in `missingHE` are written in the `BUGS` language
using the software [JAGS](https://mcmc-jags.sourceforge.io/), which
needs to be installed from its own repository and instructions for
installations under different OS can be found online. Once installed,
the software is called in `missingHE` via the `R` package
[R2jags](https://cran.r-project.org/package=R2jags). Note that the
`missingHE` package is currently under active development and therefore
it is advisable to re-install the package directly from GitHub before
each use to ensure that you are using the most updated version.

## Further details

More details on `missingHE` are available in our
[*vignettes*](https://cran.r-project.org/web/packages/missingHE/vignettes/).
Also, details about the package, including some references and links to
a pdf presentation and some posts on my own blog) are given
[here](https://angabrio.github.io/software/missingHE/).

## License

[![License: GPL
v2](https://img.shields.io/badge/License-GPLv2-blue.svg)](https://www.gnu.org/licenses/gpl-2.0)

## Contributing

Please submit contributions through `Pull Requests`, following the
[contributing
guidelines](https://github.com/AnGabrio/missingHE/blob/master/.github/CONTRIBUTING.md).
To report issues and/or seek support, please file a new ticket in the
[issue](https://github.com/AnGabrio/missingHE/issues) tracker.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/AnGabrio/missingHE/blob/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.
