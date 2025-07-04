% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pic.R
\name{pic}
\alias{pic}
\title{Predictive information criteria for Bayesian models fitted in \code{JAGS} using the funciton \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}}
\usage{
pic(x, criterion = "dic", module = "total")
}
\arguments{
\item{x}{A \code{missingHE} object containing the results of a Bayesian model fitted in cost-effectiveness analysis using the function \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}.}

\item{criterion}{type of information criteria to be produced. Available choices are \code{'dic'} for the Deviance Information Criterion, 
\code{'waic'} for the Widely Applicable Information Criterion, and \code{'looic'} for the Leave-One-Out Information Criterion.}

\item{module}{The modules with respect to which the information criteria should be computed. Available choices are \code{'total'} for the whole model, 
\code{'e'} for the effectiveness variables only, \code{'c'} for the cost variables only, and \code{'both'} for both outcome variables.}
}
\value{
A named list containing different predictive information criteria results and quantities according to the value of \code{criterion}. In all cases, the measures are 
computed on the observed data for the specific modules of the model selected in \code{module}.
\describe{
  \item{d_bar}{Posterior mean deviance (only if \code{criterion} is \code{'dic'}).}
  \item{pD}{Effective number of parameters calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}.
  \item{dic}{Deviance Information Criterion calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}. 
  \item{d_hat}{Deviance evaluated at the posterior mean of the parameters and calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}
  \item{elpd, elpd_se}{Expected log pointwise predictive density and standard error calculated on the observed data for the model nodes indicated in \code{module}
   (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
  \item{p, p_se}{Effective number of parameters and standard error calculated on the observed data for the model nodes indicated in \code{module}
   (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
  \item{looic, looic_se}{The leave-one-out information criterion and standard error calculated on the observed data for the model nodes indicated in \code{module}
   (only if \code{criterion} is \code{'loo'}).}
  \item{waic, waic_se}{The widely applicable information criterion and standard error calculated on the observed data for the model nodes indicated in \code{module}
   (only if \code{criterion} is \code{'waic'}).}
  \item{pointwise}{A matrix containing the pointwise contributions of each of the above measures calculated on the observed data for the model nodes indicated in \code{module}
   (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
  \item{pareto_k}{A vector containing the estimates of the shape parameter \eqn{k} for the generalised Pareto fit to the importance ratios for each leave-one-out distribution 
   calculated on the observed data for the model nodes indicated in \code{module} (only if \code{criterion} is \code{'loo'}). 
   See \code{\link[loo]{loo}} for details about interpreting \eqn{k}.}
}
}
\description{
Efficient approximate leave-one-out cross validation (LOO), deviance information criterion (DIC) and widely applicable information criterion (WAIC) for Bayesian models, 
calculated on the observed data.
}
\details{
The Deviance Information Criterion (DIC), Leave-One-Out Information Criterion (LOOIC) and the Widely Applicable Information Criterion (WAIC) are methods for estimating 
out-of-sample predictive accuracy from a Bayesian model using the log-likelihood evaluated at the posterior simulations of the parameters. 
DIC is computationally simple to calculate but it is known to have some problems, arising in part from it not being fully Bayesian in that it is based on a point estimate.
LOOIC can be computationally expensive but can be easily approximated using importance weights that are smoothed by fitting a generalised Pareto distribution to the upper tail 
of the distribution of the importance weights.
WAIC is fully Bayesian and closely approximates Bayesian cross-validation. Unlike DIC, WAIC is invariant to parameterisation and also works for singular models. 
In finite cases, WAIC and LOO give similar esitmates, but for influential observations WAIC underestimates the effect of leaving out one observation.
}
\examples{
 
#For examples see the function selection, pattern or hurdle 
# 
# 
}
\references{
Plummer, M. \emph{JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling.} (2003).

Vehtari, A. Gelman, A. Gabry, J. (2016a) Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. 
\emph{Statistics and Computing}. Advance online publication.

Vehtari, A. Gelman, A. Gabry, J. (2016b) Pareto smoothed importance sampling. \emph{ArXiv} preprint.

Gelman, A. Hwang, J. Vehtari, A. (2014) Understanding predictive information criteria for Bayesian models. 
\emph{Statistics and Computing} 24, 997-1016. 

Watanable, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information 
criterion in singular learning theory. \emph{Journal of Machine Learning Research} 11, 3571-3594.
}
\seealso{
\code{\link[R2jags]{jags}}, \code{\link[loo]{loo}}, \code{\link[loo]{waic}}
}
\author{
Andrea Gabrio
}
\keyword{JAGS}
\keyword{dic}
\keyword{loo}
\keyword{waic}
