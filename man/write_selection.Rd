% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_selection.R
\name{write_selection}
\alias{write_selection}
\title{An internal function to select which type of selection model to execute. 
Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
type of missingness mechanism assumed and independence or joint modelling
This function selects which type of model to execute.}
\usage{
write_selection(
  dist_e,
  dist_c,
  type,
  pe_fixed,
  pc_fixed,
  ze_fixed,
  zc_fixed,
  ind_fixed,
  pe_random,
  pc_random,
  ze_random,
  zc_random,
  ind_random,
  model_e_random,
  model_c_random,
  model_me_random,
  model_mc_random
)
}
\arguments{
\item{dist_e}{Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')}

\item{dist_c}{Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')}

\item{type}{Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)}

\item{pe_fixed}{Number of fixed effects for the effectiveness model}

\item{pc_fixed}{Number of fixed effects for the cost model}

\item{ze_fixed}{Number of fixed effects or the missingness indicators model for the effectiveness}

\item{zc_fixed}{Number of fixed effects or the missingness indicators model for the costs}

\item{ind_fixed}{Logical; if TRUE independence between effectiveness and costs is assumed, else correlation is accounted for}

\item{pe_random}{Number of random effects for the effectiveness model}

\item{pc_random}{Number of random effects for the cost model}

\item{ze_random}{Number of random effects or the missingness indicators model for the effectiveness}

\item{zc_random}{Number of random effects or the missingness indicators model for the costs}

\item{ind_random}{Logical; if TRUE independence at the level of the random effects between effectiveness and costs is assumed, else correlation is accounted for}

\item{model_e_random}{Random effects formula for the effectiveness model}

\item{model_c_random}{Random effects formula for the costs model}

\item{model_me_random}{Random effects formula for the missingness indicators model for the effectiveness}

\item{model_mc_random}{Random effects formula for the missingness indicators model for the costs}
}
\description{
An internal function to select which type of selection model to execute. 
Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
type of missingness mechanism assumed and independence or joint modelling
This function selects which type of model to execute.
}
\examples{
#Internal function only
#No examples
#
#
}
\keyword{JAGS}
\keyword{Selection}
\keyword{models}
