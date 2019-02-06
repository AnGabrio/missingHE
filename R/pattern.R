#' Full Bayesian Models to handle missingness in Economic Evaluations (Pattern Mixture Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes under different missing data 
#' mechanism assumptions, using a variatey of alternative parametric distributions for the effect and cost variables and 
#' using a pattern mixture approach to identify the model. The analysis is performed using the \code{BUGS} language, 
#' which is implemented in \code{JAGS} using the functions \code{\link[R2jags]{jags}} The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs).
#' Among these, effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#' effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#' any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#' By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#' cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#' If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#' parameter of the distribution through a linear model. A joint bivariate distribution for effects and costs can be specified by
#' including 'e' in the model for the costs.
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param Delta_e Range of values for the prior on the sensitivity parameters used to identify the mean of the effects under MNAR. The value must be set to 0 under MAR. 
#' @param Delta_c Range of values for the prior on the sensitivity parameters used to identify the mean of the costs under MNAR. The value must be set to 0 under MAR.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#' CI sample quantiles to be calculated and returned for the imputed values.
#' @param n.chains Number of chains.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of warmup iterations.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{JAGS}
#' will generate initial values for all the model parameters.
#' @param n.thin Thinning interval.
#' @param save_model Logical. If \code{save_model} is \code{TRUE} a \code{txt} file containing the model code is printed
#' in the current working directory.
#' @param prior A list containing the hyperprior values provided by the user. Each element of this list must be a vector of length two
#' containing the user-provided hyperprior values and must be named with the name of the corresponding parameter. For example, the hyperprior
#' values for the mean effect parameter can be provided using \code{prior = list('mu.prior.e' = c(0, 1))}.
#' For more information about how to provide prior hypervalues for different types of parameters and models see details. 
#' If \code{prior} is set to 'default', the default values will be used.  
#' @param ... additional input parameters provided by the user. Examples are \code{center = TRUE} to center all the covariates in the model 
#' or the additional arguments that can be provided to the function \code{\link[BCEA]{bcea}} to summarise the health economic evaluation results. 
#' @return An object of the class 'missingHE' containing the following elements
#' \describe{
#'   \item{data_set}{A list containing the original data set provided in \code{data} (see Arguments), the number of observed and missing individuals 
#'   , the total number of individuals by treatment arm and the indicator vectors for the missing values}
#'   \item{model_output}{A list containing the output of a \code{JAGS} model generated from the functions \code{\link[R2jags]{jags}}, and 
#'   the posterior samples for the main parameters of the model and the imputed values}
#'   \item{cea}{A list containing the output of the economic evaluation performed using the function \code{\link[BCEA]{bcea}}}
#'   \item{type}{A character variable that indicate which type of missingness assumption has been used to run the model, 
#'   either \code{MAR} or \code{MNAR} (see details)}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA JAGS missing data Pattern Mixture Models
#' @importFrom stats model.frame 
#' @details Depending on the distributional assumptions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} and the type of missingness mechanism assumed in the argument \code{type}, different types of pattern mixture models
#' are built and run in the background by the function \code{pattern}. The model for the outcomes is fitted in each missingness pattern 
#' and the parameters indexing the missing data distributions are identified using: the corresponding parameters identified from the observed data
#' in other patterns (under 'MAR'); or a combination of the parameters identified by the observed data and some sensitivity parameters (under 'MNAR'). 
#' A simple example can be used to show how Pattern mixture models are specified. 
#' Consider a data set comprising a response variable \eqn{y} and a set of centered covariate \eqn{X_j}. We denote with \eqn{d_i} the patterns' indicator variable for each 
#' subject in the trial \eqn{i = 1, ..., n} such that: \eqn{d_i = 1} indicates the completers (both e and c observed), \eqn{d_i = 2} and \eqn{d_i = 3} indicate that 
#' only the costs or effects are observed, respectively, while \eqn{d_i = 4} indicates that neither of the two outcomes is observed. In general, a different number of patterns 
#' can be observed between the treatment groups and \code{missingHE} accounts for this possibility by modelling a different patterns' indicator variables for each arm. 
#' For simplicity, in this example, we assume that the same number of patterns is observed in both groups. \eqn{d_i} is assigned a multinomial distribution, 
#' which probabilities are modelled using a Dirichlet prior (by default giving to each pattern the same weight). Next, the model specified in \code{dist_e} 
#' and \code{dist_c} is fitted in each pattern. The parameters that cannot be identified by the observed data in each pattern (d = 2, 3, 4), e.g. the means 
#' \eqn{mu_e[d]} and \code{mu_c[d]}, are then identified as: 
#' \deqn{mu_e[2] = \mu_e[4] = \mu_e[1] + \Delta_e}
#' \deqn{mu_c[3] = \mu_c[4] = \mu_c[1] + \Delta_c}
#' where
#' \itemize{
#' \item \eqn{\mu_e[1]} is the effects mean for the completers.
#' \item \eqn{\mu_c[1]} is the costs mean for the completers.
#' \item \eqn{\Delta_e} is the sensitivity parameters associated with the marginal effects mean.
#' \item \eqn{\Delta_c} is the sensitivity parameters associated with the marginal costs mean.
#' }
#' When \eqn{\Delta_e = 0} and \eqn{\Delta_c = 0} the model assumes a 'MAR' mechanism. When \eqn{\Delta_e != 0} and/or \eqn{\Delta_c != 0} 'MNAR' departues for the 
#' effects and/or costs are explored assuming a Uniform prior distributions for the sensitivity parameters. The range of values for these priors is defined based on the
#' boundaries specified in \code{Delta_e} and \code{Delta_c} (see Arguments), which must be provided by the user. 
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{pattern}, the elements of this list (see Arguments)
#' must be vectors of length \code{2} containing the user-provided hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names for the parameters indexing the outcome model accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item location parameters \eqn{\alpha_0} and \eqn{\beta_0}: "mean.prior.e"(effects) and/or "mean.prior.c"(costs)
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j} and \eqn{\beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' } 
#' The only exception is the missingness patterns' probability \eqn{\pi}, denoted with "patterns.prior", which hyperprior values must be provided as a list
#' formed by two elements. These must be vectors of the same length as the number of patterns in the control (first element) and intervention (second element) group.
#' 
#' @author Andrea Gabrio
#' @references  
#' Daniels, MJ. Hogan, JW. \emph{Missing Data in Longitudinal Studies: strategies for Bayesian modelling and sensitivity analysis}, CRC/Chapman Hall.
#' 
#' Baio, G.(2012). \emph{Bayesian Methods in Health Economics}. CRC/Chapman Hall, London.
#' 
#' Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). \emph{Bayesian Data Analysis, 2nd edition}, CRC Press.
#' 
#' Plummer, M. \emph{JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling.} (2003).
#' @export
#'
#' @examples
#' \dontrun{
#'#Simple example to simulate and analyse a data set
#'#Define the number of individuals fer treatment arm
#'N1 <- 150
#'N2 <- 100
#'
#'#Create the missingness indicators totally random (MCAR mechanism)
#'m_eff1 <- m_cost1 <- rbinom(N1,1,0.25)
#'m_eff2 <- m_cost2 <- rbinom(N2,1,0.25)
#'
#'#Simulate data from normal distributions for both arms
#'eff1 <- rnorm(N1, 0.5, 0.5)
#'eff2 <- rnorm(N2, 0.5, 0.5)
#'cost1 <- rnorm(N1, 90, 20)
#'cost2 <- rnorm(N2, 90, 20)
#'
#'#Set value missing if indicator is 1
#'eff1[m_eff1 == 1] <- NA
#'eff2[m_eff2 == 1] <- NA
#'cost1[m_cost1 == 1] <- NA
#'cost2[m_cost2 == 1] <- NA
#'
#'#Create treatment arm indicators
#'t1 <- rep(1, length(eff1))
#'t2 <- rep(2, length(eff2))
#'
#'#Combine variables and define a data set
#'e <- c(eff1, eff2)
#'c <- c(cost1, cost2)
#'m_eff <- c(m_eff1, m_eff2)
#'m_cost <- c(m_cost1, m_cost2)
#'t <- c(t1, t2)
#'data <- data.frame(e, c, t)
#'
#'#Run the model using pattern with JAGS assuming a MCAR missingness mechanism
#'x <- pattern(data = data, model.eff = e ~ 1, model.cost = c ~ 1,
#'Delta_e = 0, Delta_c = 0, dist_e = "norm", dist_c = "norm", type = "MAR", center = FALSE)
#'
#'#print the results of the JAGS model
#'print(x)
#'#
#'
#'#use information criteria to assess model fit
#'pic <- pic(x, criterion = "dic", module = "total")
#'#
#'
#'#Assess model convergence using graphical tools
#'#Produce histograms of the posterior samples for the mean effect
#'#parameters in the two treatment arms. 
#'diagnostic(x, type = "histogram", param = "mu.e")
#'
#'#
#'#Compare observed outcome data with imputations from the model
#'# (posteiror means and credible intervals)
#'plot(x, class = "scatter", outcome = "all")
#'
#'#
#'#Summarise the CEA information from model results
#'summary(x)
#' }
#'#
#'#


pattern <- function(data, model.eff, model.cost, dist_e, dist_c, Delta_e, Delta_c, type, prob = c(0.05, 0.95), n.chains = 2, n.iter = 20000, 
                    n.burnin = floor(n.iter / 2), inits = NULL, n.thin = 1, save_model = FALSE, prior = "default", ...) {
  filein <- NULL
  if(is.data.frame(data) == FALSE) {
    stop("data must be in data frame format")
  }
  if(!all(c("e", "c", "t") %in% names(data)) == TRUE) {
    stop("Please rename or provide variables in the data as 'e', 'c' and 't' for the effectiveness, cost and treatment indicator")
  }
  if(any(names(data) == "e") == TRUE & any(names(data) == "c") == TRUE) {
    e <- as.name("e")
    c <- as.name("c")
  }
  if(is.numeric(data$e) == FALSE | is.numeric(data$c) == FALSE) {
    stop("Effectiveness and cost data must be numeric")
  }
  cov_matrix <- subset(data, select = -c(e, c))
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  if(any(is.na(cov_matrix)) == TRUE) {
    stop("no missing covariate or treatment indicator is allowed")
  }
  if(!all(levels(as.factor(cov_matrix$t)) %in% c("1", "2")) == TRUE) {
    stop("A two arm indicator variable must be provided with '1' for the control and '2' for the other intervention")
  }
  if(is.character(type) == FALSE | is.character(dist_e) == FALSE | is.character(dist_c) == FALSE) {
    stop("you must provide character names for the objects 'type', 'dist_e' and 'dist_c'")
  }
  if(is.numeric(Delta_e) == FALSE | is.numeric(Delta_c) == FALSE) {
    stop("Delta parameters values or ranges must be numeric")
  }
  dist_e <- tolower(dist_e)
  dist_c <- tolower(dist_c)
  if(dist_e == "normal") { dist_e <- "norm" }
  if(dist_c == "normal") { dist_c <- "norm" }
  if(dist_c == "lognormal") { dist_c <- "lnorm" }
  if(!dist_e %in% c("norm", "beta") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta' for the effects and 'norm', 'gamma', 'lnorm' for the costs")
  }
  type <- toupper(type)
  if(!type %in% c("MAR", "MNAR")) {
    stop("Types available for use are 'MAR' and 'MNAR'")
  }
  if(length(prob) != 2 | is.numeric(prob) == FALSE | any(prob < 0) != FALSE | any(prob > 1) != FALSE) {
    stop("You must provide valid lower/upper quantiles for the imputed data distribution")
  }
  if(is.logical(save_model) == FALSE) {
    stop("save_model should be either TRUE or FALSE")
  }
  exArgs <- list(...)
  if(exists("center", where = exArgs)) {
    if(is.logical(exArgs$center) == FALSE) { stop("center must be either TRUE or FALSE") }
    center = exArgs$center 
  } else {center = FALSE }
  data_read <- data_read_pattern(data = data, model.eff = model.eff, model.cost = model.cost, type = type, center = center)
  if(is.vector(Delta_e) == FALSE & is.matrix(Delta_e) == FALSE) {
    stop("Delta parameters values or ranges must be provided as vectors or matrices")
  }
  if(is.vector(Delta_c) == FALSE & is.matrix(Delta_c) == FALSE) {
    stop("Delta parameters values or ranges must be provided as vectors or matrices")
  }
  if(type == "MAR") {
    if(is.vector(Delta_e) == FALSE | length(Delta_e) != 1 | Delta_e != 0 | is.vector(Delta_c) == FALSE | length(Delta_c) != 1 | Delta_c != 0){
      stop("Under MAR the Delta parameters for both outcomes must be set to zero")
   }
  }
   if(type == "MNAR") {
    if(is.matrix(Delta_e) == TRUE) {
      if(dim(Delta_e)[1] != 2 | dim(Delta_e)[2] != 2) { stop("Delta parameters ranges must be provided as 2x2 matrices")}
      if(Delta_e[1, 1] > Delta_e[1, 2] | Delta_e[2, 1] > Delta_e[2, 2]) { stop("Invalid lower/upper bounds for the ranges of the Delta parameters")}
   }
   if(is.matrix(Delta_c) == TRUE) {
      if(dim(Delta_c)[1] != 2 | dim(Delta_c)[2] != 2) { stop("Delta parameters ranges must be provided as 2x2 matrices")}
      if(Delta_c[1, 1] > Delta_c[1, 2] | Delta_c[2, 1] > Delta_c[2, 2]) { stop("Invalid lower/upper bounds for the ranges of the Delta parameters")}
   }
   if(is.vector(Delta_e) == TRUE & is.vector(Delta_c) == TRUE){
      stop("Under MNAR the range of values for at least one Delta parameters must be provided in a 2x2 matrix form")
   } else if(is.matrix(Delta_e) == TRUE & is.vector(Delta_c) == TRUE){
      type = "MNAR_eff"
   } else if(is.vector(Delta_e) == TRUE & is.matrix(Delta_c) == TRUE) {
      type = "MNAR_cost"
   } else if(is.matrix(Delta_e) == TRUE & is.matrix(Delta_c) == TRUE) {
      type = "MNAR"
   }
  } 
  N1 <- data_read$arm_lengths[1]
  N2 <- data_read$arm_lengths[2]
  pe <- ncol(data_read$covariates_effects$Intervention)
  pc <- ncol(data_read$covariates_costs$Intervention)
  n_patterns1 <- data_read$patterns_list$n_patterns[1]
  n_patterns2 <- data_read$patterns_list$n_patterns[2]
  d1_list <- data_read$patterns_list$d1
  d2_list <- data_read$patterns_list$d2
  d1 <- data_read$patterns$Control
  d2 <- data_read$patterns$Intervention
  range_e <- Delta_e
  range_c <- Delta_c
  m_eff1 <- data_read$missing_effects$Control
  m_eff2 <- data_read$missing_effects$Intervention
  m_cost1 <- data_read$missing_costs$Control
  m_cost2 <- data_read$missing_costs$Intervention
  eff1 <- data_read$raw_effects$Control
  eff2 <- data_read$raw_effects$Intervention
  cost1 <- data_read$raw_costs$Control
  cost2 <- data_read$raw_costs$Intervention
  if(length(which(is.na(c(eff1, eff2)))) == 0 & length(which(is.na(c(cost1, cost2)))) == 0) {
    stop("At leat one missing value is required in either the effects or costs variables")
  }
  N1_cc <- data_read$arm_lengths_cc[, 1]
  N2_cc <- data_read$arm_lengths_cc[, 2]
  N1_mis <- data_read$arm_missing_data[, 1]
  N2_mis <- data_read$arm_missing_data[, 2]
  X1_e <- as.matrix(data_read$covariates_effects$Control)
  X2_e <- as.matrix(data_read$covariates_effects$Intervention)
  X1_c <- as.matrix(data_read$covariates_costs$Control)
  X2_c <- as.matrix(data_read$covariates_costs$Intervention)
  if(pe == 1) {
    X1_e <- as.vector(X1_e)
    X2_e <- as.vector(X2_e)
  }
  if(pc == 1) {
    X1_c <- as.vector(X1_c)
    X2_c <- as.vector(X2_c)
  }
  mean_cov_e1 <- as.vector(data_read$mean_cov_effects$Control)
  mean_cov_e2 <- as.vector(data_read$mean_cov_effects$Intervention)
  mean_cov_c1 <- as.vector(data_read$mean_cov_costs$Control)
  mean_cov_c2 <- as.vector(data_read$mean_cov_costs$Intervention)
  corr_assumption <- model.frame(formula = model.cost, data = data)
  if("e" %in% names(corr_assumption)) {
    ind = FALSE  
  } else {ind = TRUE}
  if(anyDuplicated(names(prior)) > 0) {
    stop("you cannot provide multiple priors with the same name") 
  }
  if(any(prior == "default") == TRUE) {
    prior <- list(default = "default")
  } else if(any(prior == "default") == FALSE) {
    list_check_vector <- lapply(prior, is.vector)
    if(all(as.logical(list_check_vector)) == FALSE) {
      stop("all user-supplied priors should be in vector format")
    }
    par_prior <- c("alpha0.prior", "beta0.prior", "sigma.prior.e", "sigma.prior.c", "patterns.prior", "alpha.prior", "beta.prior", "beta_f.prior")
    stop_mes <- "priors can be assigned only using specific character names depending on the type of model assumed. Type ''help(pattern)'' for more details"
    if(!all(names(list_check_vector) %in% par_prior == TRUE)) {stop(stop_mes) }
    if(is.vector(X1_e) == TRUE & identical(X1_e,rep(1,N1))) {
      if("alpha.prior" %in% names(list_check_vector)) {stop(stop_mes) }
    }
    if(is.vector(X1_c) == TRUE & identical(X1_c,rep(1,N1))) {
      if("beta.prior" %in% names(list_check_vector)) {stop(stop_mes) }
    }
    if(ind == TRUE) {
      if("beta_f.prior" %in% names(list_check_vector)) {stop(stop_mes) } 
    }
  }
  if(exists("sigma.prior.e", where=prior)) {sigma.prior.e = prior$sigma.prior.e} else {sigma.prior.e = NULL }
  if(exists("sigma.prior.c", where=prior)) {sigma.prior.c = prior$sigma.prior.c} else {sigma.prior.c = NULL }
  if(exists("alpha0.prior", where=prior)) {alpha0.prior = prior$alpha0.prior} else {alpha0.prior = NULL }
  if(exists("beta0.prior", where=prior)) {beta0.prior = prior$beta0.prior} else {beta0.prior = NULL }
  if(exists("alpha.prior", where=prior)) {alpha.prior = prior$alpha.prior} else {alpha.prior = NULL }
  if(exists("beta.prior", where=prior)) {beta.prior = prior$beta.prior} else {beta.prior = NULL }
  if(exists("patterns.prior", where=prior)) {patterns.prior = prior$patterns.prior} else {patterns.prior = NULL }
  if(exists("beta_f.prior", where = prior)) {beta_f.prior = prior$beta_f.prior} else {beta_f.prior = NULL }
  if(n_patterns1 < 2 | n_patterns2 < 2) { 
    stop("at least two patterns are required in each group to fit the model") }
  if(any(d1 == 1, na.rm = TRUE) == FALSE | any(d2 == 1, na.rm = TRUE) == FALSE) {
    stop("some completers must be observed in both treatment groups to fit the model") }
  if(all(d1 %in% c(1,3) == TRUE) & is.matrix(Delta_e) == TRUE) {stop("Cannot introduce sensitvity parameters for effects when all effects are observed in one arm") }
  if(all(d2 %in% c(1,3) == TRUE) & is.matrix(Delta_e) == TRUE) {stop("Cannot introduce sensitvity parameters for effects when all effects are observed in one arm") }
  if(all(d1 %in% c(1,2) == TRUE) & is.matrix(Delta_c) == TRUE) {stop("Cannot introduce sensitvity parameters for costs when all costs are observed in one arm") }
  if(all(d2 %in% c(1,2) == TRUE) & is.matrix(Delta_c) == TRUE) {stop("Cannot introduce sensitvity parameters for costs when all costs are observed in one arm") }
  d_list <- list()
  d_list[[1]] <- c(n_patterns1, n_patterns2)
  d_list[[2]] <- d1_list
  d_list[[3]] <- d2_list
  names(d_list) <- c("n_patterns", "d1", "d2")
  data_set <- list("effects" = data_read$raw_effects, "costs" = data_read$raw_costs, "N in reference arm" = N1, "N in comparator arm" = N2, 
                   "N observed in reference arm" = N1_cc, "N observed in comparator arm" = N2_cc, "N missing in reference arm" = N1_mis, "N missing in comparator arm"= N2_mis, 
                   "patterns in comparator arm" = data_read$patterns$Control, "patterns in reference arm" = data_read$patterns$Intervention, "covariates_effects" = data_read$covariates_effects, 
                   "covariates_costs" = data_read$covariates_costs, "missing_effects" = data_read$missing_effects, "missing_costs" = data_read$missing_costs)
  model_output <- run_pattern(type = type, dist_e = dist_e, dist_c = dist_c, inits = inits, d_list = d_list, d1 = d1, d2 = d2)
  if(save_model == FALSE) {
    unlink(filein)
  }
    if(exists("ref", where = exArgs)) {ref = exArgs$ref } else {ref = 2 }
    if(exists("interventions", where = exArgs)) {interventions = exArgs$interventions } else {interventions = NULL }
    if(exists("Kmax", where = exArgs)) {Kmax = exArgs$Kmax } else {Kmax = 50000 }
    if(exists("wtp", where = exArgs)) {wtp = exArgs$wtp } else {wtp = NULL }
    if(exists("plot", where = exArgs)) {plot = exArgs$plot } else {plot = FALSE }
    cea <- BCEA::bcea(e = model_output$mean_effects, c = model_output$mean_costs, ref = ref, interventions = interventions, Kmax = Kmax, wtp = wtp, plot = plot)
    res <- list(data_set = data_set, model_output = model_output, cea = cea, type = type)
  class(res) <- "missingHE"
  return(res)
}
