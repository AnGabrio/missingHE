#' Full Bayesian Models to handle missingness in Economic Evaluations (Pattern Mixture Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes under different missingness 
#' mechanism assumptions, using alternative parametric distributions for the effect and cost variables and a pattern mixture approach to identify the model. 
#' The analysis is performed using the \code{BUGS} language, which is implemented in the software \code{JAGS} using the function \code{\link[R2jags]{jags}}.
#' The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find the variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs).
#' Among these, effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't', respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' effectiveness outcome ('e') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' cost outcome ('c') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. 
#' By default, covariates are placed on the "location" parameter of the distribution through a linear model. A joint bivariate distribution for effects and costs can be specified by including 'e' on the right-hand side of the formula for the costs model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param Delta_e Range of values for the prior on the sensitivity parameters used to identify the mean of the effects under MNAR. The value must be set to 0 under MAR. 
#' @param Delta_c Range of values for the prior on the sensitivity parameters used to identify the mean of the costs under MNAR. The value must be set to 0 under MAR.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param restriction type of identifying restriction to be imposed to identify the distributions of the missing data in each pattern. 
#' Available choices are: complete case restrcition ('CC') - default - or available case restriction ('AC'). 
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#' CI sample quantiles to be calculated and returned for the imputed values.
#' @param n.chains Number of chains.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of warmup iterations.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{JAGS}
#' will generate initial values for all the model parameters.
#' @param n.thin Thinning interval.
#' @param ppc Logical. If \code{ppc} is \code{TRUE}, the estimates of the parameters that can be used to generate replications from the model are saved.
#' @param save_model Logical. If \code{save_model} is \code{TRUE} a \code{txt} file containing the model code is printed
#' in the current working directory.
#' @param prior A list containing the hyperprior values provided by the user. Each element of this list must be a vector of length two
#' containing the user-provided hyperprior values and must be named with the name of the corresponding parameter. For example, the hyperprior
#' values for the standard deviation effect parameters can be provided using the list \code{prior = list('sigma.prior.e' = c(0, 100))}.
#' For more information about how to provide prior hypervalues for different types of parameters and models see details. 
#' If \code{prior} is set to 'default', the default values will be used.  
#' @param ... Additional arguments that can be provided by the user. Examples are \code{center = TRUE} to center all the covariates in the model 
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
#'   \item{data_format}{A character variable that indicate which type of analysis was conducted, either using a \code{wide} or \code{long} dataset}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA JAGS missing data Pattern Mixture Models
#' @importFrom stats model.frame terms
#' @details Depending on the distributions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} and the type of missingness mechanism specified in the argument \code{type}, different pattern mixture models
#' are built and run in the background by the function \code{pattern}. The model for the outcomes is fitted in each missingness pattern 
#' and the parameters indexing the missing data distributions are identified using: the corresponding parameters identified from the observed data
#' in other patterns (under 'MAR'); or a combination of the parameters identified by the observed data and some sensitivity parameters (under 'MNAR'). 
#' A simple example can be used to show how pattern mixture models are specified. 
#' Consider a data set comprising a response variable \eqn{y} and a set of centered covariate \eqn{X_j}. We denote with \eqn{d_i} the patterns' indicator variable for each 
#' subject in the trial \eqn{i = 1, ..., n} such that: \eqn{d_i = 1} indicates the completers (both e and c observed), \eqn{d_i = 2} and \eqn{d_i = 3} indicate that 
#' only the costs or effects are observed, respectively, while \eqn{d_i = 4} indicates that neither of the two outcomes is observed. In general, a different number of patterns 
#' can be observed between the treatment groups and \code{missingHE} accounts for this possibility by modelling a different patterns' indicator variables for each arm. 
#' For simplicity, in this example, we assume that the same number of patterns is observed in both groups. \eqn{d_i} is assigned a multinomial distribution, 
#' which probabilities are modelled using a Dirichlet prior (by default giving to each pattern the same weight). Next, the model specified in \code{dist_e} 
#' and \code{dist_c} is fitted in each pattern. The parameters that cannot be identified by the observed data in each pattern (d = 2, 3, 4), e.g. the means.
#' \eqn{mu_e[d]} and \code{mu_c[d]}, can be identified using the parameters estimated from other patterns. Two choices are currently available: the complete cases ('CC') or available cases ('AC').
#' For example, using the 'CC' restriction, the parameters indexing the distributions of the missing data are identified as: 
#' \deqn{mu_e[2] = \mu_e[4] = \mu_e[1] + \Delta_e}
#' \deqn{mu_c[3] = \mu_c[4] = \mu_c[1] + \Delta_c}
#' where
#' \itemize{
#' \item \eqn{\mu_e[1]} is the effects mean for the completers.
#' \item \eqn{\mu_c[1]} is the costs mean for the completers.
#' \item \eqn{\Delta_e} is the sensitivity parameters associated with the marginal effects mean.
#' \item \eqn{\Delta_c} is the sensitivity parameters associated with the marginal costs mean.
#' }
#' If the 'AC' restriction is chosen, only the parameters estimated from the observed data in pattern 2 (costs) and pattern 3 (effects) are used to identify those in the other patterns.  
#' When \eqn{\Delta_e = 0} and \eqn{\Delta_c = 0} the model assumes a 'MAR' mechanism. When \eqn{\Delta_e != 0} and/or \eqn{\Delta_c != 0} 'MNAR' departues for the 
#' effects and/or costs are explored assuming a Uniform prior distributions for the sensitivity parameters. The range of values for these priors is defined based on the
#' boundaries specified in \code{Delta_e} and \code{Delta_c} (see Arguments), which must be provided by the user. 
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{pattern}, the elements of this list (see Arguments)
#' must be vectors of length two containing the user-provided hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names for the parameters indexing the model which are accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item location parameters \eqn{\alpha_0} and \eqn{\beta_0}: "mean.prior.e"(effects) and/or "mean.prior.c"(costs)
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j} and \eqn{\beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' } 
#' The only exception is the missingness patterns' probability \eqn{\pi}, denoted with "patterns.prior", whose hyperprior values must be provided as a list
#' formed by two elements. These must be vectors of the same length equal to the number of patterns in the control (first element) and intervention (second element) group.
#' 
#' For each model, random effects can also be specified for each parameter by adding the term + (x | z) to each model formula, 
#' where x is the fixed regression coefficient for which also the random effects are desired and z is the clustering variable across which 
#' the random effects are specified (must be the name of a factor variable in the dataset). Multiple random effects can be specified using the 
#' notation + (x1 + x2 | site) for each covariate that was included in the fixed effects formula. Random intercepts are included by default in the models
#' if a random effects are specified but they can be removed by adding the term 0 within the random effects formula, e.g. + (0 + x | z). 
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
#' # Quck example to run using subset of MenSS dataset
#' MenSS.subset <- MenSS[50:100, ]
#' 
#' # Run the model using the pattern function assuming a SCAR mechanism
#' # Use only 100 iterations to run a quick check
#' model.pattern <- pattern(data = MenSS.subset,model.eff = e~1,model.cost = c~1,
#'    dist_e = "norm", dist_c = "norm",type = "MAR", Delta_e = 0, Delta_c = 0, 
#'    n.chains = 2, n.iter = 100, ppc = FALSE)
#' 
#' # Print the results of the JAGS model
#' print(model.pattern)
#' #
#'
#' # Use dic information criterion to assess model fit
#' pic.dic <- pic(model.pattern, criterion = "dic", module = "total")
#' pic.dic
#' #
#' 
#' # Extract regression coefficient estimates
#' coef(model.pattern)
#' #
#' 
#' \dontshow{
#' # Use waic information criterion to assess model fit
#' pic.waic <- pic(model.pattern, criterion = "waic", module = "total")
#' pic.waic
#' }
#'
#' # Assess model convergence using graphical tools
#' # Produce histograms of the posterior samples for the mean effects
#' diag.hist <- diagnostic(model.pattern, type = "histogram", param = "mu.e")
#' #
#'
#' # Compare observed effect data with imputations from the model
#' # using plots (posteiror means and credible intervals)
#' p1 <- plot(model.pattern, class = "scatter", outcome = "effects")
#' #
#'
#' # Summarise the CEA information from the model
#' summary(model.pattern)
#' 
#' \donttest{
#' # Further examples which take longer to run
#' model.pattern <- pattern(data = MenSS, model.eff = e ~ u.0,model.cost = c ~ e,
#'    Delta_e = 0, Delta_c = 0, dist_e = "norm", dist_c = "norm",
#'    type = "MAR", n.chains = 2, n.iter = 500, ppc = FALSE)
#' #
#' # Print results for all imputed values
#' print(model.pattern, value.mis = TRUE)
#' 
#' # Use looic to assess model fit
#' pic.looic<-pic(model.pattern, criterion = "looic", module = "total")
#' pic.looic
#' 
#' # Show density plots for all parameters
#' diag.hist <- diagnostic(model.pattern, type = "denplot", param = "all")
#' 
#' # Plots of imputations for all data
#' p1 <- plot(model.pattern, class = "scatter", outcome = "all")
#' 
#' # Summarise the CEA results
#' summary(model.pattern)
#' 
#' }
#' #
#' #


pattern <- function(data, model.eff, model.cost, dist_e, dist_c, Delta_e, Delta_c, type, restriction = "CC", prob = c(0.025, 0.975), n.chains = 2, n.iter = 20000, 
                    n.burnin = floor(n.iter / 2), inits = NULL, n.thin = 1, ppc = FALSE, save_model = FALSE, prior = "default", ...) {
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
  if(!restriction %in% c("CC", "AC")){
    stop("Only 'CC' or 'AC' types of restriction are allowed")
  }
  dist_e <- tolower(dist_e)
  dist_c <- tolower(dist_c)
  if(dist_e == "normal") { dist_e <- "norm" }
  if(dist_e == "exponential") { dist_e <- "exp" }
  if(dist_e == "logistic") { dist_e <- "logis" }
  if(dist_e == "bernoulli") { dist_e <- "bern" }
  if(dist_e == "poisson") { dist_e <- "pois" }
  if(dist_e == "negative binomial") { dist_e <- "nbinom" }
  if(dist_c == "normal") { dist_c <- "norm" }
  if(dist_c == "lognormal") { dist_c <- "lnorm" }
  if(!dist_e %in% c("norm", "beta", "exp", "weibull", "logis", "bern", "pois", "nbinom", "gamma") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta', 'gamma', 'logis', 'exp', 'weibull', 'nbinom', 'pois', 'bern'  for the effects and 'norm', 'gamma', 'lnorm' for the costs")
  }
  type <- toupper(type)
  if(!type %in% c("MAR", "MNAR")) {
    stop("Types available for use are 'MAR' and 'MNAR'")
  }
  if(length(prob) != 2 | is.numeric(prob) == FALSE | any(prob < 0) != FALSE | any(prob > 1) != FALSE) {
    stop("You must provide valid lower/upper quantiles for the imputed data distribution")
  }
  if(is.logical(save_model) == FALSE | is.logical(ppc) == FALSE) {
    stop("save_model and ppc are logical arguments and should be either TRUE or FALSE")
  }
  exArgs <- list(...)
  if(exists("center", where = exArgs)) {
    if(is.logical(exArgs$center) == FALSE) { stop("center must be either TRUE or FALSE") }
    center = exArgs$center 
  } else {center = FALSE }
  data_read <- data_read_pattern(data = data, model.eff = model.eff, model.cost = model.cost, type = type, center = center)
  model_e_fixed <- labels(terms(data_read$model_formula$mf_model.e_fixed))
  model_c_fixed <- labels(terms(data_read$model_formula$mf_model.c_fixed))
  if(as.character(data_read$model_formula$mf_model.e_random)[3] == "1") {
    model_e_random <- c("1")
  } else { model_e_random <- labels(terms(data_read$model_formula$mf_model.e_random)) }
  if(as.character(data_read$model_formula$mf_model.c_random)[3] == "1") {
    model_c_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.c_random)[3] == "1 + e") {
    model_c_random <- c("1", "e")
  } else { model_c_random <- labels(terms(data_read$model_formula$mf_model.c_random)) }
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
  N1 <- data_read$data_raw$arm_lengths[1]
  N2 <- data_read$data_raw$arm_lengths[2]
  pe_fixed <- ncol(data_read$data_raw$covariates_effects_fixed$Intervention)
  pc_fixed <- ncol(data_read$data_raw$covariates_costs_fixed$Intervention)
  pe_random <- ncol(data_read$data_raw$covariates_effects_random$Intervention)
  if(length(pe_random) == 0) {pe_random <- 0 }
  pc_random <- ncol(data_read$data_raw$covariates_costs_random$Intervention)
  if(length(pc_random) == 0) {pc_random <- 0 }
  n_patterns1 <- data_read$data_raw$patterns_list$n_patterns[1]
  n_patterns2 <- data_read$data_raw$patterns_list$n_patterns[2]
  d1_list <- data_read$data_raw$patterns_list$d1
  d2_list <- data_read$data_raw$patterns_list$d2
  d1 <- data_read$data_raw$patterns$Control
  d2 <- data_read$data_raw$patterns$Intervention
  range_e <- Delta_e
  range_c <- Delta_c
  m_eff1 <- data_read$data_raw$missing_effects$Control
  m_eff2 <- data_read$data_raw$missing_effects$Intervention
  m_cost1 <- data_read$data_raw$missing_costs$Control
  m_cost2 <- data_read$data_raw$missing_costs$Intervention
  eff1 <- data_read$data_raw$raw_effects$Control
  eff2 <- data_read$data_raw$raw_effects$Intervention
  cost1 <- data_read$data_raw$raw_costs$Control
  cost2 <- data_read$data_raw$raw_costs$Intervention
  clus1_e <- data_read$data_raw$clus_e$Control
  clus2_e <- data_read$data_raw$clus_e$Intervention
  clus1_c <- data_read$data_raw$clus_c$Control
  clus2_c <- data_read$data_raw$clus_c$Intervention
  n1_clus_e <- length(unique(clus1_e))
  n2_clus_e <- length(unique(clus2_e))
  n1_clus_c <- length(unique(clus1_c))
  n2_clus_c <- length(unique(clus2_c))
  if(length(which(is.na(c(eff1, eff2)))) == 0 & length(which(is.na(c(cost1, cost2)))) == 0) {
    stop("At leat one missing value is required in either the effects or costs variables")
  }
  N1_cc <- data_read$data_raw$arm_lengths_cc[, 1]
  N2_cc <- data_read$data_raw$arm_lengths_cc[, 2]
  N1_mis <- data_read$data_raw$arm_missing_data[, 1]
  N2_mis <- data_read$data_raw$arm_missing_data[, 2]
  X1_e_fixed <- as.matrix(data_read$data_raw$covariates_effects_fixed$Control)
  X2_e_fixed <- as.matrix(data_read$data_raw$covariates_effects_fixed$Intervention)
  X1_c_fixed <- as.matrix(data_read$data_raw$covariates_costs_fixed$Control)
  X2_c_fixed <- as.matrix(data_read$data_raw$covariates_costs_fixed$Intervention)
  if(pe_fixed == 1) {
    X1_e_fixed <- as.vector(X1_e_fixed)
    X2_e_fixed <- as.vector(X2_e_fixed)
  }
  if(pc_fixed == 1) {
    X1_c_fixed <- as.vector(X1_c_fixed)
    X2_c_fixed <- as.vector(X2_c_fixed)
  }
  mean_cov_e1_fixed <- as.vector(data_read$data_raw$mean_cov_effects_fixed$Control)
  mean_cov_e2_fixed <- as.vector(data_read$data_raw$mean_cov_effects_fixed$Intervention)
  mean_cov_c1_fixed <- as.vector(data_read$data_raw$mean_cov_costs_fixed$Control)
  mean_cov_c2_fixed <- as.vector(data_read$data_raw$mean_cov_costs_fixed$Intervention)
  if(length(model_e_random) != 0 & pe_random != 0) { 
    X1_e_random <- as.matrix(data_read$data_raw$covariates_effects_random$Control)
    X2_e_random <- as.matrix(data_read$data_raw$covariates_effects_random$Intervention)
    if(pe_random == 1) {
      X1_e_random <- as.vector(X1_e_random)
      X2_e_random <- as.vector(X2_e_random)
    }
    mean_cov_e1_random <- as.vector(data_read$data_raw$mean_cov_effects_random$Control)
    mean_cov_e2_random <- as.vector(data_read$data_raw$mean_cov_effects_random$Intervention)
  } else { 
    X1_e_random <- X2_e_random <- NULL
    mean_cov_e1_random <- mean_cov_e2_random <- NULL }
  if(length(model_c_random) != 0 & pc_random != 0) { 
    X1_c_random <- as.matrix(data_read$data_raw$covariates_costs_random$Control)
    X2_c_random <- as.matrix(data_read$data_raw$covariates_costs_random$Intervention)
    if(pc_random == 1) {
      X1_c_random <- as.vector(X1_c_random)
      X2_c_random <- as.vector(X2_c_random)
    }
    mean_cov_c1_random <- as.vector(data_read$data_raw$mean_cov_costs_random$Control)
    mean_cov_c2_random <- as.vector(data_read$data_raw$mean_cov_costs_random$Intervention)
  } else { 
    X1_c_random <- X2_c_random <- NULL 
    mean_cov_c1_random <- mean_cov_c2_random <- NULL }
  corr_assumption_fixed <- model.frame(formula = data_read$model_formula$mf_model.c_fixed, data = data)
  if("e" %in% names(corr_assumption_fixed)) {
    ind_fixed = FALSE  
  } else {
    ind_fixed = TRUE 
    ind_random = TRUE 
  }
  if(ind_fixed == FALSE & "e" %in% model_c_random) {
    ind_random = FALSE
  } else if(ind_fixed == FALSE & !("e" %in% model_c_random)) {
    ind_random = TRUE
  }
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
    par_prior_fixed <- c("alpha0.prior", "beta0.prior", "sigma.prior.e", "sigma.prior.c", "patterns.prior", 
                         "alpha.prior", "beta.prior", "beta_f.prior")
    par_prior_random <- c("mu.a0.prior", "mu.b0.prior", "mu.a.prior", "mu.b.prior", "mu.b_f.prior",
                          "s.a0.prior", "s.b0.prior", "s.a.prior", "s.b.prior", "s.b_f.prior")
    stop_mes <- "priors can be assigned only using specific character names depending on the type of model assumed. Type ''help(pattern)'' for more details"
    if(is.vector(X1_e_fixed) == TRUE & identical(X1_e_fixed, rep(1, N1))) {
      if("alpha.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(is.vector(X1_c_fixed) == TRUE & identical(X1_c_fixed, rep(1, N1))) {
      if("beta.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(length(model_e_random) == 0 | pe_random == 0) {
      if("mu.a0.prior" %in% names(list_check_vector) | "s.a0.prior" %in% names(list_check_vector) |
         "mu.a.prior" %in% names(list_check_vector) | "s.a.prior" %in% names(list_check_vector)) { stop(stop_mes)}
    } else if(length(model_e_random) != 0 & pe_random == 1) { 
      if("mu.a.prior" %in% names(list_check_vector) | "s.a.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(length(model_c_random) == 0 | pc_random == 0) {
      if("mu.b0.prior" %in% names(list_check_vector) | "s.b0.prior" %in% names(list_check_vector) |
         "mu.b.prior" %in% names(list_check_vector) | "s.b.prior" %in% names(list_check_vector)) { stop(stop_mes)}
    } else if(length(model_c_random) != 0 & pc_random == 1) { 
      if("mu.b.prior" %in% names(list_check_vector) | "s.b.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(ind_fixed == TRUE) {
      if("beta_f.prior" %in% names(list_check_vector)) { stop(stop_mes) } 
      if("b_f.prior" %in% names(list_check_vector)) { stop(stop_mes) } 
    }
    if(ind_fixed == FALSE & ind_random == TRUE) {
      if("mu.b_f.prior" %in% names(list_check_vector) | "s.b_f.prior" %in% names(list_check_vector)) { stop(stop_mes) } 
    }
  }
  if(length(model_c_random) == 1) {
    if(model_c_random == "e") {is_c_random_c <- TRUE} else {is_c_random_c <- FALSE}
  } else {is_c_random_c <- FALSE }
  if(length(model_c_random) == 2) {
    if(all(model_c_random == c("1", "e")) == TRUE) {is_int_c_random_c <- TRUE} else {is_int_c_random_c <- FALSE}
  } else {is_int_c_random_c <- FALSE }
  if(exists("sigma.prior.e", where = prior)) {sigma.prior.e = prior$sigma.prior.e} else {sigma.prior.e = NULL }
  if(exists("sigma.prior.c", where = prior)) {sigma.prior.c = prior$sigma.prior.c} else {sigma.prior.c = NULL }
  if(exists("alpha0.prior", where = prior)) {alpha0.prior = prior$alpha0.prior} else {alpha0.prior = NULL }
  if(exists("beta0.prior", where = prior)) {beta0.prior = prior$beta0.prior} else {beta0.prior = NULL }
  if(exists("alpha.prior", where = prior)) {alpha.prior = prior$alpha.prior} else {alpha.prior = NULL }
  if(exists("beta.prior", where = prior)) {beta.prior = prior$beta.prior} else {beta.prior = NULL }
  if(exists("patterns.prior", where = prior)) {patterns.prior = prior$patterns.prior} else {patterns.prior = NULL }
  if(exists("beta_f.prior", where = prior)) {beta_f.prior = prior$beta_f.prior} else {beta_f.prior = NULL }
  if(exists("mu.a0.prior", where = prior)) {mu.a0.prior = prior$mu.a0.prior} else {mu.a0.prior = NULL }
  if(exists("s.a0.prior", where = prior)) {s.a0.prior = prior$s.a0.prior} else {s.a0.prior = NULL }
  if(exists("mu.b0.prior", where = prior)) {mu.b0.prior = prior$mu.b0.prior} else {mu.b0.prior = NULL }
  if(exists("s.b0.prior", where = prior)) {s.b0.prior = prior$s.b0.prior} else {s.b0.prior = NULL }
  if(exists("mu.a.prior", where = prior)) {mu.a.prior = prior$mu.a.prior} else {mu.a.prior = NULL }
  if(exists("s.a.prior", where = prior)) {s.a.prior = prior$s.a.prior} else {s.a.prior = NULL }
  if(exists("mu.b.prior", where = prior)) {mu.b.prior = prior$mu.b.prior} else {mu.b.prior = NULL }
  if(exists("s.b.prior", where = prior)) {s.b.prior = prior$s.b.prior} else {s.b.prior = NULL }
  if(exists("mu.b_f.prior", where = prior)) {mu.b_f.prior = prior$mu.b_f.prior} else {mu.b_f.prior = NULL }
  if(exists("s.b_f.prior", where = prior)) {s.b_f.prior = prior$s.b_f.prior} else {s.b_f.prior = NULL }
  if(n_patterns1 < 2 | n_patterns2 < 2) { 
    stop("at least two patterns are required in each group to fit the model") }
  if(restriction == "CC"){
    if(any(d1 == 1, na.rm = TRUE) == FALSE | any(d2 == 1, na.rm = TRUE) == FALSE) {
      stop("some completers must be observed in both treatment groups to fit the model") }
  }
  if(restriction == "AC"){
    if(any(d1 == 3, na.rm = TRUE) == FALSE | any(d2 == 3, na.rm = TRUE) == FALSE |
       any(d1 == 2, na.rm = TRUE) == FALSE | any(d2 == 2, na.rm = TRUE) == FALSE) {
      stop("some non-completers for both outcomes must be observed in both treatment groups to fit the model") }
    if(ind_fixed == FALSE){
      if(any(d1 == 1, na.rm = TRUE) == FALSE | any(d2 == 1, na.rm = TRUE) == FALSE) {
        stop("some completers must be observed in both treatment groups when 'e' is included in the model of 'c' under AC restriction") }
      if(n_patterns1 == 2 | n_patterns2 == 2){
        stop("at least three patterns in each arm are required to fit the model when 'e' is included in the model of 'c' under AC restriction") }
    }
  }
  if(all(d1 %in% c(1,3) == TRUE) & is.matrix(Delta_e) == TRUE) {stop("Cannot introduce sensitvity parameters for effects when all effects are observed in one arm") }
  if(all(d2 %in% c(1,3) == TRUE) & is.matrix(Delta_e) == TRUE) {stop("Cannot introduce sensitvity parameters for effects when all effects are observed in one arm") }
  if(all(d1 %in% c(1,2) == TRUE) & is.matrix(Delta_c) == TRUE) {stop("Cannot introduce sensitvity parameters for costs when all costs are observed in one arm") }
  if(all(d2 %in% c(1,2) == TRUE) & is.matrix(Delta_c) == TRUE) {stop("Cannot introduce sensitvity parameters for costs when all costs are observed in one arm") }
  d_list <- list()
  d_list[[1]] <- c(n_patterns1, n_patterns2)
  d_list[[2]] <- d1_list
  d_list[[3]] <- d2_list
  names(d_list) <- c("n_patterns", "d1", "d2")
  data_set <- list("effects" = data_read$data_raw$raw_effects, "costs" = data_read$data_raw$raw_costs, "N in reference arm" = N1, "N in comparator arm" = N2, 
                   "N observed in reference arm" = N1_cc, "N observed in comparator arm" = N2_cc, "N missing in reference arm" = N1_mis, "N missing in comparator arm"= N2_mis, 
                   "patterns in comparator arm" = data_read$data_raw$patterns$Control, "patterns in reference arm" = data_read$data_raw$patterns$Intervention, 
                   "covariates_effects_fixed" = data_read$data_raw$covariates_effects_fixed, "covariates_costs_fixed" = data_read$data_raw$covariates_costs_fixed,
                   "covariates_effects_random" = data_read$data_raw$covariates_effects_random, "covariates_costs_random" = data_read$data_raw$covariates_costs_random, 
                   "missing_effects" = data_read$data_raw$missing_effects, "missing_costs" = data_read$data_raw$missing_costs, "clus_effects" = data_read$data_raw$clus_e, "clus_costs" = data_read$data_raw$clus_c)
  model_output <- run_pattern(type = type, dist_e = dist_e, dist_c = dist_c, inits = inits, d_list = d_list, d1 = d1, d2 = d2, restriction = restriction, ppc = ppc)
  if(save_model == FALSE) {
    unlink(filein)
  }
  if(exists("ref", where = exArgs)) {ref = exArgs$ref } else {ref = 2 }
  if(exists("interventions", where = exArgs)) {interventions = exArgs$interventions } else {interventions = NULL }
  if(exists("Kmax", where = exArgs)) {Kmax = exArgs$Kmax } else {Kmax = 50000 }
  if(exists("wtp", where = exArgs)) {wtp = exArgs$wtp } else {wtp = NULL }
  if(exists("plot", where = exArgs)) {plot = exArgs$plot } else {plot = FALSE }
  cea <- BCEA::bcea(e = model_output$mean_effects, c = model_output$mean_costs, ref = ref, interventions = interventions, Kmax = Kmax, k = wtp, plot = plot)
  format <- "wide"
  res <- list(data_set = data_set, model_output = model_output, cea = cea, type = type, data_format = format)
  class(res) <- "missingHE"
  return(res)
}