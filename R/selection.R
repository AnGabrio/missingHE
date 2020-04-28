#' Full Bayesian Models to handle missingness in Economic Evaluations (Selection Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes under different missing data 
#' mechanism assumptions, using alternative parametric distributions for the effect and cost variables and 
#' using a selection model approach to identify the model. The analysis is performed using the \code{BUGS} language, 
#' which is implemented in the software \code{JAGS} using the function \code{\link[R2jags]{jags}} The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find the variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.me}, \code{model.mc} (model formulas for the missing effect and cost models). Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't', respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' effectiveness outcome ('e') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' cost outcome ('c') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model. 
#' A joint bivariate distribution for effects and costs can be specified by including 'e' on the right-hand side of the formula for the costs model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.me A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the 
#' term 'me'(missing effects) and any covariates must be provided on the right-hand side of the formula. If there are no covariates, \code{1} should be specified on the right hand side of the formula. 
#' By default, covariates are placed on the "probability" parameter for the missing effects through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.mc A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the term 'mc'(missing costs) and any covariates must be provided on the right-hand side of the formula. 
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "probability" parameter for the missing costs through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the imputed values.
#' @param n.chains Number of chains.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of warmup iterations.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{JAGS}
#' will generate initial values for all the model parameters.
#' @param n.thin Thinning interval.
#' @param ppc Logical. If \code{ppc} is \code{TRUE}, the estimates of the parameters that can be used to generate replications from the model are saved.
#' @param save_model Logical. If \code{save_model} is \code{TRUE}, a \code{txt} file containing the model code is printed
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
#'   \item{type}{A character variable that indicate which type of missingness mechanism has been used to run the model, 
#'   either \code{MAR} or \code{MNAR} (see details)}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA JAGS missing data Selection Models
#' @importFrom stats model.frame terms
#' @details Depending on the distributions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} and the type of missingness mechanism specified in the argument \code{type}, different selection models
#' are built and run in the background by the function \code{selection}. These models consist in logistic regressions that are used to estimate
#' the probability of missingness in one or both the outcomes. A simple example can be used to show how selection models are specified. 
#' Consider a data set comprising a response variable \eqn{y} and a set of centered covariate \eqn{X_j}. For each subject in the trial \eqn{i = 1, ..., n}
#' we define an indicator variable \eqn{m_i} taking value \code{1} if the \eqn{i}-th individual is associated with a missing value and \code{0} otherwise.
#' This is modelled as:
#' \deqn{m_i ~ Bernoulli(\pi_i)}
#' \deqn{logit(\pi_i) = \gamma_0 + \sum\gamma_j X_j + \delta(y)}
#' where
#' \itemize{
#' \item \eqn{\pi_i} is the individual probability of a missing value in \eqn{y}
#' \item \eqn{\gamma_0} represents the marginal probability of a missing value in \eqn{y} on the logit scale.
#' \item \eqn{\gamma_j} represents the impact on the probability of a missing value in \eqn{y} of the centered covariates \eqn{X_j}.
#' \item \eqn{\delta} represents the impact on the probability of a missing value in \eqn{y} of the missing value itself.
#' }
#' When \eqn{\delta = 0} the model assumes a 'MAR' mechanism, while when \eqn{\delta != 0} the mechanism is 'MNAR'. For the parameters indexing the missingness model, 
#' the default prior distributions assumed are the following:
#' \itemize{
#' \item \eqn{\gamma_0 ~ Logisitc(0, 1)}
#' \item \eqn{\gamma_j ~ Normal(0, 0.01)}
#' \item \eqn{\delta ~ Normal(0, 1)}
#' }
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{selection}, the elements of this list (see Arguments)
#' must be vectors of length two containing the user-provided hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names for the parameters indexing the model which are accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item location parameters \eqn{\alpha_0} and \eqn{\beta_0}: "mean.prior.e"(effects) and/or "mean.prior.c"(costs)
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j} and \eqn{\beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' \item marginal probability of missing values \eqn{\gamma_0}: "p.prior.e"(effects) and/or "p.prior.c"(costs)
#' \item covariate parameters in the missingness model \eqn{\gamma_j} (if covariate data provided): "gamma.prior.e"(effects) and/or "gamma.prior.c"(costs)
#' \item mnar parameter \eqn{\delta}: "delta.prior.e"(effects) and/or "delta.prior.c"(costs)
#' } 
#' For simplicity, here we have assumed that the set of covariates \eqn{X_j} used in the models for the effects/costs and in the 
#' model of the missing effect/cost values is the same. However, it is possible to specify different sets of covariates for each model
#' using the arguments in the function \code{selection} (see Arguments).
#' 
#' For each model, random effects can also be specified for each parameter by adding the term + (x | z) to each model formula, 
#' where x is the fixed regression coefficient for which also the random effects are desired and z is the clustering variable across which 
#' the random effects are specified (must be the name of a factor variable in the dataset). Multiple random effects can be specified using the 
#' notation + (x1 + x2 | site) for each covariate that was included in the fixed effects formula. Random intercepts are included by default in the models
#' if a random effects are specified but they can be removed by adding the term 0 within the random effects formula, e.g. + (0 + x | z). 
#' 
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
#' # Run the model using the selection function assuming a SCAR mechanism
#' # Use only 100 iterations to run a quick check
#' model.selection <- selection(data = MenSS.subset, model.eff = e ~ 1,model.cost = c ~ 1,
#'    model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
#'    type = "MAR", n.chains = 2, n.iter = 100, ppc = TRUE)
#' 
#' # Print the results of the JAGS model
#' print(model.selection)
#' #
#'
#' # Use dic information criterion to assess model fit
#' pic.dic <- pic(model.selection, criterion = "dic", module = "total")
#' pic.dic
#' #
#' 
#' # Extract regression coefficient estimates
#' coef(model.selection)
#' #
#' 
#' \dontshow{
#' # Use waic information criterion to assess model fit
#' pic.waic <- pic(model.selection, criterion = "waic", module = "total")
#' pic.waic
#' }
#'
#' # Assess model convergence using graphical tools
#' # Produce histograms of the posterior samples for the mean effects
#' diag.hist <- diagnostic(model.selection, type = "histogram", param = "mu.e")
#' #
#'
#' # Compare observed effect data with imputations from the model
#' # using plots (posteiror means and credible intervals)
#' p1 <- plot(model.selection, class = "scatter", outcome = "effects")
#' #
#'
#' # Summarise the CEA information from the model
#' summary(model.selection)
#' 
#' \donttest{
#' # Further examples which take longer to run
#' model.selection <- selection(data = MenSS, model.eff = e ~ u.0,model.cost = c ~ e,
#'    model.se = me ~ u.0, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
#'    type = "MAR", n.chains = 2, n.iter = 500, ppc = FALSE)
#' #
#' # Print results for all imputed values
#' print(model.selection, value.mis = TRUE)
#' 
#' # Use looic to assess model fit
#' pic.looic<-pic(model.selection, criterion = "looic", module = "total")
#' pic.looic
#' 
#' # Show density plots for all parameters
#' diag.hist <- diagnostic(model.selection, type = "denplot", param = "all")
#' 
#' # Plots of imputations for all data
#' p1 <- plot(model.selection, class = "scatter", outcome = "all")
#' 
#' # Summarise the CEA results
#' summary(model.selection)
#' 
#' }
#' #
#' #


selection <- function(data, model.eff, model.cost, model.me = me ~ 1, model.mc = mc ~ 1, dist_e, dist_c, type, prob = c(0.025, 0.975), 
                      n.chains = 2, n.iter = 20000, n.burnin = floor(n.iter / 2), inits = NULL, n.thin = 1, ppc = FALSE, save_model = FALSE, prior = "default", ...) {
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
  if(!all(levels(as.factor(cov_matrix$t)) %in% c("1", "2")) == TRUE) {
    stop("A two arm indicator variable must be provided with '1' for the control and '2' for the other intervention")
  }
  if(is.character(type) == FALSE | is.character(dist_e) == FALSE | is.character(dist_c) == FALSE) {
    stop("you must provide character names for the objects 'type', 'dist_e' and 'dist_c'")
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
  data_read <- data_read_selection(data = data, model.eff = model.eff, model.cost = model.cost, 
                                   model.me = model.me, model.mc = model.mc, type = type, center = center)
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
  if(as.character(data_read$model_formula$mf_model.me_random)[3] == "1") {
    model_me_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.me_random)[3] == "1 + e") {
    model_me_random <- c("1", "e")
  } else { model_me_random <- labels(terms(data_read$model_formula$mf_model.me_random)) }
  if(as.character(data_read$model_formula$mf_model.mc_random)[3] == "1") {
    model_mc_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.mc_random)[3] == "1 + c") {
    model_mc_random <- c("1", "c")
  } else { model_mc_random <- labels(terms(data_read$model_formula$mf_model.mc_random)) }
  miss_eff_assumption <- model.frame(formula = data_read$model_formula$mf_model.me_fixed, data = data_read$data_raw$data_ind)
  miss_cost_assumption <- model.frame(formula = data_read$model_formula$mf_model.mc_fixed, data = data_read$data_raw$data_ind)
  if("e" %in% names(miss_eff_assumption) == TRUE & "c" %in% names(miss_cost_assumption) == FALSE) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'mode.mc' if 'MAR' type selected") }
    type = "MNAR_eff"
  } else if("e" %in% names(miss_eff_assumption) == FALSE & "c" %in% names(miss_cost_assumption) == TRUE) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'mode.mc' if 'MAR' type selected") }
    type = "MNAR_cost"
  } else if("e" %in% names(miss_eff_assumption) == TRUE & "c" %in% names(miss_cost_assumption) == TRUE) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'mode.mc' if 'MAR' type selected") }
    type = "MNAR"
  } else if("e" %in% names(miss_eff_assumption) == FALSE & "c" %in% names(miss_cost_assumption) == FALSE) {
    if(type == "MNAR") { stop("Please add 'e' and/or 'c' to 'model.me' and/or 'mode.mc' if 'MNAR' type selected") }
    type = "MAR"
  }
  N1 <- data_read$data_raw$arm_lengths[1]
  N2 <- data_read$data_raw$arm_lengths[2]
  pe_fixed <- ncol(data_read$data_raw$covariates_effects_fixed$Intervention)
  pc_fixed <- ncol(data_read$data_raw$covariates_costs_fixed$Intervention)
  ze_fixed <- ncol(data_read$data_raw$covariates_missing_effects_fixed$Intervention)
  zc_fixed <- ncol(data_read$data_raw$covariates_missing_costs_fixed$Intervention)
  pe_random <- ncol(data_read$data_raw$covariates_effects_random$Intervention)
  if(length(pe_random) == 0) {pe_random <- 0 }
  pc_random <- ncol(data_read$data_raw$covariates_costs_random$Intervention)
  if(length(pc_random) == 0) {pc_random <- 0 }
  ze_random <- ncol(data_read$data_raw$covariates_missing_effects_random$Intervention)
  if(length(ze_random) == 0) {ze_random <- 0 }
  zc_random <- ncol(data_read$data_raw$covariates_missing_costs_random$Intervention)
  if(length(zc_random) == 0) {zc_random <- 0 }
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
  clus1_me <- data_read$data_raw$clus_me$Control
  clus2_me <- data_read$data_raw$clus_me$Intervention
  clus1_mc <- data_read$data_raw$clus_mc$Control
  clus2_mc <- data_read$data_raw$clus_mc$Intervention
  n1_clus_e <- length(unique(clus1_e))
  n2_clus_e <- length(unique(clus2_e))
  n1_clus_c <- length(unique(clus1_c))
  n2_clus_c <- length(unique(clus2_c))
  n1_clus_me <- length(unique(clus1_me))
  n2_clus_me <- length(unique(clus2_me))
  n1_clus_mc <- length(unique(clus1_mc))
  n2_clus_mc <- length(unique(clus2_mc))
  formula.me_fixed <- all.vars(data_read$model_formula$mf_model.me_fixed)
  formula.me.length_fixed <- length(formula.me_fixed)
  formula.mc_fixed <- all.vars(data_read$model_formula$mf_model.mc_fixed)
  formula.mc.length_fixed <- length(formula.mc_fixed)
  if(length(which(is.na(c(eff1, eff2)))) == 0 & length(which(is.na(c(cost1, cost2)))) == 0) {
    stop("At leat one missing value is required in either the effects or costs variables")
  }
  if(length(which(is.na(c(eff1, eff2)))) == 0) {
    if(formula.me.length_fixed != 1 | formula.me_fixed[1] != "me") {
      stop("At least one missing effect value is required to specify the formula of the missingness effect model")
    }
  }
  if(length(which(is.na(c(cost1, cost2)))) == 0) {
    if(formula.mc.length_fixed != 1 | formula.mc_fixed[1] != "mc") {
      stop("At least one missing cost value is required to specify the formula of the missingness cost model")
    }
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
  Z1_e_fixed <- as.matrix(data_read$data_raw$covariates_missing_effects_fixed$Control)
  Z2_e_fixed <- as.matrix(data_read$data_raw$covariates_missing_effects_fixed$Intervention)
  if(ze_fixed == 1) {
    Z1_e_fixed <- as.vector(Z1_e_fixed)
    Z2_e_fixed <- as.vector(Z2_e_fixed)
  }
  mean_z_e1_fixed <- as.vector(data_read$data_raw$mean_cov_missing_effects_fixed$Control)
  mean_z_e2_fixed <- as.vector(data_read$data_raw$mean_cov_missing_effects_fixed$Intervention)
  Z1_c_fixed <- as.matrix(data_read$data_raw$covariates_missing_costs_fixed$Control)
  Z2_c_fixed <- as.matrix(data_read$data_raw$covariates_missing_costs_fixed$Intervention)
  if(zc_fixed == 1) {
    Z1_c_fixed <- as.vector(Z1_c_fixed)
    Z2_c_fixed <- as.vector(Z2_c_fixed)
  }
  mean_z_c1_fixed <- as.vector(data_read$data_raw$mean_cov_missing_costs_fixed$Control)
  mean_z_c2_fixed <- as.vector(data_read$data_raw$mean_cov_missing_costs_fixed$Intervention)
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
  if(length(model_me_random) != 0 & ze_random != 0) {
    Z1_e_random <- as.matrix(data_read$data_raw$covariates_missing_effects_random$Control)
    Z2_e_random <- as.matrix(data_read$data_raw$covariates_missing_effects_random$Intervention)
    if(ze_random == 1) {
      Z1_e_random <- as.vector(Z1_e_random)
      Z2_e_random <- as.vector(Z2_e_random)
    }
    mean_z_e1_random <- as.vector(data_read$data_raw$mean_cov_missing_effects_random$Control)
    mean_z_e2_random <- as.vector(data_read$data_raw$mean_cov_missing_effects_random$Intervention)
  } else {
    Z1_e_random <- Z2_e_random <- NULL
    mean_z_e1_random <- mean_z_e2_random <- NULL
  }
  if(length(model_mc_random) != 0 & zc_random != 0) {
    Z1_c_random <- as.matrix(data_read$data_raw$covariates_missing_costs_random$Control)
    Z2_c_random <- as.matrix(data_read$data_raw$covariates_missing_costs_random$Intervention)
    if(zc_random == 1) {
      Z1_c_random <- as.vector(Z1_c_random)
      Z2_c_random <- as.vector(Z2_c_random)
    }
    mean_z_c1_random <- as.vector(data_read$data_raw$mean_cov_missing_costs_random$Control)
    mean_z_c2_random <- as.vector(data_read$data_raw$mean_cov_missing_costs_random$Intervention)
  } else {
    Z1_c_random <- Z2_c_random <- NULL
    mean_z_c1_random <- mean_z_c2_random <- NULL
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
    par_prior_fixed <- c("alpha0.prior", "beta0.prior", "sigma.prior.e", "sigma.prior.c", "gamma.prior.e", "gamma.prior.c", 
                         "alpha.prior", "beta.prior", "gamma0.prior.e", "gamma0.prior.c", "delta.prior.e", "delta.prior.c", "beta_f.prior")
    par_prior_random <- c("mu.a0.prior", "mu.b0.prior", "mu.g.prior.e", "mu.g.prior.c", "mu.a.prior", "mu.b.prior", "mu.g0.prior.e", "mu.g0.prior.c", "mu.d.prior.e", "mu.d.prior.c", "mu.b_f.prior",
                          "s.a0.prior", "s.b0.prior", "s.g.prior.e", "s.g.prior.c", "s.a.prior", "s.b.prior", "s.g0.prior.e", "s.g0.prior.c", "s.d.prior.e", "s.d.prior.c", "s.b_f.prior")
    stop_mes <- "priors can be assigned only using specific character names depending on the type of model assumed. Type ''help(selection)'' for more details"
    if(!all(names(list_check_vector) %in% c(par_prior_fixed, par_prior_random) == TRUE)) { stop(stop_mes) }
    if(is.vector(X1_e_fixed) == TRUE & identical(X1_e_fixed, rep(1, N1))) {
      if("alpha.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(is.vector(X1_c_fixed) == TRUE & identical(X1_c_fixed, rep(1, N1))) {
      if("beta.prior" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(is.vector(Z1_e_fixed) == TRUE & identical(Z1_e_fixed, rep(1, N1))) {
      if("gamma.prior.e" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(is.vector(Z1_c_fixed) == TRUE & identical(Z1_c_fixed, rep(1, N1))) {
      if("gamma.prior.c" %in% names(list_check_vector)) { stop(stop_mes) }
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
    if(length(model_me_random) == 0 | ze_random == 0) {
      if("mu.g0.prior.e" %in% names(list_check_vector) | "mu.g.prior.e" %in% names(list_check_vector) |
         "s.g0.prior.e" %in% names(list_check_vector) | "s.g.prior.e" %in% names(list_check_vector)) { stop(stop_mes)}
    } else if(length(model_me_random) != 0 & ze_random == 1) { 
      if("mu.g.prior.e" %in% names(list_check_vector) | "s.g.prior.e" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(length(model_mc_random) == 0 | zc_random == 0) {
      if("mu.g0.prior.c" %in% names(list_check_vector) | "s.g0.prior.c" %in% names(list_check_vector) | 
         "mu.g.prior.c" %in% names(list_check_vector) | "s.g.prior.c" %in% names(list_check_vector)) { stop(stop_mes)}
    } else if(length(model_mc_random) != 0 & zc_random == 1) { 
      if("mu.g.prior.c" %in% names(list_check_vector) | "s.g.prior.c" %in% names(list_check_vector)) { stop(stop_mes) }
    }
    if(type == "MAR") {
      if("delta.prior.c" %in% names(list_check_vector) | "delta.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
      if("mu.d.prior.c" %in% names(list_check_vector) | "s.d.prior.c" %in% names(list_check_vector) | 
         "mu.d.prior.e" %in% names(list_check_vector) | "s.d.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
    } 
    if(type == "MNAR_eff") {
      if("delta.prior.c" %in% names(list_check_vector)) {stop(stop_mes) }
      if("mu.d.prior.c" %in% names(list_check_vector) | "s.d.prior.c" %in% names(list_check_vector)) {stop(stop_mes) }
      if(length(model_me_random) == 0 & "mu.d.prior.e" %in% names(list_check_vector) | 
         length(model_me_random) == 0 & "s.d.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
      if(length(model_me_random) != 0 & !("e" %in% model_me_random) & "mu.d.prior.e" %in% names(list_check_vector) | 
         length(model_me_random) != 0 & !("e" %in% model_me_random) & "s.d.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
    } 
    if(type == "MNAR_cost") {
      if("delta.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
      if("mu.d.prior.e" %in% names(list_check_vector) | "s.d.prior.e" %in% names(list_check_vector)) {stop(stop_mes) }
      if(length(model_mc_random) == 0 & "mu.d.prior.c" %in% names(list_check_vector) |
         length(model_mc_random) == 0 & "s.d.prior.c" %in% names(list_check_vector)) {stop(stop_mes) }
      if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & "mu.d.prior.c" %in% names(list_check_vector) | 
         length(model_mc_random) != 0 & !("c" %in% model_mc_random) & "s.d.prior.c" %in% names(list_check_vector)) {stop(stop_mes) }
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
  if(length(model_mc_random) == 1) {
    if(model_mc_random == "c") {is_mc_random_c <- TRUE} else {is_mc_random_c <- FALSE}
  } else {is_mc_random_c <- FALSE }
  if(length(model_mc_random) == 2) {
    if(all(model_mc_random == c("1", "c")) == TRUE) {is_int_mc_random_c <- TRUE} else {is_int_mc_random_c <- FALSE}
  } else {is_int_mc_random_c <- FALSE }
  if(length(model_me_random) == 1) {
    if(model_me_random == "e") {is_me_random_e <- TRUE} else {is_me_random_e <- FALSE}
  } else {is_me_random_e <- FALSE }
  if(length(model_me_random) == 2) {
    if(all(model_me_random == c("1", "e")) == TRUE) {is_int_me_random_e <- TRUE} else {is_int_me_random_e <- FALSE}
  } else {is_int_me_random_e <- FALSE }
  if(exists("sigma.prior.e", where = prior)) {sigma.prior.e = prior$sigma.prior.e} else {sigma.prior.e = NULL }
  if(exists("sigma.prior.c", where = prior)) {sigma.prior.c = prior$sigma.prior.c} else {sigma.prior.c = NULL }
  if(exists("alpha0.prior", where = prior)) {alpha0.prior = prior$alpha0.prior} else {alpha0.prior = NULL }
  if(exists("beta0.prior", where = prior)) {beta0.prior = prior$beta0.prior} else {beta0.prior = NULL }
  if(exists("alpha.prior", where = prior)) {alpha.prior = prior$alpha.prior} else {alpha.prior = NULL }
  if(exists("beta.prior", where = prior)) {beta.prior = prior$beta.prior} else {beta.prior = NULL }
  if(exists("gamma.prior.e", where = prior)) {gamma.prior.e = prior$gamma.prior.e} else {gamma.prior.e = NULL }
  if(exists("gamma.prior.c", where = prior)) {gamma.prior.c = prior$gamma.prior.c} else {gamma.prior.c = NULL }
  if(exists("gamma0.prior.e", where = prior)) {gamma0.prior.e = prior$gamma0.prior.e} else {gamma0.prior.e = NULL }
  if(exists("gamma0.prior.c", where = prior)) {gamma0.prior.c = prior$gamma0.prior.c} else {gamma0.prior.c = NULL }
  if(exists("delta.prior.e", where = prior)) {delta.prior.e = prior$delta.prior.e} else {delta.prior.e = NULL }
  if(exists("delta.prior.c", where = prior)) {delta.prior.c = prior$delta.prior.c} else {delta.prior.c = NULL }
  if(exists("beta_f.prior", where = prior)) {beta_f.prior = prior$beta_f.prior} else {beta_f.prior = NULL }
  if(exists("mu.a0.prior", where = prior)) {mu.a0.prior = prior$mu.a0.prior} else {mu.a0.prior = NULL }
  if(exists("s.a0.prior", where = prior)) {s.a0.prior = prior$s.a0.prior} else {s.a0.prior = NULL }
  if(exists("mu.b0.prior", where = prior)) {mu.b0.prior = prior$mu.b0.prior} else {mu.b0.prior = NULL }
  if(exists("s.b0.prior", where = prior)) {s.b0.prior = prior$s.b0.prior} else {s.b0.prior = NULL }
  if(exists("mu.a.prior", where = prior)) {mu.a.prior = prior$mu.a.prior} else {mu.a.prior = NULL }
  if(exists("s.a.prior", where = prior)) {s.a.prior = prior$s.a.prior} else {s.a.prior = NULL }
  if(exists("mu.b.prior", where = prior)) {mu.b.prior = prior$mu.b.prior} else {mu.b.prior = NULL }
  if(exists("s.b.prior", where = prior)) {s.b.prior = prior$s.b.prior} else {s.b.prior = NULL }
  if(exists("mu.g.prior.e", where = prior)) {mu.g.prior.e = prior$mu.g.prior.e} else {mu.g.prior.e = NULL }
  if(exists("s.g.prior.e", where = prior)) {s.g.prior.e = prior$s.g.prior.e} else {s.g.prior.e = NULL }
  if(exists("mu.g0.prior.e", where = prior)) {mu.g0.prior.e = prior$mu.g0.prior.e} else {mu.g0.prior.e = NULL }
  if(exists("s.g0.prior.e", where = prior)) {s.g0.prior.e = prior$s.g0.prior.e} else {s.g0.prior.e = NULL }
  if(exists("mu.g0.prior.c", where = prior)) {mu.g0.prior.c = prior$mu.g0.prior.c} else {mu.g0.prior.c = NULL }
  if(exists("s.g0.prior.c", where = prior)) {s.g0.prior.c = prior$s.g0.prior.c} else {s.g0.prior.c = NULL }
  if(exists("mu.g.prior.c", where = prior)) {mu.g.prior.c = prior$mu.g.prior.c} else {mu.g.prior.c = NULL }
  if(exists("s.g.prior.c", where = prior)) {s.g.prior.c = prior$s.g.prior.c} else {s.g.prior.c = NULL }
  if(exists("mu.d.prior.e", where = prior)) {mu.d.prior.e = prior$mu.d.prior.e} else {mu.d.prior.e = NULL }
  if(exists("s.d.prior.e", where = prior)) {s.d.prior.e = prior$s.d.prior.e} else {s.d.prior.e = NULL }
  if(exists("mu.d.prior.c", where = prior)) {mu.d.prior.c = prior$mu.d.prior.c} else {mu.d.prior.c = NULL }
  if(exists("s.d.prior.c", where = prior)) {s.d.prior.c = prior$s.d.prior.c} else {s.d.prior.c = NULL }
  if(exists("mu.b_f.prior", where = prior)) {mu.b_f.prior = prior$mu.b_f.prior} else {mu.b_f.prior = NULL }
  if(exists("s.b_f.prior", where = prior)) {s.b_f.prior = prior$s.b_f.prior} else {s.b_f.prior = NULL }
  data_set <- list("effects" = data_read$data_raw$raw_effects, "costs" = data_read$data_raw$raw_costs, "N in reference arm" = N1, "N in comparator arm" = N2, 
                   "N observed in reference arm" = N1_cc, "N observed in comparator arm" = N2_cc, "N missing in reference arm" = N1_mis, "N missing in comparator arm"=N2_mis, 
                   "covariates_effects_fixed" = data_read$data_raw$covariates_effects_fixed, "covariates_costs_fixed" = data_read$data_raw$covariates_costs_fixed, 
                   "covariates_missing_effects_fixed" = data_read$data_raw$covariates_missing_effects_fixed, 
                   "covariates_effects_random" = data_read$data_raw$covariates_effects_random, "covariates_costs_random" = data_read$data_raw$covariates_costs_random, 
                   "covariates_missing_effects_random" = data_read$data_raw$covariates_missing_effects_random, "missing_effects" = data_read$data_raw$missing_effects, 
                   "covariates_missing_costs_fixed" = data_read$data_raw$covariates_missing_costs_fixed, "covariates_missing_costs_random" = data_read$data_raw$covariates_missing_costs_random, 
                   "missing_costs" = data_read$data_raw$missing_costs, "clus_effects" = data_read$data_raw$clus_e, "clus_costs" = data_read$data_raw$clus_c, "clus_missing_effects" = data_read$data_raw$clus_me, "clus_missing_costs" = data_read$data_raw$clus_mc)
  model_output <- run_selection(type = type, dist_e = dist_e, dist_c = dist_c, inits = inits, ppc = ppc)
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