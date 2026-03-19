#' Full Bayesian Models to handle missingness in Economic Evaluations (Hurdle Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes using Hurdle models
#' under a variatey of alternative parametric distributions for the effect and cost variables. Alternative
#' assumptions about the mechanisms of the structural values are implemented using a hurdle approach. The analysis is performed using the \code{BUGS} language, 
#' which is implemented in the software \code{JAGS} using the function \code{\link[R2jags]{jags}}. The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find the variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.se}, \code{model.sc} (model formulas for the structural effect and cost models). Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 'trt', respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' effectiveness outcome ('e') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula. 
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' cost outcome ('c') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" 
#' parameter of the distribution through a linear model. A joint bivariate distribution for effects and costs can be specified by including 'e' on the right-hand side of the formula for the costs model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.se A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the 
#' term 'se'(structural effects). Any covariates in the model must be provided on the right-hand side of the formula. 
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "probability" parameter for the structural effects through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.sc A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the 
#' term 'sc'(structural costs). Any covariates in the model must be provided on the right-hand side of the formula. 
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "probability" parameter for the structural costs through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param se Structural value to be found in the effect variables defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost variables defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#' and Structural At Random (SAR).
#' @param prob A numeric vector of probabilities within the range (0, 1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the imputed values.
#' @param n.chains Number of chains.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of warmup iterations.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{JAGS}
#' will generate initial values for all the model parameters.
#' @param n.thin Thinning interval.
#' @param save_model Logical. If \code{save_model} is \code{TRUE} a \code{txt} file containing the model code is printed
#' in the current working directory.
#' @param center Logical. If \code{center} is \code{TRUE}, all covariates in both the effect and cost models will be centred.
#' @param prior A list containing the hyperprior values provided by the user. Each element of this list must be a vector
#' containing the user-provided prior distribution and parameter values and must be named with the name of the corresponding parameter. For example, the hyperprior
#' values for the standard deviation parameter for the effects can be provided using the list \code{prior = list('sigma.prior.e' = c("unif", 0, 10))}. For more information about how 
#' to provide prior hypervalues for different types of parameters and models see details. If \code{prior} is set to 'default', the default values will be used.  
#' @param ... Additional arguments that can be provided by the user. Examples are the additional arguments that can be provided to the function \code{\link[BCEA]{bcea}} to summarise the health economic evaluation results. 
#' @return An object of the class 'missingHE' containing the following elements
#' \describe{
#'   \item{data_set}{A list containing the original data set provided in \code{data} (see Arguments). Additional information is also included about, among others, 
#'   the number of observed and missing individuals, the total number of individuals by treatment arm and the indicator vectors for the structural values}
#'   \item{model_output}{A list containing the output of a \code{JAGS} model generated from the functions \code{\link[R2jags]{jags}}, and 
#'   the posterior samples for the main parameters of the model}
#'   \item{cea}{A list containing the output of the economic evaluation performed using the function \code{\link[BCEA]{bcea}}}
#'   \item{type}{A character variable that indicate which type of structural value mechanism used to run the model, either \code{SCAR} or \code{SAR} (see details)}
#'   \item{data_format}{A character variable that indicates which type of analysis was conducted, either using a \code{wide} or \code{long} dataset}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA JAGS missing data Hurdle Models 
#' @importFrom stats model.frame terms
#' @details Depending on the distributions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} and the type of structural value mechanism specified in the argument \code{type}, different hurdle models
#' are built and run in the background by the function \code{hurdle}. These are mixture models defined by two components: the first one
#' is a mass distribution at the spike, while the second is a parametric model applied to the natural range of the relevant variable.
#' Usually, a logistic regression is used to estimate the probability of incurring a "structural" value (e.g. 0 for the costs, or 1 for the
#' effects); this is then used to weigh the mean of the "non-structural" values estimated in the second component. 
#' A simple example can be used to show how hurdle models are specified. 
#' Consider a data set comprising a response variable \eqn{y} and a set of centered covariates \eqn{X_j}, for \eqn{i = j, ..., J}.Specifically, for each subject in the trial \eqn{i = 1, ..., n}
#' we define an indicator variable \eqn{s_i} taking value \code{1} if the \eqn{i}-th individual is associated with a structural value and \code{0} otherwise.
#' This is modelled as:
#' \deqn{s_i ~ Bernoulli(\pi_i)}
#' \deqn{logit(\pi_i) = \sum\gamma_j X_j}
#' where
#' \itemize{
#' \item \eqn{\pi_i} is the individual probability of a structural value in \eqn{y}.
#' \item \eqn{\gamma_j} represents the impact on the probability of a structural value in \eqn{y} of the covariate \eqn{X_j}.
#' }
#' When \eqn{\gamma_j = 0}, the model assumes a 'SCAR' mechanism, while when \eqn{\gamma_j != 0} the mechanism is 'SAR'.
#' For the parameters indexing the structural value model, the default prior distributions assumed are:
#' \itemize{
#' \item \eqn{\gamma_j ~ Normal(0, 0.01)}
#' }
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{hurdle}, the elements of this list (see Arguments)
#' must be vectors containing the user-provided distribution name and hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j, \beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' \item covariate parameters in the model of the structural values \eqn{\gamma_j} (if covariate data provided): "gamma.prior.e"(effects) and/or "gamma.prior.c"(costs)
#' } 
#' For simplicity, here we assumed that the set of covariates \eqn{X_j} used in the models for the effects/costs and in the 
#' model of the structural effect/cost values is the same. However, it is possible to specify different sets of covariates for each model
#' using the arguments in the function \code{hurdle} (see Arguments).
#' 
#' For each model, random effects can also be specified for each parameter by adding the term + (x | z) to each model formula, 
#' where x is the fixed regression coefficient for which also the random effects are desired and z is the clustering variable across which 
#' the random effects are specified (must be the name of a factor variable in the dataset). Multiple random effects can be specified using the 
#' notation + (x1 + x2 | site) for each covariate that was included in the fixed effects formula.
#' 
#' @author Andrea Gabrio
#' @references
#' Ntzoufras I. (2009). \emph{Bayesian Modelling Using WinBUGS}, John Wiley and Sons.
#' 
#' Daniels, MJ. Hogan, JW. (2008). \emph{Missing Data in Longitudinal Studies: strategies for Bayesian modelling and sensitivity analysis}, CRC/Chapman Hall.
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
#' # Run the model using the hurdle function assuming a SCAR mechanism
#' # Use only 100 iterations to run a quick check
#' model.hurdle <- hurdle(data = MenSS.subset, model.eff = e ~ trt, model.cost = c ~ trt,
#'    model.se = se ~ 1, model.sc = sc ~ 1, se = 1, sc = 0, dist_e = "norm", dist_c = "norm",
#'    type = "SCAR", n.chains = 2, n.iter = 100)
#' 
#' # Print the results of the JAGS model
#' print(model.hurdle)
#' #
#'
#' # Use dic information criterion to assess model fit
#' pic.dic <- pic(model.hurdle, criterion = "dic", cases = "cc")
#' pic.dic
#' #
#' 
#' # Extract regression coefficient estimates
#' coef(model.hurdle)
#' #
#' 
#' \dontshow{
#' # Use waic information criterion to assess model fit
#' pic.waic <- pic(model.hurdle, criterion = "waic", cases = "cc")
#' pic.waic
#' }
#'
#' # Assess model convergence using graphical tools
#' # Produce histograms of the posterior samples for the mean effects
#' diag.hist <- diagnostic(model.hurdle, type = "histogram", param = "mu.e")
#' #
#'
#' # Compare observed effect data with imputations from the model
#' # using plots (posteiror means and credible intervals)
#' p1 <- plot(model.hurdle, class = "scatter", outcome = "effects")
#' #
#'
#' # Summarise the CEA information from the model
#' summary(model.hurdle)
#' 
#' \donttest{
#' # Further examples which take longer to run
#' model.hurdle <- hurdle(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
#'    model.se = se ~ age, model.sc = sc ~ age, se = 1, sc = 0, dist_e = "norm", dist_c = "norm",
#'    type = "SAR", n.chains = 2, n.iter = 500)
#' #
#' # Print results for all imputed values
#' print(model.hurdle)
#' 
#' # Use looic to assess model fit
#' pic.looic <- pic(model.hurdle, criterion = "looic", cases = "cc")
#' pic.looic
#' 
#' # Show density plots for mean costs parameters
#' diag.den <- diagnostic(model.hurdle, type = "denplot", param = "mu.c")
#' 
#' # Plots of imputations for all data
#' p1 <- plot(model.hurdle, class = "scatter", outcome = "both")
#' 
#' # Summarise the CEA results
#' summary(model.hurdle)
#' 
#' }
#' #
#' #


hurdle <- function(data, model.eff, model.cost, 
                   model.se = se ~ 1, model.sc = sc ~ 1, se = 1, sc = 0, 
                   dist_e, dist_c, type, prob = c(0.025, 0.975), 
                   n.chains = 2, n.iter = 10000, n.burnin = floor(n.iter / 2), 
                   inits = NULL, n.thin = 1, save_model = FALSE, 
                   prior = "default", center = FALSE, ...) {
  if(is.data.frame(data) == FALSE) { stop("data must be in data frame format")}
  if(!all(c("e", "c", "trt") %in% names(data)) == TRUE) {
    stop("Please rename or provide names in the data as 'e', 
         'c' and 'trt' for the effectiveness, cost and treatment variables")}
  if(!is.numeric(data$e) | 
     !is.numeric(data$c)) { stop("Effectiveness and cost data must be numeric")}
  if(!is.factor(data$trt)) { stop("Please provide treatment `trt` as a factor variable")}
  if(!is.character(type) | !is.character(dist_e) | !is.character(dist_c)) {
    stop("Please provide valid names for 'type', 'dist_e' and 'dist_c'")}  
  dist_e <- tolower(dist_e)
  dist_c <- tolower(dist_c)
  if(dist_e %in% c("normal","gaussian")) { dist_e <- "norm"}
  if(dist_e %in% c("exponential")) { dist_e <- "exp"}
  if(dist_e %in% c("weibull")) { dist_e <- "weib"}
  if(dist_e %in% c("logistic")) { dist_e <- "logis"}
  if(dist_e %in% c("bernoulli")) { dist_e <- "bern"}
  if(dist_e %in% c("poisson")) { dist_e <- "pois"}
  if(dist_e %in% c("negative binomial", "negbinom", "nbinom", "negbinomial")) { dist_e <- "negbin"}
  if(dist_c %in% c("normal", "gaussian")) { dist_c <- "norm"}
  if(dist_c %in% c("lognormal", "lognorm", "lnormal")) { dist_c <- "lnorm"}
  if(!dist_e %in% c("norm", "beta", "exp", "weib", "logis", "bern", "pois", "negbin", "gamma") | 
     !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available are 'norm', 'beta', 'gamma', 'logis', 'exp', 'weib', 'negbin', 'pois', 'bern'  for the effects 
         and 'norm', 'gamma', 'lnorm' for the costs")}
  type <- toupper(type)
  if(!type %in% c("SCAR", "SAR")) { stop("Types available are 'SCAR' and 'SAR'")} 
  if(!is.vector(prob) | length(prob) != 2 | !is.numeric(prob) | any(prob <= 0 | prob >= 1)) {
    stop("Please provide a lower and an upper quantile for the imputed data distribution")}
  if(!any(c(length(n.chains), length(n.iter), length(n.burnin), length(n.thin)) == 1)) {
    stop("Please provide valid values for 'n.chains', 'n.iter', 'n.burnin', 'n.thin'")}
  if(!any(is.numeric(n.chains), is.numeric(n.iter), is.numeric(n.burnin), is.numeric(n.thin))) {
    stop("Please provide numeric values for 'n.chains', 'n.iter', 'n.burnin', 'n.thin'")}
  if(!any(c(n.chains, n.iter, n.burnin, n.thin) > 0) | !any(c(n.chains, n.iter, n.burnin, n.thin) %% 1 == 0)) {
    stop("Please provide valid integer values for 'n.chains', 'n.iter', 'n.burnin', 'n.thin'")}
  if(!is.logical(save_model) | !is.logical(center)) {
    stop("Please provide 'save_model', 'center' as logical values")}
  e <- as.name("e")
  c <- as.name("c")
  trt <- as.name("trt")
  cov_matrix <- subset(data, select = -c(e, c))
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  if(!inherits(model.eff, "formula") | !inherits(model.cost, "formula")) {
    stop("`model.eff` and/or `model.cost` must be formula objects")}
  fixed_e <- nobars_(model.eff)
  fixed_c <- nobars_(model.cost)
  random_e <- fb(model.eff)
  random_c <- fb(model.cost)
  clusn_e <- clusn_c <- NULL
  clusn_se <- clusn_sc <- NULL  
  if(!is.null(random_e) & length(random_e) > 1 | !is.null(random_c) & length(random_c) > 1) {
    stop("Random effects can only be included through a single expression within brackets")}
  if(!all(names(model.frame(fixed_e, data = data)) %in% c("e", names(cov_matrix))) | 
     !all(names(model.frame(fixed_c, data = data)) %in% c("c", "e", names(cov_matrix)))) {
    stop("Missing covariates cannot be included in the model")}
  if(!all(names(model.frame(fixed_e, data = data)) %in% names(data)) | 
     !all(names(model.frame(fixed_c, data = data)) %in% names(data))) {
    stop("Please provide names in the formula that correspond to those in the data")}
  if("e" %in% labels(terms(fixed_e)) | "c" %in% labels(terms(fixed_c))) {
    stop("Please remove 'e' from the right-hand side of 'model.eff' and/or 'c' from the right-hand side of 'model.cost'")}
  if(names(model.frame(fixed_e, data = data)[1]) != "e") {
    stop("Please set 'e' as the response in 'model.eff'")}
  if("c" %in% names(model.frame(fixed_e, data = data))) {
    stop("Dependence allowed only through the cost model; please remove 'c' from 'model.eff'")}
  if(names(model.frame(fixed_c, data = data)[1]) != "c") {
    stop("Please set 'c' as the response in 'model.cost'")}
  if("e" %in% labels(terms(fixed_c))) {
    if(length(grep(":e", labels(terms(fixed_c)))) != 0 | length(grep("e:", labels(terms(fixed_c)))) != 0) {
      stop("No interaction effects for 'e' is allowed")}
  }
  if(!"trt" %in% names(model.frame(fixed_c, data = data)) | !"trt" %in% names(model.frame(fixed_e, data = data))) {
    stop("Treatment indicator must be provided as covariate in 'model.eff' and/or 'model.cost'")}
  trt_lev <- levels(cov_matrix$trt)
  trt_pos <- which(trt_lev %in% levels(cov_matrix$trt))  
  data_read <- data_read_hurdle(data = data, model.eff = model.eff, 
                                   model.cost = model.cost, model.se = model.se,
                                   model.sc = model.sc, se = se, sc = sc, 
                                   cov_matrix = cov_matrix, type = type, center = center, 
                                   fixed_e = fixed_e, fixed_c = fixed_c,
                                   random_e = random_e, random_c = random_c,
                                   trt_lev = trt_lev, trt_pos = trt_pos)  
  model_e_fixed <- labels(terms(data_read$model_formula$mf_model.e_fixed))
  model_c_fixed <- labels(terms(data_read$model_formula$mf_model.c_fixed))
  if(as.character(data_read$model_formula$mf_model.e_random)[3] == "1") {
    model_e_random <- c("1")
  } else { model_e_random <- labels(terms(data_read$model_formula$mf_model.e_random))}
  if(as.character(data_read$model_formula$mf_model.c_random)[3] == "1") {
    model_c_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.c_random)[3] == "1 + e") {
    model_c_random <- c("1", "e")
  } else { model_c_random <- labels(terms(data_read$model_formula$mf_model.c_random))}
  if(as.character(data_read$model_formula$mf_model.se_random)[3] == "1") {
    model_se_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.se_random)[3] == "1 + e") {
    model_se_random <- c("1", "e")
  } else { model_se_random <- labels(terms(data_read$model_formula$mf_model.se_random))}
  if(as.character(data_read$model_formula$mf_model.sc_random)[3] == "1") {
    model_sc_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.sc_random)[3] == "1 + c") {
    model_sc_random <- c("1", "c")
  } else { model_sc_random <- labels(terms(data_read$model_formula$mf_model.sc_random))}  
  if(!is.null(se)) {
    if(as.character(data_read$model_formula$mf_model.se_random)[3] == "1") {
      model_se_random <- c("1")
    } else if(as.character(data_read$model_formula$mf_model.se_random)[3] == "0") {
      model_se_random <- NULL
    } else { model_se_random <- labels(terms(data_read$model_formula$mf_model.se_random))}
    data_read$data_raw$data$se[is.na(data_read$data_raw$data$se)] <- -999999
    str_eff_assumption <- model.frame(formula = data_read$model_formula$mf_model.se_fixed, data = data_read$data_raw$data)
    data_read$data_raw$data$se <- data_read$data_raw$se
  } else { 
    str_eff_assumption <- NULL
    model_se_random <- NULL
    }
  if(!is.null(sc)) {
    if(as.character(data_read$model_formula$mf_model.sc_random)[3] == "1") {
      model_sc_random <- c("1")
    } else if(as.character(data_read$model_formula$mf_model.sc_random)[3] == "0") {
      model_sc_random <- NULL
    } else { model_sc_random <- labels(terms(data_read$model_formula$mf_model.sc_random))}
    data_read$data_raw$data$sc[is.na(data_read$data_raw$data$sc)] <- -999999
    str_cost_assumption <- model.frame(formula = data_read$model_formula$mf_model.sc_fixed, data = data_read$data_raw$data)
    data_read$data_raw$data$sc <- data_read$data_raw$sc
  } else { 
    str_cost_assumption <- NULL
    model_sc_random <- NULL
  }
  if(length(names(str_eff_assumption)) == 0) { type_se <- "none"}
  if(length(names(str_eff_assumption)) == 1) { type_se <- "SCAR"}
  if(length(names(str_eff_assumption)) > 1) { type_se <- "SAR"}
  if(length(names(str_cost_assumption)) == 0) { type_sc <- "none"}
  if(length(names(str_cost_assumption)) == 1) { type_sc <- "SCAR"}
  if(length(names(str_cost_assumption)) > 1) { type_sc <- "SAR"}  
  if(type == "SCAR") {
    if(!type %in% c(type_se, type_sc)) {
      stop("Please remove covariates from 'model.se' and/or 'mode.sc' if 'SCAR' type selected.")}
  }
  if(type == "SAR") {
    if(!type %in% c(type_se, type_sc)) {
      stop("Please add covariates to 'model.se' and/or 'mode.sc' if 'SAR' type selected.")}
  if(type == type_se & is.null(se) | type == type_sc & is.null(sc)) {
      stop("Cannot specify a mechanism if structural effect or cost values are not provided")}
  }
  n <- sum(data_read$data_raw$n)
  pe_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_e, ncol)[1])
  pc_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_c, ncol)[1])
  if(is.list(data_read$data_raw$cov_random_e)) {
    pe_random <- as.numeric(lapply(data_read$data_raw$cov_random_e, ncol)[1])
  } else { pe_random <- 0}
  if(is.list(data_read$data_raw$cov_random_c)) {
    pc_random <- as.numeric(lapply(data_read$data_raw$cov_random_c, ncol)[1])
  } else { pc_random <- 0}
  if(!is.null(se)) {
    ze_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_se, ncol)[1])
    if(is.list(data_read$data_raw$cov_random_se)) {
      ze_random <- as.numeric(lapply(data_read$data_raw$cov_random_se, ncol)[1])
    } else { ze_random <- 0}
  } else { ze_fixed <- 0; ze_random <- 0}
  if(!is.null(sc)) {
    zc_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_sc, ncol)[1])
    if(is.list(data_read$data_raw$cov_random_sc)) {
      zc_random <- as.numeric(lapply(data_read$data_raw$cov_random_sc, ncol)[1])
    } else { zc_random <- 0}
  } else { zc_fixed <- 0; zc_random <- 0}
  eff <- data_read$data_raw$data$e
  cost <- data_read$data_raw$data$c
  m_eff <- data_read$data_raw$data$me
  m_cost <- data_read$data_raw$data$mc
  s_eff <- data_read$data_raw$se
  s_cost <- data_read$data_raw$sc
  if(!is.null(data_read$data_raw$clus_e)){ 
    clus_e <- data_read$data_raw$clusn_e
    n_clus_e <- length(data_read$data_raw$n_clus_e)
  } else { 
    clus_e <- NULL
    n_clus_e <- NULL}
  if(!is.null(data_read$data_raw$clus_c)){ 
    clus_c <- data_read$data_raw$clusn_c
    n_clus_c <- length(data_read$data_raw$n_clus_c)
  } else { 
    clus_c <- NULL
    n_clus_c <- NULL}
  if(!is.null(data_read$data_raw$clus_se)){ 
    clus_se <- data_read$data_raw$clusn_se
    n_clus_se <- length(data_read$data_raw$n_clus_se)
  } else { 
    clus_se <- NULL
    n_clus_se <- NULL}
  if(!is.null(data_read$data_raw$clus_sc)){ 
    clus_sc <- data_read$data_raw$clusn_sc
    n_clus_sc <- length(data_read$data_raw$n_clus_sc)
  } else { 
    clus_sc <- NULL
    n_clus_sc <- NULL}
  if(!is.null(se)) {
  formula.se_fixed <- all.vars(data_read$model_formula$mf_model.se_fixed)
  formula.se.length_fixed <- length(formula.se_fixed)
  } else { formula.se_fixed <- formula.se.length_fixed <- 0}
  if(!is.null(sc)) {
  formula.sc_fixed <- all.vars(data_read$model_formula$mf_model.sc_fixed)
  formula.sc.length_fixed <- length(formula.sc_fixed)
  } else { formula.sc_fixed <- formula.sc.length_fixed <- 0}
  if(!any(is.na(eff)) & !any(is.na(cost))) {
    stop("At leat one missing value is required in either the effect or cost variables")}
  if(!any(is.na(eff)) & !is.null(se)) {
    if(formula.se.length_fixed != 1 | formula.se_fixed[1] != "se") {
      stop("At least one missing effect value is required to specify `model.se`")}
  }
  if(!any(is.na(cost)) & !is.null(sc)) {
    if(formula.sc.length_fixed != 1 | formula.sc_fixed[1] != "sc") {
      stop("At least one missing cost value is required to specify `model.sc`")}
  }
  X_e_fixed <- data_read$data_raw$x_e_fixed
  X_c_fixed <- data_read$data_raw$x_c_fixed
  mean_cov_e_fixed <- apply(X_e_fixed, 2, mean, na.rm = TRUE)
  mean_cov_c_fixed <- apply(X_c_fixed, 2, mean, na.rm = TRUE)
  if(length(model_e_random) != 0 & pe_random != 0) { 
    X_e_random <- data_read$data_raw$x_e_random
    if(pe_random == 1) { 
      X_e_random <- as.vector(X_e_random)
      mean_cov_e_random <- mean(X_e_random, na.rm = TRUE)
      } else if(pe_random > 1) {
    mean_cov_e_random <- apply(X_e_random, 2, mean, na.rm = TRUE)}
    } else { X_e_random <- mean_cov_e_random <- NULL}
  if(length(model_c_random) != 0 & pc_random != 0) { 
    X_c_random <- data_read$data_raw$x_c_random
    if(pc_random == 1) { 
      X_c_random <- as.vector(X_c_random)
      mean_cov_c_random <- mean(X_c_random, na.rm = TRUE)
      } else if(pc_random > 1) {
    mean_cov_c_random <- apply(X_c_random, 2, mean, na.rm = TRUE)}
    } else { X_c_random <- mean_cov_c_random <- NULL}
  if(!is.null(se)) {
    Z_e_fixed <- data_read$data_raw$z_e_fixed
    if(ze_fixed == 1) { 
      Z_e_fixed <- as.vector(Z_e_fixed)
      mean_z_e_fixed <- mean(Z_e_fixed, na.rm = TRUE)
    } else if(ze_fixed > 1) {
      mean_z_e_fixed <- apply(Z_e_fixed, 2, mean, na.rm = TRUE)}
  } else { Z_e_fixed <- mean_z_e_fixed <- NULL}
  if(!is.null(sc)) {
    Z_c_fixed <- data_read$data_raw$z_c_fixed
    if(zc_fixed == 1) { 
      Z_c_fixed <- as.vector(Z_c_fixed)
      mean_z_c_fixed <- mean(Z_c_fixed, na.rm = TRUE)
    } else if(zc_fixed > 1) {
      mean_z_c_fixed <- apply(Z_c_fixed, 2, mean, na.rm = TRUE)}
  } else { Z_c_fixed <- mean_z_c_fixed <- NULL}
  corr_assumption_fixed <- model.frame(formula = data_read$model_formula$mf_model.c_fixed, data = data)
  if("e" %in% names(corr_assumption_fixed)) { 
    ind_fixed <- FALSE  
  } else {
    ind_fixed <- TRUE 
    ind_random <- TRUE}
  if(!ind_fixed & "e" %in% model_c_random) {
    ind_random <- FALSE
  } else if(!ind_fixed & !("e" %in% model_c_random)) {
    ind_random <- TRUE}
  if(length(model_se_random) != 0 & ze_random != 0) {
    Z_e_random <- data_read$data_raw$z_e_random
    if(ze_random == 1) { 
      Z_e_random <- as.vector(Z_e_random)
      mean_z_e_random <- mean(Z_e_random, na.rm = TRUE)
    } else if(ze_random > 1) {
      mean_z_e_random <- apply(Z_e_random, 2, mean, na.rm = TRUE)}
  } else { Z_e_random <- mean_z_e_random <- NULL}
  if(length(model_sc_random) != 0 & zc_random != 0) {
    Z_c_random <- data_read$data_raw$z_c_random
    if(zc_random == 1) { 
      Z_c_random <- as.vector(Z_c_random)
      mean_z_c_random <- mean(Z_c_random, na.rm = TRUE)
    } else if(zc_random > 1) {
      mean_z_c_random <- apply(Z_c_random, 2, mean, na.rm = TRUE)}
  } else { Z_c_random <- mean_z_c_random <- NULL}  
  if(anyDuplicated(names(prior)) > 0) {
    stop("You cannot provide multiple priors with the same name")}  
  if(any(prior == "default")) {
    prior <- list(default = "default")
  } else if(!any(prior == "default")) {
    list_check_vector <- lapply(prior, is.vector)
    if(!all(as.logical(list_check_vector))) {
      stop("All user-supplied priors should be in list format")}
    par_prior_fixed <- c("sigma.prior.e", "sigma.prior.c", "gamma.prior.e", "gamma.prior.c", 
                         "alpha.prior", "beta.prior", "beta_f.prior")
    par_prior_random <- c("mu.g.prior.e", "mu.g.prior.c", "mu.a.prior", "mu.b.prior", 
                          "mu.b_f.prior", "s.g.prior.e", "s.g.prior.c", 
                          "s.a.prior", "s.b.prior", "s.b_f.prior")
    stop_mes <- "Please specify priors using required names/values and according to assumed model structure. Type ''help(selection)'' for more details."
    if(!all(names(list_check_vector) %in% c(par_prior_fixed, par_prior_random))) { stop(stop_mes)}
    if(length(model_e_random) == 0 | pe_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.a.prior", "s.a.prior"))) { stop(stop_mes)}
    }
    if(length(model_c_random) == 0 | pc_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.b.prior", "s.b.prior"))) { stop(stop_mes)}
    }
    if(length(Z_e_fixed) == 0 | ze_fixed == 0) {
      if(any(names(list_check_vector) %in% c("gamma.prior.e"))) { stop(stop_mes)}
    }
    if(length(Z_c_fixed) == 0 | zc_fixed == 0) {
      if(any(names(list_check_vector) %in% c("gamma.prior.c"))) { stop(stop_mes)}
    }
    if(length(model_se_random) == 0 | ze_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.g.prior.e", "s.g.prior.e"))) { stop(stop_mes)}
    }
    if(length(model_sc_random) == 0 | zc_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.g.prior.c", "s.g.prior.c"))) { stop(stop_mes)}
    }
    if(ind_fixed) {
      if(any(names(list_check_vector) %in% c("beta_f.prior", "b_f.prior"))) { stop(stop_mes)}
    }
    if(!ind_fixed & ind_random) {
      if(any(names(list_check_vector) %in% c("mu.b_f.prior", "s.b_f.prior"))) { stop(stop_mes)}
    }
  }
  if(length(model_c_random) == 1) {
    if(model_c_random == "e") { 
      is_c_random_c <- TRUE} else { is_c_random_c <- FALSE}
  } else { is_c_random_c <- FALSE}
  if(length(model_c_random) == 2) {
    if(all(model_c_random == c("1", "e"))) { is_int_c_random_c <- TRUE} else { is_int_c_random_c <- FALSE}
  } else { is_int_c_random_c <- FALSE}
  if(exists("sigma.prior.e", where = prior)) { sigma.prior.e = prior$sigma.prior.e} else { sigma.prior.e = NULL}
  if(exists("sigma.prior.c", where = prior)) { sigma.prior.c = prior$sigma.prior.c} else { sigma.prior.c = NULL}
  if(exists("alpha.prior", where = prior)) { alpha.prior = prior$alpha.prior} else { alpha.prior = NULL}
  if(exists("beta.prior", where = prior)) { beta.prior = prior$beta.prior} else { beta.prior = NULL}
  if(exists("gamma.prior.e", where = prior)) { gamma.prior.e = prior$gamma.prior.e} else { gamma.prior.e = NULL}
  if(exists("gamma.prior.c", where = prior)) { gamma.prior.c = prior$gamma.prior.c} else { gamma.prior.c = NULL}
  if(exists("beta_f.prior", where = prior)) { beta_f.prior = prior$beta_f.prior} else { beta_f.prior = NULL}
  if(exists("mu.a.prior", where = prior)) { mu.a.prior = prior$mu.a.prior} else { mu.a.prior = NULL}
  if(exists("s.a.prior", where = prior)) { s.a.prior = prior$s.a.prior} else { s.a.prior = NULL}
  if(exists("mu.b.prior", where = prior)) { mu.b.prior = prior$mu.b.prior} else { mu.b.prior = NULL}
  if(exists("s.b.prior", where = prior)) { s.b.prior = prior$s.b.prior} else { s.b.prior = NULL}
  if(exists("mu.g.prior.e", where = prior)) { mu.g.prior.e = prior$mu.g.prior.e} else { mu.g.prior.e = NULL}
  if(exists("s.g.prior.e", where = prior)) { s.g.prior.e = prior$s.g.prior.e} else { s.g.prior.e = NULL}
  if(exists("mu.g.prior.c", where = prior)) { mu.g.prior.c = prior$mu.g.prior.c} else { mu.g.prior.c = NULL}
  if(exists("s.g.prior.c", where = prior)) { s.g.prior.c = prior$s.g.prior.c} else { s.g.prior.c = NULL}
  if(exists("mu.b_f.prior", where = prior)) { mu.b_f.prior = prior$mu.b_f.prior} else { mu.b_f.prior = NULL}
  if(exists("s.b_f.prior", where = prior)) { s.b_f.prior = prior$s.b_f.prior} else { s.b_f.prior = NULL}
  exArgs <- list(...)
  if(exists("se.prior", where = exArgs)) { se.prior = exArgs$se.prior} else { se.prior = 0.0000001 }
  if(exists("sc.prior", where = exArgs)) { sc.prior = exArgs$sc.prior} else { sc.prior = 0.0000001 }
  sde <- se.prior
  sdc <- sc.prior
  if(length(sde) != 1 | length(sdc) != 1) { stop("single value priors on std for structural values must be provided")}
  if(exists("s_e", where = exArgs)) { s_e = as.vector(exArgs$s_e)} else { s_e = NULL}
  if(!is.null(s_e)) {
    if(length(s_e) != length(data$e)) { stop("Please provide a valid structural value indicator vector") }
    s_eff = s_e
    if(is.null(se)) { stop("No structural values provided in the model formula")}
    data_read$data_raw$se <- s_eff
    for(i in trt_lev) {
      data_read$data_raw$s_efft[[i]] <- s_eff[data_read$data_raw$trt_index[[i]]]
      data_read$data_raw$n_ns_eff[i] <- length(na.omit(data_read$data_raw$efft[[i]])[na.omit(data_read$data_raw$efft[[i]]) != se])
    }
    data_read$data_raw$n_s_eff <- data_read$data_raw$n_obs_e - data_read$data_raw$n_ns_eff
  }
  if(exists("s_c", where = exArgs)) { s_c = as.vector(exArgs$s_c)} else { s_c = NULL}
  if(!is.null(s_c)) {
    if(length(s_c) != length(data$c)) { stop("Please provide a valid structural value indicator vector") }
    s_cost = s_c
    if(is.null(sc)) { stop("No structural values provided in the model formula")}
    data_read$data_raw$sc <- s_cost
    for(i in trt_lev) {
      data_read$data_raw$s_costt[[i]] <- s_cost[data_read$data_raw$trt_index[[i]]]
      data_read$data_raw$n_ns_cost[i] <- length(na.omit(data_read$data_raw$costt[[i]])[na.omit(data_read$data_raw$costt[[i]]) != sc])
    }
    data_read$data_raw$n_s_cost <- data_read$data_raw$n_obs_c - data_read$data_raw$n_ns_cost
  }
  data_model <- list("n" = n, "eff" = eff, "cost" = cost, "m_eff" = m_eff, "m_cost" = m_cost, 
                     "s_eff" = s_eff, "s_cost" = s_cost, "se" = se, "sc" = sc, "sde" = sde, "sdc" = sdc,
                     "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "Z_e_fixed" = Z_e_fixed, "Z_c_fixed" = Z_c_fixed,
                     "mean_cov_e_fixed" = mean_cov_e_fixed, "mean_cov_c_fixed" = mean_cov_c_fixed, 
                     "mean_z_e_fixed" = mean_z_e_fixed, "mean_z_c_fixed" = mean_z_c_fixed, 
                     "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "ze_fixed" = ze_fixed, "zc_fixed" = zc_fixed,
                     "X_e_random" = X_e_random, "X_c_random" = X_c_random, "mean_cov_e_random" = mean_cov_e_random,
                     "mean_cov_c_random" = mean_cov_c_random, "pe_random" = pe_random, "pc_random" = pc_random,
                     "clus_e" = clus_e, "clus_c" = clus_c, "n_clus_e" = n_clus_e, "n_clus_c" = n_clus_c, 
                     "Z_e_random" = Z_e_random, "Z_c_random" = Z_c_random, "mean_z_e_random" = mean_z_e_random, 
                     "mean_z_c_random" = mean_z_c_random, "ze_random" = ze_random, "zc_random" = zc_random,
                     "clus_se" = clus_se, "clus_sc" = clus_sc, "n_clus_se" = n_clus_se, "n_clus_sc" = n_clus_sc,
                     "trt_pos_e" = data_read$data_raw$trt_pos_e, "trt_pos_c" = data_read$data_raw$trt_pos_c,
                     "trt_pos_se" = data_read$data_raw$trt_pos_se, "trt_pos_sc" = data_read$data_raw$trt_pos_sc,
                     "n_trt" = data_read$data_raw$n, "trt_index" = data_read$data_raw$trt_index, "trt_lev" = trt_lev,
                     "efft" = data_read$data_raw$efft, "costt" = data_read$data_raw$costt, 
                     "m_efft" = data_read$data_raw$m_efft, "m_costt" = data_read$data_raw$m_costt,
                     "s_efft" = data_read$data_raw$s_efft, "s_costt" = data_read$data_raw$s_costt,
                     "clus_e_lev" = data_read$data_raw$clus_e_lev, "clus_c_lev" = data_read$data_raw$clus_c_lev,
                     "clus_se_lev" = data_read$data_raw$clus_se_lev, "clus_sc_lev" = data_read$data_raw$clus_sc_lev)
  if(exists("dic", where = exArgs)) { dic = exArgs$dic} else { dic = TRUE}
  if(exists("pd", where = exArgs)) { pd = exArgs$pd} else { pd = FALSE}
  if(exists("n.iter.pd", where = exArgs)) { n.iter.pd = exArgs$n.iter.pd} else { n.iter.pd = 1000}
  if(exists("n.adapt", where = exArgs)) { n.adapt = exArgs$n.adapt} else { n.adapt = 100}
  if(exists("n.mci", where = exArgs)) { n.mci = exArgs$n.mci} else { n.mci = n.iter}
  if(!is.logical(dic)) { stop("Please provide logical value for dic argument")}  
  model_info <- list("is_c_random_c" = is_c_random_c, "is_int_c_random_c" = is_int_c_random_c, 
                     "ind_random" = ind_random, "ind_fixed" = ind_fixed,
                     "model_e_fixed" = model_e_fixed, "model_c_fixed" = model_c_fixed,
                     "model_e_random" = model_e_random, "model_c_random" = model_c_random, 
                     "model_se_random" = model_se_random, "model_sc_random" = model_sc_random, 
                     "dic" = dic, "pd" = pd,  "n.iter.pd" = n.iter.pd, "n.adapt" = n.adapt, 
                     "inits" = inits, "n.chains" = n.chains, "n.iter" = n.iter, 
                     "n.burnin" = n.burnin, "n.thin" = n.thin, "prior" = prior, 
                     "n.mci" = n.mci, "prob" = prob)  
  model_output <- run_hurdle(data_model = data_model, type = type, 
                                dist_e = dist_e, dist_c = dist_c, model_info = model_info)
  if(!save_model) { unlink(model_output$filein)}
  if(exists("ref", where = exArgs)) {
    if(length(exArgs$ref) != 1 | !is.numeric(exArgs$ref)) { stop("Please provide a single numeric indicator for the reference treatment")}
    if(exArgs$ref <= 0 | !exArgs$ref %% 1 == 0) { stop("Please provide a valid indicator value for the reference treatment")}
    ref = exArgs$ref } else { ref = 1}
  if(exists("interventions", where = exArgs)) { interventions = exArgs$interventions} else { interventions = NULL}
  if(exists("Kmax", where = exArgs)) { Kmax = exArgs$Kmax} else { Kmax = 50000}
  if(exists("wtp", where = exArgs)) { wtp = exArgs$wtp} else { wtp = NULL}
  if(exists("plot", where = exArgs)) { plot = exArgs$plot} else { plot = FALSE}
  cea <- BCEA::bcea(e = model_output$mean_effects, c = model_output$mean_costs, ref = ref, 
                    interventions = interventions, Kmax = Kmax, k = wtp, plot = plot)
  format <- "wide"
  res <- list(data_set = data_read, model_output = model_output, 
              cea = cea, type = type, data_format = format)
  class(res) <- "missingHE"  
  return(res)
}