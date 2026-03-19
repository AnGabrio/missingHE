#' Full Bayesian Models to handle missingness in Economic Evaluations (Longitudinal Missing Data Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in longitudinal outcomes under different missing data 
#' mechanism assumptions, using alternative parametric distributions for the effect and cost variables. The analysis is performed using the \code{BUGS} language, 
#' which is implemented in the software \code{JAGS} using the function \code{\link[R2jags]{jags}} The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find the longitudinal variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.me}, \code{model.mc} (model formulas for the missing effect and cost models). Among these,
#' effectiveness, cost, time and treatment indicator (only two arms) variables must always be provided and named 'e', 'c', 'time' and 'trt', respectively. 
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
#' term 'me' (missing effects) and any covariates must be provided on the right-hand side of the formula. If there are no covariates, \code{1} should be specified on the right hand side of the formula. 
#' By default, covariates are placed on the "probability" parameter for the missing effects through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.mc A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the term 'mc' (missing costs) and any covariates must be provided on the right-hand side of the formula. 
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "probability" parameter for the missing costs through a logistic-linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param prob A numeric vector of probabilities within the range (0, 1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the imputed values.
#' @param time_dep Type of dependence structure assumed between effectiveness and cost outcomes. 
#' Current choices include: autoregressive structure of order one ('AR1') - default, bivariate at each time ('biv') and independence ('none').  
#' @param n.chains Number of chains.
#' @param n.iter Number of iterations.
#' @param n.burnin Number of warmup iterations.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{JAGS}
#' will generate initial values for all the model parameters.
#' @param n.thin Thinning interval.
#' @param save_model Logical. If \code{save_model} is \code{TRUE}, a \code{txt} file containing the model code is printed
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
#'   the number of observed and missing individuals, the total number of individuals by treatment arm and the indicator vectors for the missing values}
#'   \item{model_output}{A list containing the output of a \code{JAGS} model generated from the functions \code{\link[R2jags]{jags}}, and 
#'   the posterior samples for the main parameters of the model}
#'   \item{cea}{A list containing the output of the economic evaluation performed using the function \code{\link[BCEA]{bcea}}}
#'   \item{type}{A character variable that indicate which type of missing value mechanism used to run the model, either \code{MAR} or \code{MNAR} (see details)}
#'   \item{data_format}{A character variable that indicates which type of analysis was conducted, either using a \code{wide} or \code{long} dataset}
#'   \item{time_dep}{A character variable that indicate which type of time dependence assumption was made, either \code{none} or \code{AR1}}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA JAGS missing data longitudinal models
#' @importFrom stats model.frame terms
#' @details Depending on the distributions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} and the type of missingness mechanism specified in the argument \code{type}, different models
#' are built and run in the background by the function \code{lmdm}. These models fit multinomial logistic regressions to estimate
#' the probability of a given missing data pattern \code{k} (1 = completers, 2 = intermittent, 3 = dropout) in one or both longitudinal outcomes. A simple example can be used 
#' to show how longitudinal missing data models are specified. Consider a longitudinal data set comprising a response variable \eqn{y} measured at \code{S} occasions 
#' and a set of centered covariate \eqn{X_j}, for \eqn{i = j, ..., J}. For each subject in the trial \eqn{i = 1, ..., n} and time \eqn{s = 1, ..., S} we define 
#' an indicator variable \eqn{m_i} taking value \code{k = 1} if the \eqn{i}-th individual is associated with no missing value (completer), a value \code{k = 2} 
#' for intermittent missingness over the study period, and a value \code{k = 3} for dropout.
#' This is modelled as:
#' \deqn{m_i ~ Multinomial(\pi^k_i)}
#' \deqn{\pi^k_i = \phi^k_i/\sum\phi_i}
#' \deqn{log(\phi^k_i) = \sum\gamma^k_j X_j + \delta^k y_i}
#' where
#' \itemize{
#' \item \eqn{\pi^k_i} is the individual probability of a missing value in \eqn{y} for pattern \eqn{k} at a given time.
#' \item \eqn{\gamma^k_j} represents the impact on the missing data probability in \eqn{y} of the covariate \eqn{X_j} for pattern \eqn{k} at a given time.
#' \item \eqn{\delta^k} represents the impact on the missing data probability in \eqn{y} for pattern \eqn{k} of missingness itself at a given time.
#' }
#' When \eqn{\delta = 0} the model assumes a 'MAR' mechanism, while when \eqn{\delta != 0} the mechanism is 'MNAR'. For the parameters indexing the missingness model, 
#' the default prior distributions assumed are:
#' \itemize{
#' \item \eqn{\gamma^k_j ~ Normal(0, 0.01)}
#' \item \eqn{\delta^k ~ Normal(0, 1)}
#' }
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{lmdm}, the elements of this list (see Arguments)
#' must be vectors containing the user-provided distribution name and hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names for the parameters indexing the model which are accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j} and \eqn{\beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' \item covariate parameters in the missingness model \eqn{\gamma_j} (if covariate data provided): "gamma.prior.e"(effects) and/or "gamma.prior.c"(costs)
#' \item mnar parameter \eqn{\delta}: "delta.prior.e"(effects) and/or "delta.prior.c"(costs)
#' }
#' For simplicity, here we have assumed that the set of covariates \eqn{X_j} used in the models for the effects/costs and in the 
#' model of the missing effect/cost values is the same. However, it is possible to specify different sets of covariates for each model
#' using the arguments in the function \code{lmdm} (see Arguments).
#' 
#' For each model, random effects can also be specified for each parameter by adding the term + (x | z) to each model formula, 
#' where x is the fixed regression coefficient for which also the random effects are desired and z is the clustering variable across which 
#' the random effects are specified (must be the name of a factor variable in the dataset). 
#' 
#' @author Andrea Gabrio
#' @references 
#' Mason, AJ. Gomes, M. Carpenter, J. Grieve, R. (2021). \emph{Flexible Bayesian longitudinal models for cost‐effectiveness analyses with informative missing data}. Health economics, 30(12), 3138-3158.
#' 
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
#' # Quick example to run using subset of PBS dataset
#' 
#' # Load longitudinal dataset
#' 
#' PBS.long <- PBS
#' 
#' \donttest{
#' # Run the model using the long_miss function assuming a MAR mechanism
#' # Use only 100 iterations to run a quick check
#' model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt,
#'    model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
#'    type = "MAR", n.chains = 2, n.iter = 100, time_dep = "none")
#' 
#' # Extract regression coefficient estimates
#' coef(model.long)
#' #
#'
#' # Summarise the CEA information from the model
#' summary(model.long)
#' 
#' # Further examples which take longer to run
#' model.long <- lmdm(data = PBS.long, model.eff = e ~ trt, model.cost = c ~ trt + age,
#'    model.me = me ~ 1, model.mc = mc ~ 1, dist_e = "norm", dist_c = "norm",
#'    type = "MAR", n.chains = 2, n.iter = 500, time_dep = "none")
#' 
#' # Use looic to assess model fit
#' pic.looic <- pic(model.long, criterion = "looic")
#' pic.looic
#' 
#' # Show density plots for all parameters
#' diag.den <- diagnostic(model.long, type = "denplot", param = "alpha")
#' 
#' # Plots of imputations for all effect data
#' p1 <- plot(model.long, class = "scatter", outcome = "effects", time.plot = "all")
#' 
#' # Summarise the CEA results
#' summary(model.long)
#' 
#' }
#' #
#' #


lmdm <- function(data, model.eff, model.cost, 
                 model.me = me ~ 1, model.mc = mc ~ 1, 
                 dist_e, dist_c, type, time_dep = "AR1", prob = c(0.025, 0.975), 
                 n.chains = 2, n.iter = 10000, n.burnin = floor(n.iter / 2), 
                 inits = NULL, n.thin = 1, save_model = FALSE, 
                 prior = "default", center = FALSE, ...) {
  if(is.data.frame(data) == FALSE) { stop("data must be in data frame format")}
  if(!all(c("e", "c", "trt", "time") %in% names(data)) == TRUE) {
    stop("Please rename or provide names in the data as 'e', 
         'c', 'trt' and 'time' for the effectiveness, cost, treatment and time variables")}
  if(!is.numeric(data$e) | !is.numeric(data$c) | !is.numeric(data$time)) { stop("Effectiveness, cost and time data must be numeric")}
  max_time <- max(data$time, na.rm = TRUE); min_time <- min(data$time, na.rm = TRUE)
  if(min_time != 1 | max_time < 2) { stop("Please provide time as a numeric variable with consequent values starting from 1 and ending at least at 2")}
  time.int.values <- order(unique(as.integer(data$time)))
  n_times <- rep(1:max_time)
  if(!any(time.int.values == n_times)) { stop("Please provide time as a numeric variable starting at 1 with no unit gaps up to the maximum time (minimum 2)")}
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
  if(!type %in% c("MAR", "MNAR")) { stop("Types available are 'MAR' and 'MNAR'")} 
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
  if(!time_dep %in% c("AR1", "biv", "none")) { 
    stop("Time dependence structures available are 'AR1' and 'none'")}  
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
  clusn_me <- clusn_mc <- NULL  
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
  data_read <- data_read_lmdm(data = data, model.eff = model.eff, 
                                   model.cost = model.cost, model.me = model.me,
                                   model.mc = model.mc, cov_matrix = cov_matrix, 
                                   type = type, center = center, 
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
  if(as.character(data_read$model_formula$mf_model.me_random)[3] == "1") {
    model_me_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.me_random)[3] == "1 + e") {
    model_me_random <- c("1", "e")
  } else { model_me_random <- labels(terms(data_read$model_formula$mf_model.me_random))}
  if(as.character(data_read$model_formula$mf_model.mc_random)[3] == "1") {
    model_mc_random <- c("1")
  } else if(as.character(data_read$model_formula$mf_model.mc_random)[3] == "1 + c") {
    model_mc_random <- c("1", "c")
  } else { model_mc_random <- labels(terms(data_read$model_formula$mf_model.mc_random))}
  miss_eff_assumption <- model.frame(formula = data_read$model_formula$mf_model.me_fixed, data = data_read$data_raw$data)
  miss_cost_assumption <- model.frame(formula = data_read$model_formula$mf_model.mc_fixed, data = data_read$data_raw$data)
  if("e" %in% names(miss_eff_assumption) & !"c" %in% names(miss_cost_assumption)) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'model.mc' if 'MAR' type selected")}
    type <- "MNAR_eff"
  } else if(!"e" %in% names(miss_eff_assumption) & "c" %in% names(miss_cost_assumption)) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'model.mc' if 'MAR' type selected")}
    type <- "MNAR_cost"
  } else if("e" %in% names(miss_eff_assumption) & "c" %in% names(miss_cost_assumption)) {
    if(type == "MAR") { stop("Please remove 'e' and/or 'c' from 'model.me' and/or 'model.mc' if 'MAR' type selected")}
    type <- "MNAR"
  } else if(!"e" %in% names(miss_eff_assumption) & !"c" %in% names(miss_cost_assumption)) {
    if(type == "MNAR") { stop("Please add 'e' and/or 'c' to 'model.me' and/or 'model.mc' if 'MNAR' type selected")}
    type <- "MAR"}
  n <- sum(data_read$data_raw$n)
  n_long <- sum(data_read$data_raw$n_long)
  pe_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_e, ncol)[1])
  pc_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_c, ncol)[1])
  ze_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_me, ncol)[1])
  zc_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_mc, ncol)[1])
  if(is.list(data_read$data_raw$cov_random_e)) {
    pe_random <- as.numeric(lapply(data_read$data_raw$cov_random_e, ncol)[1])
  } else { pe_random <- 0}
  if(is.list(data_read$data_raw$cov_random_c)) {
    pc_random <- as.numeric(lapply(data_read$data_raw$cov_random_c, ncol)[1])
  } else { pc_random <- 0}
  if(is.list(data_read$data_raw$cov_random_me)) {
    ze_random <- as.numeric(lapply(data_read$data_raw$cov_random_me, ncol)[1])
  } else { ze_random <- 0}
  if(is.list(data_read$data_raw$cov_random_mc)) {
    zc_random <- as.numeric(lapply(data_read$data_raw$cov_random_mc, ncol)[1])
  } else { zc_random <- 0}
  eff <- data_read$data_raw$e
  cost <- data_read$data_raw$c
  eff_long <- data_read$data_raw$e_long
  cost_long <- data_read$data_raw$c_long
  m_eff <- data_read$data_raw$count_me
  m_cost <- data_read$data_raw$count_mc
  m_eff_long <- data_read$data_raw$count_me_long
  m_cost_long <- data_read$data_raw$count_mc_long
  time <- data_read$data_raw$time_long
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
  if(!is.null(data_read$data_raw$clus_me)){ 
    clus_me <- data_read$data_raw$clusn_me
    n_clus_me <- length(data_read$data_raw$n_clus_me)
  } else { 
    clus_me <- NULL
    n_clus_me <- NULL}
  if(!is.null(data_read$data_raw$clus_mc)){ 
    clus_mc <- data_read$data_raw$clusn_mc
    n_clus_mc <- length(data_read$data_raw$n_clus_mc)
  } else { 
    clus_mc <- NULL
    n_clus_mc <- NULL}
  formula.me_fixed <- all.vars(data_read$model_formula$mf_model.me_fixed)
  formula.me.length_fixed <- length(formula.me_fixed)
  formula.mc_fixed <- all.vars(data_read$model_formula$mf_model.mc_fixed)
  formula.mc.length_fixed <- length(formula.mc_fixed)
  if(!any(is.na(eff_long)) & !any(is.na(cost_long))) {
    stop("At leat one missing value is required in either the effect or cost variables")}
  if(!any(is.na(eff_long))) {
    if(formula.me.length_fixed != 1 | formula.me_fixed[1] != "me") {
      stop("At least one missing effect value is required to specify `model.me`")}
  }
  if(!any(is.na(cost_long))) {
    if(formula.mc.length_fixed != 1 | formula.mc_fixed[1] != "mc") {
      stop("At least one missing cost value is required to specify `model.mc`")}
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
  Z_e_fixed <- data_read$data_raw$z_e_fixed
  Z_c_fixed <- data_read$data_raw$z_c_fixed
  if(ze_fixed == 1) { 
    Z_e_fixed <- as.vector(Z_e_fixed)
    mean_z_e_fixed <- mean(Z_e_fixed, na.rm = TRUE)
  } else if(ze_fixed > 1) {
    mean_z_e_fixed <- apply(Z_e_fixed, 2, mean, na.rm = TRUE)
  }
  if(zc_fixed == 1) { 
    Z_c_fixed <- as.vector(Z_c_fixed)
    mean_z_c_fixed <- mean(Z_c_fixed, na.rm = TRUE)
  } else if(zc_fixed > 1) {
    mean_z_c_fixed <- apply(Z_c_fixed, 2, mean, na.rm = TRUE)
  }
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
  if(time_dep %in% c("AR1", "biv")) { 
    ind_time_fixed = FALSE
  } else if(time_dep == "none") {
    ind_time_fixed = TRUE}
  if(ind_fixed) {
    ind_time_fixed = TRUE}
  if(!ind_fixed & "e" %in% model_c_random) {
    if(time_dep %in% c("AR1", "biv")) { ind_time_fixed = FALSE}
    if(time_dep == "none") { stop("Selection of 'none' dependence does not allow inclusion of the effects in the cost model.")}
    }
  if(time_dep %in% c("AR1", "biv") & ind_fixed) {
    stop("Exclusion of the effects from the cost model does not allow time dependence.")
  }
  if(time_dep %in% c("AR1", "biv") & "e" %in% model_c_random & length(model_e_random) == 0) {
    stop("Exclusion of the random effects from the effect model does not allow random effects time dependence.")
  }
  if(length(model_me_random) != 0 & ze_random != 0) {
    Z_e_random <- data_read$data_raw$z_e_random
    if(ze_random == 1) { 
      Z_e_random <- as.vector(Z_e_random)
      mean_z_e_random <- mean(Z_e_random, na.rm = TRUE)
    } else if(ze_random > 1) {
      mean_z_e_random <- apply(Z_e_random, 2, mean, na.rm = TRUE)}
  } else { Z_e_random <- mean_z_e_random <- NULL}
  if(length(model_mc_random) != 0 & zc_random != 0) {
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
                         "alpha.prior", "beta.prior", "delta.prior.e", "delta.prior.c", "beta_f.prior",
                         "alpha_te.prior", "alpha_tc.prior", "beta_te.prior", "beta_tc.prior")
    par_prior_random <- c("mu.g.prior.e", "mu.g.prior.c", "mu.a.prior", "mu.b.prior", 
                          "mu.d.prior.e", "mu.d.prior.c", "mu.b_f.prior", "s.g.prior.e", "s.g.prior.c", 
                          "s.a.prior", "s.b.prior", "s.d.prior.e", "s.d.prior.c", "s.b_f.prior",
                          "mu.a_te.prior", "s.a_te.prior", "mu.a_tc.prior", "s.a_tc.prior", "mu.b_te.prior", 
                          "s.b_te.prior", "mu.b_tc.prior", "s.b_tc.prior")
    stop_mes <- "Please specify priors using required names/values and according to assumed model structure. Type ''help(selection)'' for more details."
    if(!all(names(list_check_vector) %in% c(par_prior_fixed, par_prior_random))) { stop(stop_mes)}
    if(length(model_e_random) == 0 | pe_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.a.prior", "s.a.prior"))) { stop(stop_mes)}
    }
    if(length(model_c_random) == 0 | pc_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.b.prior", "s.b.prior"))) { stop(stop_mes)}
    }
    if(length(model_me_random) == 0 | ze_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.g.prior.e", "s.g.prior.e"))) { stop(stop_mes)}
    }
    if(length(model_mc_random) == 0 | zc_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.g.prior.c", "s.g.prior.c"))) { stop(stop_mes)}
    }
    if(type == "MAR") {
      if(any(names(list_check_vector) %in% c("delta.prior.c", "delta.prior.e", "mu.d.prior.c", "s.d.prior.c", 
                                             "mu.d.prior.e", "s.d.prior.e"))) { stop(stop_mes)}
    } 
    if(type == "MNAR_eff") {
      if(any(names(list_check_vector) %in% c("delta.prior.c", "mu.d.prior.c", "s.d.prior.c"))) { stop(stop_mes)}
      if(length(model_me_random) == 0 & 
         any(names(list_check_vector) %in% c("mu.d.prior.e", "s.d.prior.e"))) { stop(stop_mes)}
      if(length(model_me_random) != 0 & !("e" %in% model_me_random) & 
         any(names(list_check_vector) %in% c("mu.d.prior.e", "s.d.prior.e"))) { stop(stop_mes)}
    } 
    if(type == "MNAR_cost") {
      if(any(names(list_check_vector) %in% c("delta.prior.e", "mu.d.prior.e", "s.d.prior.e"))) { stop(stop_mes)}
      if(length(model_mc_random) == 0 & 
         any(names(list_check_vector) %in% c("mu.d.prior.c", "s.d.prior.c"))) { stop(stop_mes)}
      if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & 
         any(names(list_check_vector) %in% c("mu.d.prior.c", "s.d.prior.c"))) { stop(stop_mes)}
    } 
    if(ind_fixed) {
      if(any(names(list_check_vector) %in% c("beta_f.prior", "b_f.prior"))) { stop(stop_mes)}
    }
    if(!ind_fixed & ind_random) {
      if(any(names(list_check_vector) %in% c("mu.b_f.prior", "s.b_f.prior"))) { stop(stop_mes)}
    }
    if(ind_time_fixed | time_dep == "biv") {
      if(any(c("beta_te.prior", "beta_tc.prior", "alpha_te.prior", "alpha_tc.prior",
               "mu.b_te.prior", "s.b_te.prior", "mu.b_tc.prior", "s.b_tc.prior",
               "mu.a_te.prior", "s.a_te.prior", "mu.a_tc.prior", "s.a_tc.prior") %in% names(list_check_vector))) { stop(stop_mes) } 
    }
    if(!ind_time_fixed & ind_random) {
      if(any(c("mu.b_te.prior", "s.b_te.prior", "mu.b_tc.prior", "s.b_tc.prior",
               "mu.a_te.prior", "s.a_te.prior", "mu.a_tc.prior", "s.a_tc.prior") %in% names(list_check_vector))) { stop(stop_mes) } 
    }
  }
  if(length(model_c_random) == 1) {
    if(model_c_random == "e") { 
      is_c_random_c <- TRUE} else { is_c_random_c <- FALSE}
  } else { is_c_random_c <- FALSE}
  if(length(model_c_random) == 2) {
    if(all(model_c_random == c("1", "e"))) { is_int_c_random_c <- TRUE} else { is_int_c_random_c <- FALSE}
  } else { is_int_c_random_c <- FALSE}
  if(length(model_mc_random) == 1) {
    if(model_mc_random == "c") { is_mc_random_c <- TRUE} else { is_mc_random_c <- FALSE}
  } else { is_mc_random_c <- FALSE}
  if(length(model_mc_random) == 2) {
    if(all(model_mc_random == c("1", "c"))) { is_int_mc_random_c <- TRUE} else { is_int_mc_random_c <- FALSE}
  } else { is_int_mc_random_c <- FALSE}
  if(length(model_me_random) == 1) {
    if(model_me_random == "e") { is_me_random_e <- TRUE} else { is_me_random_e <- FALSE}
  } else { is_me_random_e <- FALSE}
  if(length(model_me_random) == 2) {
    if(all(model_me_random == c("1", "e"))) { is_int_me_random_e <- TRUE} else { is_int_me_random_e <- FALSE}
  } else { is_int_me_random_e <- FALSE}
  if("e" %in% names(miss_eff_assumption) & "0 " %in% unlist(c(strsplit(as.character(data_read$model_formula$mf_model.me_fixed)[3], "+", fixed = TRUE)))) {
    stop("MNAR model specification does not allow non-intercept models.")}
  if("e" %in% names(miss_eff_assumption) & "0" %in% unlist(c(strsplit(as.character(fb(model.me)), " " , fixed = TRUE)))) {
    stop("MNAR model specification does not allow non-intercept models.")}
  if("c" %in% names(miss_eff_assumption) & "0 " %in% unlist(c(strsplit(as.character(data_read$model_formula$mf_model.mc_fixed)[3], "+", fixed = TRUE)))) {
    stop("MNAR model specification does not allow non-intercept models.")}
  if("c" %in% names(miss_eff_assumption) & "0" %in% unlist(c(strsplit(as.character(fb(model.mc)), " " , fixed = TRUE)))) {
    stop("MNAR model specification does not allow non-intercept models.")}
  if(exists("sigma.prior.e", where = prior)) { sigma.prior.e = prior$sigma.prior.e} else { sigma.prior.e = NULL}
  if(exists("sigma.prior.c", where = prior)) { sigma.prior.c = prior$sigma.prior.c} else { sigma.prior.c = NULL}
  if(exists("alpha.prior", where = prior)) { alpha.prior = prior$alpha.prior} else { alpha.prior = NULL}
  if(exists("beta.prior", where = prior)) { beta.prior = prior$beta.prior} else { beta.prior = NULL}
  if(exists("gamma.prior.e", where = prior)) { gamma.prior.e = prior$gamma.prior.e} else { gamma.prior.e = NULL}
  if(exists("gamma.prior.c", where = prior)) { gamma.prior.c = prior$gamma.prior.c} else { gamma.prior.c = NULL}
  if(exists("delta.prior.e", where = prior)) { delta.prior.e = prior$delta.prior.e} else { delta.prior.e = NULL}
  if(exists("delta.prior.c", where = prior)) { delta.prior.c = prior$delta.prior.c} else { delta.prior.c = NULL}
  if(exists("beta_f.prior", where = prior)) { beta_f.prior = prior$beta_f.prior} else { beta_f.prior = NULL}
  if(exists("mu.a.prior", where = prior)) { mu.a.prior = prior$mu.a.prior} else { mu.a.prior = NULL}
  if(exists("s.a.prior", where = prior)) { s.a.prior = prior$s.a.prior} else { s.a.prior = NULL}
  if(exists("mu.b.prior", where = prior)) { mu.b.prior = prior$mu.b.prior} else { mu.b.prior = NULL}
  if(exists("s.b.prior", where = prior)) { s.b.prior = prior$s.b.prior} else { s.b.prior = NULL}
  if(exists("mu.g.prior.e", where = prior)) { mu.g.prior.e = prior$mu.g.prior.e} else { mu.g.prior.e = NULL}
  if(exists("s.g.prior.e", where = prior)) { s.g.prior.e = prior$s.g.prior.e} else { s.g.prior.e = NULL}
  if(exists("mu.g.prior.c", where = prior)) { mu.g.prior.c = prior$mu.g.prior.c} else { mu.g.prior.c = NULL}
  if(exists("s.g.prior.c", where = prior)) { s.g.prior.c = prior$s.g.prior.c} else { s.g.prior.c = NULL}
  if(exists("mu.d.prior.e", where = prior)) { mu.d.prior.e = prior$mu.d.prior.e} else { mu.d.prior.e = NULL}
  if(exists("s.d.prior.e", where = prior)) { s.d.prior.e = prior$s.d.prior.e} else { s.d.prior.e = NULL}
  if(exists("mu.d.prior.c", where = prior)) { mu.d.prior.c = prior$mu.d.prior.c} else { mu.d.prior.c = NULL}
  if(exists("s.d.prior.c", where = prior)) { s.d.prior.c = prior$s.d.prior.c} else { s.d.prior.c = NULL}
  if(exists("mu.b_f.prior", where = prior)) { mu.b_f.prior = prior$mu.b_f.prior} else { mu.b_f.prior = NULL}
  if(exists("s.b_f.prior", where = prior)) { s.b_f.prior = prior$s.b_f.prior} else { s.b_f.prior = NULL}
  if(exists("beta_te.prior", where = prior)) {beta_te.prior = prior$beta_te.prior} else {beta_te.prior = NULL}
  if(exists("beta_tc.prior", where = prior)) {beta_tc.prior = prior$beta_tc.prior} else {beta_tc.prior = NULL}
  if(exists("alpha_te.prior", where = prior)) {alpha_te.prior = prior$alpha_te.prior} else {alpha_te.prior = NULL}
  if(exists("alpha_tc.prior", where = prior)) {alpha_tc.prior = prior$alpha_tc.prior} else {alpha_tc.prior = NULL}
  if(exists("mu.b_te.prior", where = prior)) {mu.b_te.prior = prior$mu.b_te.prior} else {mu.b_te.prior = NULL}
  if(exists("s.b_te.prior", where = prior)) {s.b_te.prior = prior$s.b_te.prior} else {s.b_te.prior = NULL}
  if(exists("mu.b_tc.prior", where = prior)) {mu.b_tc.prior = prior$mu.b_tc.prior} else {mu.b_tc.prior = NULL}
  if(exists("s.b_tc.prior", where = prior)) {s.b_tc.prior = prior$s.b_tc.prior} else {s.b_tc.prior = NULL}
  if(exists("mu.a_te.prior", where = prior)) {mu.a_te.prior = prior$mu.a_te.prior} else {mu.a_te.prior = NULL}
  if(exists("s.a_te.prior", where = prior)) {s.a_te.prior = prior$s.a_te.prior} else {s.a_te.prior = NULL}
  if(exists("mu.a_tc.prior", where = prior)) {mu.a_tc.prior = prior$mu.a_tc.prior} else {mu.a_tc.prior = NULL}
  if(exists("s.a_tc.prior", where = prior)) {s.a_tc.prior = prior$s.a_tc.prior} else {s.a_tc.prior = NULL}
  exArgs <- list(...)
  data_model <- list("n" = n, "n_long" = n_long, "eff" = eff, "cost" = cost, "eff_long" = eff_long, "cost_long" = cost_long,
                     "time" = time, "m_eff" = m_eff, "m_cost" = m_cost, "m_eff_long" = m_eff_long, "m_cost_long" = m_cost_long,
                     "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "Z_e_fixed" = Z_e_fixed, "Z_c_fixed" = Z_c_fixed,
                     "mean_cov_e_fixed" = mean_cov_e_fixed, "mean_cov_c_fixed" = mean_cov_c_fixed, 
                     "mean_z_e_fixed" = mean_z_e_fixed, "mean_z_c_fixed" = mean_z_c_fixed, 
                     "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "ze_fixed" = ze_fixed, "zc_fixed" = zc_fixed,
                     "X_e_random" = X_e_random, "X_c_random" = X_c_random, "mean_cov_e_random" = mean_cov_e_random,
                     "mean_cov_c_random" = mean_cov_c_random, "pe_random" = pe_random, "pc_random" = pc_random,
                     "clus_e" = clus_e, "clus_c" = clus_c, "n_clus_e" = n_clus_e, "n_clus_c" = n_clus_c, 
                     "Z_e_random" = Z_e_random, "Z_c_random" = Z_c_random, "mean_z_e_random" = mean_z_e_random, 
                     "mean_z_c_random" = mean_z_c_random, "ze_random" = ze_random, "zc_random" = zc_random,
                     "clus_me" = clus_me, "clus_mc" = clus_mc, "n_clus_me" = n_clus_me, "n_clus_mc" = n_clus_mc,
                     "trt_pos_e" = data_read$data_raw$trt_pos_e, "trt_pos_c" = data_read$data_raw$trt_pos_c,
                     "trt_pos_me" = data_read$data_raw$trt_pos_me, "trt_pos_mc" = data_read$data_raw$trt_pos_mc,
                     "n_trt" = data_read$data_raw$n, "n_trt_long" = data_read$data_raw$n_long, 
                     "trt_index" = data_read$data_raw$trt_index, "trt_index_long"= data_read$data_raw$trt_index_long, 
                     "trt_lev" = trt_lev, "efft" = data_read$data_raw$efft, "costt" = data_read$data_raw$costt, 
                     "timet" = data_read$data_raw$timet, "m_efft" = data_read$data_raw$m_efft, "m_costt" = data_read$data_raw$m_costt,
                     "count_met" = data_read$data_raw$count_met, "count_mct" = data_read$data_raw$count_mct,
                     "clus_e_lev" = data_read$data_raw$clus_e_lev, "clus_c_lev" = data_read$data_raw$clus_c_lev,
                     "clus_me_lev" = data_read$data_raw$clus_me_lev, "clus_mc_lev" = data_read$data_raw$clus_mc_lev)
  if(exists("dic", where = exArgs)) { dic = exArgs$dic} else { dic = TRUE}
  if(exists("pd", where = exArgs)) { pd = exArgs$pd} else { pd = FALSE}
  if(exists("n.iter.pd", where = exArgs)) { n.iter.pd = exArgs$n.iter.pd} else { n.iter.pd = 1000}
  if(exists("n.adapt", where = exArgs)) { n.adapt = exArgs$n.adapt} else { n.adapt = 100}
  if(exists("n.mci", where = exArgs)) { n.mci = exArgs$n.mci} else { n.mci = n.iter}
  if(!is.logical(dic)) { stop("Please provide logical value for dic argument")}  
  model_info <- list("is_c_random_c" = is_c_random_c, "is_int_c_random_c" = is_int_c_random_c,
                     "is_me_random_e" = is_me_random_e, "is_mc_random_c" = is_mc_random_c, 
                     "is_int_me_random_e" = is_int_me_random_e, "is_int_mc_random_c" = is_int_mc_random_c,
                     "ind_random" = ind_random, "ind_fixed" = ind_fixed, "ind_time_fixed" = ind_time_fixed,
                     "model_e_fixed" = model_e_fixed, "model_c_fixed" = model_c_fixed, "time_dep" = time_dep,
                     "model_e_random" = model_e_random, "model_c_random" = model_c_random, 
                     "model_me_random" = model_me_random, "model_mc_random" = model_mc_random, 
                     "dic" = dic, "pd" = pd,  "n.iter.pd" = n.iter.pd, "n.adapt" = n.adapt, 
                     "inits" = inits, "n.chains" = n.chains, "n.iter" = n.iter, 
                     "n.burnin" = n.burnin, "n.thin" = n.thin, "prior" = prior, 
                     "n.mci" = n.mci, "prob" = prob)
  model_output <- run_lmdm(data_model = data_model, type = type, 
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
  if(exists("w_e", where = exArgs)) { w_e = exArgs$w_e} else { w_e = 0.25}
  if(exists("w_c", where = exArgs)) { w_c = exArgs$w_c} else { w_c = 1}
  if(length(w_e) == 1) { e_calc <- rep(w_e, (max_time - 2) * 2 + 2)
  } else if(length(w_e) == (max_time - 2) * 2 + 2 | length(w_e) == max_time) { e_calc <- w_e
  } else { stop("Please provide valid effect weight values for computing aggregated quantities.")}
  if(length(w_c) == 1) { c_calc <- rep(w_c, max_time)
  } else if(length(w_c) == (max_time - 2) * 2 + 2 | length(w_c) == max_time) { c_calc <- w_c
  } else { stop("Please provide valid cost weight values for computing aggregated quantities.")}
  if(length(e_calc) == (max_time - 2) * 2 + 2) {
    e_mean_adj <- array(NA, dim = c(n.iter, length(trt_lev), (max_time - 2) * 2 + 2), dimnames = list(NULL, trt_lev, 1:c((max_time - 2) * 2 + 2))) 
    e_mean_adj[, , c(1, (max_time - 2) * 2 + 2)] <- model_output$mean_effects[, , c(1, max_time)] * e_calc[c(1, (max_time - 2) * 2 + 2)]
    e_mean_adj[, , seq(from = 1 + 1, to = (max_time - 2) * 2 + 2 - 1)] <- rep(model_output$mean_effects[, , seq(from = 1 + 1, to = max_time - 1)] * e_calc[seq(from = 1 + 1, to = (max_time - 2) * 2 + 2 - 1)], 2)
  } else if(length(e_calc) == max_time) {
    e_mean_adj <- array(NA, dim = c(n.iter, length(trt_lev), max_time), dimnames = list(NULL, trt_lev, 1:max_time)) 
    e_mean_adj[, , 1:max_time] <- model_output$mean_effects[, , 1:max_time] * e_calc[1:max_time] 
  }
  if(length(c_calc) == (max_time - 2) * 2 + 2) {
    c_mean_adj <- array(NA, dim = c(n.iter, length(trt_lev), (max_time - 2) * 2 + 2), dimnames = list(NULL, trt_lev, 1:c((max_time - 2) * 2 + 2))) 
    c_mean_adj[, , c(1, (max_time - 2) * 2 + 2)] <- model_output$mean_costs[, , c(1, max_time)] * c_calc[c(1, (max_time - 2) * 2 + 2)]
    c_mean_adj[, , seq(from = 1 + 1, to = (max_time - 2) * 2 + 2 - 1)] <- rep(model_output$mean_costs[, , seq(from = 1 + 1, to = max_time - 1)] * c_calc[seq(from = 1 + 1, to = (max_time - 2) * 2 + 2 - 1)], 2)
  } else if(length(c_calc) == max_time) {
    c_mean_adj <- array(NA, dim = c(n.iter, length(trt_lev), max_time), dimnames = list(NULL, trt_lev, 1:max_time)) 
    c_mean_adj[, , 1:max_time] <- model_output$mean_costs[, , 1:max_time] * c_calc[1:max_time] 
  }
  e_mean <- apply(e_mean_adj, 1:2, sum)
  c_mean <- apply(c_mean_adj, 1:2, sum)
  cea <- BCEA::bcea(e = e_mean, c = c_mean, ref = ref, interventions = interventions, 
                    Kmax = Kmax, k = wtp, plot = plot)
  format <- "long"
  res <- list(data_set = data_read, model_output = model_output, 
              cea = cea, type = type, data_format = format, time_dep = time_dep)
  class(res) <- "missingHE"
  return(res)
}