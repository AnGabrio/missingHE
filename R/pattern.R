#' Full Bayesian Models to handle missingness in Economic Evaluations (Pattern Mixture Models)
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes under different missing data 
#' mechanism assumptions, using alternative parametric distributions for the effect and cost variables and 
#' using a pattern mixture model approach to identify the model. The analysis is performed using the \code{BUGS} language, 
#' which is implemented in the software \code{JAGS} using the function \code{\link[R2jags]{jags}} The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find the variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs). Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 'trt', respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' effectiveness outcome ('e') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economic
#' cost outcome ('c') whose name must correspond to that used in \code{data}. Any covariates in the model must be provided on the right-hand side of the formula.
#' If there are no covariates, \code{1} should be specified on the right hand side of the formula. By default, covariates are placed on the "location" parameter of the distribution through a linear model. 
#' A joint bivariate distribution for effects and costs can be specified by including 'e' on the right-hand side of the formula for the costs model.
#' Random effects can also be specified for each model parameter. See details for how these can be specified.
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param restriction type of identifying restriction to be imposed to identify the distributions of the missing data in each pattern. 
#' Available choices are: complete case restrcition ('CC') - default - or available case restriction ('AC'). 
#' @param prob A numeric vector of probabilities within the range (0, 1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the imputed values.
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
#' Consider a data set comprising a response variable \eqn{y} and a set of centered covariates \eqn{X_j}. We denote with \eqn{d_i} the patterns' indicator variable for each 
#' subject in the trial \eqn{i = 1, ..., n} such that: \eqn{d_i = 1} indicates the completers (both e and c observed), \eqn{d_i = 2} and \eqn{d_i = 3} indicate that 
#' only the costs or effects are observed, respectively, while \eqn{d_i = 4} indicates that neither of the two outcomes is observed. \eqn{d_i} is assigned a multinomial distribution, 
#' which probabilities are modelled using a Dirichlet prior. Next, the outcomes model is fitted in each pattern and parameters that cannot be identified in each pattern (d = 2, 3, 4), 
#' e.g. \eqn{mu_e[d]} and \code{mu_c[d]}, are identified using some restrictions based on the parameters estimated from other patterns. Two choices are currently available: the complete cases ('CC') or available cases ('AC').
#' For example, using the 'CC' restriction, the parameters indexing the distributions of the missing data are identified as: 
#' \deqn{mu_e[2] = \mu_e[4] = \mu_e[1] + \delta_e}
#' \deqn{mu_c[3] = \mu_c[4] = \mu_c[1] + \delta_c}
#' where
#' \itemize{
#' \item \eqn{\mu_e[1]} is the effects mean for the completers.
#' \item \eqn{\mu_c[1]} is the costs mean for the completers.
#' \item \eqn{\delta_e} is the sensitivity parameters associated with the marginal effects mean.
#' \item \eqn{\delta_c} is the sensitivity parameters associated with the marginal costs mean.
#' }
#' If the 'AC' restriction is chosen, only the parameters estimated from the observed data in pattern 2 (costs) and pattern 3 (effects) are used to identify those in the other patterns.  
#' When \eqn{\delta_e = 0} and \eqn{\delta_c = 0} the model assumes a 'MAR' mechanism. When \eqn{\delta_e != 0} and/or \eqn{\delta_c != 0} 'MNAR' departures for the 
#' effects and/or costs are explored with priors on sensitivity parameters that must be provided by the user (see Arguments). 
#' 
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{pattern}, the elements of this list (see Arguments)
#' must be vectors containing the user-provided distribution name and hyperprior values and must take specific names according to the parameters they are associated with. 
#' Specifically, the names for the parameters indexing the model which are accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item auxiliary parameters \eqn{\sigma}: "sigma.prior.e"(effects) and/or "sigma.prior.c"(costs)
#' \item covariate parameters \eqn{\alpha_j} and \eqn{\beta_j}: "alpha.prior"(effects) and/or "beta.prior"(costs)
#' \item covariate parameters in the missingness model \eqn{\gamma_j} (if covariate data provided): "gamma.prior.e"(effects) and/or "gamma.prior.c"(costs)
#' \item sensitivity parameters \eqn{\delta}: "delta.prior.e"(effects) and/or "delta.prior.c"(costs)
#' \item missingness patterns' probabilities \eqn{\pi}: "patterns.prior"
#' }
#' 
#' For each outcome model, random effects can also be specified for each parameter by adding the term + (x | z) to each model formula, 
#' where x is the fixed regression coefficient for which also the random effects are desired and z is the clustering variable across which 
#' the random effects are specified (must be the name of a factor variable in the dataset). Multiple random effects can be specified using the 
#' notation + (x1 + x2 | site) for each covariate that was included in the fixed effects formula.
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
#' # Run the model using the pattern function assuming a MAR mechanism
#' # Use only 100 iterations to run a quick check
#' model.pattern <- pattern(data = MenSS.subset,model.eff = e ~ trt, model.cost = c ~ trt,
#'    dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 100)
#' 
#' # Print the results of the JAGS model
#' print(model.pattern)
#' #
#'
#' # Use dic information criterion to assess model fit
#' pic.dic <- pic(model.pattern, criterion = "dic", cases = "cc")
#' pic.dic
#' #
#' 
#' # Extract regression coefficient estimates
#' coef(model.pattern)
#' #
#' 
#' \dontshow{
#' # Use waic information criterion to assess model fit
#' pic.waic <- pic(model.pattern, criterion = "waic", cases = "cc")
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
#' model.pattern <- pattern(data = MenSS, model.eff = e ~ trt, model.cost = c ~ trt + e,
#'    dist_e = "norm", dist_c = "norm", type = "MAR", n.chains = 2, n.iter = 500)
#' #
#' # Print results for all imputed values
#' print(model.pattern)
#' 
#' # Use looic to assess model fit
#' pic.looic <- pic(model.pattern, criterion = "looic", cases = "cc")
#' pic.looic
#' 
#' # Show density plots for mean costs parameters
#' diag.den <- diagnostic(model.pattern, type = "denplot", param = "mu.c")
#' 
#' # Plots of imputations for all data
#' p1 <- plot(model.pattern, class = "scatter", outcome = "both")
#' 
#' # Summarise the CEA results
#' summary(model.pattern)
#' 
#' }
#' #
#' #


pattern <- function(data, model.eff, model.cost,
                    dist_e, dist_c, type, restriction = "CC", prob = c(0.025, 0.975), 
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
  if(!restriction %in% c("CC", "AC")){
    stop("Only 'CC' or 'AC' types of restriction are allowed")}
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
  data_read <- data_read_pattern(data = data, model.eff = model.eff, 
                                   model.cost = model.cost, cov_matrix = cov_matrix, 
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
  n <- sum(data_read$data_raw$n)
  pe_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_e, ncol)[1])
  pc_fixed <- as.numeric(lapply(data_read$data_raw$cov_fixed_c, ncol)[1])  
  if(is.list(data_read$data_raw$cov_random_e)) {
    pe_random <- as.numeric(lapply(data_read$data_raw$cov_random_e, ncol)[1])
  } else { pe_random <- 0}
  if(is.list(data_read$data_raw$cov_random_c)) {
    pc_random <- as.numeric(lapply(data_read$data_raw$cov_random_c, ncol)[1])
  } else { pc_random <- 0}  
  eff <- data_read$data_raw$data$e
  cost <- data_read$data_raw$data$c
  n_patterns <- max(data_read$data_raw$n_obs_pat)
  m_eff <- data_read$data_raw$data$me
  m_cost <- data_read$data_raw$data$mc
  d <- data_read$data_raw$d
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
  if(!any(is.na(eff)) & !any(is.na(cost))) {
    stop("At leat one missing value is required in either the effect or cost variables")}
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
  if(anyDuplicated(names(prior)) > 0) {
    stop("You cannot provide multiple priors with the same name")}
  if(any(prior == "default")) {
    prior <- list(default = "default")
  } else if(!any(prior == "default")) {  
    list_check_vector <- lapply(prior, is.vector)
    if(!all(as.logical(list_check_vector))) {
      stop("All user-supplied priors should be in list format")}
    par_prior_fixed <- c("sigma.prior.e", "sigma.prior.c", "alpha.prior", "beta.prior", 
                         "delta.prior.e", "delta.prior.c", "beta_f.prior", "patterns.prior")
    par_prior_random <- c("mu.a.prior", "mu.b.prior", "mu.b_f.prior", "s.a.prior", "s.b.prior", "s.b_f.prior")
    stop_mes <- "Please specify priors using required names/values and according to assumed model structure. Type ''help(selection)'' for more details."
    if(!all(names(list_check_vector) %in% c(par_prior_fixed, par_prior_random))) { stop(stop_mes)}
    if(length(model_e_random) == 0 | pe_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.a.prior", "s.a.prior"))) { stop(stop_mes)}
    }
    if(length(model_c_random) == 0 | pc_random == 0) {
      if(any(names(list_check_vector) %in% c("mu.b.prior", "s.b.prior"))) { stop(stop_mes)}
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
  if(exists("delta.prior.e", where = prior)) { delta.prior.e = prior$delta.prior.e} else { delta.prior.e = NULL}
  if(exists("delta.prior.c", where = prior)) { delta.prior.c = prior$delta.prior.c} else { delta.prior.c = NULL}
  if(exists("beta_f.prior", where = prior)) { beta_f.prior = prior$beta_f.prior} else { beta_f.prior = NULL}
  if(exists("mu.a.prior", where = prior)) { mu.a.prior = prior$mu.a.prior} else { mu.a.prior = NULL}
  if(exists("s.a.prior", where = prior)) { s.a.prior = prior$s.a.prior} else { s.a.prior = NULL}
  if(exists("mu.b.prior", where = prior)) { mu.b.prior = prior$mu.b.prior} else { mu.b.prior = NULL}
  if(exists("s.b.prior", where = prior)) { s.b.prior = prior$s.b.prior} else { s.b.prior = NULL}
  if(exists("patterns.prior", where = prior)) {patterns.prior = prior$patterns.prior} else {patterns.prior = NULL }
  if(type == "MAR" & !is.null(delta.prior.e) & !is.null(delta.prior.c)) {
    stop("Under MAR no prior for any sensitivity parameters should be specified.")} 
  if(type == "MNAR" & is.null(delta.prior.e) & is.null(delta.prior.c)) {
    stop("Under MNAR the prior for at least one sensitivity parameter should be specified.")} 
  if(type == "MNAR") {
    if(!is.null(delta.prior.e) & is.null(delta.prior.c)) { type <- "MNAR_eff"}
    if(is.null(delta.prior.e) & !is.null(delta.prior.c)) { type <- "MNAR_cost"}
  }
  if(n_patterns < 2) { 
    stop("At least two missingness patterns are required to fit the model.")}
  if(restriction == "CC"){
    if(!any(d == 1, na.rm = TRUE)) {
      stop("At least some complete cases must be observed to fit the model under 'CC' restrictions.")}
  }
  if(restriction == "AC"){
    if(!any(d == 3, na.rm = TRUE) | !any(d == 2, na.rm = TRUE)) {
      stop("At least some non-complete cases in both outcomes must be observed to fit the model under 'AC' restrictions.")}
    if(!ind_fixed){
      if(!any(d == 1, na.rm = TRUE)) {
        stop("At least some complete cases must be observed when 'e' is included in the model of 'c' under 'AC' restrictions.")}
      if(n_patterns == 2){
        stop("At least three missigness patterns must be observed to fit the model when 'e' is included in the model of 'c' under 'AC' restrictions.")}
    }
  }  
  if(all(d %in% c(1, 3)) & !is.null(delta.prior.e)) { stop("Cannot introduce sensitivity parameters for effects when all effects are observed.")}
  if(all(d %in% c(1, 2)) & !is.null(delta.prior.c)) {stop("Cannot introduce sensitvity parameters for costs when all costs are observed.")}
  exArgs <- list(...)
  data_model <- list("n" = n, "eff" = eff, "cost" = cost, "m_eff" = m_eff, "m_cost" = m_cost, "d" = d, "n_patterns" = n_patterns,
                     "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "mean_cov_e_fixed" = mean_cov_e_fixed, 
                     "mean_cov_c_fixed" = mean_cov_c_fixed, "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed,
                     "X_e_random" = X_e_random, "X_c_random" = X_c_random, "mean_cov_e_random" = mean_cov_e_random,
                     "mean_cov_c_random" = mean_cov_c_random, "pe_random" = pe_random, "pc_random" = pc_random,
                     "clus_e" = clus_e, "clus_c" = clus_c, "n_clus_e" = n_clus_e, "n_clus_c" = n_clus_c, 
                     "trt_pos_e" = data_read$data_raw$trt_pos_e, "trt_pos_c" = data_read$data_raw$trt_pos_c,
                     "n_trt" = data_read$data_raw$n, "trt_index" = data_read$data_raw$trt_index, "trt_lev" = trt_lev,
                     "efft" = data_read$data_raw$efft, "costt" = data_read$data_raw$costt, 
                     "m_efft" = data_read$data_raw$m_efft, "m_costt" = data_read$data_raw$m_costt, 
                     "dt" = data_read$data_raw$dt, "n_patternst" = data_read$data_raw$n_obs_pat, 
                     "clus_e_lev" = data_read$data_raw$clus_e_lev, "clus_c_lev" = data_read$data_raw$clus_c_lev)  
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
                     "dic" = dic, "pd" = pd,  "n.iter.pd" = n.iter.pd, "n.adapt" = n.adapt, 
                     "inits" = inits, "n.chains" = n.chains, "n.iter" = n.iter, 
                     "n.burnin" = n.burnin, "n.thin" = n.thin, "prior" = prior, "patterns.prior" = patterns.prior,
                     "n.mci" = n.mci, "prob" = prob, "restriction" = restriction)
  model_output <- run_pattern(data_model = data_model, type = type, 
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