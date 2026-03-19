#' Diagnostic checks for assessing MCMC convergence of Bayesian models fitted in \code{JAGS} using the function \code{\link{selection}}, \code{\link{pattern}}, \code{\link{hurdle}} or \code{\link{lmdm}}.
#'
#' The focus is restricted to full Bayesian models in cost-effectiveness analyses based on the function \code{\link{selection}}, \code{\link{pattern}}, 
#' \code{\link{hurdle}} and \code{\link{lmdm}}, with convergence of the MCMC chains that is assessed through graphical checks of the posterior distribution of the parameters of interest,
#' Examples are density plots, trace plots, autocorrelation plots, etc. Other types of posterior checks are related to some summary MCMC statistics 
#' that are able to detect possible issues in the convergence of the algorithm, such as the potential scale reduction factor or the effective sample size.
#' Different types of diagnostic tools and statistics are used to assess model convergence using functions contained in the package \strong{ggmcmc}. 
#' Graphics and plots are managed using functions contained in the package \strong{ggplot2} and \strong{ggthemes}.
#' @keywords diagnostics MCMC convergence checks
#' @param x An object of class "missingHE" containing the posterior results of a full Bayesian model implemented using the function \code{\link{selection}}, 
#' \code{\link{pattern}}, \code{\link{hurdle}} or \code{\link{lmdm}}.
#' @param type Type of diagnostic check to be plotted for the model parameter selected. Available choices include: 'histogram' for histogram plots,
#' 'denplot' for density plots, 'traceplot' for trace plots, 'acf' for autocorrelation plots, 'running' for running mean plots,
#' 'compare' for comparing the distribution of the whole chain with only its last part, 'cross' for cross correlation plots, 'Rhat' for the potential scale reduction factor, 'geweke' for the geweke diagnostic,
#' 'pairs' for posterior correlation among the parameters,'caterpillar' for caterpillar plots.
#' @param param Name of the family of parameters to process, as given by a regular expression. For example the mean parameters 
#' for the effect and cost variables can be specified using 'mu.e' and 'mu.c', respectively. Different types of
#' models may have different parameters depending on the assumed distributions and missing data assumptions. 
#' To see a complete list of all possible parameters by types of models assumed see details.
#' @param theme Type of ggplot theme among some pre-defined themes, mostly taken from the package \strong{ggthemes}. For a full list of available themes see details.
#' @param ... Additional parameters that can be provided to manage the graphical output of \code{diagnostic}.
#' @return A \strong{ggplot} object containing the plots specified in the argument \code{type}
#' @seealso \code{\link[ggmcmc]{ggs}} \code{\link{selection}} \code{\link{pattern}} \code{\link{hurdle}} \code{\link{lmdm}}.
#' @details Depending on the types of plots specified in the argument \code{type}, the output of \code{diagnostic} can produce
#' different combinations of MCMC visual posterior checks for the family of parameters indicated in the argument \code{param}.
#' For a full list of the available plots see the description of the argument \code{type} or see the corresponding plots in the package \strong{ggmcmc}.
#' 
#' The parameters that can be assessed through \code{diagnostic} are only those included in the object \code{x} (see Arguments). Specific character names
#' must be specified in the argument \code{param} according to the specific model implemented. The available names and the parameters associated with them are:
#' \itemize{
#' \item "mu.e" the mean of the effects across treatment arms.
#' \item "mu.c" the mean of the costs across treatment arms.
#' \item "sd.e" the standard deviation of the effects.
#' \item "sd.c" the standard deviation of the costs.
#' \item "alpha" the regression coefficients for the effects.
#' \item "beta" the regression coefficients for the costs.
#' \item "beta.f" the regression coefficients for the costs related to the effects predictor.
#' \item "alpha.time" the autoregressive coefficients for the effects (only with the function \code{lmdm}).
#' \item "beta.time" the autoregressive coefficients for the costs (only with the function \code{lmdm}).
#' \item "random.alpha" the regression random effects coefficients for the effects.
#' \item "random.beta" the regression random effects coefficients for the costs.
#' \item "random.alpha.time" the autoregressive random effects coefficients for the effects (only with the function \code{lmdm}).
#' \item "random.beta.time" the autoregressive random effects coefficients for the costs (only with the function \code{lmdm}).
#' \item "p.e" the probability of missingness or structural values for the effects (only with the function \code{selection}, \code{hurdle} or \code{lmdm}).
#' \item "p.c" the probability of missingness or structural values for the costs (only with the function \code{selection}, \code{hurdle} or \code{lmdm}).
#' \item "gamma.e" the regression coefficients of missingness or structural values for the effects (only with the function \code{selection}, \code{hurdle}).
#' \item "gamma.c" the regression coefficientd of missingness or structural values for the costs (only with the function \code{selection}, \code{hurdle}).
#' \item "random.gamma.e" the random effects regression coefficients of missingness or structural values for the effects (only with the function \code{selection}, \code{hurdle} or \code{lmdm}).
#' \item "random.gamma.c" the random effects regression coefficients of missingness or structural values for the costs (only with the function \code{selection}, \code{hurdle} or \code{lmdm}).
#' \item "pattern" the probabilities of the missingness patterns (only with the function \code{pattern}).
#' \item "delta.e" the mnar parameter for the effects (only with the function \code{selection}, \code{pattern} or \code{lmdm}).
#' \item "delta.c" the mnar parameters for the costs (only with the function \code{selection}, \code{pattern} or \code{lmdm}).
#' \item "random.delta.e" the random effects mnar parameters for the effects (only with the function \code{selection} or \code{lmdm}).
#' \item "random.delta.c" the random effects mnar parameters for the costs (only with the function \code{selection} or \code{lmdm}).
#' \item "all" all available parameters stored in \code{x}.
#' }
#' When the object \code{x} is created using the function \code{pattern}, pattern-specific standard deviation ("sd.e", "sd.c") and regression coefficient 
#' parameters ("alpha", "beta") for both outcomes can be visualised. The parameters associated with a missingness mechanism can be accessed only when \code{x}
#' is created using the function \code{selection}, \code{pattern} or \code{lmdm}, while the parameters associated with the model for the structural values mechanism
#' can be accessed only when \code{x} is created using the function \code{hurdle}.
#' 
#' The argument \code{theme} allows to customise the graphical output of the plots generated by \code{diagnostic} and
#' allows to choose among a set of possible pre-defined themes taken form the package \strong{ggtheme}. For a complete list of the available character names
#' for each theme, see \strong{ggthemes}.
#' 
#' @author Andrea Gabrio
#' @references 
#' Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). \emph{Bayesian Data Analysis, 2nd edition}, CRC Press.
#'
#' Brooks, S. Gelman, A. Jones, JL. Meng, XL. (2011). \emph{Handbook of Markov Chain Monte Carlo}, CRC/Chapman and Hall.
#' @import ggplot2 mcmcr
#' @importFrom stats quantile
#' @importFrom coda varnames varnames<- as.mcmc
#' @export 
#' @examples 
#' # For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
#' # \code{\link{hurdle}} or \code{\link{lmdm}}
#' #
#' #

diagnostic <- function(x, type = "denplot", param = "all", theme = NULL, ...) {
  if(!inherits(x, "missingHE")) {
    stop("Only objects of class 'missingHE' can be used")}
  if(!isTRUE(requireNamespace("ggmcmc")) | !isTRUE(requireNamespace("mcmcr")) | !isTRUE(requireNamespace("coda"))) {
    stop("You need to install the R packages 'ggmcmc', 'coda' and 'mcmcr'. 
         Please run in your R terminal:\n install.packages('ggmcmc', 'coda', 'mcmcr')")}
  if(is.null(theme)) { theme = "classic"}
  if(!theme %in% c("classic", "base", "calc", "economist", "excel", "few", "538", "gdocs", 
                   "hc", "par", "pander", "solarized", "stata", "tufte", "wsj")) {
    stop("Please provide a valid theme type among those available. Type ''help(ggthemes)'' for more details.")}
  if(length(param) != 1) {
    stop("You can only visualise diagnostic checks for one family of parameters at a time 
         or for all parameters together by setting param = 'all'")}
  if(!type %in% c("histogram", "running", "denplot", "compare", "traceplot", 
                  "acf", "cross", "Rhat", "geweke", "caterpillar", "pairs")) {
    stop("Types of diagnostics available for use are: 'histogram', 'running', 'denplot', 'compare', 
         'traceplot', 'acf', 'cross', 'Rhat', 'geweke', 'caterpillar', 'pairs'")}
  par_sel_mar <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c",
                     "random.alpha", "random.beta", "random.beta.f", "random.gamma.e", "random.gamma.c")
  par_sel_mnar_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e" ,
                       "random.alpha", "random.beta", "random.beta.f", "random.gamma.e", "random.gamma.c", "random.delta.e")
  par_sel_mnar_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c", "delta.c",
                       "random.alpha", "random.beta", "random.beta.f", "random.gamma.e", "random.gamma.c", "random.delta.c")
  par_sel_mnar_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e", "delta.c",
                        "random.alpha", "random.beta", "random.beta.f", "random.gamma.e", "random.gamma.c", "random.delta.e", "random.delta.c")
  par_hur_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "gamma.e", "random.alpha", "random.beta", "random.beta.f", "random.gamma.e")
  par_hur_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.c", "gamma.c", "random.alpha", "random.beta", "random.beta.f", "random.gamma.c")
  par_hur_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c" , 
                  "random.alpha", "random.beta", "random.beta.f", "random.gamma.e", "random.gamma.c")
  par_pat_mar <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.p", "random.alpha", "random.beta", "random.beta.f")
  par_pat_mnar_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.p", "delta.e", "random.alpha", "random.beta", "random.beta.f")
  par_pat_mnar_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.p", "delta.c", "random.alpha", "random.beta", "random.beta.f")
  par_pat_mnar_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "p.p", "delta.e", "delta.c", "random.alpha", "random.beta", "random.beta.f")
  par_lmdm_mar <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "alpha.time", "beta.time", "p.e", "p.c", "gamma.e", "gamma.c",
                    "random.alpha", "random.beta", "random.beta.f", "random.alpha.time", "random.beta.time", "random.gamma.e", "random.gamma.c")
  par_lmdm_mnar_e <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "alpha.time", "beta.time", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e" ,
                       "random.alpha", "random.beta", "random.beta.f", "random.alpha.time", "random.beta.time", "random.gamma.e", "random.gamma.c", "random.delta.e")
  par_lmdm_mnar_c <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "alpha.time", "beta.time", "beta", "beta.f", "p.e", "p.c", "gamma.e", "gamma.c", "delta.c",
                       "random.alpha", "random.beta", "random.beta.f", "random.alpha.time", "random.beta.time", "random.gamma.e", "random.gamma.c", "random.delta.c")
  par_lmdm_mnar_ec <- c("all", "mu.e", "mu.c", "sd.e", "sd.c", "alpha", "beta", "beta.f", "alpha.time", "beta.time", "p.e", "p.c", "gamma.e", "gamma.c", "delta.e", "delta.c",
                        "random.alpha", "random.beta", "random.beta.f", "random.alpha.time", "random.beta.time", "random.gamma.e", "random.gamma.c", "random.delta.e", "random.delta.c")
  if(x$model_output$method == "SELECTION") {
    if(x$type == "MAR") {
      if(!param %in% par_sel_mar) { stop("Please provide valid parameter names based on the fitted model.")}
      }
    if(x$type == "MNAR_eff") {
      if(!param %in% par_sel_mnar_e) { stop("Please provide valid parameter names based on the fitted model.")}
      }
    if(x$type == "MNAR_cost") {
      if(!param %in% par_sel_mnar_c) { stop("Please provide valid parameter names based on the fitted model.")}
      }
    if(x$type == "MNAR") {
      if(!param %in% par_sel_mnar_ec) { stop("Please provide valid parameter names based on the fitted model.")}
    }
  }
  if(x$model_output$method == "HURDLE") {
    if(!is.null(x$model_output$mean_nostr_effects) & !is.null(x$model_output$mean_nostr_costs)) {
      if(!param %in% par_hur_ec) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(!is.null(x$model_output$mean_nostr_effects) & is.null(x$model_output$mean_nostr_costs)) {
      if(!param %in% par_hur_e) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(is.null(x$model_output$mean_nostr_effects) & !is.null(x$model_output$mean_nostr_costs)) {
      if(!param %in% par_hur_c) { stop("Please provide valid parameter names based on the fitted model.")}
    }    
  }
  if(x$model_output$method == "PATTERN") {
    if(x$type == "MAR") {
      if(!param %in% par_pat_mar) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR_eff") {
      if(!param %in% par_pat_mnar_e) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR_cost") {
      if(!param %in% par_pat_mnar_c) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR") {
      if(!param %in% par_pat_mnar_ec) { stop("Please provide valid parameter names based on the fitted model.")}
    }
  }  
  if(x$model_output$method == "LMDM") {
    if(x$type == "MAR") {
      if(!param %in% par_lmdm_mar) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR_eff") {
      if(!param %in% par_lmdm_mnar_e) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR_cost") {
      if(!param %in% par_lmdm_mnar_c) { stop("Please provide valid parameter names based on the fitted model.")}
    }
    if(x$type == "MNAR") {
      if(!param %in% par_lmdm_mnar_ec) { stop("Please provide valid parameter names based on the fitted model.")}
    }
  }
  labs <- param
  labs[pmatch("mu.e", labs)] <- "tmu_e"
  labs[pmatch("mu.c", labs)] <- "tmu_c"
  labs[pmatch("sd.e", labs)] <- "s_e"
  labs[pmatch("sd.c", labs)] <- "s_c"
  labs[pmatch("alpha", labs)] <- "alpha"
  labs[pmatch("beta", labs)] <- "beta"
  if(param == "beta.f") { labs <- "beta_f"}
  labs[pmatch("p.e", labs)] <- "p_e"
  labs[pmatch("p.c", labs)] <- "p_c"
  labs[pmatch("gamma.e", labs)] <- "gamma_e"
  labs[pmatch("gamma.c", labs)] <- "gamma_c"
  labs[pmatch("delta.e", labs)] <- "delta_e"
  labs[pmatch("delta.c", labs)] <- "delta_c"
  labs[pmatch("random.alpha", labs)] <- "a"
  labs[pmatch("random.beta", labs)] <- "b"
  labs[pmatch("random.beta.f", labs)] <- "b_f"
  labs[pmatch("random.gamma.e", labs)] <- "g_e"
  labs[pmatch("random.gamma.c", labs)] <- "g_c"
  labs[pmatch("random.delta.e", labs)] <- "d_e"
  labs[pmatch("random.delta.c", labs)] <- "d_c"
  labs[pmatch("p.p", labs)] <- "p_prob"
  if(param == "alpha.time") { labs <- "alpha_t"}
  if(param == "beta.time") { labs <- "beta_t"}
  if(param == "random.alpha.time") { labs <- "a_t"}
  if(param == "random.beta.time") { labs <- "b_t"}
  mcmc_object <- as.mcmc(x$model_output$model)
  mcmc_nchains <- mcmcr::nchains(mcmc_object)
  v_name <- varnames(mcmc_object[, , drop = FALSE])
  check_name_eff <- grepl("eff", v_name)
  check_index_eff <- which(check_name_eff, TRUE)
  check_name_cost <- grepl("cost", v_name)
  check_index_cost <- which(check_name_cost, TRUE)
  check_name_loglik <- grepl("loglik", v_name)
  check_index_loglik <- which(check_name_loglik, TRUE)
  check_name_deviance <- grepl("dev", v_name)
  check_index_deviance <- which(check_name_deviance, TRUE)
  check_name_params_cmu_e <- grepl("cmu_e", v_name)
  check_name_params_fcmu_e <- rep(FALSE, length(v_name))
  check_name_params_cmu_c <- grepl("cmu_c", v_name)
  check_name_params_fcmu_c <- rep(FALSE, length(v_name))
  if(x$model_output$dist_e %in% c("norm", "logis", "negbin")) {
    check_name_params_tau_e <- grepl("tau_e", v_name)
    check_index_params_e <- c(which(check_name_params_cmu_e, TRUE), which(check_name_params_fcmu_e, TRUE), which(check_name_params_tau_e, TRUE))
  } 
  if(x$model_output$dist_e %in% c("exp", "bern", "pois")) {
    check_index_params_e <- c(which(check_name_params_cmu_e, TRUE), which(check_name_params_fcmu_e, TRUE))
  }
  if(x$model_output$dist_e %in% c("beta", "gamma", "weib")) {
    check_name_params_tau_e <- grepl("tau_e", v_name)
    check_name_params_ftau_e <- rep(FALSE, length(v_name))    
    check_index_params_e <- c(which(check_name_params_cmu_e, TRUE), which(check_name_params_fcmu_e, TRUE), which(check_name_params_tau_e, TRUE), which(check_name_params_ftau_e, TRUE))
  }
  if(x$model_output$dist_c %in% c("norm")) {
    check_name_params_tau_c <- grepl("tau_c", v_name)
    check_index_params_c <- c(which(check_name_params_cmu_c, TRUE), which(check_name_params_fcmu_c, TRUE), which(check_name_params_tau_c, TRUE))
  }
  if(x$model_output$dist_c %in% c("norm")) {
    check_name_params_tau_c <- grepl("tau_c", v_name)
    check_index_params_c <- c(which(check_name_params_cmu_c, TRUE), which(check_name_params_fcmu_c, TRUE), which(check_name_params_tau_c, TRUE))
  }
  if(x$model_output$dist_c %in% c("gamma")) {
    check_name_params_tau_c <- grepl("tau_c", v_name)
    check_name_params_ftau_c <- rep(FALSE, length(v_name))    
    check_index_params_c <- c(which(check_name_params_cmu_c, TRUE), which(check_name_params_fcmu_c, TRUE), which(check_name_params_tau_c, TRUE), which(check_name_params_ftau_c, TRUE))
  }
  if(x$model_output$dist_c %in% c("lnorm")) {
    check_name_params_cmu_c <- grepl("lcmu_c", v_name)
    check_name_params_fcmu_c <- rep(FALSE, length(v_name))
    check_name_params_tau_c <- grepl("ltau_c", v_name)
    check_index_params_c <- c(which(check_name_params_cmu_c, TRUE), which(check_name_params_fcmu_c, TRUE), which(check_name_params_tau_c, TRUE))
  }
  check_index <- c(check_index_cost, check_index_eff, check_index_loglik, check_index_deviance)
  check_index <- c(check_index, check_index_params_e, check_index_params_c) 
  parameters <- v_name[-check_index]
  parameters <- gsub("\\[|\\]", "", parameters)
  parameters <- gsub('[[:digit:]]+', '', parameters)
  parameters <- gsub(",", '', parameters)
  parameters <- paste(unique(parameters))
  if(param == "random.alpha") {
    if(!"a" %in% parameters) { stop("no random effects for effects found")}
    if(all(c("a_te", "a_tc") %in% parameters)) {
      index_at <- which(parameters %in% c("a_te", "a_tc"))
      parameters <- parameters[-index_at]
    }
    if(all(c("alpha_te", "alpha_tc") %in% parameters)) {
      index_alphat <- which(parameters %in% c("alpha_te", "alpha_tc"))
      parameters <- parameters[-index_alphat]
    }
    if(x$model_output$method == "PATTERN") { 
      index_alpha <- which(parameters %in% c("alpha_p", "beta_p", "beta_f_p", "delta_e", "delta_c"))
      parameters <- parameters[-index_alpha]
      } else {
    index_alpha <- which(parameters %in% c("alpha", "beta", "beta_f", "gamma_e", "gamma_c", "delta_e", "delta_c"))
    parameters <- parameters[-index_alpha]}
  }
  if(param == "random.beta") {
    if(!"b" %in% parameters) { stop("no random effects for costs found")}
    if(all(c("b_te", "b_tc") %in% parameters)) {
      index_bt <- which(parameters %in% c("b_te", "b_tc"))
      parameters <- parameters[-index_bt]
    }
    if("b_f" %in% parameters) {
      index_b_f <- which(parameters == "b_f")
      parameters <- parameters[-index_b_f]
    }
    if(all(c("beta_te", "beta_tc") %in% parameters)) {
      index_betat <- which(parameters %in% c("beta_te", "beta_tc"))
      parameters <- parameters[-index_betat]
    }
    if(any(c("beta_f", "beta_f_p") %in% parameters)) {
      if(x$model_output$method == "PATTERN") { 
        index_beta <- which(parameters %in% c("beta_f_p", "p_prob"))
        parameters <- parameters[-index_beta]
      } else {
      index_betaf <- which(parameters %in% "beta_f")
      parameters <- parameters[-index_betaf]}
    }
    if(x$model_output$method == "PATTERN") { 
      index_beta <- which(parameters %in% c("beta_p", "p_prob"))
      parameters <- parameters[-index_beta]
    } else{
    index_beta <- which(parameters %in% "beta")
    parameters <- parameters[-index_beta]}
  }
  if(param == "random.beta.f") {
    if("b" %in% parameters) {
      index_b <- which(parameters %in% c("b", "b_tc", "b_te"))
      parameters <- parameters[-index_b]
    }
    if(x$model_output$method == "PATTERN") { 
      index_beta <- which(parameters %in% c("beta_p", "beta_f_p", "p_prob"))
      parameters <- parameters[-index_beta]
    } else {
    index_beta <- which(parameters %in% c("beta", "beta_f"))
    parameters <- parameters[-index_beta]}
  }
  if(param == "random.beta.f") {
    if(!"b_f" %in% parameters) { stop("no random effects for costs found related effects predictor")}
  }
  if(param == "random.gamma.e") {
    if(!"g_e" %in% parameters) { stop("no random effects for missing effects found")}
  }
  if(param == "random.gamma.c") {
    if(!"g_c" %in% parameters) { stop("no random effects for missing costs found")}
  }
  if(param == "random.delta.e") {
    if(!"d_e" %in% parameters) { stop("no random effects for mnar effects found")}
  }
  if(param == "random.delta.c") {
    if(!"d_c" %in% parameters) { stop("no random effects for mnar costs found")}
  }
  if(param == "alpha") {
    if(all(c("alpha_te", "alpha_tc") %in% parameters)) {
      index_alphat <- which(parameters %in% c("alpha_te", "alpha_tc"))
      parameters <- parameters[-index_alphat]
    }
  }
  if(param == "beta") {
    if(all(c("beta_te", "beta_tc") %in% parameters)) {
      index_betat <- which(parameters %in% c("beta_te", "beta_tc"))
      parameters <- parameters[-index_betat]
    }
    if(any(c("beta_f", "beta_f_p") %in% parameters)) {
      if(x$model_output$method == "PATTERN") { 
        index_betaf <- which(parameters %in% c("beta_f_p", "p_prob"))
        parameters <- parameters[-index_betaf]
      } else {
      index_betaf <- which(parameters %in% "beta_f")
      parameters <- parameters[-index_betaf]}
    }
  }
  if(param == "beta.f") {
    if(!any(c("beta_f", "beta_f_p") %in% parameters)) { stop("no fixed effects for costs found related effects predictor")}
  }
  if(param %in% "alpha.time") {
    if(!any(c("alpha_te", "alpha_tc") %in% parameters)) { stop("no autoregressive fixed effects for effects found")}
  }
  if(param %in% "beta.time") {
    if(!any(c("beta_te", "beta_tc") %in% parameters)) { stop("no autoregressive fixed effects for costs found")}
  }
  if(param %in% "random.alpha.time") {
    if(!any(c("a_te", "a_tc") %in% parameters)) { stop("no random effects for effects found")}
  }
  if(param == "random.beta.time") {
    if(!any(c("b_te", "b_tc") %in% parameters)) { stop("no random effects for costs found")}
  }
  mcmc_object_subset <- subset(mcmc_object, pars = parameters)
  exArgs <- list(...)
    if(param == "all") {
      varnames(mcmc_object_subset) <- paste("par", varnames(mcmc_object_subset), sep=".")
      family <- "par"
    } else { 
      family <- labs
      if(x$model_output$method == "HURDLE") {
        if(family %in% c("s_c") & !is.null(x$model_output$mean_nostr_costs)) { family <- paste(family, "\\[1\\]", sep = "")}  
        if(family %in% c("s_e") & !is.null(x$model_output$mean_nostr_effects)) { family <- paste(family, "\\[1\\]", sep = "")}  
        if(family %in% c("beta") & !is.null(x$model_output$mean_nostr_costs)) {
          if(any(varnames(mcmc_object_subset) %in% c("beta_f[1]"))) {
            if(param == "beta.f") { family <- paste(family, "_f\\[1\\]", sep = "")}
            if(!param == "beta.f") { family <- paste(family, "\\[[[:digit:]]+,1\\]", sep = "")}
          }
        }
        if(family %in% c("alpha") & !is.null(x$model_output$mean_nostr_effects)) { 
          family <- paste(family, "\\[[[:digit:]]+,1\\]", sep = "")}
      }
    }
  ggmcmc_object <- ggmcmc::ggs(mcmc_object_subset)
  if(type == "histogram") {
    if(exists("bins", where = exArgs)) { bins = exArgs$bins} else { bins = 30}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_histogram(ggmcmc_object, family = family, bins = bins, greek = greek)
  } 
  if(type == "denplot") {
    if(exists("rug", where = exArgs)) { rug = exArgs$rug} else { rug = FALSE}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_density(ggmcmc_object, family = family, rug = rug, greek = greek)
  }
  if(type == "running") {
    if(exists("original_burnin", where=exArgs)) { original_burnin = exArgs$original_burnin} else { original_burnin = TRUE}
    if(exists("original_thin", where=exArgs)) { original_thin = exArgs$original_thin} else { original_thin = TRUE}
    if(exists("greek", where=exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_running(ggmcmc_object, family = family, original_burnin = original_burnin, original_thin = original_thin, greek = greek)
  } else if(type == "compare") {
    if(exists("partial", where = exArgs)) { partial = exArgs$partial} else { partial = 0.1}
    if(exists("rug", where = exArgs)) { rug = exArgs$rug} else { rug = FALSE}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_compare_partial(ggmcmc_object, family = family, partial = partial, rug = rug, greek = greek)
  }
  if(type == "traceplot") {
    if(exists("original_burnin", where = exArgs)) { original_burnin = exArgs$original_burnin} else { original_burnin = TRUE}
    if(exists("original_thin", where = exArgs)) { original_thin = exArgs$original_thin} else { original_thin = TRUE}
    if(exists("simplify", where = exArgs)) { simplify = exArgs$simplify} else { simplify = NULL}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_traceplot(ggmcmc_object, family = family, original_burnin = original_burnin, original_thin = original_thin, 
                                        simplify = simplify, greek = greek)
  }
  if(type == "acf") {
    if(exists("nLags", where = exArgs)) { nLags = exArgs$nLags} else { nLags = 50}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_autocorrelation(ggmcmc_object, family = family, nLags = nLags, greek = greek)
  }
  if(type=="cross") {
    if(exists("absolute_scale", where = exArgs)) { absolute_scale = exArgs$absolute_scale} else { absolute_scale = TRUE}
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_crosscorrelation(ggmcmc_object, family = family, absolute_scale = absolute_scale, greek = greek)
  }
  if(type == "Rhat") {
    if(exists("scaling", where = exArgs)) { scaling = exArgs$scaling} else { scaling = 1.5}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_Rhat(ggmcmc_object, family = family, scaling = scaling, greek = greek) + ggplot2::xlab("R_hat")
  } 
  if(type == "geweke") {
    if(exists("frac1", where = exArgs)) { frac1 = exArgs$frac1} else { frac1 = 0.1}
    if(exists("frac2", where = exArgs)) { frac2 = exArgs$frac2} else { frac2 = 0.5}
    if(exists("shadow_limit", where = exArgs)) {shadow_limit = exArgs$shadow_limit} else {shadow_limit = TRUE }
    if(exists("greek", where = exArgs)) {greek = exArgs$greek} else {greek = FALSE }
    ggmcmc_out <- ggmcmc::ggs_geweke(ggmcmc_object, family = family, frac1 = frac1, frac2 = frac2, shadow_limit = shadow_limit, greek = greek)
  }
  if(type == "caterpillar") {
    if(exists("X", where = exArgs)) { X = exArgs$X} else { X = NA}
    if(exists("thick_ci", where = exArgs)) { thick_ci = exArgs$thick_ci} else { thick_ci = c(0.05, 0.95)}
    if(exists("thin_ci", where = exArgs)) { thin_ci = exArgs$thin_ci} else { thin_ci = c(0.025, 0.975)}
    if(exists("line", where = exArgs)) { line = exArgs$line} else { line = NA }
    if(exists("horizontal", where = exArgs)) { horizontal = exArgs$horizontal} else {horizontal = TRUE}
    if(exists("model_labels", where = exArgs)) { model_labels = exArgs$model_labels} else {model_labels = NULL}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out<-ggmcmc::ggs_caterpillar(ggmcmc_object, family = family, X = X, thick_ci = thick_ci, thin_ci = thin_ci, line = line, 
                                        horizontal = horizontal, model_labels = model_labels, greek = greek)
  }
  if(type == "pairs") {
    if(exists("title", where = exArgs)) { title = exArgs$title} else { title = NULL}
    if(exists("upper", where = exArgs)) { upper = exArgs$upper} else { upper = list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na")}
    if(exists("lower", where = exArgs)) { lower = exArgs$lower} else { lower = list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na")}
    if(exists("diag", where = exArgs)) { diag = exArgs$diag} else { diag = list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")}
    if(exists("xlab", where = exArgs)) { xlab = exArgs$xlab} else { xlab = NULL}
    if(exists("ylab", where = exArgs)) { ylab = exArgs$ylab} else { ylab = NULL}
    if(exists("axisLabels", where = exArgs)) { axisLabels = exArgs$axisLabels} else { axisLabels = c("show", "internal", "none")}
    if(exists("labeller", where = exArgs)) { labeller = exArgs$labeller} else { labeller = "label_value"}
    if(exists("showStrips", where = exArgs)) { showStrips = exArgs$showStrips} else { showStrips = NULL}
    if(exists("legend", where = exArgs)) { legend = exArgs$legend} else { legend = NULL}
    if(exists("greek", where = exArgs)) { greek = exArgs$greek} else { greek = FALSE}
    ggmcmc_out <- ggmcmc::ggs_pairs(ggmcmc_object, family = family,greek = greek, title = title, upper = upper, lower = lower, diag = diag, 
                                    xlab = xlab, ylab = ylab, axisLabels = axisLabels, labeller = labeller, showStrips = showStrips, legend = legend)
  }
  ggmcmc_out <- ggmcmc_out + ggplot2::theme(legend.position = "top")
  if(theme == "classic") {
    ggmcmc_out <- ggmcmc_out + ggplot2::theme(panel.grid.major = ggplot2::element_blank(), 
                                    panel.grid.minor = ggplot2::element_blank(), 
                                    panel.background = ggplot2::element_blank(), 
                                    axis.line = ggplot2::element_line(colour = "black"), 
                                    strip.background = ggplot2::element_rect(fill="white"), 
                                    strip.text.x = ggplot2::element_text(size = 12, color = "black", face = "bold"),
                                    strip.text.y = ggplot2::element_text(size = 12, color = "black", face = "bold"))
  }
  if(theme == "base") ggmcmc_out <- ggmcmc_out + ggthemes::theme_base()
  if(theme == "calc") ggmcmc_out <- ggmcmc_out + ggthemes::theme_calc()
  if(theme == "economist") ggmcmc_out <- ggmcmc_out + ggthemes::theme_economist()
  if(theme == "excel") ggmcmc_out <- ggmcmc_out + ggthemes::theme_excel()
  if(theme == "few") ggmcmc_out <- ggmcmc_out + ggthemes::theme_few()
  if(theme == "538") ggmcmc_out <- ggmcmc_out + ggthemes::theme_fivethirtyeight()
  if(theme == "gdocs") ggmcmc_out <- ggmcmc_out + ggthemes::theme_gdocs()
  if(theme == "hc") ggmcmc_out <- ggmcmc_out + ggthemes::theme_hc()
  if(theme == "par") ggmcmc_out <- ggmcmc_out + ggthemes::theme_par()
  if(theme == "solarized") ggmcmc_out <- ggmcmc_out + ggthemes::theme_solarized()
  if(theme == "pander") ggmcmc_out <- ggmcmc_out + ggthemes::theme_pander()
  if(theme == "stata") ggmcmc_out <- ggmcmc_out + ggthemes::theme_stata()
  if(theme == "tufte") ggmcmc_out <- ggmcmc_out + ggthemes::theme_tufte()
  if(theme == "wsj") ggmcmc_out <- ggmcmc_out + ggthemes::theme_wsj()
  return(print(ggmcmc_out))
}