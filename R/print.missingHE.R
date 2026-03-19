#' Print method for the posterior results contained in the objects of class \code{missingHE}
#' 
#' Prints the summary table for the model fitted, with the estimate of the parameters and/or missing values.
#' @keywords print JAGS missing data  
#' @param x A \code{missingHE} object containing the results of the Bayesian model run using the function \code{\link{selection}}, \code{\link{pattern}},
#' \code{\link{hurdle}} or \code{\link{lmdm}}.
#' @param display Set of parameters for which posterior summaries are displayed. Choices are: fixed effects ('fixed'), random effects ('random'),
#' individual log likelihood values ('loglik'), individual conditional parameter values ('conditional') and individual imputed values ('mis').
#' @param digits Integer indicating the number of decimal places to be used for rounding (default = 3).
#' @param ... additional arguments affecting the output produced. For example: \code{param = } JAGS name assigned to the specific model parameters 
#' to be shown in the printed table. 
#' @author Andrea Gabrio
#' @export
#' @examples  
#' # For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
#' # \code{\link{hurdle}} or \code{\link{lmdm}}
#' # 
#' # 


print.missingHE <- function(x, display = "fixed", digits = 3, ...) {
  exArgs <- list(...)
  if(!display %in% c("fixed", "random", "loglik", "conditional", "mis")) {
    stop("Please select a valid input for 'display' argument")}
  if(!inherits(x, "missingHE")) { 
    stop("Only objects of class 'missingHE' can be used")}
  if(x$model_output$model$BUGSoutput$n.chains == 1){
    stop("Output is not available when 'n.chain' = 1")}
  stat_display_name <- c("mean", "sd", "2.5%", "50%", "97.5%", "Rhat", "n.eff")
  if(x$model_output$method == "SELECTION") { params_print <- c("alpha", "beta", "beta_f", 
                                                               "gamma_e", "gamma_c", "delta_e", "delta_c", "s_c", "s_e",
                                                               "p_e", "p_c", "tau_e", "tau_c", "tmu_e", "tmu_c", "ls_c", "ltau_c")}
  if(x$model_output$method == "HURDLE") { params_print <- c("alpha", "beta", "beta_f", 
                                                            "gamma_e", "gamma_c", "p_e", "p_c", "s_c", "s_e",
                                                            "tau_e", "tau_c", "tmu_e", "tmu_c", "ls_c", "ltau_c")}
  if(x$model_output$method == "PATTERN") { params_print <- c("alpha_p", "bet_pa", "beta_f_p", 
                                                            "delta_e", "delta_c", "s_c_p", "s_e_p",
                                                            "p_prob", "tau_e_p", "tau_c_p", "tmu_e", "tmu_c", "ls_c_p", "ltau_c_p")}
  if(x$model_output$method == "LMDM") { params_print <- c("alpha", "beta", "beta_f", "alpha_te", "alpha_tc", "beta_te", "beta_tc",
                                                          "gamma_e", "gamma_c", "delta_e", "delta_c", "s_c", "s_e",
                                                          "p_e", "p_c", "tau_e", "tau_c", "tmu_e", "tmu_c", "ls_c", "ltau_c")}
  print_sum_fixed <- jagsresults(x = x$model_output$model, params = params_print, invert = FALSE)
  if(x$model_output$method == "SELECTION") {
  if(!is.null(x$model_output$effects_random) | !is.null(x$model_output$costs_random) |
     !is.null(x$model_output$mis_effects_random) | !is.null(x$model_output$mis_costs_random)) {
    params_print_re <- c("a", "b", "b_f", "g_e", "g_c", "d_e", "d_c")
    print_sum_random <- jagsresults(x = x$model_output$model, params = params_print_re, invert = FALSE)
    } else {print_sum_random <- print_sum_fixed}
  }
  if(x$model_output$method == "HURDLE") {
    if(!is.null(x$model_output$effects_random) | !is.null(x$model_output$costs_random) |
       !is.null(x$model_output$str_effects_random) | !is.null(x$model_output$str_costs_random)) {
      params_print_re <- c("a", "b", "b_f", "g_e", "g_c")
      print_sum_random <- jagsresults(x = x$model_output$model, params = params_print_re, invert = FALSE)
    } else {print_sum_random <- print_sum_fixed}
  }  
  if(x$model_output$method == "PATTERN") {
    if(!is.null(x$model_output$effects_random) | !is.null(x$model_output$costs_random)) {
      params_print_re <- c("a", "b", "b_f")
      print_sum_random <- jagsresults(x = x$model_output$model, params = params_print_re, invert = FALSE)
    } else {print_sum_random <- print_sum_fixed}
  }  
  if(x$model_output$method == "LMDM") {
    if(!is.null(x$model_output$effects_random) | !is.null(x$model_output$costs_random) |
       !is.null(x$model_output$mis_effects_random) | !is.null(x$model_output$mis_costs_random)) {
      params_print_re <- c("a", "b", "b_f", "a_te", "a_tc", "b_te", "b_tc", "g_e", "g_c", "d_e", "d_c")
      print_sum_random <- jagsresults(x = x$model_output$model, params = params_print_re, invert = FALSE)
    } else {print_sum_random <- print_sum_fixed}
  }
  if(x$model_output$method == "SELECTION") { params_print_logl <- c("loglik_e", "loglik_c", "loglik_me", "loglik_mc")}
  if(x$model_output$method == "HURDLE") { params_print_logl <- c("loglik_e", "loglik_c", "loglik_se", "loglik_sc")}
  if(x$model_output$method == "PATTERN") { params_print_logl <- c("loglik_e", "loglik_c", "loglik_d")}
  if(x$model_output$method == "LMDM") { params_print_logl <- c("loglik_e", "loglik_c", "loglik_me", "loglik_mc")}
  print_sum_loglik <- jagsresults(x = x$model_output$model, params = params_print_logl, invert = FALSE)
  print_sum_conditional <- jagsresults(x = x$model_output$model, params = c("cmu_e", "cmu_c", "ctau_e", "ctau_c", "clmu_c"), invert = FALSE)
  if(x$data_format == "wide") {
    mis_index_e <- paste("eff[", x$model_output$imputed$effects$index, "]", sep = "")
    mis_index_c <- paste("cost[", x$model_output$imputed$costs$index, "]", sep = "")
    }
  if(x$data_format == "long") {
    mis_index_e_arr <- which(is.na(x$data_set$data_raw$e), arr.ind = TRUE)
    mis_index_c_arr <- which(is.na(x$data_set$data_raw$c), arr.ind = TRUE)
    mis_index_e <- paste("eff[", mis_index_e_arr[, 1], ",", mis_index_e_arr[, 2], "]", sep = "")
    mis_index_c <- paste("cost[", mis_index_c_arr[, 1], ",", mis_index_c_arr[, 2], "]", sep = "")
  }
  print_sum_mis <- jagsresults(x = x$model_output$model, params = c(mis_index_e, mis_index_c), invert = FALSE)
  if(display == "fixed") { print_sum <- print_sum_fixed}
  if(display == "random") { print_sum <- print_sum_random}
  if(display == "loglik") { print_sum <- print_sum_loglik}
  if(display == "conditional") { print_sum <- print_sum_conditional}
  if(display == "mis") { print_sum <- print_sum_mis}
  if(exists("param", where = exArgs)) { param = exArgs$param } else { param = NULL}
  if(!is.null(param)) {
    print_sum <- jagsresults(x = x$model_output$model, params = param, invert = FALSE)
  }
  print_sum <- print_sum[, stat_display_name]
  print(print_sum, digits = digits)
}
