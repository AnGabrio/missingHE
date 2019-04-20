#' An internal function to execute a JAGS pattern mixture model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} function and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm') or Beta ('beta').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters
#' @param d_list a list of the number and types of patterns in the data
#' @param d1 Patterns in the control
#' @param d2 Patterns in the intervention
#' @keywords JAGS Bayesian pattern mixture models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_pattern <- function(type, dist_e, dist_c, inits, d_list, d1, d2) eval.parent(substitute( {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  if(!dist_e %in% c("norm", "beta") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta' for the effects and 'norm', 'gamma', 'lnorm' for the costs")
  }
  if(!type %in% c("MAR", "MNAR", "MNAR_eff", "MNAR_cost")) {
    stop("Types available for use are 'MAR', 'MNAR_eff', 'MNAR_cost'and 'MNAR'")
  }
  if(is.null(inits) == FALSE) {inits = inits }
  n_patterns1 <- length(unique(d1))
  n_patterns2 <- length(unique(d2))
  if(is.null(patterns.prior) == TRUE) {
  pp1 <- rep(1, d_list$n_patterns[1])
  pp2 <- rep(1, d_list$n_patterns[2]) }
  if(is.null(patterns.prior) == FALSE) {
    if(is.list(patterns.prior) == FALSE | length(patterns.prior) != 2) {
      stop("priors for patterns in both arms must be provided as a list object") }
    if(length(patterns.prior[[1]]) != n_patterns1 | is.vector(patterns.prior[[1]]) == FALSE | any(is.numeric(patterns.prior[[1]])) == FALSE) {
      stop("provide correct hyper prior values") }
    if(length(patterns.prior[[2]]) != n_patterns2 | is.vector(patterns.prior[[2]]) == FALSE | any(is.numeric(patterns.prior[[2]])) == FALSE) {
      stop("provide correct hyper prior values") }
    pp1 <- patterns.prior[[1]]
    pp2 <- patterns.prior[[2]] }
  model <- write_pattern(type = type , dist_e = dist_e, dist_c = dist_c, pe = pe, pc = pc, ind = ind, d_list = d_list, d1 = d1, d2 = d2)
  filein <- model
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "n_patterns1", "n_patterns2", 
                   "X1_e", "X2_e", "X1_c", "X2_c", "mean_cov_e1", "mean_cov_e2", "mean_cov_c1", 
                   "mean_cov_c2", "pe", "pc", "pp1", "pp2", "d1", "d2", "range_e", "range_c")
  if(pe == 1) {pe_index <- match("pe", datalist)
  datalist <- datalist[-pe_index] }
  if(pc == 1) {pc_index <- match("pc", datalist)
  datalist <- datalist[-pc_index] }
  if(type == "MAR" | type == "MNAR_cost") {range_e_index <- match("range_e", datalist)
  datalist <- datalist[-range_e_index] }
  if(type == "MAR" | type == "MNAR_eff") {range_c_index <- match("range_c", datalist)
  datalist <- datalist[-range_c_index] }
  DIC <- TRUE
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "mu_c_p1", "mu_c_p2", "mu_e_p1", "mu_e_p2",
              "s_c_p1", "s_c_p2", "s_e_p1", "s_e_p2", "p_prob1", "p_prob2",
              "beta_p1", "beta_p2","alpha_p1", "alpha_p2", "Delta_e", "Delta_c")
  params <- c(params, "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_d1", "loglik_d2")
  if(ind == FALSE) {params <- c(params, "beta_f_p1", "beta_f_p2") }
  if(type == "MNAR_cost" | type == "MAR") {deltae_index <- match("Delta_e", params)
  params <- params[-deltae_index] }
  if(type == "MNAR_eff" | type == "MAR") {deltac_index <- match("Delta_c", params)
  params <- params[-deltac_index] }
  modelN1 <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, model.file = filein, n.chains = n.chains, 
                          n.iter = n.iter, n.burnin = n.burnin, DIC = DIC, n.thin = n.thin)
    mu_e <- modelN1$BUGSoutput$sims.list$mu_e
    mu_c <- modelN1$BUGSoutput$sims.list$mu_c
    mu_c_p1 <- modelN1$BUGSoutput$sims.list$mu_c_p1
    mu_c_p2 <- modelN1$BUGSoutput$sims.list$mu_c_p2
    mu_e_p1 <- modelN1$BUGSoutput$sims.list$mu_e_p1
    mu_e_p2 <- modelN1$BUGSoutput$sims.list$mu_e_p2
    s_c_p1 <- modelN1$BUGSoutput$sims.list$s_c_p1
    s_c_p2 <- modelN1$BUGSoutput$sims.list$s_c_p2
    s_e_p1 <- modelN1$BUGSoutput$sims.list$s_e_p1
    s_e_p2 <- modelN1$BUGSoutput$sims.list$s_e_p2
    alpha_p1 <- modelN1$BUGSoutput$sims.list$alpha_p1
    alpha_p2 <- modelN1$BUGSoutput$sims.list$alpha_p2
    beta_p1 <- modelN1$BUGSoutput$sims.list$beta_p1
    beta_p2 <- modelN1$BUGSoutput$sims.list$beta_p2
    p_prob1 <- modelN1$BUGSoutput$sims.list$p_prob1
    p_prob2 <- modelN1$BUGSoutput$sims.list$p_prob2
    eff1_pos <- matrix(eff1, N1 ,3)
    cost1_pos <- matrix(cost1, N1, 3)
    eff2_pos <- matrix(eff2, N2, 3)
    cost2_pos <- matrix(cost2, N2, 3)
    eff1_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, mean)
    eff1_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, quantile, probs = prob[1])
    eff1_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, quantile, probs = prob[2])
    eff2_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, mean)
    eff2_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, quantile, probs = prob[1])
    eff2_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, quantile, probs = prob[2])
    cost1_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, mean)
    cost1_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, quantile, probs = prob[1])
    cost1_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, quantile, probs = prob[2])
    cost2_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, mean)
    cost2_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, quantile, probs = prob[1])
    cost2_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, quantile, probs = prob[2])
    loglik_e1 <- modelN1$BUGSoutput$sims.list$loglik_e1
    loglik_e2 <- modelN1$BUGSoutput$sims.list$loglik_e2
    loglik_c1 <- modelN1$BUGSoutput$sims.list$loglik_c1
    loglik_c2 <- modelN1$BUGSoutput$sims.list$loglik_c2
    loglik_d1 <- modelN1$BUGSoutput$sims.list$loglik_d1
    loglik_d2 <- modelN1$BUGSoutput$sims.list$loglik_d2
  if(ind == FALSE & dist_e == "norm" & dist_c == "norm") {
    beta_f_p1 <- modelN1$BUGSoutput$sims.list$beta_f_p1
    beta_f_p2 <- modelN1$BUGSoutput$sims.list$beta_f_p2
    beta_p1 <- list("beta_p1" = beta_p1, "beta_f_p1" = beta_f_p1)
    beta_p2 <- list("beta_p2" = beta_p2, "beta_f_p2" = beta_f_p2)
  }
  if(type == "MNAR") {
    Delta_e <- modelN1$BUGSoutput$sims.list$Delta_e
    Delta_c <- modelN1$BUGSoutput$sims.list$Delta_c
  } else if(type == "MNAR_eff") {
    Delta_e <- modelN1$BUGSoutput$sims.list$Delta_e
  } else if(type == "MNAR_cost") {
    Delta_c <- modelN1$BUGSoutput$sims.list$Delta_c
  }
   if(n.chains > 1) {
     model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_e1', 'loglik_e2',
                                                            'loglik_c1', 'loglik_c2', 'loglik_d1', 'loglik_d2'), invert = TRUE), digits = 3)
   } else{model_sum <- NULL }
    loglik_e <- list("control" = loglik_e1, "intervention" = loglik_e2)
    loglik_c <- list("control" = loglik_c1, "intervention" = loglik_c2)
    loglik_d <- list("control" = loglik_d1, "intervention" = loglik_d2)
    loglik <- list("effects" = loglik_e, "costs" = loglik_c, "pattern indicators" = loglik_d)
    mu_e_p <- list("control" = mu_e_p1, "intervention" = mu_e_p2)
    mu_c_p <- list("control" = mu_c_p1, "intervention" = mu_c_p2)
    s_e_p <- list("control" = s_e_p1, "intervention" = s_e_p2)
    s_c_p <- list("control" = s_c_p1, "intervention" = s_c_p2)
    alpha_p <- list("control" = alpha_p1, "intervention" = alpha_p2)
    beta_p <- list("control" = beta_p1, "intervention" = beta_p2)
    prob_p <- list("control" = p_prob1, "intervention" = p_prob2)
    colnames(eff1_pos) <- c("mean", "LB", "UB")
    colnames(eff2_pos) <- c("mean", "LB", "UB")
    colnames(cost1_pos) <- c("mean", "LB", "UB")
    colnames(cost2_pos) <- c("mean", "LB", "UB")
    imputed <- list("effects1" = eff1_pos, "effects2" = eff2_pos, "costs1" = cost1_pos, "costs2" = cost2_pos)
    if(type == "MAR") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_pattern" = alpha_p, "covariate_parameter_costs_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN", "ind" = ind)
    } else if(type == "MNAR") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_pattern" = alpha_p, "covariate_parameter_costs_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "mnar_parameter_effects" = Delta_e, "mnar_parameter_costs" = Delta_c, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_ec", "ind" = ind)
    } else if(type == "MNAR_eff") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_pattern" = alpha_p, "covariate_parameter_costs_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "mnar_parameter_effects" = Delta_e, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_e", "ind" = ind)
    } else if(type == "MNAR_cost") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_pattern" = alpha_p, "covariate_parameter_costs_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "mnar_parameter_costs" = Delta_c, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_c", "ind" = ind)
    }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))