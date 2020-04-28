#' An internal function to execute a JAGS pattern mixture model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR).
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters.
#' @param d_list a list of the number and types of patterns in the data.
#' @param d1 Patterns in the control.
#' @param d2 Patterns in the intervention.
#' @param restriction type of identifying restriction to be imposed.
#' @param ppc Logical. If \code{ppc} is \code{TRUE}, the estimates of the parameters that can be used to generate replications from the model are saved.
#' @keywords JAGS Bayesian pattern mixture models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_pattern <- function(type, dist_e, dist_c, inits, d_list, d1, d2, restriction, ppc) eval.parent(substitute( {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  if(!dist_e %in% c("norm", "beta", "exp", "bern", "nbinom", "weibull", "gamma", "logis", "pois") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta', 'gamma', 'weibull', 'logis', 'exp', 'bern', 'pois', 'nbinom' for the effects and 'norm', 'gamma', 'lnorm' for the costs")
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
    pp2 <- patterns.prior[[2]] 
  }
  model <- write_pattern(type = type , dist_e = dist_e, dist_c = dist_c, pe_fixed = pe_fixed, pc_fixed = pc_fixed, ind_fixed = ind_fixed,
                         pe_random = pe_random, pc_random = pc_random, ind_random = ind_random, model_e_random = model_e_random, model_c_random = model_c_random, 
                         d_list = d_list, d1 = d1, d2 = d2, restriction = restriction)
  filein <- model
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "n_patterns1", "n_patterns2", "X1_e_fixed", "X2_e_fixed", "X1_c_fixed", "X2_c_fixed", 
                   "mean_cov_e1_fixed", "mean_cov_e2_fixed", "mean_cov_c1_fixed", "mean_cov_c2_fixed", "pe_fixed", "pc_fixed", 
                   "X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e", 
                   "X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c", 
                   "pp1", "pp2", "d1", "d2", "range_e", "range_c")
  e_random_list <- c("X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e")
  c_random_list <- c("X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c")
  if(pe_fixed == 1) {pe_fixed_index <- match("pe_fixed", datalist)
  datalist <- datalist[-pe_fixed_index] }
  if(pc_fixed == 1) {pc_fixed_index <- match("pc_fixed", datalist)
  datalist <- datalist[-pc_fixed_index] }
  if(length(model_e_random) != 0) {
    if(pe_random == 1) {pe_random_index <- match("pe_random", datalist)
    datalist <- datalist[-pe_random_index] }
  } else if(length(model_e_random) == 0) { e_random_index <- match(e_random_list, datalist)
  datalist <- datalist[-e_random_index] }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE) {
    if(pc_random == 1 | length(model_c_random) != 0 & is_c_random_c == TRUE) {
      pc_random_index <- match("pc_random", datalist)
      datalist <- datalist[-pc_random_index] }
  } else if(length(model_c_random) != 0 & is_c_random_c == TRUE) {
    c_random_index <- match(c_random_list[1:5], datalist)
    datalist <- datalist[-c_random_index] 
  } else if(length(model_c_random) == 0) { 
    c_random_index <- match(c_random_list, datalist)
    datalist <- datalist[-c_random_index] 
  }
  if(type == "MAR" | type == "MNAR_cost") {range_e_index <- match("range_e", datalist)
  datalist <- datalist[-range_e_index] }
  if(type == "MAR" | type == "MNAR_eff") {range_c_index <- match("range_c", datalist)
  datalist <- datalist[-range_c_index] }
  DIC <- TRUE
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "mu_c_p1", "mu_c_p2", "mu_e_p1", "mu_e_p2",
              "s_c_p1", "s_c_p2", "s_e_p1", "s_e_p2", "p_prob1", "p_prob2",
              "beta_p1", "beta_p2","alpha_p1", "alpha_p2", "Delta_e", "Delta_c")
  params <- c(params, "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_d1", "loglik_d2")
  if(ind_fixed == FALSE) {params <- c(params, "beta_f_p1", "beta_f_p2") }
  if(length(model_e_random) != 0){params <- c(params, "a1", "a2") }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE){params <- c(params, "b1", "b2") }
  if(length(model_c_random) != 0 & ind_random == FALSE) {params <- c(params, "b1_f", "b2_f") }
  if(type == "MNAR_cost" | type == "MAR") {deltae_index <- match("Delta_e", params)
  params <- params[-deltae_index] }
  if(type == "MNAR_eff" | type == "MAR") {deltac_index <- match("Delta_c", params)
  params <- params[-deltac_index] }
    if(ppc == TRUE) { 
     if(dist_e %in% c("norm", "nbinom", "logis")) {
       ppc_e_params <- c("mu_e1", "mu_e2", "tau_e_p1", "tau_e_p2") 
     } else if(dist_e %in% c("beta", "gamma", "weibull")) {
       ppc_e_params <- c("mu_e1", "tau_e1", "mu_e2", "tau_e2")
     } else if(dist_e %in% c("exp", "bern", "pois")) {
       ppc_e_params <- c("mu_e1", "mu_e2")
     } 
     if(dist_c == "norm") {
       ppc_c_params <- c("mu_c1", "mu_c2", "tau_c_p1", "tau_c_p2") 
     } else if(dist_c == "gamma") {
       ppc_c_params <- c("mu_c1", "tau_c1", "mu_c2", "tau_c2")
     } else if(dist_c == "lnorm") {
       ppc_c_params <- c("lmu_c1", "lmu_c2", "ltau_c_p1", "ltau_c_p2")
      } 
      params <- c(params, ppc_e_params, ppc_c_params)
    }
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
    a <- b <- NULL
    if(length(model_e_random) != 0) { 
      a1 <- modelN1$BUGSoutput$sims.list$a1
      a2 <- modelN1$BUGSoutput$sims.list$a2 
      a <- list("a1" = a1, "a2" = a2) 
    }
    if(length(model_c_random) != 0 & is_c_random_c == FALSE) { 
      b1 <- modelN1$BUGSoutput$sims.list$b1
      b2 <- modelN1$BUGSoutput$sims.list$b2 
      b <- list("b1" = b1, "b2" = b2) 
    }
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
    if(ind_fixed == FALSE) {
      beta_f_p1 <- modelN1$BUGSoutput$sims.list$beta_f_p1
      beta_f_p2 <- modelN1$BUGSoutput$sims.list$beta_f_p2
      beta_p1 <- list("beta_p1" = beta_p1, "beta_f_p1" = beta_f_p1)
      beta_p2 <- list("beta_p2" = beta_p2, "beta_f_p2" = beta_f_p2)
    }
    if(length(model_c_random) != 0 & ind_random == FALSE) {
      b1_f <- modelN1$BUGSoutput$sims.list$b1_f
      b2_f <- modelN1$BUGSoutput$sims.list$b2_f
      if(length(model_c_random) != 0 & "e" %in% model_c_random & is_c_random_c == TRUE) {
        b <- list("b1_f" = b1_f, "b2_f" = b2_f) 
      } else if(length(model_c_random) != 0 & "e" %in% model_c_random & is_c_random_c == FALSE) {
        b <- list("b1" = b1, "b1_f" = b1_f, "b2" = b2, "b2_f" = b2_f)
      }
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
                                "covariate_parameter_effects_fixed_pattern" = alpha_p, "covariate_parameter_costs_fixed_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                                "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    } else if(type == "MNAR") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_fixed_pattern" = alpha_p, "covariate_parameter_costs_fixed_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                                "mnar_parameter_effects" = Delta_e, "mnar_parameter_costs" = Delta_c, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_ec", 
                                "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    } else if(type == "MNAR_eff") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_fixed_pattern" = alpha_p, "covariate_parameter_costs_fixed_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                                "mnar_parameter_effects" = Delta_e, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_e", "ind_fixed" = ind_fixed, "ind_random" = ind_random,
                                "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    } else if(type == "MNAR_cost") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "mean_effects_pattern" = mu_e_p,
                                "mean_costs_pattern" = mu_c_p, "sd_effects_pattern" = s_e_p, "sd_costs_pattern" = s_c_p, 
                                "covariate_parameter_effects_fixed_pattern" = alpha_p, "covariate_parameter_costs_fixed_pattern" = beta_p, "pattern_probability" = prob_p, 
                                "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                                "mnar_parameter_costs" = Delta_c, "imputed" = imputed, "loglik" = loglik,"type" = "PATTERN_c", "ind_fixed" = ind_fixed, "ind_random" = ind_random,
                                "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))