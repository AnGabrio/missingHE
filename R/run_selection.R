#' An internal function to execute a JAGS selection model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR).
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters.
#' @param ppc Logical. If \code{ppc} is \code{TRUE}, the estimates of the parameters that can be used to generate replications from the model are saved.
#' @keywords JAGS Bayesian selection models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_selection <- function(type, dist_e, dist_c, inits, ppc) eval.parent(substitute( {
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
  model <- write_selection(type = type , dist_e = dist_e, dist_c = dist_c, pe_fixed = pe_fixed, pc_fixed = pc_fixed, ze_fixed = ze_fixed, zc_fixed = zc_fixed, ind_fixed = ind_fixed,
                           pe_random = pe_random, pc_random = pc_random, ze_random = ze_random, zc_random = zc_random, ind_random = ind_random, 
                           model_e_random = model_e_random, model_c_random = model_c_random, model_me_random = model_me_random, model_mc_random = model_mc_random)
  filein <- model
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "m_eff1", "m_eff2", "m_cost1", "m_cost2", 
                   "X1_e_fixed", "X2_e_fixed", "X1_c_fixed", "X2_c_fixed", "Z1_e_fixed", "Z2_e_fixed", "Z1_c_fixed", "Z2_c_fixed", 
                   "mean_cov_e1_fixed", "mean_cov_e2_fixed", "mean_cov_c1_fixed", "mean_cov_c2_fixed", "mean_z_e1_fixed", "mean_z_e2_fixed", "mean_z_c1_fixed", "mean_z_c2_fixed", 
                   "pe_fixed", "pc_fixed", "ze_fixed", "zc_fixed", "X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", 
                   "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e", "X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", 
                   "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c", "Z1_e_random", "Z2_e_random","mean_z_e1_random", "mean_z_e2_random", "ze_random", 
                   "clus1_me", "clus2_me", "n1_clus_me", "n2_clus_me", "Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", 
                   "clus1_mc", "clus2_mc", "n1_clus_mc", "n2_clus_mc")
  e_random_list <- c("X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e")
  c_random_list <- c("X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c")
  me_random_list <- c("Z1_e_random", "Z2_e_random","mean_z_e1_random", "mean_z_e2_random", "ze_random", "clus1_me", "clus2_me", "n1_clus_me", "n2_clus_me")
  mc_random_list <- c("Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", "clus1_mc", "clus2_mc", "n1_clus_mc", "n2_clus_mc")
  if(pe_fixed == 1) {pe_fixed_index <- match("pe_fixed", datalist)
  datalist <- datalist[-pe_fixed_index] }
  if(pc_fixed == 1) {pc_fixed_index <- match("pc_fixed", datalist)
  datalist <- datalist[-pc_fixed_index] }
  if(ze_fixed == 1) {ze_fixed_index <- match("ze_fixed", datalist)
  datalist <- datalist[-ze_fixed_index] }
  if(zc_fixed == 1) {zc_fixed_index <- match("zc_fixed", datalist)
  datalist <- datalist[-zc_fixed_index] }
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
  if(length(model_me_random) != 0 & is_me_random_e == FALSE) {
    if(ze_random == 1 | length(model_me_random) != 0 & is_me_random_e == TRUE) {
    ze_random_index <- match("ze_random", datalist)
    datalist <- datalist[-ze_random_index] }
  } else if(length(model_me_random) != 0 & is_me_random_e == TRUE) {
  me_random_index <- match(me_random_list[1:5], datalist)
  datalist <- datalist[-me_random_index] 
  } else if(length(model_me_random) == 0) { 
    me_random_index <- match(me_random_list, datalist)
    datalist <- datalist[-me_random_index] 
  }
  if(length(model_mc_random) != 0 & is_mc_random_c == FALSE) {
    if(zc_random == 1 | length(model_mc_random) != 0 & is_mc_random_c == TRUE) {
      zc_random_index <- match("zc_random", datalist)
      datalist <- datalist[-zc_random_index] }
  } else if(length(model_mc_random) != 0 & is_mc_random_c == TRUE) {
    mc_random_index <- match(mc_random_list[1:5], datalist)
    datalist <- datalist[-mc_random_index] 
  } else if(length(model_mc_random) == 0) { 
    mc_random_index <- match(mc_random_list, datalist)
    datalist <- datalist[-mc_random_index] 
  }
  DIC <- TRUE
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "s_e", "s_c", "p_e", "p_c", "beta", "alpha", "gamma_e", "gamma_c", "delta_e", "delta_c")
  params <- c(params, "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_me1", "loglik_me2", "loglik_mc1", "loglik_mc2")
  if(ind_fixed == FALSE) {params <- c(params, "beta_f") }
  if(length(model_e_random) != 0){params <- c(params, "a1", "a2") }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE){params <- c(params, "b1", "b2") }
  if(length(model_me_random) != 0 & is_me_random_e == FALSE){params <- c(params, "g1_e", "g2_e") }
  if(length(model_mc_random) != 0 & is_mc_random_c == FALSE){params <- c(params, "g1_c", "g2_c") }
  if(length(model_c_random) != 0 & ind_random == FALSE) {params <- c(params, "b1_f", "b2_f") }
  if(type == "MNAR_cost" | type == "MAR") {
    deltae_index <- match("delta_e", params)
    params <- params[-deltae_index] 
  }
  if(type == "MNAR_eff" | type == "MAR") {
    deltac_index <- match("delta_c", params)
    params <- params[-deltac_index]
  }
  if(type == "MNAR" | type == "MNAR_eff") {
    if(length(model_me_random) != 0 & "e" %in% model_me_random){
      params <- c(params, "d1_e", "d2_e") 
    }
  }
  if(type == "MNAR" | type == "MNAR_cost") {
    if(length(model_mc_random) != 0 & "c" %in% model_mc_random){
      params <- c(params, "d1_c", "d2_c") 
    }
  }
  if(ppc == TRUE) { 
    if(dist_e %in% c("norm", "nbinom", "logis")) {
      ppc_e_params <- c("mu_e1", "mu_e2", "tau_e") 
    } else if(dist_e %in% c("beta", "gamma", "weibull")) {
      ppc_e_params <- c("mu_e1", "tau_e1", "mu_e2", "tau_e2")
    } else if(dist_e %in% c("exp", "bern", "pois")) {
      ppc_e_params <- c("mu_e1", "mu_e2")
    } 
    if(dist_c == "norm") {
      ppc_c_params <- c("mu_c1", "mu_c2", "tau_c") 
    } else if(dist_c == "gamma") {
      ppc_c_params <- c("mu_c1", "tau_c1", "mu_c2", "tau_c2")
    } else if(dist_c == "lnorm") {
      ppc_c_params <- c("lmu_c1", "lmu_c2", "ltau_c")
    } 
    params <- c(params, ppc_e_params, ppc_c_params)
  }
  modelN1 <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, model.file = filein, n.chains = n.chains, 
                          n.iter = n.iter, n.burnin = n.burnin, DIC = DIC, n.thin = n.thin)
  mu_e <- modelN1$BUGSoutput$sims.list$mu_e
  mu_c <- modelN1$BUGSoutput$sims.list$mu_c
  s_e <- modelN1$BUGSoutput$sims.list$s_e
  s_c <- modelN1$BUGSoutput$sims.list$s_c
  alpha <- modelN1$BUGSoutput$sims.list$alpha
  beta <- modelN1$BUGSoutput$sims.list$beta
  p_e <- modelN1$BUGSoutput$sims.list$p_e
  p_c <- modelN1$BUGSoutput$sims.list$p_c
  gamma_e <- modelN1$BUGSoutput$sims.list$gamma_e
  gamma_c <- modelN1$BUGSoutput$sims.list$gamma_c
  a <- b <- g_e <- g_c <- d_e <- d_c <- NULL
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
  if(length(model_me_random) != 0 & is_me_random_e == FALSE) { 
    g1_e <- modelN1$BUGSoutput$sims.list$g1_e
    g2_e <- modelN1$BUGSoutput$sims.list$g2_e
    g_e <- list("g1_e" = g1_e, "g2_e" = g2_e) 
  }
  if(length(model_mc_random) != 0 & is_mc_random_c == FALSE) { 
    g1_c <- modelN1$BUGSoutput$sims.list$g1_c
    g2_c <- modelN1$BUGSoutput$sims.list$g2_c
    g_c <- list("g1_c" = g1_c, "g2_c" = g2_c) 
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
  loglik_me1 <- modelN1$BUGSoutput$sims.list$loglik_me1
  loglik_me2 <- modelN1$BUGSoutput$sims.list$loglik_me2
  loglik_mc1 <- modelN1$BUGSoutput$sims.list$loglik_mc1
  loglik_mc2 <- modelN1$BUGSoutput$sims.list$loglik_mc2
  if(ind_fixed == FALSE) {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    beta <- list("beta" = beta, "beta_f" = beta_f)
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
    delta_e <- modelN1$BUGSoutput$sims.list$delta_e
    delta_c <- modelN1$BUGSoutput$sims.list$delta_c
    if(length(model_me_random) != 0 & "e" %in% model_me_random) {
      d1_e <- modelN1$BUGSoutput$sims.list$d1_e
      d2_e <- modelN1$BUGSoutput$sims.list$d2_e
      d_e <- list("d1_e" = d1_e, "d2_e" = d2_e) 
    }
    if(length(model_mc_random) != 0 & "c" %in% model_mc_random) {
      d1_c <- modelN1$BUGSoutput$sims.list$d1_c
      d2_c <- modelN1$BUGSoutput$sims.list$d2_c
      d_c <- list("d1_c" = d1_c, "d2_c" = d2_c) 
    }
  } else if(type == "MNAR_eff") {
    delta_e <- modelN1$BUGSoutput$sims.list$delta_e
    if(length(model_me_random) != 0 & "e" %in% model_me_random) {
      d1_e <- modelN1$BUGSoutput$sims.list$d1_e
      d2_e <- modelN1$BUGSoutput$sims.list$d2_e
      d_e <- list("d1_e" = d1_e, "d2_e" = d2_e) 
    }
  } else if(type == "MNAR_cost") {
    delta_c <- modelN1$BUGSoutput$sims.list$delta_c
    if(length(model_mc_random) != 0 & "c" %in% model_mc_random) {
      d1_c <- modelN1$BUGSoutput$sims.list$d1_c
      d2_c <- modelN1$BUGSoutput$sims.list$d2_c
      d_c <- list("d1_c" = d1_c, "d2_c" = d2_c) 
    }
  }
  if(n.chains > 1) {
    model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_e1', 'loglik_e2',
                                                           'loglik_c1', 'loglik_c2', 'loglik_me1', 'loglik_me2',
                                                           'loglik_mc1', 'loglik_mc2'), invert = TRUE), digits = 3)
  } else{model_sum <- NULL }
  loglik_e <- list("control" = loglik_e1, "intervention" = loglik_e2)
  loglik_c <- list("control" = loglik_c1, "intervention" = loglik_c2)
  loglik_me <- list("control" = loglik_me1, "intervention" = loglik_me2)
  loglik_mc <- list("control" = loglik_mc1, "intervention" = loglik_mc2)
  loglik <- list("effects" = loglik_e, "costs" = loglik_c, "missing indicators effects" = loglik_me, "missing indicators costs" = loglik_mc)
  colnames(eff1_pos) <- c("mean", "LB", "UB")
  colnames(eff2_pos) <- c("mean", "LB", "UB")
  colnames(cost1_pos) <- c("mean", "LB", "UB")
  colnames(cost2_pos) <- c("mean", "LB", "UB")
  imputed <- list("effects1" = eff1_pos, "effects2" = eff2_pos, "costs1" = cost1_pos, "costs2" = cost2_pos)
  if(type == "MAR") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_e, 
                              "missingness_parameter_effects_fixed" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, "missingness_parameter_effects_random" = g_e, "missingness_parameter_costs_random" = g_c, 
                              "imputed" = imputed, "loglik" = loglik,"type" = "SELECTION", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
  } else if(type == "MNAR") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_e, 
                              "missingness_parameter_effects_fixed" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_effects_fixed" = delta_e, "mnar_parameter_costs_fixed" = delta_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_e, "missingness_parameter_costs_random" = g_c, "mnar_parameter_effects_random" = d_e, "mnar_parameter_costs_random" = d_c, 
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_ec", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
  } else if(type == "MNAR_eff") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_e, 
                              "missingness_parameter_effects_fixed" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_effects_fixed" = delta_e, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_e, "missingness_parameter_costs_random" = g_c, "mnar_parameter_effects_random" = d_e, 
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_e", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
  } else if(type == "MNAR_cost") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_e, 
                              "missingness_parameter_effects_fixed" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_costs_fixed" = delta_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_e, "missingness_parameter_costs_random" = g_c, "mnar_parameter_costs_random" = d_c,
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_c", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
  }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))