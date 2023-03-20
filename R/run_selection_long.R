#' An internal function to execute a JAGS selection model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR).
#' @param dist_u distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
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


run_selection_long <- function(type, dist_u, dist_c, inits, ppc) eval.parent(substitute( {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  if(!dist_u %in% c("norm", "beta", "exp", "bern", "nbinom", "weibull", "gamma", "logis", "pois") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta', 'gamma', 'weibull', 'logis', 'exp', 'bern', 'pois', 'nbinom' for the effects and 'norm', 'gamma', 'lnorm' for the costs")
  }
  if(!type %in% c("MAR", "MNAR", "MNAR_eff", "MNAR_cost")) {
    stop("Types available for use are 'MAR', 'MNAR_eff', 'MNAR_cost'and 'MNAR'")
  }
  if(is.null(inits) == FALSE) {inits = inits }
  model <- write_selection_long(type = type , dist_u = dist_u, dist_c = dist_c, pu_fixed = pu_fixed, pc_fixed = pc_fixed, zu_fixed = zu_fixed, zc_fixed = zc_fixed, ind_fixed = ind_fixed,
                                ind_time_fixed = ind_time_fixed, pu_random = pu_random, pc_random = pc_random, zu_random = zu_random, zc_random = zc_random, ind_random = ind_random, 
                           model_u_random = model_u_random, model_c_random = model_c_random, model_mu_random = model_mu_random, model_mc_random = model_mc_random)
  filein <- model
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "m_eff1", "m_eff2", "m_cost1", "m_cost2", "max_time",
                   "X1_u_fixed", "X2_u_fixed", "X1_c_fixed", "X2_c_fixed", "Z1_u_fixed", "Z2_u_fixed", "Z1_c_fixed", "Z2_c_fixed", 
                   "mean_cov_u1_fixed", "mean_cov_u2_fixed", "mean_cov_c1_fixed", "mean_cov_c2_fixed", "mean_z_u1_fixed", "mean_z_u2_fixed", "mean_z_c1_fixed", "mean_z_c2_fixed", 
                   "pu_fixed", "pc_fixed", "zu_fixed", "zc_fixed", "X1_u_random", "X2_u_random","mean_cov_u1_random", "mean_cov_u2_random", "pu_random", 
                   "clus1_u", "clus2_u", "n1_clus_u", "n2_clus_u", "X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", 
                   "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c", "Z1_u_random", "Z2_u_random","mean_z_u1_random", "mean_z_u2_random", "zu_random", 
                   "clus1_mu", "clus2_mu", "n1_clus_mu", "n2_clus_mu", "Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", 
                   "clus1_mc", "clus2_mc", "n1_clus_mc", "n2_clus_mc")
  u_random_list <- c("X1_u_random", "X2_u_random","mean_cov_u1_random", "mean_cov_u2_random", "pu_random", "clus1_u", "clus2_u", "n1_clus_u", "n2_clus_u")
  c_random_list <- c("X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c")
  mu_random_list <- c("Z1_u_random", "Z2_u_random","mean_z_u1_random", "mean_z_u2_random", "zu_random", "clus1_mu", "clus2_mu", "n1_clus_mu", "n2_clus_mu")
  mc_random_list <- c("Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", "clus1_mc", "clus2_mc", "n1_clus_mc", "n2_clus_mc")
  if(pu_fixed == 1) {pu_fixed_index <- match("pu_fixed", datalist)
  datalist <- datalist[-pu_fixed_index] }
  if(pc_fixed == 1) {pc_fixed_index <- match("pc_fixed", datalist)
  datalist <- datalist[-pc_fixed_index] }
  if(pc_fixed == 0 & ind_fixed == FALSE) {
    pc_fixed_index <- match("pc_fixed", datalist)
    X_c_fixed_index <- match(c("X1_c_fixed", "X2_c_fixed", "mean_cov_c1_fixed", "mean_cov_c2_fixed"), datalist)
    datalist <- datalist[-c(pc_fixed_index, X_c_fixed_index)] 
  }
  if(zu_fixed == 1) {zu_fixed_index <- match("zu_fixed", datalist)
  datalist <- datalist[-zu_fixed_index] }
  if(zc_fixed == 1) {zc_fixed_index <- match("zc_fixed", datalist)
  datalist <- datalist[-zc_fixed_index] }
  if(length(model_u_random) != 0) {
    if(pu_random == 1) {pu_random_index <- match("pu_random", datalist)
    datalist <- datalist[-pu_random_index] }
  } else if(length(model_u_random) == 0) { u_random_index <- match(u_random_list, datalist)
  datalist <- datalist[-u_random_index] }
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
  if(length(model_mu_random) != 0 & is_mu_random_u == FALSE) {
    if(zu_random == 1 | length(model_mu_random) != 0 & is_mu_random_u == TRUE) {
    zu_random_index <- match("zu_random", datalist)
    datalist <- datalist[-zu_random_index] }
  } else if(length(model_mu_random) != 0 & is_mu_random_u == TRUE) {
  mu_random_index <- match(mu_random_list[1:5], datalist)
  datalist <- datalist[-mu_random_index] 
  } else if(length(model_mu_random) == 0) { 
    mu_random_index <- match(mu_random_list, datalist)
    datalist <- datalist[-mu_random_index] 
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
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_u", "mu_c", "s_u", "s_c", "p_u", "p_c", "beta", "alpha", "gamma_u", "gamma_c", "delta_u", "delta_c")
  params <- c(params, "loglik_u1", "loglik_u2", "loglik_c1", "loglik_c2", "loglik_mu1", "loglik_mu2", "loglik_mc1", "loglik_mc2")
  if(pc_fixed == 0 & ind_fixed == FALSE) {
    beta_index <- match("beta", params)
    params <- params[-beta_index]
  }
  if(ind_fixed == FALSE) {params <- c(params, "beta_f") }
  if(time_dep == "AR1") {params <- c(params, "beta_tu", "beta_tc", "alpha_tu", "alpha_tc") }
  if(length(model_u_random) != 0){params <- c(params, "a1", "a2") }
  if(length(model_u_random) != 0 & ind_time_fixed == FALSE & time_dep == "AR1"){params <- c(params, "a1_tu", "a2_tu", "a1_tc", "a2_tc") }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE){params <- c(params, "b1", "b2") }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE & ind_time_fixed == FALSE & time_dep == "AR1"){
    params <- c(params, "b1_tu", "b1_tc", "b2_tu", "b2_tc") 
  }
  if(length(model_c_random) != 0 & is_c_random_c == TRUE & ind_time_fixed == FALSE & time_dep == "AR1"){
    if(length(model_c_random) == 1 & "u" %in% model_c_random)
    params <- c(params, "b1_tu", "b1_tc", "b2_tu", "b2_tc") 
  }
  if(length(model_mu_random) != 0 & is_mu_random_u == FALSE){params <- c(params, "g1_u", "g2_u") }
  if(length(model_mc_random) != 0 & is_mc_random_c == FALSE){params <- c(params, "g1_c", "g2_c") }
  if(length(model_c_random) != 0 & ind_random == FALSE) {params <- c(params, "b1_f", "b2_f") }
  if(ind_random == TRUE | time_dep == "none"){
    if("a1_tu" %in% params | "a1_tc" %in% params | "a2_tu" %in% params | "a2_tc" %in% params) {
      a_corr_index <- match(c("a1_tu", "a1_tc", "a2_tu", "a2_tc"), params)
      params <- params[-a_corr_index]
    }
    if("b1_tu" %in% params | "b1_tc" %in% params | "b2_tu" %in% params | "b2_tc" %in% params) {
      b_corr_index <- match(c("b1_tu", "b1_tc", "b2_tu", "b2_tc"), params)
      params <- params[-b_corr_index]
    }
  }
  if(type == "MNAR_cost" | type == "MAR") {
    deltau_index <- match("delta_u", params)
    params <- params[-deltau_index] 
  }
  if(type == "MNAR_eff" | type == "MAR") {
    deltac_index <- match("delta_c", params)
    params <- params[-deltac_index]
  }
  if(type == "MNAR" | type == "MNAR_eff") {
    if(length(model_mu_random) != 0 & "u" %in% model_mu_random){
      params <- c(params, "d1_u", "d2_u") 
    }
  }
  if(type == "MNAR" | type == "MNAR_cost") {
    if(length(model_mc_random) != 0 & "c" %in% model_mc_random){
      params <- c(params, "d1_c", "d2_c") 
    }
  }
  if(ppc == TRUE) { 
    if(dist_u %in% c("norm", "nbinom", "logis")) {
      ppc_u_params <- c("mu_u1", "mu_u2", "tau_u") 
    } else if(dist_u %in% c("beta", "gamma", "weibull")) {
      ppc_u_params <- c("mu_u1", "tau_u1", "mu_u2", "tau_u2")
    } else if(dist_u %in% c("exp", "bern", "pois")) {
      ppc_u_params <- c("mu_u1", "mu_u2")
    } 
    if(dist_c == "norm") {
      ppc_c_params <- c("mu_c1", "mu_c2", "tau_c") 
    } else if(dist_c == "gamma") {
      ppc_c_params <- c("mu_c1", "tau_c1", "mu_c2", "tau_c2")
    } else if(dist_c == "lnorm") {
      ppc_c_params <- c("lmu_c1", "lmu_c2", "ltau_c")
    } 
    params <- c(params, ppc_u_params, ppc_c_params)
  }
  modelN1 <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, model.file = filein, n.chains = n.chains, 
                          n.iter = n.iter, n.burnin = n.burnin, DIC = DIC, n.thin = n.thin)
  mu_u <- modelN1$BUGSoutput$sims.list$mu_u
  mu_c <- modelN1$BUGSoutput$sims.list$mu_c
  s_u <- modelN1$BUGSoutput$sims.list$s_u
  s_c <- modelN1$BUGSoutput$sims.list$s_c
  alpha <- modelN1$BUGSoutput$sims.list$alpha
  if(pc_fixed == 0 & ind_fixed == FALSE) { 
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
  } else {
    beta <- modelN1$BUGSoutput$sims.list$beta
  }
  p_u <- modelN1$BUGSoutput$sims.list$p_u
  p_c <- modelN1$BUGSoutput$sims.list$p_c
  gamma_u <- modelN1$BUGSoutput$sims.list$gamma_u
  gamma_c <- modelN1$BUGSoutput$sims.list$gamma_c
  a <- b <- g_u <- g_c <- d_u <- d_c <- NULL
  if(length(model_u_random) != 0) { 
    a1 <- modelN1$BUGSoutput$sims.list$a1
    a2 <- modelN1$BUGSoutput$sims.list$a2 
    a <- list("a1" = a1, "a2" = a2) 
  }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE) { 
    b1 <- modelN1$BUGSoutput$sims.list$b1
    b2 <- modelN1$BUGSoutput$sims.list$b2 
    b <- list("b1" = b1, "b2" = b2) 
  }
  if(length(model_mu_random) != 0 & is_mu_random_u == FALSE) { 
    g1_u <- modelN1$BUGSoutput$sims.list$g1_u
    g2_u <- modelN1$BUGSoutput$sims.list$g2_u
    g_u <- list("g1_u" = g1_u, "g2_u" = g2_u) 
  }
  if(length(model_mc_random) != 0 & is_mc_random_c == FALSE) { 
    g1_c <- modelN1$BUGSoutput$sims.list$g1_c
    g2_c <- modelN1$BUGSoutput$sims.list$g2_c
    g_c <- list("g1_c" = g1_c, "g2_c" = g2_c) 
  }
  eff1_pos <- array(NA, dim = c(N1, max_time, 3))
  cost1_pos <- array(NA, dim = c(N1, max_time, 3))
  eff2_pos <- array(NA, dim = c(N2, max_time, 3))
  cost2_pos <- array(NA, dim = c(N2, max_time, 3))
  for(time in 1:max_time){
    eff1_pos[, time, 1] <- apply(modelN1$BUGSoutput$sims.list$eff1[, , time], 2, mean)
    eff1_pos[, time, 2] <- apply(modelN1$BUGSoutput$sims.list$eff1[, , time], 2, quantile, probs = prob[1])
    eff1_pos[, time, 3] <- apply(modelN1$BUGSoutput$sims.list$eff1[, , time], 2, quantile, probs = prob[2])
    eff2_pos[, time, 1] <- apply(modelN1$BUGSoutput$sims.list$eff2[, , time], 2, mean)
    eff2_pos[, time, 2] <- apply(modelN1$BUGSoutput$sims.list$eff2[, , time], 2, quantile, probs = prob[1])
    eff2_pos[, time, 3] <- apply(modelN1$BUGSoutput$sims.list$eff2[, , time], 2, quantile, probs = prob[2])
    cost1_pos[, time, 1] <- apply(modelN1$BUGSoutput$sims.list$cost1[, , time], 2, mean)
    cost1_pos[, time, 2] <- apply(modelN1$BUGSoutput$sims.list$cost1[, , time], 2, quantile, probs = prob[1])
    cost1_pos[, time, 3] <- apply(modelN1$BUGSoutput$sims.list$cost1[, , time], 2, quantile, probs = prob[2])
    cost2_pos[, time, 1] <- apply(modelN1$BUGSoutput$sims.list$cost2[, , time], 2, mean)
    cost2_pos[, time, 2] <- apply(modelN1$BUGSoutput$sims.list$cost2[, , time], 2, quantile, probs = prob[1])
    cost2_pos[, time, 3] <- apply(modelN1$BUGSoutput$sims.list$cost2[, , time], 2, quantile, probs = prob[2])
  }
  loglik_u1 <- modelN1$BUGSoutput$sims.list$loglik_u1
  loglik_u2 <- modelN1$BUGSoutput$sims.list$loglik_u2
  loglik_c1 <- modelN1$BUGSoutput$sims.list$loglik_c1
  loglik_c2 <- modelN1$BUGSoutput$sims.list$loglik_c2
  loglik_mu1 <- modelN1$BUGSoutput$sims.list$loglik_mu1
  loglik_mu2 <- modelN1$BUGSoutput$sims.list$loglik_mu2
  loglik_mc1 <- modelN1$BUGSoutput$sims.list$loglik_mc1
  loglik_mc2 <- modelN1$BUGSoutput$sims.list$loglik_mc2
  if(ind_fixed == FALSE & ind_time_fixed == FALSE & pc_fixed != 0) {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    if(time_dep == "AR1") {
    beta_tu <- modelN1$BUGSoutput$sims.list$beta_tu
    beta_tc <- modelN1$BUGSoutput$sims.list$beta_tc
    alpha_tu <- modelN1$BUGSoutput$sims.list$alpha_tu
    alpha_tc <- modelN1$BUGSoutput$sims.list$alpha_tc
    beta <- list("beta" = beta, "beta_f" = beta_f, "beta_tu" = beta_tu, "beta_tc" = beta_tc)
    alpha <- list("alpha" = alpha, "alpha_tu" = alpha_tu, "alpha_tc" = alpha_tc)
    } else if(time_dep == "none") {
      beta <- list("beta" = beta, "beta_f" = beta_f)
    }
  } else if(ind_fixed == FALSE & ind_time_fixed == FALSE & pc_fixed == 0) {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    if(time_dep == "AR1") {
      beta_tu <- modelN1$BUGSoutput$sims.list$beta_tu
      beta_tc <- modelN1$BUGSoutput$sims.list$beta_tc
      alpha_tu <- modelN1$BUGSoutput$sims.list$alpha_tu
      alpha_tc <- modelN1$BUGSoutput$sims.list$alpha_tc
      beta <- list("beta_f" = beta_f, "beta_tu" = beta_tu, "beta_tc" = beta_tc)
      alpha <- list("alpha" = alpha, "alpha_tu" = alpha_tu, "alpha_tc" = alpha_tc)
    } else if(time_dep == "none") {
      beta <- list("beta_f" = beta_f)
    }
  }
  if(ind_fixed == TRUE & ind_time_fixed == FALSE) {
    if(time_dep == "AR1") {
    beta_tu <- modelN1$BUGSoutput$sims.list$beta_tu
    beta_tc <- modelN1$BUGSoutput$sims.list$beta_tc
    alpha_tu <- modelN1$BUGSoutput$sims.list$alpha_tu
    alpha_tc <- modelN1$BUGSoutput$sims.list$alpha_tc
    beta <- list("beta" = beta, "beta_tu" = beta_tu, "beta_tc" = beta_tc)
    alpha <- list("alpha" = alpha, "alpha_tu" = alpha_tu, "alpha_tc" = alpha_tc)
    }
  }
  if(ind_fixed == FALSE & ind_time_fixed == TRUE & pc_fixed != 0) {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    beta <- list("beta" = beta, "beta_f" = beta_f)
  } else if(ind_fixed == FALSE & ind_time_fixed == TRUE & pc_fixed == 0) {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    beta <- list("beta_f" = beta_f)
  }
  if(length(model_c_random) != 0 & ind_random == FALSE & ind_time_fixed == FALSE) {
    b1_f <- modelN1$BUGSoutput$sims.list$b1_f
    b2_f <- modelN1$BUGSoutput$sims.list$b2_f
    b1_tu <- modelN1$BUGSoutput$sims.list$b1_tu
    b2_tu <- modelN1$BUGSoutput$sims.list$b2_tu
    b1_tc <- modelN1$BUGSoutput$sims.list$b1_tc
    b2_tc <- modelN1$BUGSoutput$sims.list$b2_tc
    a1_tu <- modelN1$BUGSoutput$sims.list$a1_tu
    a2_tu <- modelN1$BUGSoutput$sims.list$a2_tu
    a1_tc <- modelN1$BUGSoutput$sims.list$a1_tc
    a2_tc <- modelN1$BUGSoutput$sims.list$a2_tc
    if(length(model_c_random) != 0 & "u" %in% model_c_random & is_c_random_c == TRUE) {
    b <- list("b1_f" = b1_f, "b2_f" = b2_f, "b1_tu" = b1_tu, "b2_tu" = b2_tu, "b1_tc" = b1_tc, "b2_tc" = b2_tc) 
    } else if(length(model_c_random) != 0 & "u" %in% model_c_random & is_c_random_c == FALSE) {
    b <- list("b1" = b1, "b1_f" = b1_f, "b2" = b2, "b2_f" = b2_f, "b1_tu" = b1_tu, "b2_tu" = b2_tu, "b1_tc" = b1_tc, "b2_tc" = b2_tc)
    }
    if(length(model_u_random) != 0) {
      a <- list("a1" = a1, "a2" = a2, "a1_tu" = a1_tu, "a2_tu" = a2_tu, "a1_tc" = a1_tc, "a2_tc" = a2_tc) 
    }
  }
  if(length(model_c_random) != 0 & ind_random == FALSE & ind_time_fixed == TRUE) {
    b1_f <- modelN1$BUGSoutput$sims.list$b1_f
    b2_f <- modelN1$BUGSoutput$sims.list$b2_f
    if(length(model_c_random) != 0 & "u" %in% model_c_random & is_c_random_c == TRUE) {
      b <- list("b1_f" = b1_f, "b2_f" = b2_f) 
    } else if(length(model_c_random) != 0 & "u" %in% model_c_random & is_c_random_c == FALSE) {
      b <- list("b1" = b1, "b1_f" = b1_f, "b2" = b2, "b2_f" = b2_f)
    }
  }
  if(type == "MNAR") {
    delta_u <- modelN1$BUGSoutput$sims.list$delta_u
    delta_c <- modelN1$BUGSoutput$sims.list$delta_c
    if(length(model_mu_random) != 0 & "u" %in% model_mu_random) {
      d1_u <- modelN1$BUGSoutput$sims.list$d1_u
      d2_u <- modelN1$BUGSoutput$sims.list$d2_u
      d_u <- list("d1_u" = d1_u, "d2_u" = d2_u) 
    }
    if(length(model_mc_random) != 0 & "c" %in% model_mc_random) {
      d1_c <- modelN1$BUGSoutput$sims.list$d1_c
      d2_c <- modelN1$BUGSoutput$sims.list$d2_c
      d_c <- list("d1_c" = d1_c, "d2_c" = d2_c) 
    }
  } else if(type == "MNAR_eff") {
    delta_u <- modelN1$BUGSoutput$sims.list$delta_u
    if(length(model_mu_random) != 0 & "u" %in% model_mu_random) {
      d1_u <- modelN1$BUGSoutput$sims.list$d1_u
      d2_u <- modelN1$BUGSoutput$sims.list$d2_u
      d_u <- list("d1_u" = d1_u, "d2_u" = d2_u) 
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
    model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_u1', 'loglik_u2',
                                                           'loglik_c1', 'loglik_c2', 'loglik_mu1', 'loglik_mu2',
                                                           'loglik_mc1', 'loglik_mc2'), invert = TRUE), digits = 3)
  } else{model_sum <- NULL }
  loglik_u <- list("control" = loglik_u1, "intervention" = loglik_u2)
  loglik_c <- list("control" = loglik_c1, "intervention" = loglik_c2)
  loglik_mu <- list("control" = loglik_mu1, "intervention" = loglik_mu2)
  loglik_mc <- list("control" = loglik_mc1, "intervention" = loglik_mc2)
  loglik <- list("effects" = loglik_u, "costs" = loglik_c, "missing indicators effects" = loglik_mu, "missing indicators costs" = loglik_mc)
  colnames(eff1_pos) <- c("mean", "LB", "UB")
  colnames(eff2_pos) <- c("mean", "LB", "UB")
  colnames(cost1_pos) <- c("mean", "LB", "UB")
  colnames(cost2_pos) <- c("mean", "LB", "UB")
  imputed <- list("effects1" = eff1_pos, "effects2" = eff2_pos, "costs1" = cost1_pos, "costs2" = cost2_pos)
  if(type == "MAR") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_u, "mean_costs" = mu_c, "sd_effects" = s_u, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_u, 
                              "missingness_parameter_effects_fixed" = gamma_u, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, "missingness_parameter_effects_random" = g_u, "missingness_parameter_costs_random" = g_c, 
                              "imputed" = imputed, "loglik" = loglik,"type" = "SELECTION", "ind_fixed" = ind_fixed, "ind_time_fixed" = ind_time_fixed,"ind_random" = ind_random, "ppc" = ppc, "dist_u" = dist_u, "dist_c" = dist_c)
  } else if(type == "MNAR") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_u, "mean_costs" = mu_c, "sd_effects" = s_u, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_u, 
                              "missingness_parameter_effects_fixed" = gamma_u, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_effects_fixed" = delta_u, "mnar_parameter_costs_fixed" = delta_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_u, "missingness_parameter_costs_random" = g_c, "mnar_parameter_effects_random" = d_u, "mnar_parameter_costs_random" = d_c, 
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_uc", "ind_fixed" = ind_fixed, "ind_time_fixed" = ind_time_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_u" = dist_u, "dist_c" = dist_c)
  } else if(type == "MNAR_eff") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_u, "mean_costs" = mu_c, "sd_effects" = s_u, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_u, 
                              "missingness_parameter_effects_fixed" = gamma_u, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_effects_fixed" = delta_u, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_u, "missingness_parameter_costs_random" = g_c, "mnar_parameter_effects_random" = d_u, 
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_u", "ind_fixed" = ind_fixed, "ind_time_fixed" = ind_time_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_u" = dist_u, "dist_c" = dist_c)
  } else if(type == "MNAR_cost") {
    model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_u, "mean_costs" = mu_c, "sd_effects" = s_u, "sd_costs" = s_c, 
                              "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "missingness_probability_effects" = p_u, 
                              "missingness_parameter_effects_fixed" = gamma_u, "missingness_probability_costs" = p_c, "missingness_parameter_costs_fixed" = gamma_c, 
                              "mnar_parameter_costs_fixed" = delta_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b, 
                              "missingness_parameter_effects_random" = g_u, "missingness_parameter_costs_random" = g_c, "mnar_parameter_costs_random" = d_c,
                              "imputed" = imputed, "loglik" = loglik, "type" = "SELECTION_c", "ind_fixed" = ind_fixed, "ind_time_fixed" = ind_time_fixed, "ind_random" = ind_random, "ppc" = ppc, "dist_u" = dist_u, "dist_c" = dist_c)
  }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))