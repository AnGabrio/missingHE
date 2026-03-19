#'An internal function to select which type of hurdle model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of structural value mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Hurdle models
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR) and Structural At Random (SAR).
#' @param model_txt_info list containing model specification information used to write the txt file of the JAGS model.
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_hurdle <- function(dist_e , dist_c, type, model_txt_info) {
  model_string_jags <- "
  model{

  #control
  for(i in 1:n) {
  #costs and effects model
  cost[i] ~ dnorm(cmu_c[i], tau_c[s_cost[i] + 1])
  eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])
  
  #derive mean and std effects
  #derive mean and std costs

  #mean regression
  cmu_c[i] <- inprod(X_c_fixed[i, ], beta[, s_cost[i] + 1]) + beta_f[s_cost[i] + 1] * (eff[i] - tmu_e) + inprod(X_c_random[i, ], b[, clus_c[i]]) * ((s_cost[i] - 1) * -1) + b_f[clus_c[i]] * (eff[i] - tmu_e) * ((s_cost[i] - 1) * -1)
  cmu_e[i] <- inprod(X_e_fixed[i, ], alpha[, s_eff[i] + 1]) + inprod(X_e_random[i, ], a[, clus_e[i]]) * ((s_eff[i] - 1) * -1)

  #structural values mechanism
  s_eff[i] ~ dbern(pq[i])
  logit(pq[i]) <- inprod(Z_e_fixed[i, ], gamma_e[]) + inprod(Z_e_random[i, ], g_e[, clus_se[i]])
  s_cost[i] ~ dbern(pc[i])
  logit(pc[i]) <- inprod(Z_c_fixed[i, ], gamma_c[]) + inprod(Z_c_random[i, ], g_c[, clus_sc[i]])

  #loglikelihood
  loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])
  loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c[s_cost[i] + 1])
  loglik_se[i] <- logdensity.bern(s_eff[i], pq[i])
  loglik_sc[i] <- logdensity.bern(s_cost[i], pc[i])
  }
  
  #transformation of parameters
  for (s in 1:2) {#begin transformation costs
  tau_c[s] <- 1 / ss_c[s]
  ss_c[s] <- s_c[s] * s_c[s]
  }#end transformation costs
  #std for lnorm 
  #mean for lnorm
  for (s in 1:2) {#begin transformation effects
  tau_e[s] <- 1 / ss_e[s]
  ss_e[s] <- s_e[s] * s_e[s]
  }#end transformation effects

  #transformation of random effects parameters
  # begin transformation random effects
  for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]
  ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }
  for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]
  ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }
  for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]
  ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }
  for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]
  ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }
  #end transformation of random effects 

  #calculate means at mean of covariates for non-structural values
  tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1]) + inprod(mean_cov_c_random[], mu_b_hat[]) 
  tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]) 
  
  #weighted mean parameters
  tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e
  tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c
  
  #structural values probability
  p_c <- ilogit(inprod(mean_z_c_fixed[], gamma_c[]) + inprod(mean_z_c_random[], mu_g_c_hat[]))
  p_e <- ilogit(inprod(mean_z_e_fixed[], gamma_e[]) + inprod(mean_z_e_random[], mu_g_e_hat[]))

  #priors
  
  #priors for mean regression coefficients
  for (j in 2:pe_fixed) {#begin alpha priors effects
  alpha[j, 1] ~ dnorm(0, 0.0000001)
  alpha[j, 2] <- 0
  }#end alpha priors effects
  alpha[1, 1] ~ dnorm(0, 0.0000001)
  alpha[1, 2] <- se

  for (j in 2:pc_fixed) {#begin beta priors costs
  beta[j, 1] ~ dnorm(0, 0.0000001)
  beta[j, 2] <- 0
  }#end beta priors costs
  beta[1, 1] ~ dnorm(0, 0.0000001)
  beta[1, 2] <- sc

  #priors for mean regression random coefficients
  for (j in 1:pe_random) {#begin a priors effects
  for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }
  }#end a priors effects

  for (j in 1:pc_random) {#begin b priors costs
  for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }
  }#end b priors costs
  
  #standard deviation priors
  s_c[1] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_c[2] <- sdc
  s_e[2] <- sde

  #correlation
  beta_f[1] ~ dnorm(0, 0.0000001)
  beta_f[2] <- 0
  
  # mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)
  s_b_hat[j] ~ dunif(0, 100) }
  for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)
  s_a_hat[j] ~ dunif(0, 100) }
  for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)
  s_g_c_hat[j] ~ dunif(0, 100) }
  for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)
  s_g_e_hat[j] ~ dunif(0, 100) }

  # correlation random effects
  for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }
  mu_b_f_hat ~ dnorm(0, 0.001)
  tau_b_f_hat <- 1 / ss_b_f_hat
  ss_b_f_hat <- s_b_f_hat * s_b_f_hat
  s_b_f_hat ~ dunif(0, 100)
  
  #priors on structural values mechanism
  for (j in 1:ze_fixed) {#begin gamma priors effects
  gamma_e[j] ~ dnorm(0, 0.01)
  }#end gamma priors effects
  
  for (j in 1:zc_fixed) {#begin gamma priors costs
  gamma_c[j] ~ dnorm(0, 0.01)
  }#end gamma priors costs
  
  #priors on random effects structural values mechanism
  for (j in 1:ze_random) {#begin g_e priors effects
  for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }
  }#end g_e priors effects

  for (j in 1:zc_random) {#begin g_c priors costs
  for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }
  }#end g_c priors costs
  
  }
  "
if(is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dnorm(cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- inprod(X_e_fixed[i, ], alpha[, s_eff[i] + 1]) + inprod(X_e_random[i, ], a[, clus_e[i]]) * ((s_eff[i] - 1) * -1)", "cmu_e[i] <- inprod(X_e_fixed[i, ], alpha[]) + inprod(X_e_random[i, ], a[, clus_e[i]])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_eff[i] ~ dbern(pq[i])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("logit(pq[i]) <- inprod(Z_e_fixed[i, ], gamma_e[]) + inprod(Z_e_random[i, ], g_e[, clus_se[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "tau_e <- 1 / ss_e", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "ss_e <- s_e * s_e", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "#begin transformation effects", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end transformation effects", "#end transformation effects", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", 
                            "tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("p_e <- ilogit(inprod(mean_z_e_fixed[], gamma_e[]) + inprod(mean_z_e_random[], mu_g_e_hat[]))", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[j, 1] ~ dnorm(0, 0.0000001)", "alpha[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[1, 1] ~ dnorm(0, 0.0000001)", "alpha[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[1, 2] <- se", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[2] <- sde", "", model_string_jags, fixed = TRUE)  
  model_string_jags <- gsub("for (j in 1:ze_fixed) {#begin gamma priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_e[j] ~ dnorm(0, 0.01)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end gamma priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_se[i] <- logdensity.bern(s_eff[i], pq[i])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
}
if(!is.null(model_txt_info$se) & is.null(model_txt_info$sc)) {
  model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c[s_cost[i] + 1])", "cost[i] ~ dnorm(cmu_c[i], tau_c)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i] <- inprod(X_c_fixed[i, ], beta[, s_cost[i] + 1]) + beta_f[s_cost[i] + 1] * (eff[i] - tmu_e) + inprod(X_c_random[i, ], b[, clus_c[i]]) * ((s_cost[i] - 1) * -1) + b_f[clus_c[i]] * (eff[i] - tmu_e) * ((s_cost[i] - 1) * -1)", "cmu_c[i] <- inprod(X_c_fixed[i, ], beta[]) + beta_f * (eff[i] - tmu_e) + inprod(X_c_random[i, ], b[, clus_c[i]]) + b_f[clus_c[i]] * (eff[i] - tmu_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_cost[i] ~ dbern(pc[i])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("logit(pc[i]) <- inprod(Z_c_fixed[i, ], gamma_c[]) + inprod(Z_c_random[i, ], g_c[, clus_sc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[s] <- 1 / ss_c[s]", "tau_c <- 1 / ss_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[s] <- s_c[s] * s_c[s]", "ss_c <- s_c * s_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (s in 1:2) {#begin transformation costs", "#begin transformation costs", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end transformation costs", "#end transformation costs", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1])", "tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c", "tmu_c <- tmu_ns_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("p_c <- ilogit(inprod(mean_z_c_fixed[], gamma_c[]) + inprod(mean_z_c_random[], mu_g_c_hat[]))", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[j, 1] ~ dnorm(0, 0.0000001)", "beta[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[1, 1] ~ dnorm(0, 0.0000001)", "beta[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[1, 2] <- sc", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c ~ dt(0, pow(2.5, -2), 1)T(0,)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[2] <- sdc", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_fixed) {#begin gamma priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_c[j] ~ dnorm(0, 0.01)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end gamma priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta_f[1] ~ dnorm(0, 0.0000001)", "beta_f ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta_f[2] <- 0", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c[s_cost[i] + 1])", "loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_sc[i] <- logdensity.bern(s_cost[i], pc[i])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
}
if(length(model_txt_info$model_e_random) == 0) {
  if(is.null(model_txt_info$se)) { model_string_jags <- gsub(" + inprod(X_e_random[i, ], a[, clus_e[i]])", "", model_string_jags, fixed = TRUE)}
  if(!is.null(model_txt_info$se)) { model_string_jags <- gsub(" + inprod(X_e_random[i, ], a[, clus_e[i]]) * ((s_eff[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)}
  model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_e_random[], mu_a_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
}
if(length(model_txt_info$model_c_random) == 0) {
  if(is.null(model_txt_info$sc)) { model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]])", "", model_string_jags, fixed = TRUE)}
  if(!is.null(model_txt_info$sc)) { model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]]) * ((s_cost[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)}
  model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_c_random[], mu_b_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e) * ((s_cost[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE)}
} 
if(length(model_txt_info$model_se_random) == 0) {
  if(!is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
  model_string_jags <- gsub(" + inprod(Z_e_random[i, ], g_e[, clus_se[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e_random[], mu_g_e_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g2_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(!is.null(model_txt_info$se) & is.null(model_txt_info$sc)) {
  model_string_jags <- gsub(" + inprod(Z_e_random[i, ], g_e[, clus_se[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e_random[], mu_g_e_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_sc_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects structural values mechanism", "", model_string_jags, fixed = TRUE) 
    if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | 
       length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# end mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }    
  }
}
if(length(model_txt_info$model_sc_random) == 0) {
  if(!is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
  model_string_jags <- gsub(" + inprod(Z_c_random[i, ], g_c[, clus_sc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c_random[], mu_g_c_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_se_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects structural values mechanism", "", model_string_jags, fixed = TRUE) 
    if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | 
       length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# end mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   } 
  } 
  if(is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
    model_string_jags <- gsub(" + inprod(Z_c_random[i, ], g_c[, clus_sc[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(mean_z_c_random[], mu_g_c_hat[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#priors on random effects structural values mechanism", "", model_string_jags, fixed = TRUE) 
    if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | 
       length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# end mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }  
  }
}
if(model_txt_info$ind_random) {
  if(!is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
  model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e) * ((s_cost[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  }
  if(is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e) * ((s_cost[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)   
  }
  if(!is.null(model_txt_info$se) & is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)  
  }
}
if(dist_c == "norm") {
  model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
} 
if(dist_c == "gamma") {
  if(!is.null(model_txt_info$sc)) {
  model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c[s_cost[i] + 1])", "cost[i] ~ dgamma(cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs", "ctau_c[i] <- cmu_c[i] / pow(s_c[s_cost[i] + 1], 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i] <- ", "log(cmu_c[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (s in 1:2) {#begin transformation costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[s] <- 1 / ss_c[s]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[s] <- s_c[s] * s_c[s]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end transformation costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c", "tmu_c <- tmu_ns_c * (1 - p_c) + exp(sc) * p_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c[s_cost[i] + 1])", "loglik_c[i] <- logdensity.gamma(cost[i], cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_c_random) == 0) {
    model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1])", "tmu_ns_c <- exp(inprod(mean_cov_c_fixed[], beta[, 1]))", model_string_jags, fixed = TRUE)}
  if(length(model_txt_info$model_c_random) != 0) {
    model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1]) + inprod(mean_cov_c_random[], mu_b_hat[])", "tmu_ns_c <- exp(inprod(mean_cov_c_fixed[], beta[, 1]) + inprod(mean_cov_c_random[], mu_b_hat[]))", model_string_jags, fixed = TRUE)}
  }
  if(is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c)", "cost[i] ~ dgamma(cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs", "ctau_c[i] <- cmu_c[i] / pow(s_c, 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_c[i] <- ", "log(cmu_c[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c <- 1 / ss_c", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c <- s_c * s_c", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)", "loglik_c[i] <- logdensity.gamma(cost[i], cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
      model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[])", "tmu_ns_c <- exp(inprod(mean_cov_c_fixed[], beta[]))", model_string_jags, fixed = TRUE)}
    if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", "tmu_ns_c <- exp(inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[]))", model_string_jags, fixed = TRUE)}
  }
}
if(dist_c == "lnorm") {
    if(!is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c[s_cost[i] + 1])", "cost[i] ~ dlnorm(clmu_c[i], ltau_c[s_cost[i] + 1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_c[i] <- ", "clmu_c[i] <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c[s] <- 1 / ss_c[s]", "ltau_c[s] <- 1 / lss_c[s]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c[s] <- s_c[s] * s_c[s]", "lss_c[s] <- ls_c[s] * ls_c[s]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c", "tmu_c <- tmu_ns_c * (1 - p_c) + exp(sc) * p_c", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm", "tmu_ns_c <- exp(ltmu_ns_c + lss_c[1] / 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c[2] <- sdc", "ls_c[2] <- sdc", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#std for lnorm", "s_c <- sqrt(exp(2 * tmu_ns_c + lss_c[1]) * (exp(lss_c[1]) - 1))",model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c[s_cost[i] + 1])", "loglik_c[i] <- logdensity.lnorm(cost[i], clmu_c[i], ltau_c[s_cost[i] + 1])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
      model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1])", "ltmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1])", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1]) + inprod(mean_cov_c_random[], mu_b_hat[])", "ltmu_ns_c <- inprod(mean_cov_c_fixed[], beta[, 1]) + inprod(mean_cov_c_random[], mu_b_hat[])", model_string_jags, fixed = TRUE)
    }
    }
    if(is.null(model_txt_info$sc)) {
      model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c)", "cost[i] ~ dlnorm(clmu_c[i], ltau_c)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("cmu_c[i] <- ", "clmu_c[i] <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c <- 1 / ss_c", "ltau_c <- 1 / lss_c", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c <- s_c * s_c", "lss_c <- ls_c * ls_c", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tmu_c <- tmu_ns_c * (1 - p_c) + sc * p_c", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm", "tmu_c <- exp(ltmu_ns_c + lss_c / 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#std for lnorm", "s_c <- sqrt(exp(2 * tmu_ns_c + lss_c) * (exp(lss_c) - 1))",model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)", "loglik_c[i] <- logdensity.lnorm(cost[i], clmu_c[i], ltau_c)", model_string_jags, fixed = TRUE)
      if(length(model_txt_info$model_c_random) == 0) {
        model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[])", "ltmu_ns_c <- inprod(mean_cov_c_fixed[], beta[])", model_string_jags, fixed = TRUE)
      }
      if(length(model_txt_info$model_c_random) != 0) {
        model_string_jags <- gsub("tmu_ns_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", "ltmu_ns_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", model_string_jags, fixed = TRUE)
      }    
    }
}
if(dist_e == "norm") {
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
}
if(dist_e == "beta") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dbeta(cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- (cmu_e[i] * (1 - cmu_e[i]) / pow(s_e[s_eff[i] + 1], 2) - 1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] ~ dunif(0, sqrt(tmu_ns_e * (1 - tmu_ns_e)))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + ilogit(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.beta(eff[i], cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dbeta(cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- (cmu_e[i] * (1 - cmu_e[i]) / pow(s_e, 2) - 1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, sqrt(tmu_ns_e * (1 - tmu_ns_e)))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.beta(eff[i], cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }    
  }
}
if(dist_e == "gamma") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dgamma(cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- cmu_e[i] / pow(s_e[s_eff[i] + 1], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + exp(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.gamma(eff[i], cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dgamma(cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- cmu_e[i] / pow(s_e, 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.gamma(eff[i], cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }  
  }
} 
if(dist_e == "exp") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dexp(1 / cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] <- tmu_ns_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + exp(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.exp(eff[i], 1 / cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dexp(1 / cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- tmu_ns_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.exp(eff[i], 1 / cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }  
  }
}
if(dist_e == "weib") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dweib(ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- pow(s_e[s_eff[i] + 1] / cmu_e[i], - 1.086)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + exp(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.weib(eff[i], ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dweib(ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- pow(s_e / cmu_e[i], - 1.086)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.weib(eff[i], ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }  
  }
}
if(dist_e == "logis") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dlogis(cmu_e[i], tau_e[s_eff[i] + 1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "tau_e[s] <- 1 / sqrt((3 * ss_e[s]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.logis(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", model_string_jags, fixed = TRUE)
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dlogis(cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "tau_e <- 1 / sqrt((3 * ss_e) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.logis(eff[i], cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
  }
}
if(dist_e == "bern") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dbern(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] <- sqrt(tmu_ns_e * (1 - tmu_ns_e))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + ilogit(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.bern(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dbern(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- sqrt(tmu_ns_e * (1 - tmu_ns_e))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.bern(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }    
  }
}
if(dist_e == "pois") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dpois(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[1] <- sqrt(tmu_ns_e)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + exp(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.pois(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dpois(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- sqrt(tmu_ns_e)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.pois(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }  
  }
}
if(dist_e == "negbin") {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e[s_eff[i] + 1])", "eff[i] ~ dnegbin(tau_e[s_eff[i] + 1] / (tau_e[s_eff[i] + 1] + cmu_e[i]), tau_e[s_eff[i] + 1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (s in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e[s] <- 1 / ss_e[s]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e[s] <- s_e[s] * s_e[s]", "s_e[1] <- sqrt((tau_e[1] / (tau_e[1] + tmu_ns_e)) * tau_e[1]) / (1 - (tau_e[1] / (tau_e[1] + tmu_ns_e)))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e[2] <- sde", "tau_e[2] <- sde", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "tmu_e <- tmu_ns_e * (1 - p_e) + exp(se) * p_e", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e[s_eff[i] + 1])", "loglik_e[i] <- logdensity.negbin(eff[i], tau_e[s_eff[i] + 1] / (tau_e[s_eff[i] + 1] + cmu_e[i]), tau_e[s_eff[i] + 1])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[, 1]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dnegbin(tau_e / (tau_e + cmu_e[i]), tau_e)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e <- s_e * s_e", "s_e <- sqrt((tau_e / (tau_e + tmu_ns_e)) * tau_e) / (1 - (tau_e / (tau_e + tmu_ns_e)))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tmu_e <- tmu_ns_e * (1 - p_e) + se * p_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.negbin(eff[i], tau_e / (tau_e + cmu_e[i]), tau_e)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("tmu_ns_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_ns_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
}
if(dist_e %in% c("beta", "gamma", "weib", "bern", "pois", "exp")) {
  if(dist_c == "gamma"){
    model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
  }
}
if(length(model_txt_info$model_c_random) != 0 & model_txt_info$is_c_random_c & !model_txt_info$is_int_c_random_c) {
  if(!is.null(model_txt_info$sc)) {
    model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]]) * ((s_cost[i] - 1) * -1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(mean_cov_c_random[], mu_b_hat[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE)}
  }
  if(is.null(model_txt_info$sc)) {
    model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(mean_cov_c_random[], mu_b_hat[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE) }
  }
}
if(length(model_txt_info$model_e_random) != 0 & model_txt_info$pe_random == 1) {
  if(!is.null(model_txt_info$se)) {
    model_string_jags <- gsub("inprod(X_e_random[i, ], a[, clus_e[i]]) * ((s_eff[i] - 1) * -1)", "X_e_random[i] * a[clus_e[i]] * ((s_eff[i] - 1) * -1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_e_random[], mu_a_hat[])", "mean_cov_e_random * mu_a_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "#begin a priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }", "for(s in 1:n_clus_e) {a[s] ~ dnorm(mu_a_hat, tau_a_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end a priors effects", "#end a1 priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)", "mu_a_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_a_hat[j] ~ dunif(0, 100) }", "s_a_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]", "tau_a_hat <- 1 / ss_a_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }", "ss_a_hat <- s_a_hat * s_a_hat", model_string_jags, fixed = TRUE)
  }
  if(is.null(model_txt_info$se)) {
    model_string_jags <- gsub("inprod(X_e_random[i, ], a[, clus_e[i]])", "X_e_random[i] * a[clus_e[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_e_random[], mu_a_hat[])", "mean_cov_e_random * mu_a_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "#begin a priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }", "for(s in 1:n_clus_e) {a[s] ~ dnorm(mu_a_hat, tau_a_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end a priors effects", "#end a priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)", "mu_a_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_a_hat[j] ~ dunif(0, 100) }", "s_a_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]", "tau_a_hat <- 1 / ss_a_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }", "ss_a_hat <- s_a_hat * s_a_hat", model_string_jags, fixed = TRUE)
  }
}  
if(length(model_txt_info$model_c_random) != 0 & model_txt_info$pc_random == 1) {
  if(!is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i]]) * ((s_cost[i] - 1) * -1)", "X_c_random[i] * b[clus_c[i]] * ((s_cost[i] - 1) * -1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_c_random[], mu_b_hat[])", "mean_cov_c_random * mu_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "#begin b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "for(s in 1:n_clus_c) {b[s] ~ dnorm(mu_b_hat, tau_b_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "#end b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "mu_b_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "s_b_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "tau_b_hat <- 1 / ss_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "ss_b_hat <- s_b_hat * s_b_hat", model_string_jags, fixed = TRUE)
  }
  if(is.null(model_txt_info$sc)) {
    model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i]])", "X_c_random[i] * b[clus_c[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_c_random[], mu_b_hat[])", "mean_cov_c_random * mu_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "#begin b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "for(s in 1:n_clus_c) {b[s] ~ dnorm(mu_b_hat, tau_b_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "#end b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "mu_b_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "s_b_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "tau_b_hat <- 1 / ss_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "ss_b_hat <- s_b_hat * s_b_hat", model_string_jags, fixed = TRUE)
  }
}
if(!is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
    if(model_txt_info$ze_fixed == 1) {
    inprod_e <- "Z_e_fixed[i] * gamma_e"
    inprod_mean_e <- "mean_z_e_fixed * gamma_e"
    begin_prior_gamma <- "#begin gamma priors effects"
    prior_gamma_e <- "gamma_e ~ dnorm(0, 0.01)"
    end_prior_gamma <- "#end gamma priors effects"
    model_string_jags <- gsub("inprod(Z_e_fixed[i, ], gamma_e[])", inprod_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_z_e_fixed[], gamma_e[])", inprod_mean_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:ze_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("gamma_e[j] ~ dnorm(0, 0.01)", prior_gamma_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
    }
  if(model_txt_info$ze_random == 1 & length(model_txt_info$model_se_random) != 0) {
    model_string_jags <- gsub("inprod(Z_e_random[i, ], g_e[, clus_se[i]])", "Z_e_random[i] * g_e[clus_se[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "tau_g_e_hat <- 1 / ss_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "ss_g_e_hat <- s_g_e_hat * s_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[])", "mean_z_e_random * mu_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "#begin g_e priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "for(s in 1:n_clus_se) {g_e[s] ~ dnorm(mu_g_e_hat, tau_g_e_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_e priors effects", "#end g_e priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "mu_g_e_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "s_g_e_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  }  
   if(model_txt_info$zc_fixed == 1) {
      inprod_c <- "Z_c_fixed[i] * gamma_c"
      inprod_mean_c <- "mean_z_c_fixed * gamma_c"
      begin_prior_gamma <- "#begin gamma priors costs"
      prior_gamma_c <- "gamma_c ~ dnorm(0, 0.01)"
      end_prior_gamma <- "#end gamma priors costs"
      model_string_jags <- gsub("inprod(Z_c_fixed[i, ], gamma_c[])", inprod_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c_fixed[], gamma_c[])", inprod_mean_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[j] ~ dnorm(0, 0.01)", prior_gamma_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
    }
  if(model_txt_info$zc_random == 1 & length(model_txt_info$model_sc_random) != 0) {
    model_string_jags <- gsub("inprod(Z_c_random[i, ], g_c[, clus_sc[i]])", "Z_c_random[i] * g_c[clus_sc[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "tau_g_c_hat <- 1 / ss_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "ss_g_c_hat <- s_g_c_hat * s_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[])", "mean_z_c_random * mu_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "#begin g_c priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "for(s in 1:n_clus_sc) {g_c[s] ~ dnorm(mu_g_c_hat, tau_g_c_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_c priors costs", "#end g_c priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "mu_g_c_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "s_g_c_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  }
}
if(is.null(model_txt_info$se) & !is.null(model_txt_info$sc)) {
  if(model_txt_info$zc_fixed == 1) {
  inprod_c <- "Z_c_fixed[i] * gamma_c"
  inprod_mean_c <- "ilogit(mean_z_c_fixed * gamma_c)"
  begin_prior_gamma <- "#begin gamma priors costs"
  prior_gamma_c <- "gamma_c ~ dnorm(0, 0.01)"
  end_prior_gamma <- "#end gamma priors costs"
  model_string_jags <- gsub("inprod(Z_c_fixed[i, ], gamma_c[])", inprod_c, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ilogit(inprod(mean_z_c_fixed[], gamma_c[]))", inprod_mean_c, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_c[j] ~ dnorm(0, 0.01)", prior_gamma_c, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
  }
  if(model_txt_info$zc_random == 1 & length(model_txt_info$model_sc_random) != 0) {
    model_string_jags <- gsub("inprod(Z_c_random[i, ], g_c[, clus_sc[i]])", "Z_c_random[i] * g_c[clus_sc[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "tau_g_c_hat <- 1 / ss_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "ss_g_c_hat <- s_g_c_hat * s_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[])", "mean_z_c_random * mu_g_c_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "#begin g_c priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_sc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "for(s in 1:n_clus_sc) {g_c[s] ~ dnorm(mu_g_c_hat, tau_g_c_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_c priors costs", "#end g_c priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "mu_g_c_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "s_g_c_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  }
}
if(!is.null(model_txt_info$se) & is.null(model_txt_info$sc)) {
  if(model_txt_info$ze_fixed == 1) {
    inprod_e <- "Z_e_fixed[i] * gamma_e"
    inprod_mean_e <- "ilogit(mean_z_e_fixed * gamma_e)"
    begin_prior_gamma <- "#begin gamma priors effects"
    prior_gamma_e <- "gamma_e ~ dnorm(0, 0.01)"
    end_prior_gamma <- "#end gamma priors effects"
    model_string_jags <- gsub("inprod(Z_e_fixed[i, ], gamma_e[])", inprod_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ilogit(inprod(mean_z_e_fixed[], gamma_e[]))", inprod_mean_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:ze_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("gamma_e[j] ~ dnorm(0, 0.01)", prior_gamma_e, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
  }
  if(model_txt_info$ze_random == 1 & length(model_txt_info$model_se_random) != 0) {
    model_string_jags <- gsub("inprod(Z_e_random[i, ], g_e[, clus_se[i]])", "Z_e_random[i] * g_e[clus_se[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "tau_g_e_hat <- 1 / ss_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "ss_g_e_hat <- s_g_e_hat * s_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[])", "mean_z_e_random * mu_g_e_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "#begin g_e priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_se) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "for(s in 1:n_clus_se) {g_e[s] ~ dnorm(mu_g_e_hat, tau_g_e_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_e priors effects", "#end g_e priors effects", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "mu_g_e_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "s_g_e_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  }
}
if(model_txt_info$ind_fixed) {
  if(!is.null(model_txt_info$sc)) {
    model_string_jags <- gsub(" + beta_f[s_cost[i] + 1] * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[2] <- 0", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } else if(is.null(model_txt_info$sc)) {
    model_string_jags <- gsub(" + beta_f * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  }
}
model_string_jags <- prior_hurdle(type = type, dist_e = dist_e, dist_c = dist_c, 
          model_txt_info = model_txt_info, model_string_jags = model_string_jags)
writeLines(model_string_jags, "hurdle.txt")
model_string <- "hurdle.txt"
return(model_string)
}