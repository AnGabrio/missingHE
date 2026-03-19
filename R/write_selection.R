#'An internal function to select which type of selection model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Selection models
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param model_txt_info list containing model specification information used to write the txt file of the JAGS model.
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_selection <- function(dist_e , dist_c, type, model_txt_info) {
  model_string_jags<-  "
  model {

  #full sample
  for(i in 1:n) {
  #costs and effects model
  cost[i] ~ dnorm(cmu_c[i], tau_c)
  eff[i] ~ dnorm(cmu_e[i], tau_e)

  #derive mean and std effects 
  #derive mean and std costs

  #mean regression
  cmu_c[i] <- inprod(X_c_fixed[i, ], beta[]) + beta_f * (eff[i] - tmu_e) + inprod(X_c_random[i, ], b[, clus_c[i]]) + b_f[clus_c[i]] * (eff[i] - tmu_e)
  cmu_e[i] <- inprod(X_e_fixed[i, ], alpha[]) + inprod(X_e_random[i, ], a[, clus_e[i]])   

  #missing data mechanism
  m_eff[i] ~ dbern(pq[i])
  logit(pq[i]) <- inprod(Z_e_fixed[i, ], gamma_e[]) + delta_e * eff[i] + inprod(Z_e_random[i, ], g_e[, clus_me[i]]) + d_e[clus_me[i]] * eff[i]
  m_cost[i] ~ dbern(pc[i])
  logit(pc[i]) <- inprod(Z_c_fixed[i, ], gamma_c[]) + delta_c * cost[i] + inprod(Z_c_random[i, ], g_c[, clus_mc[i]]) + d_c[clus_mc[i]] * cost[i]

  #loglikelihood
  loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)
  loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)
  loglik_me[i] <- logdensity.bern(m_eff[i], pq[i])
  loglik_mc[i] <- logdensity.bern(m_cost[i], pc[i])
  }
  
  #transformation of parameters
  #begin transformation
  tau_c <- 1 / ss_c
  ss_c <- s_c * s_c
  #std for lnorm 
  #mean for lnorm
  tau_e <- 1 / ss_e
  ss_e <- s_e * s_e
  # end transformation
  
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
  tau_d_c_hat <- 1 / ss_d_c_hat
  ss_d_c_hat <- s_d_c_hat * s_d_c_hat
  tau_d_e_hat <- 1 / ss_d_e_hat
  ss_d_e_hat <- s_d_e_hat * s_d_e_hat
  #end transformation of random effects 
  
  #missingness probability
  p_c <- ilogit(inprod(mean_z_c_fixed[], gamma_c[]) + delta_c * mean(cost[]) + inprod(mean_z_c_random[], mu_g_c_hat[]) + mu_d_c_hat * mean(cost[]))
  p_e <- ilogit(inprod(mean_z_e_fixed[], gamma_e[]) + delta_e * mean(eff[]) + inprod(mean_z_e_random[], mu_g_e_hat[]) + mu_d_e_hat * mean(eff[]))

  #calculate means at mean of covariates
  tmu_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[]) 
  tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])

  #priors
  
  #priors for mean regression coefficients
  for (j in 1:pe_fixed) {#begin alpha priors effects
  alpha[j] ~ dnorm(0, 0.0000001)
  }#end alpha priors effects

  for (j in 1:pc_fixed) {#begin beta priors costs
  beta[j] ~ dnorm(0, 0.0000001)
  }#end beta priors costs

  #priors for mean regression random coefficients
  for (j in 1:pe_random) {#begin a priors effects
  for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }
  }#end a priors effects

  for (j in 1:pc_random) {#begin b priors costs
  for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }
  }#end b priors costs
  
  #standard deviation priors
  s_c ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_e ~ dt(0, pow(2.5, -2), 1)T(0,)

  #correlation
  beta_f ~ dnorm(0, 0.0000001)
  
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
  
  #priors on missing data mechanism
  for (j in 1:ze_fixed) {#begin gamma priors effects
  gamma_e[j] ~ dnorm(0, 0.01)
  }#end gamma priors effects

  for (j in 1:zc_fixed) {#begin gamma priors costs
  gamma_c[j] ~ dnorm(0, 0.01)
  }#end gamma priors costs

  #priors on random effects missing data mechanism
  for (j in 1:ze_random) {#begin g_e priors effects
  for(s in 1:n_clus_me) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }
  }#end g_e priors effects

  for (j in 1:zc_random) {#begin g_c priors costs
  for(s in 1:n_clus_mc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }
  }#end g_c priors costs
  
  #mnar parameters
  # begin mnar priors
  delta_e ~ dnorm(0, 1)
  delta_c ~ dnorm(0, 1)
  #end mnar priors

  #mnar random effects priors
  for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }
  for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }

  # mean and sd mean mnar random effects priors
  #begin mnar random effects priors
  mu_d_e_hat ~ dnorm(0, 1)
  mu_d_c_hat ~ dnorm(0, 1)
  s_d_e_hat ~ dunif(0, 1)
  s_d_c_hat ~ dunif(0, 1) 
  # end mnar random effects priors

  }
  "
 if(length(model_txt_info$model_e_random) == 0) {
  model_string_jags <- gsub(" + inprod(X_e_random[i, ], a[, clus_e[i]])", "", model_string_jags, fixed = TRUE)
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
  model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_c_random[], mu_b_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE)}
  } else if(length(model_txt_info$model_c_random) != 0 & model_txt_info$is_c_random_c & !model_txt_info$is_int_c_random_c) {
    model_string_jags <- gsub(" + inprod(X_c_random[i, ], b[, clus_c[i]])", "", model_string_jags, fixed = TRUE)
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
  if(length(model_txt_info$model_me_random) == 0) {
  model_string_jags <- gsub(" + inprod(Z_e_random[i, ], g_e[, clus_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e_random[], mu_g_e_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d_e[clus_me[i]] * eff[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_e_hat * mean(eff[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_e_hat <- 1 / ss_d_e_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_e_hat <- s_d_e_hat * s_d_e_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_e_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_e_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_mc_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_txt_info$model_me_random) != 0 & model_txt_info$is_me_random_e & !model_txt_info$is_int_me_random_e) {
  model_string_jags <- gsub(" + inprod(Z_e_random[i, ], g_e[, clus_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e_random[], mu_g_e_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_txt_info$model_mc_random) == 0) {
  model_string_jags <- gsub(" + inprod(Z_c_random[i, ], g_c[, clus_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c_random[], mu_g_c_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d_c[clus_mc[i]] * cost[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_c_hat * mean(cost[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_c_hat <- 1 / ss_d_c_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_c_hat <- s_d_c_hat * s_d_c_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_c_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_c_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_me_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    if(length(model_txt_info$model_e_random) == 0 & length(model_txt_info$model_c_random) == 0 | length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_txt_info$model_mc_random) != 0 & model_txt_info$is_mc_random_c & !model_txt_info$is_int_mc_random_c) {
  model_string_jags <- gsub(" + inprod(Z_c_random[i, ], g_c[, clus_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c_random[], mu_g_c_hat[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(type == "MAR") {
  model_string_jags <- gsub(" + delta_e * mean(eff[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e * eff[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_e ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c * mean(cost[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c * cost[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_c ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# begin mnar priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#end mnar priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mnar parameters", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_me_random) != 0) {
    model_string_jags <- gsub("+ d_e[clus_me[i]] * eff[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat * mean(eff[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat <- 1 / ss_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat <- s_d_e_hat * s_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_mc_random) != 0) {
    model_string_jags <- gsub("+ d_c[clus_mc[i]] * cost[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat * mean(cost[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat <- 1 / ss_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat <- s_d_c_hat * s_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }", "", model_string_jags, fixed = TRUE)
    }
  } else if(type == "MNAR_eff") {
  model_string_jags <- gsub(" + delta_c * mean(cost[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c * cost[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_c ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_mc_random) != 0) {
    model_string_jags <- gsub("+ d_c[clus_mc[i]] * cost[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat * mean(cost[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat <- 1 / ss_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat <- s_d_c_hat * s_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
    model_string_jags <- gsub("+ d_e[clus_me[i]] * eff[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat * mean(eff[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat <- 1 / ss_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat <- s_d_e_hat * s_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR_cost") {
  model_string_jags <- gsub(" + delta_e * mean(eff[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e * eff[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_e ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_me_random) != 0) {
    model_string_jags <- gsub("+ d_e[clus_me[i]] * eff[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat * mean(eff[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat <- 1 / ss_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat <- s_d_e_hat * s_d_e_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random)) {
    model_string_jags <- gsub("+ d_c[clus_mc[i]] * cost[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat * mean(cost[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat <- 1 / ss_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat <- s_d_c_hat * s_d_c_hat", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR") {
    if(length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("+ d_e[clus_me[i]] * eff[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat * mean(eff[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_e_hat <- 1 / ss_d_e_hat", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_e_hat <- s_d_e_hat * s_d_e_hat", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_e_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_e_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s] ~ dnorm(mu_d_e_hat, tau_d_e_hat) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random)) {
      model_string_jags <- gsub("+ d_c[clus_mc[i]] * cost[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat * mean(cost[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat <- 1 / ss_d_c_hat", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat <- s_d_c_hat * s_d_c_hat", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s] ~ dnorm(mu_d_c_hat, tau_d_c_hat) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  }
  if(model_txt_info$ind_random) {
  model_string_jags <- gsub("+ b_f[clus_c[i]] * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_f_hat ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat <- 1 / ss_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat <- s_b_f_hat * s_b_f_hat", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  } 
  if(model_txt_info$ind_fixed) {
  model_string_jags <- gsub(" + beta_f * (eff[i] - tmu_e)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta_f ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
  if(dist_c == "norm") {
  model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
  model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c)", "cost[i] ~ dgamma(cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs", "ctau_c[i] <- cmu_c[i] / pow(s_c, 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i] <- ", "log(cmu_c[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c <- 1 / ss_c", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c <- s_c * s_c", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)", "loglik_c[i] <- logdensity.gamma(cost[i], cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
    model_string_jags <- gsub("tmu_c <- inprod(mean_cov_c_fixed[], beta[])", "tmu_c <- exp(inprod(mean_cov_c_fixed[], beta[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) != 0) {
    model_string_jags <- gsub("tmu_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", "tmu_c <- exp(inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_c == "lnorm") {
  model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c)", "cost[i] ~ dlnorm(clmu_c[i], ltau_c)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i] <- ", "clmu_c[i] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c <- 1 / ss_c", "ltau_c <- 1 / lss_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c <- s_c * s_c", "lss_c <- ls_c * ls_c", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "s_c <- sqrt(exp(2 * tlmu_c + lss_c) * (exp(lss_c) - 1))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c ~ dunif(0, 10)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "tmu_c <- exp(tlmu_c + lss_c / 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c)", "loglik_c[i] <- logdensity.lnorm(cost[i], clmu_c[i], ltau_c)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
    model_string_jags <- gsub("tmu_c <- inprod(mean_cov_c_fixed[], beta[])", "tlmu_c <- inprod(mean_cov_c_fixed[], beta[])", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) != 0) {
    model_string_jags <- gsub("tmu_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", "tlmu_c <- inprod(mean_cov_c_fixed[], beta[]) + inprod(mean_cov_c_random[], mu_b_hat[])", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e == "norm") {
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "beta") {
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dbeta(cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- (cmu_e[i] * (1 - cmu_e[i]) / pow(s_e, 2) - 1)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, sqrt(tmu_e * (1 - tmu_e)))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.beta(eff[i], cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "gamma"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dgamma(cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- cmu_e[i] / pow(s_e, 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.gamma(eff[i], cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "exp"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dexp(1 / cmu_e[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- tmu_e", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.exp(eff[i], 1/cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "weib"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dweib(ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- pow(s_e / cmu_e[i], - 1.086)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.weib(eff[i], ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "logis"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dlogis(cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "tau_e <- 1 / sqrt((3 * ss_e) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.logis(eff[i], cmu_e[i], tau_e)", model_string_jags, fixed = TRUE)
  } else if(dist_e == "bern"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dbern(cmu_e[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- sqrt(tmu_e * (1 - tmu_e))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.bern(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- ilogit(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "pois"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dpois(cmu_e[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e <- sqrt(tmu_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.pois(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "negbin"){
  model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e)", "eff[i] ~ dnegbin(tau_e / (tau_e + cmu_e[i]), tau_e)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e <- 1 / ss_e", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e <- s_e * s_e", "s_e <- sqrt((tau_e / (tau_e + tmu_e)) * tau_e) / (1 - (tau_e / (tau_e + tmu_e)))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e)", "loglik_e[i] <- logdensity.negbin(eff[i], tau_e / (tau_e + cmu_e[i]), tau_e)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e <- inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[])", "tmu_e <- exp(inprod(mean_cov_e_fixed[], alpha[]) + inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e %in% c("beta", "gamma", "weib", "bern", "pois", "exp")) {
     if(dist_c == "gamma"){
     model_string_jags <- gsub("#begin transformation", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("# end transformation", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
   }
   if(length(model_txt_info$model_e_random) != 0 & model_txt_info$pe_random == 1) {
       model_string_jags <- gsub("inprod(X_e_random[i, ], a[, clus_e[i]])", "X_e_random[i] * a[clus_e[i]]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e_random[], mu_a_hat[])", "mean_cov_e_random * mu_a_hat", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "#begin a priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }", 
                                 "for(s in 1:n_clus_e) {a[s] ~ dnorm(mu_a_hat, tau_a_hat) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a priors effects", "#end a priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)", "mu_a_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_a_hat[j] ~ dunif(0, 100) }", "s_a_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]", "tau_a_hat <- 1 / ss_a_hat", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }", "ss_a_hat <- s_a_hat * s_a_hat", model_string_jags, fixed = TRUE)
       }
   if(length(model_txt_info$model_c_random) != 0 & model_txt_info$pc_random == 1) {
   model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i]])", "X_c_random[i] * b[clus_c[i]]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c_random[], mu_b_hat[])", "mean_cov_c_random * mu_b_hat", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "#begin b priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", 
                             "for(s in 1:n_clus_c) {b[s] ~ dnorm(mu_b_hat, tau_b_hat) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b priors costs", "#end b priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "mu_b_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "s_b_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "tau_b_hat <- 1 / ss_b_hat", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "ss_b_hat <- s_b_hat * s_b_hat", model_string_jags, fixed = TRUE)
   }
   if(model_txt_info$ze_fixed == 1) {
   inprod_e <- "Z_e_fixed[i] * gamma_e"
     if(type %in% c("MAR", "MNAR_cost")) {
     inprod_mean_e <- "mean_z_e_fixed * gamma_e"
     model_string_jags <- gsub("inprod(mean_z_e_fixed[], gamma_e[])", 
                               inprod_mean_e, model_string_jags, fixed = TRUE)
     }
     if(type %in% c("MNAR_eff", "MNAR")) {
     inprod_mean_e <- "mean_z_e_fixed * gamma_e + delta_e * mean(eff[])"
     model_string_jags <- gsub("inprod(mean_z_e_fixed[], gamma_e[]) + delta_e * mean(eff[])", 
                               inprod_mean_e, model_string_jags, fixed = TRUE)
     }
   begin_prior_gamma <- "#begin gamma priors effects"
   prior_gamma_e <- "gamma_e ~ dnorm(0, 0.01)"
   end_prior_gamma <- "#end gamma priors effects"
   model_string_jags <- gsub("inprod(Z_e_fixed[i, ], gamma_e[])", inprod_e, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:ze_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_e[j] ~ dnorm(0, 0.01)", prior_gamma_e, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
   } 
   if(model_txt_info$ze_random == 1 & length(model_txt_info$model_me_random) != 0) {
      if(type %in% c("MAR", "MNAR_cost")) {
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[])", "mean_z_e_random * mu_g_e_hat", 
                                model_string_jags, fixed = TRUE)
      }
      if(type %in% c("MNAR_eff", "MNAR")) {
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[]) + mu_d_e_hat * mean(eff[])", "mean_z_e_random * mu_g_e_hat + mu_d_e_hat * mean(eff[])", 
                                model_string_jags, fixed = TRUE)
      }
      model_string_jags <- gsub("inprod(Z_e_random[i, ], g_e[, clus_me[i]])", "Z_e_random[i] * g_e[clus_me[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j] <- 1 / ss_g_e_hat[j]", "tau_g_e_hat <- 1 / ss_g_e_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_e_hat[j] <- s_g_e_hat[j] * s_g_e_hat[j] }", "ss_g_e_hat <- s_g_e_hat * s_g_e_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[])", "mean_z_e_random * mu_g_e_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "#begin g_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s] ~ dnorm(mu_g_e_hat[j], tau_g_e_hat[j]) }", "for(s in 1:n_clus_me) {g_e[s] ~ dnorm(mu_g_e_hat, tau_g_e_hat) }", 
                                model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g_e priors effects", "#end g_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j] ~ dnorm(0, 0.001)", "mu_g_e_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_e_hat[j] ~ dunif(0, 100) }", "s_g_e_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
   if(model_txt_info$zc_fixed == 1) {
      inprod_c <- "Z_c_fixed[i] * gamma_c"
          if(type %in% c("MAR", "MNAR_eff")) {
          inprod_mean_c <- "mean_z_c_fixed * gamma_c"
          model_string_jags <- gsub("inprod(mean_z_c_fixed[], gamma_c[])", inprod_mean_c, model_string_jags, fixed = TRUE)
          }
          if(type %in% c("MNAR_cost", "MNAR")) {
          inprod_mean_c <- "mean_z_c_fixed * gamma_c + delta_c * mean(cost[])"
          model_string_jags <- gsub("inprod(mean_z_c_fixed[], gamma_c[]) + delta_c * mean(cost[])", inprod_mean_c, model_string_jags, fixed = TRUE)
          }
      begin_prior_gamma <- "#begin gamma priors costs"
      prior_gamma_c <- "gamma_c ~ dnorm(0, 0.01)"
      end_prior_gamma <- "#end gamma priors costs"
      model_string_jags <- gsub("inprod(Z_c_fixed[i, ], gamma_c[])", inprod_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[j] ~ dnorm(0, 0.01)", prior_gamma_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
      }
      if(model_txt_info$zc_random == 1 & length(model_txt_info$model_mc_random) != 0) {
         if(type %in% c("MAR", "MNAR_eff")) {
         model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[])", "mean_z_c_random * mu_g_c_hat", model_string_jags, fixed = TRUE)
         }
         if(type %in% c("MNAR_cost", "MNAR")) {
         model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[]) + mu_d_c_hat * mean(cost[])", "mean_z_c_random * mu_g_c_hat + mu_d_c_hat * mean(cost[])", 
                                   model_string_jags, fixed = TRUE)
         }
      model_string_jags <- gsub("inprod(Z_c_random[i, ], g_c[, clus_mc[i]])", "Z_c_random[i] * g_c[clus_mc[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j] <- 1 / ss_g_c_hat[j]", "tau_g_c_hat <- 1 / ss_g_c_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_c_hat[j] <- s_g_c_hat[j] * s_g_c_hat[j] }", "ss_g_c_hat <- s_g_c_hat * s_g_c_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[])", "mean_z_c_random * mu_g_c_hat", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "#begin g_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s] ~ dnorm(mu_g_c_hat[j], tau_g_c_hat[j]) }", "for(s in 1:n_clus_mc) {g_c[s] ~ dnorm(mu_g_c_hat, tau_g_c_hat) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g_c priors costs", "#end g_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j] ~ dnorm(0, 0.001)", "mu_g_c_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_c_hat[j] ~ dunif(0, 100) }", "s_g_c_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
      model_string_jags <- prior_selection(type = type, dist_e = dist_e, dist_c = dist_c, 
                                           model_txt_info = model_txt_info, model_string_jags = model_string_jags)
   writeLines(model_string_jags, "selection.txt")
   model_string <- "selection.txt"
   return(model_string)
}