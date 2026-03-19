#'An internal function to select which type of longitudinal missing data model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Longitudinal missing data models
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


write_lmdm <- function(dist_e , dist_c, type, model_txt_info) {
  model_string_jags<-  "
  model {

  #full sample
  for(i in 1:n) {
  for(time in 1:max_time) {# loop through time
  #costs and effects model by time point
  cost[i, time] ~ dnorm(cmu_c[i, time], tau_c[time])
  eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])
  }#end loop time

  #derive mean and std effects
  #derive mean and std costs

  #mean regression at time one
  cmu_c[i, 1] <- inprod(X_c_fixed[i, ], beta[, 1]) + beta_f[1] * (eff[i, 1] - tmu_e[1]) + inprod(X_c_random[i, ], b[, clus_c[i], 1]) + b_f[clus_c[i], 1] * (eff[i, 1] - tmu_e[1]) 
  cmu_e[i, 1] <- inprod(X_e_fixed[i, ], alpha[, 1]) + inprod(X_e_random[i, ], a[, clus_e[i], 1])   

  #mean regression at later times
  for(time in 2:max_time) {# loop through time
  cmu_c[i, time] <- inprod(X_c_fixed[i, ], beta[, time]) + beta_f[time] * (eff[i, time] - tmu_e[time]) + inprod(X_c_random[i, ], b[, clus_c[i], time]) + b_f[clus_c[i], time] * (eff[i, time] - tmu_e[time]) + beta_te[time] * (eff[i, time - 1] - tmu_e[time - 1]) + beta_tc[time] * (cost[i, time - 1] - tmu_c[time - 1]) + b_te[clus_c[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + b_tc[clus_c[i], time] * (cost[i, time - 1] - tmu_c[time - 1]) 
  cmu_e[i, time] <- inprod(X_e_fixed[i, ], alpha[, time]) + inprod(X_e_random[i, ], a[, clus_e[i], time]) + alpha_te[time] * (eff[i, time - 1] - tmu_e[time - 1]) + alpha_tc[time] * (cost[i, time - 1] - tmu_c[time - 1]) + a_te[clus_e[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + a_tc[clus_e[i], time] * (cost[i, time - 1] - tmu_c[time - 1])
  }#end loop time


  #missing data mechanism
  for(time in 1:max_time) {# loop through time
  m_eff[i, time, 1:n_patt_e] ~ dmulti(pq[i, time, 1:n_patt_e], 1)
  m_cost[i, time, 1:n_patt_c] ~ dmulti(pc[i, time, 1:n_patt_c], 1)
  for(mdrop in 1:n_patt_e){# loop through dropoute
  pq[i, time, mdrop] <- phiq[i, time, mdrop]/sum(phiq[i, time, ])
  }#end loop dropoute
  for(mdrop in 1:n_patt_c){# loop through dropoutc
  pc[i, time, mdrop] <- phic[i, time, mdrop]/sum(phic[i, time, ])
  }#end loop dropoutc
  }#end loop time
  
  #log missing regression at time one
  for(mdrop in 1:n_patt_e){# loop through dropoute
  log(phiq[i, 1, mdrop]) <- inprod(Z_e_fixed[i, ], gamma_e[, 1, mdrop]) + delta_e[1, mdrop] * eff[i, 1] + inprod(Z_e_random[i, ], g_e[, clus_me[i], 1]) + d_e[clus_me[i], 1] * eff[i, 1]
  }#end loop dropoute
  for(mdrop in 1:n_patt_c){# loop through dropoutc
  log(phic[i, 1, mdrop]) <- inprod(Z_c_fixed[i, ], gamma_c[, 1, mdrop]) + delta_c[1, mdrop] * cost[i, 1] + inprod(Z_c_random[i, ], g_c[, clus_mc[i], 1]) + d_c[clus_mc[i], 1] * cost[i, 1]
  }#end loop dropoutc

  #log missing regression at later times
  for(time in 2:max_time) {# loop through time
  for(mdrop in 1:n_patt_e){# loop through dropoute
  log(phiq[i, time, mdrop]) <- inprod(Z_e_fixed[i, ], gamma_e[, time, mdrop]) + delta_e[time, mdrop] * eff[i, time] + inprod(Z_e_random[i, ], g_e[, clus_me[i], time]) + d_e[clus_me[i], time] * eff[i, time]
  }#end loop dropoute
  
  for(mdrop in 1:n_patt_c){# loop through dropoutc
  log(phic[i, time, mdrop]) <- inprod(Z_c_fixed[i, ], gamma_c[, time, mdrop]) + delta_c[time, mdrop] * cost[i, time] + inprod(Z_c_random[i, ], g_c[, clus_mc[i], time]) + d_c[clus_mc[i], time] * cost[i, time]
  }#end loop dropoutc
  }#end loop time
  
  #loglikelihood
  for(time in 1:max_time) {# loop through time
  loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])
  loglik_c[i, time] <- logdensity.norm(cost[i, time], cmu_c[i, time], tau_c[time])
  loglik_me[i, time] <- logdensity.multi(m_eff[i, time, 1:n_patt_e], pq[i, time, 1:n_patt_e], 1)
  loglik_mc[i, time] <- logdensity.multi(m_cost[i, time, 1:n_patt_c], pc[i, time, 1:n_patt_c], 1)
  }#end loop time
  }
  
  #transformation of parameters
  #begin transformation of params
  for(time in 1:max_time) {# loop through time transformation of params
  tau_c[time] <- 1 / ss_c[time]
  ss_c[time] <- s_c[time] * s_c[time]
  #std for lnorm 
  #mean for lnorm
  tau_e[time] <- 1 / ss_e[time]
  ss_e[time] <- s_e[time] * s_e[time]
  }#end loop time transformation of params
  # end transformation of params
  
  #transformation of random effects parameters
  # begin transformation random effects
  for(time in 1:max_time) {# loop through time transformation of random effects
  for(j in 1:pc_random) {tau_b_hat[j, time] <- 1 / ss_b_hat[j, time]
  ss_b_hat[j, time] <- s_b_hat[j, time] * s_b_hat[j, time] }
  for(j in 1:pe_random) {tau_a_hat[j, time] <- 1 / ss_a_hat[j, time]
  ss_a_hat[j, time] <- s_a_hat[j, time] * s_a_hat[j, time] }
  for(j in 1:zc_random) {tau_g_c_hat[j, time] <- 1 / ss_g_c_hat[j, time]
  ss_g_c_hat[j, time] <- s_g_c_hat[j, time] * s_g_c_hat[j, time] }
  for(j in 1:ze_random) {tau_g_e_hat[j, time] <- 1 / ss_g_e_hat[j, time]
  ss_g_e_hat[j, time] <- s_g_e_hat[j, time] * s_g_e_hat[j, time] }
  tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]
  ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]
  tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]
  ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]
  }#end loop time transformation of random effects
  #end transformation of random effects 
  
  #missingness probability
  for(time in 1:max_time) {# loop through time
  for(mdrop in 1:n_patt_c){# loop through mpatternc
  p_c[time, mdrop] <- exp(inprod(mean_z_c_fixed[], gamma_c[, time, mdrop]) + delta_c[time, mdrop] * mean(cost[, time]) + inprod(mean_z_c_random[], mu_g_c_hat[, time]) + mu_d_c_hat[time] * mean(cost[, time]))
  } # end loop mpatternc
  for(mdrop in 1:n_patt_e){# loop through mpatterne
  p_e[time, mdrop] <- exp(inprod(mean_z_e_fixed[], gamma_e[, time, mdrop]) + delta_e[time, mdrop] * mean(eff[, time]) + inprod(mean_z_e_random[], mu_g_e_hat[, time]) + mu_d_e_hat[time] * mean(eff[, time]))
  } # end loop mpatterne
  }#end loop time
  
  #calculate means at mean of covariates
  for(time in 1:max_time) {# loop through time
  tmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time]) + inprod(mean_cov_c_random[], mu_b_hat[, time]) 
  tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])
  }#end loop time

  #priors
  
  #priors for mean regression coefficients
  for(time in 1:max_time) {# loop through time alpha
  for (j in 1:pe_fixed) {#begin alpha priors effects
  alpha[j, time] ~ dnorm(0, 0.0000001)
  }#end alpha priors effects
  }#end alpha loop time

  for(time in 1:max_time) {# loop through time beta
  for (j in 1:pc_fixed) {#begin beta priors costs
  beta[j, time] ~ dnorm(0, 0.0000001)
  }#end beta priors costs
  }#end beta loop time

  #priors for mean regression random coefficients
  for(time in 1:max_time) {# loop through time a 
  for (j in 1:pe_random) {#begin a priors effects
  for(s in 1:n_clus_e) {a[j, s, time] ~ dnorm(mu_a_hat[j, time], tau_a_hat[j, time]) }
  }#end a priors effects
  }#end a loop time

  for(time in 1:max_time) {# loop through time b 
  for (j in 1:pc_random) {#begin b priors costs
  for(s in 1:n_clus_c) {b[j, s, time] ~ dnorm(mu_b_hat[j, time], tau_b_hat[j, time]) }
  }#end b priors costs
  }#end b loop time

  #standard deviation priors
  for(time in 1:max_time) {# loop through time standard deviation priors
  s_c[time] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)

  #correlation
  beta_f[time] ~ dnorm(0, 0.0000001)
  
  #time dependence
  beta_te[time] ~ dnorm(0, 0.0000001)
  beta_tc[time] ~ dnorm(0, 0.0000001)
  alpha_te[time] ~ dnorm(0, 0.0000001)
  alpha_tc[time] ~ dnorm(0, 0.0000001)

  # mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j, time] ~ dnorm(0, 0.001)
  s_b_hat[j, time] ~ dunif(0, 100) }
  for(j in 1:pe_random) {mu_a_hat[j, time] ~ dnorm(0, 0.001)
  s_a_hat[j, time] ~ dunif(0, 100) }
  for(j in 1:zc_random) {mu_g_c_hat[j, time] ~ dnorm(0, 0.001)
  s_g_c_hat[j, time] ~ dunif(0, 100) }
  for(j in 1:ze_random) {mu_g_e_hat[j, time] ~ dnorm(0, 0.001)
  s_g_e_hat[j, time] ~ dunif(0, 100) }
  }#end loop through time standard deviation priors
  
  # correlation random effects
  for(time in 1:max_time) {# loop through time correlation priors random effects
  for(s in 1:n_clus_c) {b_f[s, time] ~ dnorm(mu_b_f_hat[time], tau_b_f_hat[time]) }
  mu_b_f_hat[time] ~ dnorm(0, 0.001)
  tau_b_f_hat[time] <- 1 / ss_b_f_hat[time]
  ss_b_f_hat[time] <- s_b_f_hat[time] * s_b_f_hat[time]
  s_b_f_hat[time] ~ dunif(0, 100)
  
  # correlation time random effects
  for(s in 1:n_clus_c) {b_te[s, time] ~ dnorm(mu_b_te_hat[time], tau_b_te_hat[time]) 
  b_tc[s, time] ~ dnorm(mu_b_tc_hat[time], tau_b_tc_hat[time])}
  for(s in 1:n_clus_e) {a_te[s, time] ~ dnorm(mu_a_te_hat[time], tau_a_te_hat[time]) 
  a_tc[s, time] ~ dnorm(mu_a_tc_hat[time], tau_a_tc_hat[time])}

  # correlation time priors random effects
  mu_b_te_hat[time] ~ dnorm(0, 0.001)
  tau_b_te_hat[time] <- 1 / ss_b_te_hat[time]
  ss_b_te_hat[time] <- s_b_te_hat[time] * s_b_te_hat[time]
  s_b_te_hat[time] ~ dunif(0, 100) 
  mu_b_tc_hat[time] ~ dnorm(0, 0.001)
  tau_b_tc_hat[time] <- 1 / ss_b_tc_hat[time]
  ss_b_tc_hat[time] <- s_b_tc_hat[time] * s_b_tc_hat[time]
  s_b_tc_hat[time] ~ dunif(0, 100) 
  mu_a_te_hat[time] ~ dnorm(0, 0.001)
  tau_a_te_hat[time] <- 1 / ss_a_te_hat[time]
  ss_a_te_hat[time] <- s_a_te_hat[time] * s_a_te_hat[time]
  s_a_te_hat[time] ~ dunif(0, 100) 
  mu_a_tc_hat[time] ~ dnorm(0, 0.001)
  tau_a_tc_hat[time] <- 1 / ss_a_tc_hat[time]
  ss_a_tc_hat[time] <- s_a_tc_hat[time] * s_a_tc_hat[time]
  s_a_tc_hat[time] ~ dunif(0, 100)
  #end correlation time priors random effects 
  }#end loop time correlation priors random effects

  #priors on missing data mechanism
  for(time in 1:max_time) {# loop through time gamma priors effects
  for (j in 1:ze_fixed) {#begin gamma priors effects
  for(mdrop in 1:n_patt_e){# loop through mpatterne gamma priors effects
  gamma_e[j, time, mdrop] ~ dnorm(0, 0.01)
  }#end loop mpatterne gamma priors effects
  }#end gamma priors effects
  }#end loop time gamma priors effects

  for(time in 1:max_time) {# loop through time gamma priors costs
  for (j in 1:zc_fixed) {#begin gamma priors costs
  for(mdrop in 1:n_patt_c){# loop through mpatternc gamma priors costs
  gamma_c[j, time, mdrop] ~ dnorm(0, 0.01)
  }#end loop mpatterne gamma priors costs
  }#end gamma priors costs
  }#end loop time gamma priors costs

  #priors on random effects missing data mechanism
  for(time in 1:max_time) {# loop through time g_e priors effects
  for (j in 1:ze_random) {#begin g_e priors effects
  for(s in 1:n_clus_me) {g_e[j, s, time] ~ dnorm(mu_g_e_hat[j, time], tau_g_e_hat[j, time]) }
  }#end g_e priors effects
  }#end loop time g_e priors effects

  for(time in 1:max_time) {# loop through time g_c priors costs
  for (j in 1:zc_random) {#begin g_c priors costs
  for(s in 1:n_clus_mc) {g_c[j, s, time] ~ dnorm(mu_g_c_hat[j, time], tau_g_c_hat[j, time]) }
  }#end g_c priors costs
  }#end loop time g_c priors costs

  #mnar parameters
  for(time in 1:max_time) {# loop through time mnar priors
  for(mdrop in 1:n_patt_e){# loop through mpatterne mnar priors
  delta_e[time, mdrop] ~ dnorm(0, 1)
  }#end loop mpatterne mnar priors
  for(mdrop in 1:n_patt_c){# loop through mpatternc mnar priors
  delta_c[time, mdrop] ~ dnorm(0, 1)
  }#end loop mpatternc mnar priors
  }#end loop time mnar priors

  #begin mnar random effects priors
  for(time in 1:max_time) {# loop through time mnar random effects priors
  for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }
  for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }
  }#end loop time mnar random effects priors

  # mean and sd mean mnar random effects priors
  for(time in 1:max_time) {# loop through time mnar random effects priors
  mu_d_e_hat[time] ~ dnorm(0, 1)
  mu_d_c_hat[time] ~ dnorm(0, 1)
  s_d_e_hat[time] ~ dunif(0, 1)
  s_d_c_hat[time] ~ dunif(0, 1) 
  } #end mnar random effects priors

}
 "
 if(length(model_txt_info$model_e_random) == 0) {
  model_string_jags <- gsub(" + inprod(X_e_random[i, ], a[, clus_e[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X_e_random[i, ], a[, clus_e[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j, time] <- 1 / ss_a_hat[j, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_hat[j, time] <- s_a_hat[j, time] * s_a_hat[j, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s, time] ~ dnorm(mu_a_hat[j, time], tau_a_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ a_te[clus_e[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + a_tc[clus_e[i], time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_e) {a_te[s, time] ~ dnorm(mu_a_te_hat[time], tau_a_te_hat[time]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("a_tc[s, time] ~ dnorm(mu_a_tc_hat[time], tau_a_tc_hat[time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_a_te_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_a_te_hat[time] <- 1 / ss_a_te_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_te_hat[time] <- s_a_te_hat[time] * s_a_te_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_te_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_a_tc_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_a_tc_hat[time] <- 1 / ss_a_tc_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_tc_hat[time] <- s_a_tc_hat[time] * s_a_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_tc_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
 } 
 if(length(model_txt_info$model_c_random) == 0) {
  model_string_jags <- gsub("+ inprod(X_c_random[i, ], b[, clus_c[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(X_c_random[i, ], b[, clus_c[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, time] <- 1 / ss_b_hat[j, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_hat[j, time] <- s_b_hat[j, time] * s_b_hat[j, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_c_random[], mu_b_hat[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s, time] ~ dnorm(mu_b_hat[j, time], tau_b_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b_f[clus_c[i], 1] * (eff[i, 1] - tmu_e[1]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b_f[clus_c[i], time] * (eff[i, time] - tmu_e[time]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b_te[clus_c[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + b_tc[clus_c[i], time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s, time] ~ dnorm(mu_b_f_hat[time], tau_b_f_hat[time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_f_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat[time] <- 1 / ss_b_f_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat[time] <- s_b_f_hat[time] * s_b_f_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_c) {b_te[s, time] ~ dnorm(mu_b_te_hat[time], tau_b_te_hat[time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("b_tc[s, time] ~ dnorm(mu_b_tc_hat[time], tau_b_tc_hat[time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_te_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_te_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_te_hat[time] <- 1 / ss_b_te_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_te_hat[time] <- s_b_te_hat[time] * s_b_te_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_te_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_tc_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_tc_hat[time] <- 1 / ss_b_tc_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_tc_hat[time] <- s_b_tc_hat[time] * s_b_tc_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_tc_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  } else if(length(model_txt_info$model_c_random) != 0 & model_txt_info$is_c_random_c & !model_txt_info$is_int_c_random_c) {
    model_string_jags <- gsub("+ inprod(X_c_random[i, ], b[, clus_c[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(X_c_random[i, ], b[, clus_c[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, time] <- 1 / ss_b_hat[j, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j, time] <- s_b_hat[j, time] * s_b_hat[j, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_cov_c_random[], mu_b_hat[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s, time] ~ dnorm(mu_b_hat[j, time], tau_b_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_txt_info$model_me_random) == 0) {
  model_string_jags <- gsub("+ inprod(Z_e_random[i, ], g_e[, clus_me[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(Z_e_random[i, ], g_e[, clus_me[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, time] <- 1 / ss_g_e_hat[j, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j, time] <- s_g_e_hat[j, time] * s_g_e_hat[j, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_z_e_random[], mu_g_e_hat[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s, time] ~ dnorm(mu_g_e_hat[j, time], tau_g_e_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end loop time g_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  if(length(model_txt_info$model_mc_random) == 0 & !"c" %in% model_txt_info$model_mc_random) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_txt_info$model_me_random) != 0 & model_txt_info$is_me_random_e & !model_txt_info$is_int_me_random_e) {
    model_string_jags <- gsub("+ inprod(Z_e_random[i, ], g_e[, clus_me[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z_e_random[i, ], g_e[, clus_me[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, time] <- 1 / ss_g_e_hat[j, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_e_hat[j, time] <- s_g_e_hat[j, time] * s_g_e_hat[j, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_e_random[], mu_g_e_hat[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g_e priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s, time] ~ dnorm(mu_g_e_hat[j, time], tau_g_e_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_e priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g_e priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_e_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_txt_info$model_mc_random) == 0) {
    model_string_jags <- gsub("+ inprod(Z_c_random[i, ], g_c[, clus_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z_c_random[i, ], g_c[, clus_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, time] <- 1 / ss_g_c_hat[j, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j, time] <- s_g_c_hat[j, time] * s_g_c_hat[j, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c_random[], mu_g_c_hat[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s, time] ~ dnorm(mu_g_c_hat[j, time], tau_g_c_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
  if(length(model_txt_info$model_me_random) == 0 & !"e" %in% model_txt_info$model_me_random) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0 & length(model_txt_info$model_c_random) == 0 | length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_txt_info$model_mc_random) != 0 & model_txt_info$is_mc_random_c & !model_txt_info$is_int_mc_random_c) {
    model_string_jags <- gsub("+ inprod(Z_c_random[i, ], g_c[, clus_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z_c_random[i, ], g_c[, clus_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, time] <- 1 / ss_g_c_hat[j, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j, time] <- s_g_c_hat[j, time] * s_g_c_hat[j, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c_random[], mu_g_c_hat[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_e priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s, time] ~ dnorm(mu_g_c_hat[j, time], tau_g_c_hat[j, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
  }
  if(type == "MAR") {
    model_string_jags <- gsub("+ delta_e[time, mdrop] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_e[1, mdrop] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_e[time, mdrop] * eff[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_e[time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[time, mdrop] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, mdrop] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[time, mdrop] * cost[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c[time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(mdrop in 1:n_patt_e){# loop through mpatterne mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop mpatterne mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(mdrop in 1:n_patt_c){# loop through mpatternc mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop mpatternc mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mnar parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_me_random) != 0) {
    model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_mc_random) != 0) {
      model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
    }
  } else if(type == "MNAR_eff") {
    model_string_jags <- gsub("+ delta_c[time, mdrop] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, mdrop] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[time, mdrop] * cost[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(mdrop in 1:n_patt_c){# loop through mpatternc mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop mpatternc mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c[time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
     if(length(model_txt_info$model_mc_random) != 0) {
      model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!("c" %in% model_txt_info$model_mc_random) & !"e" %in% model_txt_info$model_me_random) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(length(model_txt_info$model_mc_random) == 0 & !"e" %in% model_txt_info$model_me_random) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    }
  } else if(type == "MNAR_cost") {
    model_string_jags <- gsub("+ delta_e[time, mdrop] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_e[1, mdrop] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_e[time, mdrop] * eff[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(mdrop in 1:n_patt_e){# loop through mpatterne mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop mpatterne mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_e[time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_me_random) != 0) {
      model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!("e" %in% model_txt_info$model_me_random) & !"c" %in% model_txt_info$model_mc_random) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random)) {
      model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(length(model_txt_info$model_me_random) == 0 & !"c" %in% model_txt_info$model_mc_random) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR") {
    if(length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("+ d_e[clus_me[i], 1] * eff[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_e[clus_me[i], time] * eff[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat[time] * mean(eff[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {d_e[s, time] ~ dnorm(mu_d_e_hat[time], tau_d_e_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_e_hat[time] <- 1 / ss_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_e_hat[time] <- s_d_e_hat[time] * s_d_e_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_e_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_e_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!"c" %in% model_txt_info$model_mc_random) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random)) {
      model_string_jags <- gsub("+ d_c[clus_mc[i], 1] * cost[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d_c[clus_mc[i], time] * cost[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[time] * mean(cost[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {d_c[s, time] ~ dnorm(mu_d_c_hat[time], tau_d_c_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[time] <- 1 / ss_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[time] <- s_d_c_hat[time] * s_d_c_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!"e" %in% model_txt_info$model_me_random) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_txt_info$model_mc_random) != 0 & !("c" %in% model_txt_info$model_mc_random) & length(model_txt_info$model_me_random) != 0 & !("e" %in% model_txt_info$model_me_random)) {
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  }
  if(model_txt_info$ind_random | model_txt_info$time_dep %in% c("none", "biv")) {
    if(model_txt_info$ind_random) {
      model_string_jags <- gsub("+ b_f[clus_c[i], 1] * (eff[i, 1] - tmu_e[1]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ b_f[clus_c[i], time] * (eff[i, time] - tmu_e[time]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_c) {b_f[s, time] ~ dnorm(mu_b_f_hat[time], tau_b_f_hat[time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_f_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_f_hat[time] <- 1 / ss_b_f_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_f_hat[time] <- s_b_f_hat[time] * s_b_f_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_f_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time correlation priors random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time correlation priors random effects", "", model_string_jags, fixed = TRUE)
    }
      model_string_jags <- gsub("+ a_te[clus_e[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + a_tc[clus_e[i], time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_e) {a_te[s, time] ~ dnorm(mu_a_te_hat[time], tau_a_te_hat[time]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("a_tc[s, time] ~ dnorm(mu_a_tc_hat[time], tau_a_tc_hat[time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_a_te_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_a_te_hat[time] <- 1 / ss_a_te_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_a_te_hat[time] <- s_a_te_hat[time] * s_a_te_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_a_te_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_a_tc_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_a_tc_hat[time] <- 1 / ss_a_tc_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_a_tc_hat[time] <- s_a_tc_hat[time] * s_a_tc_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_a_tc_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)  
      model_string_jags <- gsub("+ b_te[clus_c[i], time] * (eff[i, time - 1] - tmu_e[time - 1]) + b_tc[clus_c[i], time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_c) {b_te[s, time] ~ dnorm(mu_b_te_hat[time], tau_b_te_hat[time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("b_tc[s, time] ~ dnorm(mu_b_tc_hat[time], tau_b_tc_hat[time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# correlation time priors random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_te_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_te_hat[time] <- 1 / ss_b_te_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_te_hat[time] <- s_b_te_hat[time] * s_b_te_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_te_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_tc_hat[time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_tc_hat[time] <- 1 / ss_b_tc_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_tc_hat[time] <- s_b_tc_hat[time] * s_b_tc_hat[time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_tc_hat[time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop correlation time priors random effects", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time a ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end a loop time", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end b loop time", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) == 1) {
      if(model_txt_info$model_c_random == "e") {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b ", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end b loop time", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_txt_info$model_e_random) == 0 & length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_me_random) == 0 & length(model_txt_info$model_mc_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time transformation of random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time transformation of random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# correlation time random effects", "", model_string_jags, fixed = TRUE)
    }
  }
  if(!model_txt_info$ind_random) {
   if(length(model_txt_info$model_c_random) == 1) {
     if(model_txt_info$model_c_random == "e") {
       model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b ", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end b loop time", "", model_string_jags, fixed = TRUE)
     }
   }
  }
  if(model_txt_info$ind_fixed) {
    model_string_jags <- gsub("+ beta_f[1] * (eff[i, 1] - tmu_e[1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ beta_f[time] * (eff[i, time] - tmu_e[time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
 if(model_txt_info$time_dep %in% c("none", "biv")) {
   model_string_jags <- gsub("+ beta_te[time] * (eff[i, time - 1] - tmu_e[time - 1]) + beta_tc[time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("+ alpha_te[time] * (eff[i, time - 1] - tmu_e[time - 1]) + alpha_tc[time] * (cost[i, time - 1] - tmu_c[time - 1])", "", model_string_jags, fixed = TRUE)  
   model_string_jags <- gsub("beta_te[time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta_tc[time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha_te[time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha_tc[time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#time dependence", "", model_string_jags, fixed = TRUE)
 }  
  if(dist_c == "norm") {
  model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
  model_string_jags <- gsub("cost[i, time] ~ dnorm(cmu_c[i, time], tau_c[time])", "cost[i, time] ~ dgamma(cmu_c[i, time] * ctau_c[i, time], ctau_c[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs", "for(time in 1:max_time) { ctau_c[i, time] <- cmu_c[i, time] / pow(s_c[time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i, 1] <- ", "log(cmu_c[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i, time] <- ", "log(cmu_c[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[time] <- 1 / ss_c[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[time] <- s_c[time] * s_c[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c[time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i, time] <- logdensity.norm(cost[i, time], cmu_c[i, time], tau_c[time])", "loglik_c[i, time] <- logdensity.gamma(cost[i, time], cmu_c[i, time] * ctau_c[i, time], ctau_c[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
      model_string_jags <- gsub("tmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time])", "tmu_c[time] <- exp(inprod(mean_cov_c_fixed[], beta[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("tmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time]) + inprod(mean_cov_c_random[], mu_b_hat[, time])", "tmu_c[time] <- exp(inprod(mean_cov_c_fixed[], beta[, time]) + inprod(mean_cov_c_random[], mu_b_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_c == "lnorm") {
  model_string_jags <- gsub("cost[i, time] ~ dnorm(cmu_c[i, time], tau_c[time])", "cost[i, time] ~ dlnorm(clmu_c[i, time], ltau_c[time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i, 1] <- ", "clmu_c[i, 1] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_c[i, time] <- ", "clmu_c[i, time] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[time] <- 1 / ss_c[time]", "ltau_c[time] <- 1 / lss_c[time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[time] <- s_c[time] * s_c[time]", "lss_c[time] <- ls_c[time] * ls_c[time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#std for lnorm", "s_c[time] <- sqrt(exp(2 * tlmu_c[time] + lss_c[time]) * (exp(lss_c[time]) - 1))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c[time] ~ dunif(0, 10)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "tmu_c[time] <- exp(tlmu_c[time] + lss_c[time] / 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c[i, time] <- logdensity.norm(cost[i, time], cmu_c[i, time], tau_c[time])", "loglik_c[i, time] <- logdensity.lnorm(cost[i, time], clmu_c[i, time], ltau_c[time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) == 0) {
      model_string_jags <- gsub("tmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time])", "tlmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time])", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("tmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time]) + inprod(mean_cov_c_random[], mu_b_hat[, time])", "tlmu_c[time] <- inprod(mean_cov_c_fixed[], beta[, time]) + inprod(mean_cov_c_random[], mu_b_hat[, time])", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e == "norm") {
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "beta") {
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dbeta(cmu_e[i, time] * ctau_e[i, time], (1 - cmu_e[i, time]) * ctau_e[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "for(time in 1:max_time) { ctau_e[i, time] <- (cmu_e[i, time] * (1 - cmu_e[i, time]) / pow(s_e[time], 2) - 1) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "logit(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "logit(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] ~ dunif(0, sqrt(tmu_e[time] * (1 - tmu_e[time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.beta(eff[i, time], cmu_e[i, time] * ctau_e[i, time], (1 - cmu_e[i, time]) * ctau_e[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- ilogit(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- ilogit(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])  + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "gamma"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dgamma(cmu_e[i, time] * ctau_e[i, time], ctau_e[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "for(time in 1:max_time) { ctau_e[i, time] <- cmu_e[i, time] / pow(s_e[time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "log(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "log(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.gamma(eff[i, time], cmu_e[i, time] * ctau_e[i, time], ctau_e[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "exp"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dexp(1 / cmu_e[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "log(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "log(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] <- tmu_e[time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.exp(eff[i, time], 1/cmu_e[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "weib"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dweib(ctau_e[i, time], cmu_e[i, time] / exp(loggam(1 + 1/ctau_e[i, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "for(time in 1:max_time) { ctau_e[i, time] <- pow(s_e[time] / cmu_e[i, time], - 1.086) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "log(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "log(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.weib(eff[i, time], ctau_e[i, time], cmu_e[i, time] / exp(loggam(1 + 1/ctau_e[i, time])))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "logis"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dlogis(cmu_e[i, time], tau_e[time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "tau_e[time] <- 1 / sqrt((3 * ss_e[time]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.logis(eff[i, time], cmu_e[i, time], tau_e[time])", model_string_jags, fixed = TRUE)
  } else if(dist_e == "bern"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dbern(cmu_e[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "logit(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "logit(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] <- sqrt(tmu_e[time] * (1 - tmu_e[time]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.bern(eff[i, time], cmu_e[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- ilogit(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- ilogit(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "pois"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dpois(cmu_e[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "log(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "log(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e[time] <- sqrt(tmu_e[time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.pois(eff[i, time], cmu_e[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "negbin"){
  model_string_jags <- gsub("eff[i, time] ~ dnorm(cmu_e[i, time], tau_e[time])", "eff[i, time] ~ dnegbin(tau_e[time] / (tau_e[time] + cmu_e[i, time]), tau_e[time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, 1] <- ", "log(cmu_e[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cmu_e[i, time] <- ", "log(cmu_e[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[time] <- 1 / ss_e[time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[time] <- s_e[time] * s_e[time]", "s_e[time] <- sqrt((tau_e[time] / (tau_e[time] + tmu_e[time])) * tau_e[time]) / (1 - (tau_e[time] / (tau_e[time] + tmu_e[time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[time] ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e[i, time] <- logdensity.norm(eff[i, time], cmu_e[i, time], tau_e[time])", "loglik_e[i, time] <- logdensity.negbin(eff[i, time], tau_e[time] / (tau_e[time] + cmu_e[i, time]), tau_e[time])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) == 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_txt_info$model_e_random) != 0) {
    model_string_jags <- gsub("tmu_e[time] <- inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time])", "tmu_e[time] <- exp(inprod(mean_cov_e_fixed[], alpha[, time]) + inprod(mean_cov_e_random[], mu_a_hat[, time]))", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e %in% c("beta", "gamma", "weib", "bern", "pois", "exp")) {
     if(dist_c == "gamma"){
     model_string_jags <- gsub("#begin transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for(time in 1:max_time) {# loop through time transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end loop time transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("# end transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
   }
   if(length(model_txt_info$model_e_random) != 0 & model_txt_info$pe_random == 1) {
       model_string_jags <- gsub("inprod(X_e_random[i, ], a[, clus_e[i], 1])", "X_e_random[i] * a[clus_e[i], 1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X_e_random[i, ], a[, clus_e[i], time])", "X_e_random[i] * a[clus_e[i], time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e_random[], mu_a_hat[, time])", "mean_cov_e_random * mu_a_hat[time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pe_random) {#begin a priors effects", "#begin a priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n_clus_e) {a[j, s, time] ~ dnorm(mu_a_hat[j, time], tau_a_hat[j, time]) }", 
                                 "for(s in 1:n_clus_e) {a[s, time] ~ dnorm(mu_a_hat[time], tau_a_hat[time]) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a priors effects", "#end a priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j, time] ~ dnorm(0, 0.001)", "mu_a_hat[time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_a_hat[j, time] ~ dunif(0, 100) }", "s_a_hat[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j, time] <- 1 / ss_a_hat[j, time]", "tau_a_hat[time] <- 1 / ss_a_hat[time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_a_hat[j, time] <- s_a_hat[j, time] * s_a_hat[j, time] }", "ss_a_hat[time] <- s_a_hat[time] * s_a_hat[time]", model_string_jags, fixed = TRUE)
   }
   if(model_txt_info$pc_fixed == 0 & !model_txt_info$ind_fixed) {
   prior_beta <- "#"
   end_prior_beta <- "#end beta priors costs"
   model_string_jags <- gsub("inprod(X_c_fixed[i, ], beta[, 1]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X_c_fixed[i, ], beta[, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c_fixed[], beta[, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c_fixed[], beta[, time])", "beta_f[time] * (mean(eff[, time]) - tmu_e[time])", model_string_jags, fixed = TRUE)
   if(model_txt_info$pc_random == 0 & !model_txt_info$ind_random) {
   model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b ", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b loop time", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c_fixed[], beta[, time])", "beta_f[time] * (mean(eff[, time]) - tmu_e[time])", model_string_jags, fixed = TRUE)
   }
   model_string_jags <- gsub("for (j in 1:pc_fixed) {#begin beta priors costs", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[j, time] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta priors costs", end_prior_beta,model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(time in 1:max_time) {# loop through time beta", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#end beta priors costs", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta loop time", "", model_string_jags, fixed = TRUE)
   }
   if(length(model_txt_info$model_c_random) != 0 & model_txt_info$pc_random == 1) {
   model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i], 1])", "X_c_random[i] * b[clus_c[i], 1]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i], time])", "X_c_random[i] * b[clus_c[i], time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c_random[], mu_b_hat[, time])", "mean_cov_c_random * mu_b_hat[time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "#begin b priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s, time] ~ dnorm(mu_b_hat[j, time], tau_b_hat[j, time]) }", "for(s in 1:n_clus_c) {b[s, time] ~ dnorm(mu_b_hat[time], tau_b_hat[time]) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b priors costs", "#end b priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, time] ~ dnorm(0, 0.001)", "mu_b_hat[time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("s_b_hat[j, time] ~ dunif(0, 100) }", "s_b_hat[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, time] <- 1 / ss_b_hat[j, time]", "tau_b_hat[time] <- 1 / ss_b_hat[time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("ss_b_hat[j, time] <- s_b_hat[j, time] * s_b_hat[j, time] }", "ss_b_hat[time] <- s_b_hat[time] * s_b_hat[time]", model_string_jags, fixed = TRUE)
   }
   if(model_txt_info$ze_fixed == 1) {
   inprod_e_base <- "Z_e_fixed[i] * gamma_e[1, mdrop]"
   inprod_e <- "Z_e_fixed[i] * gamma_e[time, mdrop]"
     if(type %in% c("MAR", "MNAR_cost")) {
     inprod_mean_e <- "mean_z_e_fixed * gamma_e[time, mdrop]"
     model_string_jags <- gsub("inprod(mean_z_e_fixed[], gamma_e[, time, mdrop])", inprod_mean_e, model_string_jags, fixed = TRUE)
     }
     if(type %in% c("MNAR_eff", "MNAR")) {
     inprod_mean_e <- "mean_z_e_fixed * gamma_e[time, mdrop] + delta_e[time, mdrop] * mean(eff[, time])"
     model_string_jags <- gsub("inprod(mean_z_e_fixed[], gamma_e[, time, mdrop]) + delta_e[time, mdrop] * mean(eff[, time])", inprod_mean_e, model_string_jags, fixed = TRUE)
     }
   begin_prior_gamma <- "#begin gamma priors effects"
   prior_gamma_e <- "gamma_e[time, mdrop] ~ dnorm(0, 0.01)"
   end_prior_gamma <- "#end gamma priors effects"
   model_string_jags <- gsub("inprod(Z_e_fixed[i, ], gamma_e[, 1, mdrop])", inprod_e_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(Z_e_fixed[i, ], gamma_e[, time, mdrop])", inprod_e, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:ze_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_e[j, time, mdrop] ~ dnorm(0, 0.01)", prior_gamma_e, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_e[1, time, mdrop] ~ dnorm(0, 0.01)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
   #   model_string_jags <- gsub("for(mdrop in 1:n_patt_e){# loop through mpatterne gamma priors effects", "", model_string_jags, fixed = TRUE)
   #   model_string_jags <- gsub("}#end loop mpatterne gamma priors effects", "", model_string_jags, fixed = TRUE)
   } 
   if(model_txt_info$ze_random == 1 & length(model_txt_info$model_me_random) != 0) {
      if(type %in% c("MAR", "MNAR_cost")) {
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[, time])", "mean_z_e_random * mu_g_e_hat[time]", model_string_jags, fixed = TRUE)
      }
      if(type %in% c("MNAR_eff", "MNAR")) {
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[, time]) + mu_d_e_hat[time] * mean(eff[, time])", "mean_z_e_random * mu_g_e_hat[time] + mu_d_e_hat[time] * mean(eff[, time])", model_string_jags, fixed = TRUE)
      }
      model_string_jags <- gsub("inprod(Z_e_random[i, ], g_e[, clus_me[i], 1])", "Z_e_random[i] * g_e[clus_me[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z_e_random[i, ], g_e[, clus_me[i], time])", "Z_e_random[i] * g_e[clus_me[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, time] <- 1 / ss_g_e_hat[j, time]", "tau_g_e_hat[time] <- 1 / ss_g_e_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_e_hat[j, time] <- s_g_e_hat[j, time] * s_g_e_hat[j, time] }", "ss_g_e_hat[time] <- s_g_e_hat[time] * s_g_e_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e_random[], mu_g_e_hat[, time])", "mean_z_e_random * mu_g_e_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:ze_random) {#begin g_e priors effects", "#begin g_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_me) {g_e[j, s, time] ~ dnorm(mu_g_e_hat[j, time], tau_g_e_hat[j, time]) }", "for(s in 1:n_clus_me) {g_e[s, time] ~ dnorm(mu_g_e_hat[time], tau_g_e_hat[time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g_e priors effects", "#end g_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, time] ~ dnorm(0, 0.001)", "mu_g_e_hat[time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_e_hat[j, time] ~ dunif(0, 100) }", "s_g_e_hat[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
   if(model_txt_info$zc_fixed == 1) {
      inprod_c_base <- "Z_c_fixed[i] * gamma_c[1, mdrop]"
      inprod_c <- "Z_c_fixed[i] * gamma_c[time, mdrop]"
          if(type %in% c("MAR", "MNAR_eff")) {
          inprod_mean_c <- "mean_z_c_fixed * gamma_c[time, mdrop]"
          model_string_jags <- gsub("inprod(mean_z_c_fixed[], gamma_c[, time, mdrop])", inprod_mean_c, model_string_jags, fixed = TRUE)
          }
          if(type %in% c("MNAR_cost", "MNAR")) {
          inprod_mean_c <- "mean_z_c_fixed * gamma_c[time, mdrop] + delta_c[time, mdrop] * mean(cost[, time])"
          model_string_jags <- gsub("inprod(mean_z_c_fixed[], gamma_c[, time, mdrop]) + delta_c[time, mdrop] * mean(cost[, time])", inprod_mean_c, model_string_jags, fixed = TRUE)
          }
      begin_prior_gamma <- "#begin gamma priors costs"
      prior_gamma_c <- "gamma_c[time, mdrop] ~ dnorm(0, 0.01)"
      end_prior_gamma <- "#end gamma priors costs"
      model_string_jags <- gsub("inprod(Z_c_fixed[i, ], gamma_c[, 1, mdrop])", inprod_c_base, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z_c_fixed[i, ], gamma_c[, time, mdrop])", inprod_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[j, time, mdrop] ~ dnorm(0, 0.01)", prior_gamma_c, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[1, time, mdrop] ~ dnorm(0, 0.01)", "#", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
      #      model_string_jags <- gsub("for(mdrop in 1:n_patt_c){# loop through mpatternc gamma priors costs", "", model_string_jags, fixed = TRUE)
      #      model_string_jags <- gsub("}#end loop mpatternc gamma priors costs", "", model_string_jags, fixed = TRUE)
      }
      if(model_txt_info$zc_random == 1 & length(model_txt_info$model_mc_random) != 0) {
         if(type %in% c("MAR", "MNAR_eff")) {
         model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[, time])", "mean_z_c_random * mu_g_c_hat[time]", model_string_jags, fixed = TRUE)
         }
         if(type %in% c("MNAR_cost", "MNAR")) {
         model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[, time]) + mu_d_c_hat[time] * mean(cost[, time])", "mean_z_c_random * mu_g_c_hat[time] + mu_d_c_hat[time] * mean(cost[, time])", model_string_jags, fixed = TRUE)
         }
      model_string_jags <- gsub("inprod(Z_c_random[i, ], g_c[, clus_mc[i], 1])", "Z_c_random[i] * g_c[clus_mc[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z_c_random[i, ], g_c[, clus_mc[i], time])", "Z_c_random[i] * g_c[clus_mc[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, time] <- 1 / ss_g_c_hat[j, time]", "tau_g_c_hat[time] <- 1 / ss_g_c_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_c_hat[j, time] <- s_g_c_hat[j, time] * s_g_c_hat[j, time] }", "ss_g_c_hat[time] <- s_g_c_hat[time] * s_g_c_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c_random[], mu_g_c_hat[, time])", "mean_z_c_random * mu_g_c_hat[time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g_c priors costs", "#begin g_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n_clus_mc) {g_c[j, s, time] ~ dnorm(mu_g_c_hat[j, time], tau_g_c_hat[j, time]) }", "for(s in 1:n_clus_mc) {g_c[s, time] ~ dnorm(mu_g_c_hat[time], tau_g_c_hat[time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g_c priors costs", "#end g_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, time] ~ dnorm(0, 0.001)", "mu_g_c_hat[time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_c_hat[j, time] ~ dunif(0, 100) }", "s_g_c_hat[time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
      model_string_jags <- prior_lmdm(type = type, dist_e = dist_e, dist_c = dist_c, 
                                     model_txt_info = model_txt_info, model_string_jags = model_string_jags)
   writeLines(model_string_jags, "lmdm.txt")
   model_string <- "lmdm.txt"
   return(model_string)
}