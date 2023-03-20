#'An internal function to select which type of selection model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Selection models
#' @param dist_u Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param ind_fixed Logical; if TRUE independence between effectiveness and costs at the same time is assumed, else correlation is accounted for
#' @param ind_time_fixed Logical; if TRUE independence between effectiveness and costs over time is assumed, else an AR1 correlation structure is accounted for
#' @param ind_random Logical; if TRUE independence at the level of the random effects between effectiveness and costs is assumed, else correlation is accounted for
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param pu_fixed Number of fixed effects for the effectiveness model
#' @param pc_fixed Number of fixed effects for the cost model
#' @param zu_fixed Number of fixed effects or the missingness indicators model for the effectiveness
#' @param zc_fixed Number of fixed effects or the missingness indicators model for the costs
#' @param pu_random Number of random effects for the effectiveness model
#' @param pc_random Number of random effects for the cost model
#' @param zu_random Number of random effects or the missingness indicators model for the effectiveness
#' @param zc_random Number of random effects or the missingness indicators model for the costs
#' @param model_u_random Random effects formula for the effectiveness model
#' @param model_c_random Random effects formula for the costs model
#' @param model_mu_random Random effects formula for the missingness indicators model for the effectiveness
#' @param model_mc_random Random effects formula for the missingness indicators model for the costs
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_selection_long <- function(dist_u , dist_c, type, pu_fixed, pc_fixed, zu_fixed, zc_fixed, ind_fixed, ind_time_fixed, pu_random, pc_random, zu_random, zc_random, ind_random, 
                            model_u_random, model_c_random, model_mu_random, model_mc_random) eval.parent(substitute( {
  model_string_jags<-  "
  model {

  #control
  for(i in 1:N1) {
  for(time in 1:max_time) {# loop through time
  #costs and effects model by time point
  cost1[i, time] ~ dnorm(mu_c1[i, time], tau_c[1, time])
  eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])
  }#end loop time

  #derive mean and std effects1 
  #derive mean and std costs1

  #mean regression at baseline
  mu_c1[i, 1] <- inprod(X1_c_fixed[i, ], beta[, 1, 1]) + beta_f[1, 1] * (eff1[i, 1] - mu_u[1, 1]) + inprod(X1_c_random[i, ], b1[, clus1_c[i], 1]) + b1_f[clus1_c[i], 1] * (eff1[i, 1] - mu_u[1, 1]) 
  mu_u1[i, 1] <- inprod(X1_u_fixed[i, ], alpha[, 1, 1]) + inprod(X1_u_random[i, ], a1[, clus1_u[i], 1])   

  #mean regression at followup
  for(time in 2:max_time) {# loop through time
  mu_c1[i, time] <- inprod(X1_c_fixed[i, ], beta[, 1, time]) + beta_f[1, time] * (eff1[i, time] - mu_u[1, time]) + inprod(X1_c_random[i, ], b1[, clus1_c[i], time]) + b1_f[clus1_c[i], time] * (eff1[i, time] - mu_u[1, time]) + beta_tu[1, time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + beta_tc[1, time] * (cost1[i, time - 1] - mu_c[1, time - 1]) + b1_tu[clus1_c[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + b1_tc[clus1_c[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1]) 
  mu_u1[i, time] <- inprod(X1_u_fixed[i, ], alpha[, 1, time]) + inprod(X1_u_random[i, ], a1[, clus1_u[i], time]) + alpha_tu[1, time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + alpha_tc[1, time] * (cost1[i, time - 1] - mu_c[1, time - 1]) + a1_tu[clus1_u[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + a1_tc[clus1_u[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1])
  }#end loop time


  #missing data mechanism
  for(time in 1:max_time) {# loop through time
  m_eff1[i, time] ~ dcat(pq_1[i, time, 1:3])
  m_cost1[i, time] ~ dcat(pc_1[i, time, 1:3])
  for(mdrop in 1:3){# loop through dropout
  pq_1[i, time, mdrop] <- phiq_1[i, time, mdrop]/sum(phiq_1[i, time, ])
  pc_1[i, time, mdrop] <- phic_1[i, time, mdrop]/sum(phic_1[i, time, ])
  }#end loop dropout
  }#end loop time
  
  #log missing regression at baseline
  for(mdrop in 1:3){# loop through dropout
  log(phiq_1[i, 1, mdrop]) <- inprod(Z1_u_fixed[i, ], gamma_u[, 1, 1, mdrop]) + delta_u[1, 1, mdrop] * eff1[i, 1] + inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], 1]) + d1_u[clus1_mu[i], 1] * eff1[i, 1]
  log(phic_1[i, 1, mdrop]) <- inprod(Z1_c_fixed[i, ], gamma_c[, 1, 1, mdrop]) + delta_c[1, 1, mdrop] * cost1[i, 1] + inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], 1]) + d1_c[clus1_mc[i], 1] * cost1[i, 1]
  }#end loop dropout

  #log missing regression at followup
  for(time in 2:max_time) {# loop through time
  for(mdrop in 1:3){# loop through dropout
  log(phiq_1[i, time, mdrop]) <- inprod(Z1_u_fixed[i, ], gamma_u[, 1, time, mdrop]) + delta_u[1, time, mdrop] * eff1[i, time] + inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], time]) + d1_u[clus1_mu[i], time] * eff1[i, time]
  
  log(phic_1[i, time, mdrop]) <- inprod(Z1_c_fixed[i, ], gamma_c[, 1, time, mdrop]) + delta_c[1, time, mdrop] * cost1[i, time] + inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], time]) + d1_c[clus1_mc[i], time] * cost1[i, time]
  }#end loop dropout
  }#end loop time
  
  #loglikelihood
  for(time in 1:max_time) {# loop through time
  loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])
  loglik_c1[i, time] <- logdensity.norm(cost1[i, time], mu_c1[i, time], tau_c[1, time])
  loglik_mu1[i, time] <- logdensity.cat(m_eff1[i, time], pq_1[i, time, ])
  loglik_mc1[i, time] <- logdensity.cat(m_cost1[i, time], pc_1[i, time, ])
  }#end loop time
  }
  

  #intervention
  for(i in 1:N2) {
  for(time in 1:max_time) {# loop through time
  #costs and effects model by time point
  cost2[i, time] ~ dnorm(mu_c2[i, time], tau_c[2, time])
  eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])
  }#end loop time

  #derive mean and std effects2 
  #derive mean and std costs2

  #mean regression at baseline
  mu_c2[i, 1] <- inprod(X2_c_fixed[i, ], beta[, 2, 1]) + beta_f[2, 1] * (eff2[i, 1] - mu_u[2, 1]) + inprod(X2_c_random[i, ], b2[, clus2_c[i], 1]) + b2_f[clus2_c[i], 1] * (eff2[i, 1] - mu_u[2, 1]) 
  mu_u2[i, 1] <- inprod(X2_u_fixed[i, ], alpha[, 2, 1]) + inprod(X2_u_random[i, ], a2[, clus2_u[i], 1])   

  #mean regression at followup
  for(time in 2:max_time) {# loop through time
  mu_c2[i, time] <- inprod(X2_c_fixed[i, ], beta[, 2, time]) + beta_f[2, time] * (eff2[i, time] - mu_u[2, time]) + inprod(X2_c_random[i, ], b2[, clus2_c[i], time]) + b2_f[clus2_c[i], time] * (eff2[i, time] - mu_u[2, time]) + beta_tu[2, time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + beta_tc[2, time] * (cost2[i, time - 1] - mu_c[2, time - 1]) + b2_tu[clus2_c[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + b2_tc[clus2_c[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1]) 
  mu_u2[i, time] <- inprod(X2_u_fixed[i, ], alpha[, 2, time]) + inprod(X2_u_random[i, ], a2[, clus2_u[i], time]) + alpha_tu[2, time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + alpha_tc[2, time] * (cost2[i, time - 1] - mu_c[2, time - 1]) + a2_tu[clus2_u[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + a2_tc[clus2_u[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1])
  }#end loop time


  #missing data mechanism
  for(time in 1:max_time) {# loop through time
  m_eff2[i, time] ~ dcat(pq_2[i, time, 1:3])
  m_cost2[i, time] ~ dcat(pc_2[i, time, 1:3])
  for(mdrop in 1:3){# loop through dropout
  pq_2[i, time, mdrop] <- phiq_2[i, time, mdrop]/sum(phiq_2[i, time, ])
  pc_2[i, time, mdrop] <- phic_2[i, time, mdrop]/sum(phic_2[i, time, ])
  }#end loop dropout
  }#end loop time
  
  #log missing regression at baseline
  for(mdrop in 1:3){# loop through dropout
  log(phiq_2[i, 1, mdrop]) <- inprod(Z2_u_fixed[i, ], gamma_u[, 2, 1, mdrop]) + delta_u[2, 1, mdrop] * eff2[i, 1] + inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], 1]) + d2_u[clus2_mu[i], 1] * eff2[i, 1]
  log(phic_2[i, 1, mdrop]) <- inprod(Z2_c_fixed[i, ], gamma_c[, 2, 1, mdrop]) + delta_c[2, 1, mdrop] * cost2[i, 1] + inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], 1]) + d2_c[clus2_mc[i], 1] * cost2[i, 1]
  }#end loop dropout

  #log missing regression at followup
  for(time in 2:max_time) {# loop through time
  for(mdrop in 1:3){# loop through dropout
  log(phiq_2[i, time, mdrop]) <- inprod(Z2_u_fixed[i, ], gamma_u[, 2, time, mdrop]) + delta_u[2, time, mdrop] * eff2[i, time] + inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], time]) + d2_u[clus2_mu[i], time] * eff2[i, time]
  log(phic_2[i, time, mdrop]) <- inprod(Z2_c_fixed[i, ], gamma_c[, 2, time, mdrop]) + delta_c[2, time, mdrop] * cost2[i, time] + inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], time]) + d2_c[clus2_mc[i], time] * cost2[i, time]
  }#end loop dropout
  }#end loop time
  
  #loglikelihood
  for(time in 1:max_time) {# loop through time
  loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])
  loglik_c2[i, time] <- logdensity.norm(cost2[i, time], mu_c2[i, time], tau_c[2, time])
  loglik_mu2[i, time] <- logdensity.cat(m_eff2[i, time], pq_2[i, time, ])
  loglik_mc2[i, time] <- logdensity.cat(m_cost2[i, time], pc_2[i, time, ])
  }#end loop time
  }
  
  #transformation of parameters
  for (t in 1:2) {#begin transformation
  for(time in 1:max_time) {# loop through time transformation of params
  tau_c[t, time] <- 1 / ss_c[t, time]
  ss_c[t, time] <- s_c[t, time] * s_c[t, time]
  s_c[t, time] <- exp(ls_c[t, time])
  #mean for lnorm
  tau_u[t, time] <- 1 / ss_u[t, time]
  ss_u[t, time] <- s_u[t, time] * s_u[t, time]
  s_u[t, time] <- exp(ls_u[t, time])
  }#end loop time transformation of params
  }# end transformation of params
  
  #transformation of random effects parameters
  for (t in 1:2) {# begin transformation random effects
  for(time in 1:max_time) {# loop through time transformation of random effects
  for(j in 1:pc_random) {tau_b_hat[j, t, time] <- 1 / ss_b_hat[j, t, time]
  ss_b_hat[j, t, time] <- s_b_hat[j, t, time] * s_b_hat[j, t, time] }
  for(j in 1:pu_random) {tau_a_hat[j, t, time] <- 1 / ss_a_hat[j, t, time]
  ss_a_hat[j, t, time] <- s_a_hat[j, t, time] * s_a_hat[j, t, time] }
  for(j in 1:zc_random) {tau_g_c_hat[j, t, time] <- 1 / ss_g_c_hat[j, t, time]
  ss_g_c_hat[j, t, time] <- s_g_c_hat[j, t, time] * s_g_c_hat[j, t, time] }
  for(j in 1:zu_random) {tau_g_u_hat[j, t, time] <- 1 / ss_g_u_hat[j, t, time]
  ss_g_u_hat[j, t, time] <- s_g_u_hat[j, t, time] * s_g_u_hat[j, t, time] }
  tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]
  ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]
  tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]
  ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]
  }#end loop time transformation of random effects
  }#end transformation of random effects 
  
  #missingness probability
  for(time in 1:max_time) {# loop through time
  for(mdrop in 1:3){# loop through mpattern
  p_c[1, time, mdrop] <- exp(inprod(mean_z_c1_fixed[], gamma_c[, 1, time, mdrop]) + delta_c[1, time, mdrop] * mean(cost1[, time]) + inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time]) + mu_d_c_hat[1, time] * mean(cost1[, time]))
  p_c[2, time, mdrop] <- exp(inprod(mean_z_c2_fixed[], gamma_c[, 2, time, mdrop]) + delta_c[2, time, mdrop] * mean(cost2[, time]) + inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time]) + mu_d_c_hat[2, time] * mean(cost2[, time]))
  p_u[1, time, mdrop] <- exp(inprod(mean_z_u1_fixed[], gamma_u[, 1, time, mdrop]) + delta_u[1, time, mdrop] * mean(eff1[, time]) + inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time]) + mu_d_u_hat[1, time] * mean(eff1[, time]))
  p_u[2, time, mdrop] <- exp(inprod(mean_z_u2_fixed[], gamma_u[, 2, time, mdrop]) + delta_u[2, time, mdrop] * mean(eff2[, time]) + inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time]) + mu_d_u_hat[2, time] * mean(eff2[, time]))
  }#end loop mpattern
  }#end loop time
  
  #calculate means at mean of covariates
  for(time in 1:max_time) {# loop through time
  mu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1, time]) 
  mu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])
  mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])
  mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])
  }#end loop time

  #priors
  
  #priors for mean regression coefficients
  for(time in 1:max_time) {# loop through time alpha
  for (j in 2:pu_fixed) {#begin alpha priors effects
  for(t in 1:2) {alpha[j, t, time] ~ dnorm(0, 0.0000001) }
  }#end alpha priors effects
  alpha[1, 1, time] ~ dnorm(0, 0.0000001)
  alpha[1, 2, time] ~ dnorm(0, 0.0000001)
  }#end alpha loop time

  for(time in 1:max_time) {# loop through time beta
  for (j in 2:pc_fixed) {#begin beta priors costs
  for(t in 1:2) {beta[j, t, time] ~ dnorm(0, 0.0000001) }
  }#end beta priors costs
  beta[1, 1, time] ~ dnorm(0, 0.0000001)
  beta[1, 2, time] ~ dnorm(0, 0.0000001)
  }#end beta loop time

  #priors for mean regression random coefficients
  for(time in 1:max_time) {# loop through time a1
  for (j in 1:pu_random) {#begin a1 priors effects
  for(s in 1:n1_clus_u) {a1[j, s, time] ~ dnorm(mu_a_hat[j, 1, time], tau_a_hat[j, 1, time]) }
  }#end a1 priors effects
  }#end a1 loop time
  for(time in 1:max_time) {# loop through time a2
  for (j in 1:pu_random) {#begin a2 priors effects
  for(s in 1:n2_clus_u) {a2[j, s, time] ~ dnorm(mu_a_hat[j, 2, time], tau_a_hat[j, 2, time]) }
  }#end a2 priors effects
  }#end a2 loop time

  for(time in 1:max_time) {# loop through time b1 
  for (j in 1:pc_random) {#begin b1 priors costs
  for(s in 1:n1_clus_c) {b1[j, s, time] ~ dnorm(mu_b_hat[j, 1, time], tau_b_hat[j, 1, time]) }
  }#end b1 priors costs
  }#end b1 loop time
  for(time in 1:max_time) {# loop through time b2 
  for (j in 1:pc_random) {#begin b2 priors costs
  for(s in 1:n2_clus_c) {b2[j, s, time] ~ dnorm(mu_b_hat[j, 2, time], tau_b_hat[j, 2, time]) }
  }#end b2 priors costs
  }#end b2 loop time

  #standard deviation priors
  for(time in 1:max_time) {# loop through time standard deviation priors
  for(t in 1:2) {# loop through trt standard deviation priors
  ls_c[t, time] ~ dunif(-5, 10)
  ls_u[t, time] ~ dunif(-5, 10)

  #correlation
  beta_f[t, time] ~ dnorm(0, 0.0000001)
  
  #time dependence
  beta_tu[t, time] ~ dnorm(0, 0.0000001)
  beta_tc[t, time] ~ dnorm(0, 0.0000001)
  alpha_tu[t, time] ~ dnorm(0, 0.0000001)
  alpha_tc[t, time] ~ dnorm(0, 0.0000001)

  # mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j, t, time] ~ dnorm(0, 0.001)
  s_b_hat[j, t, time] ~ dunif(0, 100) }
  for(j in 1:pu_random) {mu_a_hat[j, t, time] ~ dnorm(0, 0.001)
  s_a_hat[j, t, time] ~ dunif(0, 100) }
  for(j in 1:zc_random) {mu_g_c_hat[j, t, time] ~ dnorm(0, 0.001)
  s_g_c_hat[j, t, time] ~ dunif(0, 100) }
  for(j in 1:zu_random) {mu_g_u_hat[j, t, time] ~ dnorm(0, 0.001)
  s_g_u_hat[j, t, time] ~ dunif(0, 100) }
  }#end loop through trt standard deviation priors
  }#end loop through time standard deviation priors
  
  # correlation random effects
  for(time in 1:max_time) {# loop through time correlation priors random effects
  for(s in 1:n1_clus_c) {b1_f[s, time] ~ dnorm(mu_b_f_hat[1, time], tau_b_f_hat[1, time]) }
  for(s in 1:n2_clus_c) {b2_f[s, time] ~ dnorm(mu_b_f_hat[2, time], tau_b_f_hat[2, time]) }
  for(t in 1:2) {mu_b_f_hat[t, time] ~ dnorm(0, 0.001)
  tau_b_f_hat[t, time] <- 1 / ss_b_f_hat[t, time]
  ss_b_f_hat[t, time] <- s_b_f_hat[t, time] * s_b_f_hat[t, time]
  s_b_f_hat[t, time] ~ dunif(0, 100) }
  
  # correlation time random effects
  for(s in 1:n1_clus_c) {b1_tu[s, time] ~ dnorm(mu_b_tu_hat[1, time], tau_b_tu_hat[1, time]) 
  b1_tc[s, time] ~ dnorm(mu_b_tc_hat[1, time], tau_b_tc_hat[1, time])}
  for(s in 1:n2_clus_c) {b2_tu[s, time] ~ dnorm(mu_b_tu_hat[2, time], tau_b_tu_hat[2, time]) 
  b2_tc[s, time] ~ dnorm(mu_b_tc_hat[2, time], tau_b_tc_hat[2, time])}
  for(s in 1:n1_clus_u) {a1_tu[s, time] ~ dnorm(mu_a_tu_hat[1, time], tau_a_tu_hat[1, time]) 
  a1_tc[s, time] ~ dnorm(mu_a_tc_hat[1, time], tau_a_tc_hat[1, time])}
  for(s in 1:n2_clus_u) {a2_tu[s, time] ~ dnorm(mu_a_tu_hat[2, time], tau_a_tu_hat[2, time]) 
  a2_tc[s, time] ~ dnorm(mu_a_tc_hat[2, time], tau_a_tc_hat[2, time])}
  
  for(t in 1:2) {# correlation time priors random effects
  mu_b_tu_hat[t, time] ~ dnorm(0, 0.001)
  tau_b_tu_hat[t, time] <- 1 / ss_b_tu_hat[t, time]
  ss_b_tu_hat[t, time] <- s_b_tu_hat[t, time] * s_b_tu_hat[t, time]
  s_b_tu_hat[t, time] ~ dunif(0, 100) 
  mu_b_tc_hat[t, time] ~ dnorm(0, 0.001)
  tau_b_tc_hat[t, time] <- 1 / ss_b_tc_hat[t, time]
  ss_b_tc_hat[t, time] <- s_b_tc_hat[t, time] * s_b_tc_hat[t, time]
  s_b_tc_hat[t, time] ~ dunif(0, 100) 
  mu_a_tu_hat[t, time] ~ dnorm(0, 0.001)
  tau_a_tu_hat[t, time] <- 1 / ss_a_tu_hat[t, time]
  ss_a_tu_hat[t, time] <- s_a_tu_hat[t, time] * s_a_tu_hat[t, time]
  s_a_tu_hat[t, time] ~ dunif(0, 100) 
  mu_a_tc_hat[t, time] ~ dnorm(0, 0.001)
  tau_a_tc_hat[t, time] <- 1 / ss_a_tc_hat[t, time]
  ss_a_tc_hat[t, time] <- s_a_tc_hat[t, time] * s_a_tc_hat[t, time]
  s_a_tc_hat[t, time] ~ dunif(0, 100)
  }#end loop correlation time priors random effects 
  }#end loop time correlation priors random effects

  #priors on missing data mechanism
  for(time in 1:max_time) {# loop through time gamma priors effects
  for (j in 2:zu_fixed) {#begin gamma priors effects
  for(mdrop in 1:3){# loop through mpattern gamma priors effects
  for(t in 1:2) {gamma_u[j, t, time, mdrop] ~ dnorm(0, 0.01) }
  }#end loop mpattern gamma priors effects
  }#end gamma priors effects
  for(mdrop in 1:3){# loop through mpattern gamma int priors effects
  gamma_u[1, 1, time, mdrop] ~ dlogis(0, 1)
  gamma_u[1, 2, time, mdrop] ~ dlogis(0, 1)
  }#end loop mpattern gamma int priors effects
  }#end loop time gamma priors effects

  for(time in 1:max_time) {# loop through time gamma priors costs
  for (j in 2:zc_fixed) {#begin gamma priors costs
  for(mdrop in 1:3){# loop through mpattern gamma priors costs
  for(t in 1:2) {gamma_c[j, t, time, mdrop] ~ dnorm(0, 0.01) }
  }#end loop mpattern gamma priors costs
  }#end gamma priors costs
  for(mdrop in 1:3){# loop through mpattern gamma int priors costs
  gamma_c[1, 1, time, mdrop] ~ dlogis(0, 1)
  gamma_c[1, 2, time, mdrop] ~ dlogis(0, 1)
  }#end loop mpattern gamma int priors costs
  }#end loop time gamma priors costs

  #priors on random effects missing data mechanism
  for(time in 1:max_time) {# loop through time g1_u priors effects
  for (j in 1:zu_random) {#begin g1_u priors effects
  for(s in 1:n1_clus_mu) {g1_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 1, time], tau_g_u_hat[j, 1, time]) }
  }#end g1_u priors effects
  }#end loop time g1_u priors effects
  for(time in 1:max_time) {# loop through time g2_u priors effects
  for (j in 1:zu_random) {#begin g2_u priors effects
  for(s in 1:n2_clus_mu) {g2_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 2, time], tau_g_u_hat[j, 2, time]) }
  }#end g2_u priors effects
  }#end loop time g2_u priors effects

  for(time in 1:max_time) {# loop through time g1_c priors costs
  for (j in 1:zc_random) {#begin g1_c priors costs
  for(s in 1:n1_clus_mc) {g1_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 1, time], tau_g_c_hat[j, 1, time]) }
  }#end g1_c priors costs
  }#end loop time g1_c priors costs
  for(time in 1:max_time) {# loop through time g2_c priors costs
  for (j in 1:zc_random) {#begin g2_c priors costs
  for(s in 1:n2_clus_mc) {g2_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 2, time], tau_g_c_hat[j, 2, time]) }
  }#end g2_c priors costs
  }#end loop time g2_c priors costs

  #mnar parameters
  for(time in 1:max_time) {# loop through time mnar priors
  for(mdrop in 1:3){# loop through mpattern mnar priors
  for(t in 1:2) {# begin mnar priors
  delta_u[t, time, mdrop] ~ dnorm(0, 1)
  delta_c[t, time, mdrop] ~ dnorm(0, 1)
  }#end mnar priors
  }#end loop mpattern mnar priors
  }#end loop time mnar priors

  #mnar random effects priors
  for(time in 1:max_time) {# loop through time mnar random effects priors
  for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }
  for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }
  for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }
  for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }
  }#end loop time mnar random effects priors

  # mean and sd mean mnar random effects priors
  for(time in 1:max_time) {# loop through time mnar random effects priors
  for(t in 1:2) {#begin mnar random effects priors
  mu_d_u_hat[t, time] ~ dnorm(0, 1)
  mu_d_c_hat[t, time] ~ dnorm(0, 1)
  s_d_u_hat[t, time] ~ dunif(0, 1)
  s_d_c_hat[t, time] ~ dunif(0, 1) 
  }#end loop time mnar random effects priors
  } #end mnar random effects priors

}
 "
 if(length(model_u_random) == 0) {
  model_string_jags <- gsub(" + inprod(X1_u_random[i, ], a1[, clus1_u[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X1_u_random[i, ], a1[, clus1_u[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X2_u_random[i, ], a2[, clus2_u[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X2_u_random[i, ], a2[, clus2_u[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pu_random) {tau_a_hat[j, t, time] <- 1 / ss_a_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_hat[j, t, time] <- s_a_hat[j, t, time] * s_a_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pu_random) {#begin a1 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_u) {a1[j, s, time] ~ dnorm(mu_a_hat[j, 1, time], tau_a_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a1 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pu_random) {#begin a2 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_u) {a2[j, s, time] ~ dnorm(mu_a_hat[j, 2, time], tau_a_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a2 priors effects", "", model_string_jags, fixed = TRUE)
  
  model_string_jags <- gsub("for(j in 1:pu_random) {mu_a_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ a1_tu[clus1_u[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + a1_tc[clus1_u[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ a2_tu[clus2_u[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + a2_tc[clus2_u[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_u) {a1_tu[s, time] ~ dnorm(mu_a_tu_hat[1, time], tau_a_tu_hat[1, time]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("a1_tc[s, time] ~ dnorm(mu_a_tc_hat[1, time], tau_a_tc_hat[1, time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_u) {a2_tu[s, time] ~ dnorm(mu_a_tu_hat[2, time], tau_a_tu_hat[2, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("a2_tc[s, time] ~ dnorm(mu_a_tc_hat[2, time], tau_a_tc_hat[2, time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_a_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_a_tu_hat[t, time] <- 1 / ss_a_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_tu_hat[t, time] <- s_a_tu_hat[t, time] * s_a_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_tu_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_a_tc_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_a_tc_hat[t, time] <- 1 / ss_a_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_tc_hat[t, time] <- s_a_tc_hat[t, time] * s_a_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_tc_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
 } 
  if(length(model_c_random) == 0) {
  model_string_jags <- gsub("+ inprod(X1_c_random[i, ], b1[, clus1_c[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(X1_c_random[i, ], b1[, clus1_c[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(X2_c_random[i, ], b2[, clus2_c[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(X2_c_random[i, ], b2[, clus2_c[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t, time] <- 1 / ss_b_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_hat[j, t, time] <- s_b_hat[j, t, time] * s_b_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s, time] ~ dnorm(mu_b_hat[j, 1, time], tau_b_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b1 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s, time] ~ dnorm(mu_b_hat[j, 2, time], tau_b_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b2 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b1_f[clus1_c[i], 1] * (eff1[i, 1] - mu_u[1, 1]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b1_f[clus1_c[i], time] * (eff1[i, time] - mu_u[1, time]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b1_tu[clus1_c[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + b1_tc[clus1_c[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b2_f[clus2_c[i], 1] * (eff2[i, 1] - mu_u[2, 1]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b2_f[clus2_c[i], time] * (eff2[i, time] - mu_u[2, time]) ", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b2_tu[clus2_c[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + b2_tc[clus2_c[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_f[s, time] ~ dnorm(mu_b_f_hat[1, time], tau_b_f_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_f[s, time] ~ dnorm(mu_b_f_hat[2, time], tau_b_f_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {mu_b_f_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat[t, time] <- 1 / ss_b_f_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat[t, time] <- s_b_f_hat[t, time] * s_b_f_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat[t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_tu[s, time] ~ dnorm(mu_b_tu_hat[1, time], tau_b_tu_hat[1, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("b1_tc[s, time] ~ dnorm(mu_b_tc_hat[1, time], tau_b_tc_hat[1, time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_tu[s, time] ~ dnorm(mu_b_tu_hat[2, time], tau_b_tu_hat[2, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("b2_tc[s, time] ~ dnorm(mu_b_tc_hat[2, time], tau_b_tc_hat[2, time])}", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_tu_hat[t, time] <- 1 / ss_b_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_tu_hat[t, time] <- s_b_tu_hat[t, time] * s_b_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_tu_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_b_tc_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_tc_hat[t, time] <- 1 / ss_b_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_tc_hat[t, time] <- s_b_tc_hat[t, time] * s_b_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_tc_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
  } else if(length(model_c_random) != 0 & is_c_random_c == TRUE & is_int_c_random_c == FALSE) {
    model_string_jags <- gsub("+ inprod(X1_c_random[i, ], b1[, clus1_c[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(X1_c_random[i, ], b1[, clus1_c[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(X2_c_random[i, ], b2[, clus2_c[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(X2_c_random[i, ], b2[, clus2_c[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t, time] <- 1 / ss_b_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j, t, time] <- s_b_hat[j, t, time] * s_b_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s, time] ~ dnorm(mu_b_hat[j, 1, time], tau_b_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b1 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s, time] ~ dnorm(mu_b_hat[j, 2, time], tau_b_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b2 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_mu_random) == 0) {
  model_string_jags <- gsub("+ inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zu_random) {tau_g_u_hat[j, t, time] <- 1 / ss_g_u_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_u_hat[j, t, time] <- s_g_u_hat[j, t, time] * s_g_u_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g1_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zu_random) {#begin g1_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_mu) {g1_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 1, time], tau_g_u_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g1_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end loop time g1_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g2_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zu_random) {#begin g2_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_mu) {g2_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 2, time], tau_g_u_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g2_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end loop time g2_u priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zu_random) {mu_g_u_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_u_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  if(length(model_mc_random) == 0 & !"c" %in% unlist(c(strsplit(as.character(fb(model.mc)), " " , fixed = TRUE)))) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_c_random) == 0 & length(model_u_random) == 0 | length(model_u_random) == 0 & pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for (t in 1:2) {# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_mu_random) != 0 & is_mu_random_u == TRUE & is_int_mu_random_u == FALSE) {
    model_string_jags <- gsub("+ inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zu_random) {tau_g_u_hat[j, t, time] <- 1 / ss_g_u_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_u_hat[j, t, time] <- s_g_u_hat[j, t, time] * s_g_u_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g1_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zu_random) {#begin g1_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mu) {g1_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 1, time], tau_g_u_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g1_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g1_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g2_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zu_random) {#begin g2_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mu) {g2_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 2, time], tau_g_u_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g2_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g2_u priors effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zu_random) {mu_g_u_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_u_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_mc_random) == 0) {
    model_string_jags <- gsub("+ inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t, time] <- 1 / ss_g_c_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j, t, time] <- s_g_c_hat[j, t, time] * s_g_c_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 1, time], tau_g_c_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 2, time], tau_g_c_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
  if(length(model_mu_random) == 0 & !"u" %in% unlist(c(strsplit(as.character(fb(model.mu)), " " , fixed = TRUE)))) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0 & length(model_c_random) == 0 | length(model_u_random) == 0 & pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for (t in 1:2) {# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_mc_random) != 0 & is_mc_random_c == TRUE & is_int_mc_random_c == FALSE) {
    model_string_jags <- gsub("+ inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t, time] <- 1 / ss_g_c_hat[j, t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_g_c_hat[j, t, time] <- s_g_c_hat[j, t, time] * s_g_c_hat[j, t, time] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_u priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 1, time], tau_g_c_hat[j, 1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g1_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 2, time], tau_g_c_hat[j, 2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time g2_c priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_g_c_hat[j, t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
  }
  if(type == "MAR") {
    model_string_jags <- gsub("+ delta_u[1, time, mdrop] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, time, mdrop] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[1, 1, mdrop] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[1, time, mdrop] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, 1, mdrop] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, time, mdrop] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_u[t, time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, time, mdrop] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, time, mdrop] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, 1, mdrop] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, time, mdrop] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, 1, mdrop] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, time, mdrop] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c[t, time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(mdrop in 1:3){# loop through mpattern mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(t in 1:2) {# begin mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop mpattern mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mnar parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_mu_random) != 0) {
    model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_mc_random) != 0) {
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)  
    }
  } else if(type == "MNAR_eff") {
    model_string_jags <- gsub("+ delta_c[1, time, mdrop] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, time, mdrop] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, 1, mdrop] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[1, time, mdrop] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, 1, mdrop] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_c[2, time, mdrop] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c[t, time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_mc_random) != 0) {
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!("c" %in% model_mc_random) & !"u" %in% unlist(c(strsplit(as.character(fb(model.mu)), " " , fixed = TRUE)))) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_mu_random) != 0 & !("u" %in% model_mu_random)) {
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(length(model_mc_random) == 0 & !"u" %in% unlist(c(strsplit(as.character(fb(model.mu)), " " , fixed = TRUE)))) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_mu_random) != 0 & !("u" %in% model_mu_random)) {
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    }
  } else if(type == "MNAR_cost") {
    model_string_jags <- gsub("+ delta_u[1, time, mdrop] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, time, mdrop] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[1, 1, mdrop] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[1, time, mdrop] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, 1, mdrop] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ delta_u[2, time, mdrop] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_u[t, time, mdrop] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_mu_random) != 0) {
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!("u" %in% model_mu_random) & !"c" %in% unlist(c(strsplit(as.character(fb(model.mc)), " " , fixed = TRUE)))) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random)) {
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(length(model_mu_random) == 0 & !"c" %in% unlist(c(strsplit(as.character(fb(model.mc)), " " , fixed = TRUE)))) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_mu_random) != 0 & !("u" %in% model_mu_random)) {
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR") {
    if(length(model_mu_random) != 0 & !("u" %in% model_mu_random)) {
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], 1] * eff1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], 1] * eff2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_u[clus1_mu[i], time] * eff1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_u[clus2_mu[i], time] * eff2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[1, time] * mean(eff1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_u_hat[2, time] * mean(eff2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mu) {d1_u[s, time] ~ dnorm(mu_d_u_hat[1, time], tau_d_u_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mu) {d2_u[s, time] ~ dnorm(mu_d_u_hat[2, time], tau_d_u_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_u_hat[t, time] <- 1 / ss_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_u_hat[t, time] <- s_d_u_hat[t, time] * s_d_u_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!"c" %in% unlist(c(strsplit(as.character(fb(model.mc)), " " , fixed = TRUE)))) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random)) {
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], 1] * cost1[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], 1] * cost2[i, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d1_c[clus1_mc[i], time] * cost1[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i], time] * cost2[i, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[1, time] * mean(cost1[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[2, time] * mean(cost2[, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s, time] ~ dnorm(mu_d_c_hat[1, time], tau_d_c_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s, time] ~ dnorm(mu_d_c_hat[2, time], tau_d_c_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[t, time] <- 1 / ss_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[t, time] <- s_d_c_hat[t, time] * s_d_c_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      if(!"u" %in% unlist(c(strsplit(as.character(fb(model.mu)), " " , fixed = TRUE)))) {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      }
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_mu_random) != 0 & !("u" %in% model_mu_random)) {
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time mnar random effects priors", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  }
  if(ind_random == TRUE | time_dep == "none") {
    if(ind_random == TRUE) {
      model_string_jags <- gsub("+ b1_f[clus1_c[i], 1] * (eff1[i, 1] - mu_u[1, 1]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ b1_f[clus1_c[i], time] * (eff1[i, time] - mu_u[1, time]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ b2_f[clus2_c[i], 1] * (eff2[i, 1] - mu_u[2, 1]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ b2_f[clus2_c[i], time] * (eff2[i, time] - mu_u[2, time]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_f[s, time] ~ dnorm(mu_b_f_hat[1, time], tau_b_f_hat[1, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_f[s, time] ~ dnorm(mu_b_f_hat[2, time], tau_b_f_hat[2, time]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {mu_b_f_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_f_hat[t, time] <- 1 / ss_b_f_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_f_hat[t, time] <- s_b_f_hat[t, time] * s_b_f_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_f_hat[t, time] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time correlation priors random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time correlation priors random effects", "", model_string_jags, fixed = TRUE)
    }
      model_string_jags <- gsub("+ a1_tu[clus1_u[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + a1_tc[clus1_u[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ a2_tu[clus2_u[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + a2_tc[clus2_u[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_u) {a1_tu[s, time] ~ dnorm(mu_a_tu_hat[1, time], tau_a_tu_hat[1, time]) ", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("a1_tc[s, time] ~ dnorm(mu_a_tc_hat[1, time], tau_a_tc_hat[1, time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_u) {a2_tu[s, time] ~ dnorm(mu_a_tu_hat[2, time], tau_a_tu_hat[2, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("a2_tc[s, time] ~ dnorm(mu_a_tc_hat[2, time], tau_a_tc_hat[2, time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_a_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_a_tu_hat[t, time] <- 1 / ss_a_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_a_tu_hat[t, time] <- s_a_tu_hat[t, time] * s_a_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_a_tu_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_a_tc_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_a_tc_hat[t, time] <- 1 / ss_a_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_a_tc_hat[t, time] <- s_a_tc_hat[t, time] * s_a_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_a_tc_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)  
      model_string_jags <- gsub("+ b1_tu[clus1_c[i], time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + b1_tc[clus1_c[i], time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ b2_tu[clus2_c[i], time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + b2_tc[clus2_c[i], time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_tu[s, time] ~ dnorm(mu_b_tu_hat[1, time], tau_b_tu_hat[1, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("b1_tc[s, time] ~ dnorm(mu_b_tc_hat[1, time], tau_b_tc_hat[1, time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_tu[s, time] ~ dnorm(mu_b_tu_hat[2, time], tau_b_tu_hat[2, time])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("b2_tc[s, time] ~ dnorm(mu_b_tc_hat[2, time], tau_b_tc_hat[2, time])}", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {# correlation time priors random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_tu_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_tu_hat[t, time] <- 1 / ss_b_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_tu_hat[t, time] <- s_b_tu_hat[t, time] * s_b_tu_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_tu_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_b_tc_hat[t, time] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_b_tc_hat[t, time] <- 1 / ss_b_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_b_tc_hat[t, time] <- s_b_tc_hat[t, time] * s_b_tc_hat[t, time]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_b_tc_hat[t, time] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop correlation time priors random effects", "", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time a1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end a1 loop time", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time a2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end a2 loop time", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end b1 loop time", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end b2 loop time", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) == 1) {
      if(model_c_random == "u") {
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b1", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end b1 loop time", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b2", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("}#end b2 loop time", "", model_string_jags, fixed = TRUE)
      }
    }
    if(length(model_u_random) == 0 & length(model_c_random) == 0 & length(model_mu_random) == 0 & length(model_mc_random) == 0) {
      model_string_jags <- gsub("for(time in 1:max_time) {# loop through time transformation of random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop time transformation of random effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("# correlation time random effects", "", model_string_jags, fixed = TRUE)
    }
  }
  if(ind_random == FALSE) {
   if(length(model_c_random) == 1) {
     if(model_c_random == "u") {
       model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end b1 loop time", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end b2 loop time", "", model_string_jags, fixed = TRUE)
     }
   }
  }
  if(ind_fixed == TRUE) {
    model_string_jags <- gsub("+ beta_f[1, 1] * (eff1[i, 1] - mu_u[1, 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ beta_f[2, 1] * (eff2[i, 1] - mu_u[2, 1])", "", model_string_jags, fixed = TRUE)  
    model_string_jags <- gsub("+ beta_f[1, time] * (eff1[i, time] - mu_u[1, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ beta_f[2, time] * (eff2[i, time] - mu_u[2, time])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[t, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
 if(ind_time_fixed == TRUE | time_dep == "none") {
   model_string_jags <- gsub("+ beta_tu[1, time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + beta_tc[1, time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("+ alpha_tu[1, time] * (eff1[i, time - 1] - mu_u[1, time - 1]) + alpha_tc[1, time] * (cost1[i, time - 1] - mu_c[1, time - 1])", "", model_string_jags, fixed = TRUE)  
   model_string_jags <- gsub("+ beta_tu[2, time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + beta_tc[2, time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("+ alpha_tu[2, time] * (eff2[i, time - 1] - mu_u[2, time - 1]) + alpha_tc[2, time] * (cost2[i, time - 1] - mu_c[2, time - 1])", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta_tu[t, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta_tc[t, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha_tu[t, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha_tc[t, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#time dependence", "", model_string_jags, fixed = TRUE)
 }  
  if(dist_c == "norm") {
  model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
  model_string_jags <- gsub("cost1[i, time] ~ dnorm(mu_c1[i, time], tau_c[1, time])", "cost1[i, time] ~ dgamma(mu_c1[i, time] * tau_c1[i, time], tau_c1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs1", "for(time in 1:max_time) { tau_c1[i, time] <- mu_c1[i, time] / pow(s_c[1, time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i, 1] <- ", "log(mu_c1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i, time] <- ", "log(mu_c1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cost2[i, time] ~ dnorm(mu_c2[i, time], tau_c[2, time])", "cost2[i, time] ~ dgamma(mu_c2[i, time] * tau_c2[i, time], tau_c2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "for(time in 1:max_time) { tau_c2[i, time] <- mu_c2[i, time] / pow(s_c[2, time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i, 1] <- ", "log(mu_c2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i, time] <- ", "log(mu_c2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[t, time] <- 1 / ss_c[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[t, time] <- s_c[t, time] * s_c[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[t, time] <- exp(ls_c[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_c[t, time] ~ dunif(-5, 10)", "s_c[t, time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c1[i, time] <- logdensity.norm(cost1[i, time], mu_c1[i, time], tau_c[1, time])", "loglik_c1[i, time] <- logdensity.gamma(cost1[i, time], mu_c1[i, time] * tau_c1[i, time], tau_c1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c2[i, time] <- logdensity.norm(cost2[i, time], mu_c2[i, time], tau_c[2, time])", "loglik_c2[i, time] <- logdensity.gamma(cost2[i, time], mu_c2[i, time] * tau_c2[i, time], tau_c2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
    if(length(model_c_random) == 0) {
      model_string_jags <- gsub("mu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time])", "mu_c[1, time] <- exp(inprod(mean_cov_c1_fixed[], beta[, 1, time]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time])", "mu_c[2, time] <- exp(inprod(mean_cov_c2_fixed[], beta[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) != 0) {
      model_string_jags <- gsub("mu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", "mu_c[1, time] <- exp(inprod(mean_cov_c1_fixed[], beta[, 1, time]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1, time]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", "mu_c[2, time] <- exp(inprod(mean_cov_c2_fixed[], beta[, 2, time]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_c == "lnorm") {
  model_string_jags <- gsub("cost1[i, time] ~ dnorm(mu_c1[i, time], tau_c[1, time])", "cost1[i, time] ~ dlnorm(lmu_c1[i, time], ltau_c[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i, 1] <- ", "lmu_c1[i, 1] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i, time] <- ", "lmu_c1[i, time] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cost2[i, time] ~ dnorm(mu_c2[i, time], tau_c[2, time])", "cost2[i, time] ~ dlnorm(lmu_c2[i, time], ltau_c[2, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i, 1] <- ", "lmu_c2[i, 1] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i, time] <- ", "lmu_c2[i, time] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[t, time] <- 1 / ss_c[t, time]", "ltau_c[t, time] <- 1 / lss_c[t, time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[t, time] <- s_c[t, time] * s_c[t, time]", "lss_c[t, time] <- ls_c[t, time] * ls_c[t, time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[t, time] <- exp(ls_c[t, time])", "s_c[t, time] <- sqrt(exp(2 * lmu_c[t, time] + lss_c[t, time]) * (exp(lss_c[t, time]) - 1))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_c[t, time] ~ dunif(-5, 10)", "ls_c[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "mu_c[t, time] <- exp(lmu_c[t, time] + lss_c[t, time] / 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c1[i, time] <- logdensity.norm(cost1[i, time], mu_c1[i, time], tau_c[1, time])", "loglik_c1[i, time] <- logdensity.lnorm(cost1[i, time], lmu_c1[i, time], ltau_c[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c2[i, time] <- logdensity.norm(cost2[i, time], mu_c2[i, time], tau_c[2, time])", "loglik_c2[i, time] <- logdensity.lnorm(cost2[i, time], lmu_c2[i, time], ltau_c[2, time])", model_string_jags, fixed = TRUE)
    if(length(model_c_random) == 0) {
      model_string_jags <- gsub("mu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time])", "lmu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time])", "lmu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time])", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) != 0) {
      model_string_jags <- gsub("mu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", "lmu_c[1, time] <- inprod(mean_cov_c1_fixed[], beta[, 1, time]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", "lmu_c[2, time] <- inprod(mean_cov_c2_fixed[], beta[, 2, time]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_u == "norm") {
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  }
  if(dist_u == "beta") {
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dbeta(mu_u1[i, time] * tau_u1[i, time], (1 - mu_u1[i, time]) * tau_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "for(time in 1:max_time) { tau_u1[i, time] <- (mu_u1[i, time] * (1 - mu_u1[i, time]) / pow(s_u[1, time], 2) - 1) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "logit(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "logit(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dbeta(mu_u2[i, time] * tau_u2[i, time], (1 - mu_u2[i, time]) * tau_u2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "for(time in 1:max_time) { tau_u2[i, time] <- (mu_u2[i, time] * (1 - mu_u2[i, time]) / pow(s_u[2, time], 2) - 1) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "logit(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "logit(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] ~ dunif(0, sqrt(mu_u[t, time] * (1 - mu_u[t, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.beta(eff1[i, time], mu_u1[i, time] * tau_u1[i, time], (1 - mu_u1[i, time]) * tau_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.beta(eff2[i, time], mu_u2[i, time] * tau_u2[i, time], (1 - mu_u2[i, time]) * tau_u2[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- ilogit(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- ilogit(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- ilogit(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])  + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- ilogit(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])  + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "gamma"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dgamma(mu_u1[i, time] * tau_u1[i, time], tau_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "for(time in 1:max_time) { tau_u1[i, time] <- mu_u1[i, time] / pow(s_u[1, time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "log(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "log(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dgamma(mu_u2[i, time] * tau_u2[i, time], tau_u2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "for(time in 1:max_time) { tau_u2[i, time] <- mu_u2[i, time] / pow(s_u[2, time], 2) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "log(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "log(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.gamma(eff1[i, time], mu_u1[i, time] * tau_u1[i, time], tau_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.gamma(eff2[i, time], mu_u2[i, time] * tau_u2[i, time], tau_u2[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "exp"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dexp(1 / mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "log(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "log(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dexp(1 / mu_u2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "log(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "log(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] <- mu_u[t, time]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.exp(eff1[i, time], 1/mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.exp(eff2[i, time], 1/mu_u2[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "weibull"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dweib(tau_u1[i, time], mu_u1[i, time] / exp(loggam(1 + 1/tau_u1[i, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "for(time in 1:max_time) { tau_u1[i, time] <- pow(s_u[1, time] / mu_u1[i, time], - 1.086) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "log(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "log(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dweib(tau_u2[i, time], mu_u2[i, time] / exp(loggam(1 + 1/tau_u2[i, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "for(time in 1:max_time) { tau_u2[i, time] <- pow(s_u[2, time] / mu_u2[i, time], - 1.086) }", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "log(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "log(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.weib(eff1[i, time], tau_u1[i, time], mu_u1[i, time] / exp(loggam(1 + 1/tau_u1[i, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.weib(eff2[i, time], tau_u2[i, time], mu_u2[i, time] / exp(loggam(1 + 1/tau_u2[i, time])))", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "logis"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dlogis(mu_u1[i, time], tau_u[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dlogis(mu_u2[i, time], tau_u[2, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "tau_u[t, time] <- 1 / sqrt((3 * ss_u[t, time]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.logis(eff1[i, time], mu_u1[i, time], tau_u[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.logis(eff2[i, time], mu_u2[i, time], tau_u[2, time])", model_string_jags, fixed = TRUE)
  } else if(dist_u == "bern"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dbern(mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "logit(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "logit(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dbern(mu_u2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "logit(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "logit(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] <- sqrt(mu_u[t, time] * (1 - mu_u[t, time]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.bern(eff1[i, time], mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.bern(eff2[i, time], mu_u2[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- ilogit(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- ilogit(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_e[1, time] <- ilogit(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_e[2, time] <- ilogit(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "pois"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dpois(mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "log(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "log(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dpois(mu_u2[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "log(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "log(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "s_u[t, time] <- sqrt(mu_u[t, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.pois(eff1[i, time], mu_u1[i, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.pois(eff2[i, time], mu_u2[i, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_u == "nbinom"){
  model_string_jags <- gsub("eff1[i, time] ~ dnorm(mu_u1[i, time], tau_u[1, time])", "eff1[i, time] ~ dnegbin(tau_u[1, time] / (tau_u[1, time] + mu_u1[i, time]), tau_u[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, 1] <- ", "log(mu_u1[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u1[i, time] <- ", "log(mu_u1[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i, time] ~ dnorm(mu_u2[i, time], tau_u[2, time])", "eff2[i, time] ~ dnegbin(tau_u[2, time] / (tau_u[2, time] + mu_u2[i, time]), tau_u[2, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, 1] <- ", "log(mu_u2[i, 1]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_u2[i, time] <- ", "log(mu_u2[i, time]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_u[t, time] <- 1 / ss_u[t, time]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_u[t, time] <- s_u[t, time] * s_u[t, time]", "s_u[t, time] <- sqrt((tau_u[t, time] / (tau_u[t, time] + mu_u[t, time])) * tau_u[t, time]) / (1 - (tau_u[t, time] / (tau_u[t, time] + mu_u[t, time])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_u[t, time] <- exp(ls_u[t, time])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10)", "tau_u[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u1[i, time] <- logdensity.norm(eff1[i, time], mu_u1[i, time], tau_u[1, time])", "loglik_u1[i, time] <- logdensity.negbin(eff1[i, time], tau_u[1, time] / (tau_u[1, time] + mu_u1[i, time]), tau_u[1, time])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_u2[i, time] <- logdensity.norm(eff2[i, time], mu_u2[i, time], tau_u[2, time])", "loglik_u2[i, time] <- logdensity.negbin(eff2[i, time], tau_u[2, time] / (tau_u[2, time] + mu_u2[i, time]), tau_u[2, time])", model_string_jags, fixed = TRUE)
    if(length(model_u_random) == 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_u_random) != 0) {
    model_string_jags <- gsub("mu_u[1, time] <- inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mu_u[1, time] <- exp(inprod(mean_cov_u1_fixed[], alpha[, 1, time]) + inprod(mean_cov_u1_random[], mu_a_hat[, 1, time]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_u[2, time] <- inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mu_u[2, time] <- exp(inprod(mean_cov_u2_fixed[], alpha[, 2, time]) + inprod(mean_cov_u2_random[], mu_a_hat[, 2, time]))", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_u %in% c("beta", "gamma", "weibull", "bern", "pois", "exp")) {
     if(dist_c == "gamma"){
     model_string_jags <- gsub("for (t in 1:2) {#begin transformation", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for(time in 1:max_time) {# loop through time transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end loop time transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}# end transformation of params", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
   }
   if(pu_fixed == 1) {
   inprod_u1_base <- "X1_u_fixed[i] * alpha[1, 1]"
   inprod_u1 <- "X1_u_fixed[i] * alpha[1, time]"
   inprod_u2_base <- "X2_u_fixed[i] * alpha[2, 1]"
   inprod_u2 <- "X2_u_fixed[i] * alpha[2, time]"
   inprod_mean_u1 <- "mean_cov_u1_fixed * alpha[1, time]"
   inprod_mean_u2 <- "mean_cov_u2_fixed * alpha[2, time]"
   begin_prior_beta <- "#begin alpha priors effects"
   prior_beta <- "#"
   end_prior_beta <- "#end alpha priors effects"
   prior_beta_u1 <- "alpha[1, time] ~ dnorm(0, 0.0000001)"
   prior_beta_u2 <- "alpha[2, time] ~ dnorm(0, 0.0000001)"
   model_string_jags <- gsub("inprod(X1_u_fixed[i, ], alpha[, 1, 1])", inprod_u1_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X1_u_fixed[i, ], alpha[, 1, time])", inprod_u1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_u_fixed[i, ], alpha[, 2, 1])", inprod_u2_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_u_fixed[i, ], alpha[, 2, time])", inprod_u2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_u1_fixed[], alpha[, 1, time])", inprod_mean_u1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_u2_fixed[], alpha[, 2, time])", inprod_mean_u2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:pu_fixed) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {alpha[j, t, time] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha[1, 1, time] ~ dnorm(0, 0.0000001)", prior_beta_u1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha[1, 2, time] ~ dnorm(0, 0.0000001)", prior_beta_u2, model_string_jags, fixed = TRUE)
   }
   if(length(model_u_random) != 0 & pu_random == 1) {
       model_string_jags <- gsub("inprod(X1_u_random[i, ], a1[, clus1_u[i], 1])", "X1_u_random[i] * a1[clus1_u[i], 1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X1_u_random[i, ], a1[, clus1_u[i], time])", "X1_u_random[i] * a1[clus1_u[i], time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X2_u_random[i, ], a2[, clus2_u[i], 1])", "X2_u_random[i] * a2[clus2_u[i], 1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X2_u_random[i, ], a2[, clus2_u[i], time])", "X2_u_random[i] * a2[clus2_u[i], time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_u1_random[], mu_a_hat[, 1, time])", "mean_cov_u1_random * mu_a_hat[1, time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_u2_random[], mu_a_hat[, 2, time])", "mean_cov_u2_random * mu_a_hat[2, time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pu_random) {#begin a1 priors effects", "#begin a1 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n1_clus_u) {a1[j, s, time] ~ dnorm(mu_a_hat[j, 1, time], tau_a_hat[j, 1, time]) }", "for(s in 1:n1_clus_u) {a1[s, time] ~ dnorm(mu_a_hat[1, time], tau_a_hat[1, time]) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a1 priors effects", "#end a1 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pu_random) {#begin a2 priors effects", "#begin a2 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n2_clus_u) {a2[j, s, time] ~ dnorm(mu_a_hat[j, 2, time], tau_a_hat[j, 2, time]) }", "for(s in 1:n2_clus_u) {a2[s, time] ~ dnorm(mu_a_hat[2, time], tau_a_hat[2, time]) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a2 priors effects", "#end a2 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pu_random) {mu_a_hat[j, t, time] ~ dnorm(0, 0.001)", "mu_a_hat[t, time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_a_hat[j, t, time] ~ dunif(0, 100) }", "s_a_hat[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pu_random) {tau_a_hat[j, t, time] <- 1 / ss_a_hat[j, t, time]", "tau_a_hat[t, time] <- 1 / ss_a_hat[t, time]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_a_hat[j, t, time] <- s_a_hat[j, t, time] * s_a_hat[j, t, time] }", "ss_a_hat[t, time] <- s_a_hat[t, time] * s_a_hat[t, time]", model_string_jags, fixed = TRUE)
   }
   if(pc_fixed == 1) {
   inprod_c1_base <- "X1_c_fixed[i] * beta[1, 1]"
   inprod_c1 <- "X1_c_fixed[i] * beta[1, time]"
   inprod_c2_base <- "X2_c_fixed[i] * beta[2, 1]"
   inprod_c2 <- "X2_c_fixed[i] * beta[2, time]"
   inprod_mean_c1 <- "mean_cov_c1_fixed * beta[1, time]"
   inprod_mean_c2 <- "mean_cov_c2_fixed * beta[2, time]"
   begin_prior_beta <- "#begin beta priors costs"
   prior_beta <- "#"
   end_prior_beta <- "#end beta priors costs"
   prior_beta_c1 <- "beta[1, time] ~ dnorm(0, 0.0000001)"
   prior_beta_c2 <- "beta[2, time] ~ dnorm(0, 0.0000001)"
   model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta[, 1, 1])", inprod_c1_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta[, 1, time])", inprod_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta[, 2, 1])", inprod_c2_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta[, 2, time])", inprod_c2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta[, 1, time])", inprod_mean_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta[, 2, time])", inprod_mean_c2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:pc_fixed) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {beta[j, t, time] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta priors costs", end_prior_beta,model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 1, time] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 2, time] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
   }
   if(pc_fixed == 0 & ind_fixed == FALSE) {
   prior_beta <- "#"
   end_prior_beta <- "#end beta priors costs"
   model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta[, 1, 1]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta[, 1, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta[, 2, 1]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta[, 2, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta[, 1, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta[, 2, time]) +", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta[, 1, time])", "beta_f[1, time] * (mean(eff1[, time]) - mu_u[1, time])", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta[, 2, time])", "beta_f[2, time] * (mean(eff2[, time]) - mu_u[2, time])", model_string_jags, fixed = TRUE)
   if(pc_random == 0 & ind_random == FALSE) {
   model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b1", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b1 loop time", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(time in 1:max_time) {# loop through time b2", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b2 loop time", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta[, 1, time])", "beta_f[1, time] * (mean(eff1[, time]) - mu_u[1, time])", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta[, 2, time])", "beta_f[2, time] * (mean(eff2[, time]) - mu_u[2, time])", model_string_jags, fixed = TRUE)
   }
   model_string_jags <- gsub("for (j in 2:pc_fixed) {#begin beta priors costs", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {beta[j, t, time] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta priors costs", end_prior_beta,model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 1, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 2, time] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(time in 1:max_time) {# loop through time beta", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#end beta priors costs", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta loop time", "", model_string_jags, fixed = TRUE)
   }
   if(length(model_c_random) != 0 & pc_random == 1) {
   model_string_jags <- gsub("inprod(X1_c_random[i, ], b1[, clus1_c[i], 1])", "X1_c_random[i] * b1[clus1_c[i], 1]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X1_c_random[i, ], b1[, clus1_c[i], time])", "X1_c_random[i] * b1[clus1_c[i], time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_random[i, ], b2[, clus2_c[i], 1])", "X2_c_random[i] * b2[clus2_c[i], 1]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_random[i, ], b2[, clus2_c[i], time])", "X2_c_random[i] * b2[clus2_c[i], time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_random[], mu_b_hat[, 1, time])", "mean_cov_c1_random * mu_b_hat[1, time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_random[], mu_b_hat[, 2, time])", "mean_cov_c2_random * mu_b_hat[2, time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "#begin b1 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s, time] ~ dnorm(mu_b_hat[j, 1, time], tau_b_hat[j, 1, time]) }", "for(s in 1:n1_clus_c) {b1[s, time] ~ dnorm(mu_b_hat[1, time], tau_b_hat[1, time]) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b1 priors costs", "#end b1 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "#begin b2 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s, time] ~ dnorm(mu_b_hat[j, 2, time], tau_b_hat[j, 2, time]) }", "for(s in 1:n2_clus_c) {b2[s, time] ~ dnorm(mu_b_hat[2, time], tau_b_hat[2, time]) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b2 priors costs", "#end b2 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t, time] ~ dnorm(0, 0.001)", "mu_b_hat[t, time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("s_b_hat[j, t, time] ~ dunif(0, 100) }", "s_b_hat[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t, time] <- 1 / ss_b_hat[j, t, time]", "tau_b_hat[t, time] <- 1 / ss_b_hat[t, time]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("ss_b_hat[j, t, time] <- s_b_hat[j, t, time] * s_b_hat[j, t, time] }", "ss_b_hat[t, time] <- s_b_hat[t, time] * s_b_hat[t, time]", model_string_jags, fixed = TRUE)
   }
   if(zu_fixed == 1) {
   inprod_u1_base <- "Z1_u_fixed[i] * gamma_u[1, 1, mdrop]"
   inprod_u1 <- "Z1_u_fixed[i] * gamma_u[1, time, mdrop]"
   inprod_u2_base <- "Z2_u_fixed[i] * gamma_u[2, 1, mdrop]"
   inprod_u2 <- "Z2_u_fixed[i] * gamma_u[2, time, mdrop]"
     if(type == "MAR" | type == "MNAR_cost") {
     inprod_mean_u1 <- "mean_z_u1_fixed * gamma_u[1, time, mdrop]"
     inprod_mean_u2 <- "mean_z_u2_fixed * gamma_u[2, time, mdrop]"
     model_string_jags <- gsub("inprod(mean_z_u1_fixed[], gamma_u[, 1, time, mdrop])", inprod_mean_u1, model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("inprod(mean_z_u2_fixed[], gamma_u[, 2, time, mdrop])", inprod_mean_u2, model_string_jags, fixed = TRUE)
     }
     if(type == "MNAR_eff" | type == "MNAR") {
     inprod_mean_u1 <- "mean_z_u1_fixed * gamma_u[1, time, mdrop] + delta_u[1, time, mdrop] * mean(eff1[, time])"
     inprod_mean_u2 <- "mean_z_u2_fixed * gamma_u[2, time, mdrop] + delta_u[2, time, mdrop] * mean(eff2[, time])"
     model_string_jags <- gsub("inprod(mean_z_u1_fixed[], gamma_u[, 1, time, mdrop]) + delta_u[1, time, mdrop] * mean(eff1[, time])", inprod_mean_u1, model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("inprod(mean_z_u2_fixed[], gamma_u[, 2, time, mdrop]) + delta_u[2, time, mdrop] * mean(eff2[, time])", inprod_mean_u2, model_string_jags, fixed = TRUE)
     }
   begin_prior_gamma <- "#begin gamma priors effects"
   begin_prior_gamma2 <- "#"
   prior_gamma_u1 <- "gamma_u[1, time, mdrop] ~ dlogis(0, 1)"
   prior_gamma_u2 <- "gamma_u[2, time, mdrop] ~ dlogis(0, 1)"
   end_prior_gamma <- "#end gamma priors effects"
   model_string_jags <- gsub("inprod(Z1_u_fixed[i, ], gamma_u[, 1, 1, mdrop])", inprod_u1_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(Z1_u_fixed[i, ], gamma_u[, 1, time, mdrop])", inprod_u1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(Z2_u_fixed[i, ], gamma_u[, 2, 1, mdrop])", inprod_u2_base, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(Z2_u_fixed[i, ], gamma_u[, 2, time, mdrop])", inprod_u2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:zu_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {gamma_u[j, t, time, mdrop] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_u[1, 1, time, mdrop] ~ dlogis(0, 1)", prior_gamma_u1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_u[1, 2, time, mdrop] ~ dlogis(0, 1)", prior_gamma_u2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(mdrop in 1:3){# loop through mpattern gamma priors effects", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end loop mpattern gamma priors effects", "", model_string_jags, fixed = TRUE)
   } 
   if(zu_random == 1 & length(model_mu_random) != 0) {
      if(type == "MAR" | type == "MNAR_cost") {
      model_string_jags <- gsub("inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time])", "mean_z_u1_random * mu_g_u_hat[1, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time])", "mean_z_u2_random * mu_g_u_hat[2, time]", model_string_jags, fixed = TRUE)
      }
      if(type == "MNAR_eff" | type == "MNAR") {
      model_string_jags <- gsub("inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time]) + mu_d_u_hat[1, time] * mean(eff1[, time])", "mean_z_u1_random * mu_g_u_hat[1, time] + mu_d_u_hat[1, time] * mean(eff1[, time])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time]) + mu_d_u_hat[2, time] * mean(eff2[, time])", "mean_z_u2_random * mu_g_u_hat[2, time] + mu_d_u_hat[2, time] * mean(eff2[, time])", model_string_jags, fixed = TRUE)
      }
      model_string_jags <- gsub("inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], 1])", "Z1_u_random[i] * g1_u[clus1_mu[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z1_u_random[i, ], g1_u[, clus1_mu[i], time])", "Z1_u_random[i] * g1_u[clus1_mu[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], 1])", "Z2_u_random[i] * g2_u[clus2_mu[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_u_random[i, ], g2_u[, clus2_mu[i], time])", "Z2_u_random[i] * g2_u[clus2_mu[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zu_random) {tau_g_u_hat[j, t, time] <- 1 / ss_g_u_hat[j, t, time]", "tau_g_u_hat[t, time] <- 1 / ss_g_u_hat[t, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_u_hat[j, t, time] <- s_g_u_hat[j, t, time] * s_g_u_hat[j, t, time] }", "ss_g_u_hat[t, time] <- s_g_u_hat[t, time] * s_g_u_hat[t, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_u1_random[], mu_g_u_hat[, 1, time])", "mean_z_u1_random * mu_g_u_hat[1, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_u2_random[], mu_g_u_hat[, 2, time])", "mean_z_u2_random * mu_g_u_hat[2, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zu_random) {#begin g1_u priors effects", "#begin g1_u priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mu) {g1_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 1, time], tau_g_u_hat[j, 1, time]) }", "for(s in 1:n1_clus_mu) {g1_u[s, time] ~ dnorm(mu_g_u_hat[1, time], tau_g_u_hat[1, time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g1_u priors effects", "#end g1_u priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zu_random) {#begin g2_u priors effects", "#begin g2_u priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mu) {g2_u[j, s, time] ~ dnorm(mu_g_u_hat[j, 2, time], tau_g_u_hat[j, 2, time]) }", "for(s in 1:n2_clus_mu) {g2_u[s, time] ~ dnorm(mu_g_u_hat[2, time], tau_g_u_hat[2, time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g2_u priors effects", "#end g2_u priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zu_random) {mu_g_u_hat[j, t, time] ~ dnorm(0, 0.001)", "mu_g_u_hat[t, time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_u_hat[j, t, time] ~ dunif(0, 100) }", "s_g_u_hat[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
   if(zc_fixed == 1) {
      inprod_c1_base <- "Z1_c_fixed[i] * gamma_c[1, 1, mdrop]"
      inprod_c1 <- "Z1_c_fixed[i] * gamma_c[1, time, mdrop]"
      inprod_c2_base <- "Z2_c_fixed[i] * gamma_c[2, 1, mdrop]"
      inprod_c2 <- "Z2_c_fixed[i] * gamma_c[2, time, mdrop]"
          if(type == "MAR" | type == "MNAR_eff") {
          inprod_mean_c1 <- "mean_z_c1_fixed * gamma_c[1, time, mdrop]"
          inprod_mean_c2 <- "mean_z_c2_fixed * gamma_c[2, time, mdrop]"
          model_string_jags <- gsub("inprod(mean_z_c1_fixed[], gamma_c[, 1, time, mdrop])", inprod_mean_c1, model_string_jags, fixed = TRUE)
          model_string_jags <- gsub("inprod(mean_z_c2_fixed[], gamma_c[, 2, time, mdrop])", inprod_mean_c2, model_string_jags, fixed = TRUE)
          }
          if(type == "MNAR_cost" | type == "MNAR") {
          inprod_mean_c1 <- "mean_z_c1_fixed * gamma_c[1, time, mdrop] + delta_c[1, time, mdrop] * mean(cost1[, time])"
          inprod_mean_c2 <- "mean_z_c2_fixed * gamma_c[2, time, mdrop] + delta_c[2, time, mdrop] * mean(cost2[, time])"
          model_string_jags <- gsub("inprod(mean_z_c1_fixed[], gamma_c[, 1, time, mdrop]) + delta_c[1, time, mdrop] * mean(cost1[, time])", inprod_mean_c1, model_string_jags, fixed = TRUE)
          model_string_jags <- gsub("inprod(mean_z_c2_fixed[], gamma_c[, 2, time, mdrop]) + delta_c[2, time, mdrop] * mean(cost2[, time])", inprod_mean_c2, model_string_jags, fixed = TRUE)
          }
      begin_prior_gamma <- "#begin gamma priors costs"
      begin_prior_gamma2 <- "#"
      prior_gamma_c1 <- "gamma_c[1, time, mdrop] ~ dlogis(0, 1)"
      prior_gamma_c2 <- "gamma_c[2, time, mdrop] ~ dlogis(0, 1)"
      end_prior_gamma <- "#end gamma priors costs"
      model_string_jags <- gsub("inprod(Z1_c_fixed[i, ], gamma_c[, 1, 1, mdrop])", inprod_c1_base, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z1_c_fixed[i, ], gamma_c[, 1, time, mdrop])", inprod_c1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_fixed[i, ], gamma_c[, 2, 1, mdrop])", inprod_c2_base, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_fixed[i, ], gamma_c[, 2, time, mdrop])", inprod_c2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 2:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t, time, mdrop] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[1, 1, time, mdrop] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[1, 2, time, mdrop] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(mdrop in 1:3){# loop through mpattern gamma priors costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end loop mpattern gamma priors costs", "", model_string_jags, fixed = TRUE)
      }
      if(zc_random == 1 & length(model_mc_random) != 0) {
         if(type == "MAR" | type == "MNAR_eff") {
         model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time])", "mean_z_c1_random * mu_g_c_hat[1, time]", model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time])", "mean_z_c2_random * mu_g_c_hat[2, time]", model_string_jags, fixed = TRUE)
         }
         if(type == "MNAR_cost" | type == "MNAR") {
         model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time]) + mu_d_c_hat[1, time] * mean(cost1[, time])", "mean_z_c1_random * mu_g_c_hat[1, time] + mu_d_c_hat[1, time] * mean(cost1[, time])", model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time]) + mu_d_c_hat[2, time] * mean(cost2[, time])", "mean_z_c2_random * mu_g_c_hat[2, time] + mu_d_c_hat[2, time] * mean(cost2[, time])", model_string_jags, fixed = TRUE)
         }
      model_string_jags <- gsub("inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], 1])", "Z1_c_random[i] * g1_c[clus1_mc[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i], time])", "Z1_c_random[i] * g1_c[clus1_mc[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], 1])", "Z2_c_random[i] * g2_c[clus2_mc[i], 1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i], time])", "Z2_c_random[i] * g2_c[clus2_mc[i], time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t, time] <- 1 / ss_g_c_hat[j, t, time]", "tau_g_c_hat[t, time] <- 1 / ss_g_c_hat[t, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_c_hat[j, t, time] <- s_g_c_hat[j, t, time] * s_g_c_hat[j, t, time] }", "ss_g_c_hat[t, time] <- s_g_c_hat[t, time] * s_g_c_hat[t, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1, time])", "mean_z_c1_random * mu_g_c_hat[1, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2, time])", "mean_z_c2_random * mu_g_c_hat[2, time]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_c priors costs", "#begin g1_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 1, time], tau_g_c_hat[j, 1, time]) }", "for(s in 1:n1_clus_mc) {g1_c[s, time] ~ dnorm(mu_g_c_hat[1, time], tau_g_c_hat[1, time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g1_c priors costs", "#end g1_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "#begin g2_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s, time] ~ dnorm(mu_g_c_hat[j, 2, time], tau_g_c_hat[j, 2, time]) }", "for(s in 1:n2_clus_mc) {g2_c[s, time] ~ dnorm(mu_g_c_hat[2, time], tau_g_c_hat[2, time]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g2_c priors costs", "#end g2_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t, time] ~ dnorm(0, 0.001)", "mu_g_c_hat[t, time] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_c_hat[j, t, time] ~ dunif(0, 100) }", "s_g_c_hat[t, time] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
      model_string_jags <- prior_selection_long(type = type, dist_u = dist_u, dist_c = dist_c, pu_fixed = pu_fixed, pc_fixed = pc_fixed, zu_fixed = zu_fixed, zc_fixed = zc_fixed, 
                                           model_u_random = model_u_random, model_c_random = model_c_random, model_mu_random = model_mu_random, model_mc_random = model_mc_random,
                                           pu_random = pu_random, pc_random = pc_random, zu_random = zu_random, zc_random = zc_random)
   writeLines(model_string_jags, "selection_long.txt")
   model_string <- "selection_long.txt"
   return(model_string)
   }))