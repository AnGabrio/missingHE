#'An internal function to select which type of selection model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Selection models
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param ind_fixed Logical; if TRUE independence between effectiveness and costs is assumed, else correlation is accounted for
#' @param ind_random Logical; if TRUE independence at the level of the random effects between effectiveness and costs is assumed, else correlation is accounted for
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param pe_fixed Number of fixed effects for the effectiveness model
#' @param pc_fixed Number of fixed effects for the cost model
#' @param ze_fixed Number of fixed effects or the missingness indicators model for the effectiveness
#' @param zc_fixed Number of fixed effects or the missingness indicators model for the costs
#' @param pe_random Number of random effects for the effectiveness model
#' @param pc_random Number of random effects for the cost model
#' @param ze_random Number of random effects or the missingness indicators model for the effectiveness
#' @param zc_random Number of random effects or the missingness indicators model for the costs
#' @param model_e_random Random effects formula for the effectiveness model
#' @param model_c_random Random effects formula for the costs model
#' @param model_me_random Random effects formula for the missingness indicators model for the effectiveness
#' @param model_mc_random Random effects formula for the missingness indicators model for the costs
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_selection <- function(dist_e , dist_c, type, pe_fixed, pc_fixed, ze_fixed, zc_fixed, ind_fixed, pe_random, pc_random, ze_random, zc_random, ind_random, 
                            model_e_random, model_c_random, model_me_random, model_mc_random) eval.parent(substitute( {
  model_string_jags<-  "
  model {

  #control
  for(i in 1:N1) {
  #costs and effects model
  cost1[i] ~ dnorm(mu_c1[i], tau_c[1])
  eff1[i] ~ dnorm(mu_e1[i], tau_e[1])

  #derive mean and std effects1 
  #derive mean and std costs1

  #mean regression
  mu_c1[i] <- inprod(X1_c_fixed[i, ], beta[, 1]) + beta_f[1] * (eff1[i] - mu_e[1]) + inprod(X1_c_random[i, ], b1[, clus1_c[i]]) + b1_f[clus1_c[i]] * (eff1[i] - mu_e[1])
  mu_e1[i] <- inprod(X1_e_fixed[i, ], alpha[, 1]) + inprod(X1_e_random[i, ], a1[, clus1_e[i]])   

  #missing data mechanism
  m_eff1[i] ~ dbern(pq_1[i])
  logit(pq_1[i]) <- inprod(Z1_e_fixed[i, ], gamma_e[, 1]) + delta_e[1] * eff1[i] + inprod(Z1_e_random[i, ], g1_e[, clus1_me[i]]) + d1_e[clus1_me[i]] * eff1[i]
  m_cost1[i] ~ dbern(pc_1[i])
  logit(pc_1[i]) <- inprod(Z1_c_fixed[i, ], gamma_c[, 1]) + delta_c[1] * cost1[i] + inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i]]) + d1_c[clus1_mc[i]] * cost1[i]

  #loglikelihood
  loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])
  loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c[1])
  loglik_me1[i] <- logdensity.bern(m_eff1[i], pq_1[i])
  loglik_mc1[i] <- logdensity.bern(m_cost1[i], pc_1[i])
  }
  

  #intervention
  for(i in 1:N2) {
  #costs and effects model
  cost2[i] ~ dnorm(mu_c2[i], tau_c[2])
  eff2[i] ~ dnorm(mu_e2[i], tau_e[2])

  #derive mean and std effects2 
  #derive mean and std costs2

  #mean regression
  mu_c2[i] <- inprod(X2_c_fixed[i, ], beta[, 2]) + beta_f[2] * (eff2[i] - mu_e[2]) + inprod(X2_c_random[i, ], b2[, clus2_c[i]]) + b2_f[clus2_c[i]] * (eff2[i] - mu_e[2])
  mu_e2[i] <- inprod(X2_e_fixed[i, ], alpha[, 2]) + inprod(X2_e_random[i, ], a2[, clus2_e[i]])

  #missing data mechanism
  m_eff2[i] ~ dbern(pq_2[i])
  logit(pq_2[i]) <- inprod(Z2_e_fixed[i, ], gamma_e[, 2]) + delta_e[2] * eff2[i] + inprod(Z2_e_random[i, ], g2_e[, clus2_me[i]]) + d2_e[clus2_me[i]] * eff2[i]
  m_cost2[i] ~ dbern(pc_2[i])
  logit(pc_2[i]) <- inprod(Z2_c_fixed[i, ], gamma_c[, 2]) + delta_c[2] * cost2[i] + inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i]]) + d2_c[clus2_mc[i]] * cost2[i]

  #loglikelihood
  loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])
  loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c[2])
  loglik_me2[i] <- logdensity.bern(m_eff2[i], pq_2[i])
  loglik_mc2[i] <- logdensity.bern(m_cost2[i], pc_2[i])
  }
  
  #transformation of parameters
  for (t in 1:2) {#begin transformation
  tau_c[t] <- 1 / ss_c[t]
  ss_c[t] <- s_c[t] * s_c[t]
  s_c[t] <- exp(ls_c[t])
  #mean for lnorm
  tau_e[t] <- 1 / ss_e[t]
  ss_e[t] <- s_e[t] * s_e[t]
  s_e[t] <- exp(ls_e[t])
  }# end transformation
  
  #transformation of random effects parameters
  for (t in 1:2) {# begin transformation random effects
  for(j in 1:pc_random) {tau_b_hat[j, t] <- 1 / ss_b_hat[j, t]
  ss_b_hat[j, t] <- s_b_hat[j, t] * s_b_hat[j, t] }
  for(j in 1:pe_random) {tau_a_hat[j, t] <- 1 / ss_a_hat[j, t]
  ss_a_hat[j, t] <- s_a_hat[j, t] * s_a_hat[j, t] }
  for(j in 1:zc_random) {tau_g_c_hat[j, t] <- 1 / ss_g_c_hat[j, t]
  ss_g_c_hat[j, t] <- s_g_c_hat[j, t] * s_g_c_hat[j, t] }
  for(j in 1:ze_random) {tau_g_e_hat[j, t] <- 1 / ss_g_e_hat[j, t]
  ss_g_e_hat[j, t] <- s_g_e_hat[j, t] * s_g_e_hat[j, t] }
  tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]
  ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]
  tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]
  ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]
  }#end transformation of random effects 
  
  #missingness probability
  p_c[1] <- ilogit(inprod(mean_z_c1_fixed[], gamma_c[, 1]) + delta_c[1] * mean(cost1[]) + inprod(mean_z_c1_random[], mu_g_c_hat[, 1]) + mu_d_c_hat[1] * mean(cost1[]))
  p_c[2] <- ilogit(inprod(mean_z_c2_fixed[], gamma_c[, 2]) + delta_c[2] * mean(cost2[]) + inprod(mean_z_c2_random[], mu_g_c_hat[, 2]) + mu_d_c_hat[2] * mean(cost2[]))
  p_e[1] <- ilogit(inprod(mean_z_e1_fixed[], gamma_e[, 1]) + delta_e[1] * mean(eff1[]) + inprod(mean_z_e1_random[], mu_g_e_hat[, 1]) + mu_d_e_hat[1] * mean(eff1[]))
  p_e[2] <- ilogit(inprod(mean_z_e2_fixed[], gamma_e[, 2]) + delta_e[2] * mean(eff2[]) + inprod(mean_z_e2_random[], mu_g_e_hat[, 2]) + mu_d_e_hat[2] * mean(eff2[]))
  
  #calculate means at mean of covariates
  mu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1]) 
  mu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2])
  mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])
  mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])
  
  #priors
  
  #priors for mean regression coefficients
  for (j in 2:pe_fixed) {#begin alpha priors effects
  for(t in 1:2) {alpha[j, t] ~ dnorm(0, 0.0000001) }
  }#end alpha priors effects
  alpha[1, 1] ~ dnorm(0, 0.0000001)
  alpha[1, 2] ~ dnorm(0, 0.0000001)
  
  for (j in 2:pc_fixed) {#begin beta priors costs
  for(t in 1:2) {beta[j, t] ~ dnorm(0, 0.0000001) }
  }#end beta priors costs
  beta[1, 1] ~ dnorm(0, 0.0000001)
  beta[1, 2] ~ dnorm(0, 0.0000001)
  
  #priors for mean regression random coefficients
  for (j in 1:pe_random) {#begin a1 priors effects
  for(s in 1:n1_clus_e) {a1[j, s] ~ dnorm(mu_a_hat[j, 1], tau_a_hat[j, 1]) }
  }#end a1 priors effects
  for (j in 1:pe_random) {#begin a2 priors effects
  for(s in 1:n2_clus_e) {a2[j, s] ~ dnorm(mu_a_hat[j, 2], tau_a_hat[j, 2]) }
  }#end a2 priors effects

  for (j in 1:pc_random) {#begin b1 priors costs
  for(s in 1:n1_clus_c) {b1[j, s] ~ dnorm(mu_b_hat[j, 1], tau_b_hat[j, 1]) }
  }#end b1 priors costs
  for (j in 1:pc_random) {#begin b2 priors costs
  for(s in 1:n2_clus_c) {b2[j, s] ~ dnorm(mu_b_hat[j, 2], tau_b_hat[j, 2]) }
  }#end b2 priors costs
  
  #standard deviation priors
  for(t in 1:2) {
  ls_c[t] ~ dunif(-5, 10)
  ls_e[t] ~ dunif(-5, 10)

  #correlation
  beta_f[t] ~ dnorm(0, 0.0000001)
  
  # mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j, t] ~ dnorm(0, 0.001)
  s_b_hat[j, t] ~ dunif(0, 100) }
  for(j in 1:pe_random) {mu_a_hat[j, t] ~ dnorm(0, 0.001)
  s_a_hat[j, t] ~ dunif(0, 100) }
  for(j in 1:zc_random) {mu_g_c_hat[j, t] ~ dnorm(0, 0.001)
  s_g_c_hat[j, t] ~ dunif(0, 100) }
  for(j in 1:ze_random) {mu_g_e_hat[j, t] ~ dnorm(0, 0.001)
  s_g_e_hat[j, t] ~ dunif(0, 100) }
  }
  
  # correlation random effects
  for(s in 1:n1_clus_c) {b1_f[s] ~ dnorm(mu_b_f_hat[1], tau_b_f_hat[1]) }
  for(s in 1:n2_clus_c) {b2_f[s] ~ dnorm(mu_b_f_hat[2], tau_b_f_hat[2]) }
  for(t in 1:2) {mu_b_f_hat[t] ~ dnorm(0, 0.001)
  tau_b_f_hat[t] <- 1 / ss_b_f_hat[t]
  ss_b_f_hat[t] <- s_b_f_hat[t] * s_b_f_hat[t]
  s_b_f_hat[t] ~ dunif(0, 100) }
  
  #priors on missing data mechanism
  for (j in 2:ze_fixed) {#begin gamma priors effects
  for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }
  }#end gamma priors effects
  gamma_e[1, 1] ~ dlogis(0, 1)
  gamma_e[1, 2] ~ dlogis(0, 1)
  
  for (j in 2:zc_fixed) {#begin gamma priors costs
  for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }
  }#end gamma priors costs
  gamma_c[1, 1] ~ dlogis(0, 1)
  gamma_c[1, 2] ~ dlogis(0, 1)
  
  #priors on random effects missing data mechanism
  for (j in 1:ze_random) {#begin g1_e priors effects
  for(s in 1:n1_clus_me) {g1_e[j, s] ~ dnorm(mu_g_e_hat[j, 1], tau_g_e_hat[j, 1]) }
  }#end g1_e priors effects
  for (j in 1:ze_random) {#begin g2_e priors effects
  for(s in 1:n2_clus_me) {g2_e[j, s] ~ dnorm(mu_g_e_hat[j, 2], tau_g_e_hat[j, 2]) }
  }#end g2_e priors effects

  for (j in 1:zc_random) {#begin g1_c priors costs
  for(s in 1:n1_clus_mc) {g1_c[j, s] ~ dnorm(mu_g_c_hat[j, 1], tau_g_c_hat[j, 1]) }
  }#end g1_c priors costs
  for (j in 1:zc_random) {#begin g2_c priors costs
  for(s in 1:n2_clus_mc) {g2_c[j, s] ~ dnorm(mu_g_c_hat[j, 2], tau_g_c_hat[j, 2]) }
  }#end g2_c priors costs
  
  #mnar parameters
  for(t in 1:2) {# begin mnar priors
  delta_e[t] ~ dnorm(0, 1)
  delta_c[t] ~ dnorm(0, 1)
  }#end mnar priors

  #mnar random effects priors
  for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }
  for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }
  for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }
  for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }

  # mean and sd mean mnar random effects priors
  for(t in 1:2) {#begin mnar random effects priors
  mu_d_e_hat[t] ~ dnorm(0, 1)
  mu_d_c_hat[t] ~ dnorm(0, 1)
  s_d_e_hat[t] ~ dunif(0, 1)
  s_d_c_hat[t] ~ dunif(0, 1) 
  } #end mnar random effects priors

}
 "
 if(length(model_e_random) == 0) {
  model_string_jags <- gsub(" + inprod(X1_e_random[i, ], a1[, clus1_e[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X2_e_random[i, ], a2[, clus2_e[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j, t] <- 1 / ss_a_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_a_hat[j, t] <- s_a_hat[j, t] * s_a_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pe_random) {#begin a1 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_e) {a1[j, s] ~ dnorm(mu_a_hat[j, 1], tau_a_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a1 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pe_random) {#begin a2 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_e) {a2[j, s] ~ dnorm(mu_a_hat[j, 2], tau_a_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end a2 priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_a_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_c_random) == 0) {
  model_string_jags <- gsub(" + inprod(X1_c_random[i, ], b1[, clus1_c[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(X2_c_random[i, ], b2[, clus2_c[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t] <- 1 / ss_b_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_hat[j, t] <- s_b_hat[j, t] * s_b_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s] ~ dnorm(mu_b_hat[j, 1], tau_b_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b1 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s] ~ dnorm(mu_b_hat[j, 2], tau_b_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end b2 priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b1_f[clus1_c[i]] * (eff1[i] - mu_e[1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b2_f[clus2_c[i]] * (eff2[i] - mu_e[2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_f[s] ~ dnorm(mu_b_f_hat[1], tau_b_f_hat[1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_f[s] ~ dnorm(mu_b_f_hat[2], tau_b_f_hat[2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {mu_b_f_hat[t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat[t] <- 1 / ss_b_f_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat[t] <- s_b_f_hat[t] * s_b_f_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat[t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  if(length(model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE) }
  } else if(length(model_c_random) != 0 & is_c_random_c == TRUE & is_int_c_random_c == FALSE) {
    model_string_jags <- gsub(" + inprod(X1_c_random[i, ], b1[, clus1_c[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(X2_c_random[i, ], b2[, clus2_c[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t] <- 1 / ss_b_hat[j, t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j, t] <- s_b_hat[j, t] * s_b_hat[j, t] }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s] ~ dnorm(mu_b_hat[j, 1], tau_b_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b1 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s] ~ dnorm(mu_b_hat[j, 2], tau_b_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b2 priors costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE) }
  }
  if(length(model_me_random) == 0) {
  model_string_jags <- gsub(" + inprod(Z1_e_random[i, ], g1_e[, clus1_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(Z2_e_random[i, ], g2_e[, clus2_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, t] <- 1 / ss_g_e_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j, t] <- s_g_e_hat[j, t] * s_g_e_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e1_random[], mu_g_e_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e2_random[], mu_g_e_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g1_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_me) {g1_e[j, s] ~ dnorm(mu_g_e_hat[j, 1], tau_g_e_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g1_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g2_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_me) {g2_e[j, s] ~ dnorm(mu_g_e_hat[j, 2], tau_g_e_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g2_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d1_e[clus1_me[i]] * eff1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d2_e[clus2_me[i]] * eff2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_e_hat[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_e_hat[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_e_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_e_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }", "", model_string_jags, fixed = TRUE)
  if(length(model_mc_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    if(length(model_c_random) == 0 & length(model_e_random) == 0 | length(model_e_random) == 0 & pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for (t in 1:2) {# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_me_random) != 0 & is_me_random_e == TRUE & is_int_me_random_e == FALSE) {
  model_string_jags <- gsub(" + inprod(Z1_e_random[i, ], g1_e[, clus1_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(Z2_e_random[i, ], g2_e[, clus2_me[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, t] <- 1 / ss_g_e_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_e_hat[j, t] <- s_g_e_hat[j, t] * s_g_e_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e1_random[], mu_g_e_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_e2_random[], mu_g_e_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g1_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_me) {g1_e[j, s] ~ dnorm(mu_g_e_hat[j, 1], tau_g_e_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g1_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:ze_random) {#begin g2_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_me) {g2_e[j, s] ~ dnorm(mu_g_e_hat[j, 2], tau_g_e_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g2_e priors effects", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_e_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(length(model_mc_random) == 0) {
  model_string_jags <- gsub(" + inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t] <- 1 / ss_g_c_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j, t] <- s_g_c_hat[j, t] * s_g_c_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c1_random[], mu_g_c_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c2_random[], mu_g_c_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s] ~ dnorm(mu_g_c_hat[j, 1], tau_g_c_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g1_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s] ~ dnorm(mu_g_c_hat[j, 2], tau_g_c_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g2_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d1_c[clus1_mc[i]] * cost1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ d2_c[clus2_mc[i]] * cost2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_c_hat[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ mu_d_c_hat[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_d_c_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_d_c_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }", "", model_string_jags, fixed = TRUE)
  if(length(model_me_random) == 0) { 
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    if(length(model_e_random) == 0 & length(model_c_random) == 0 | length(model_e_random) == 0 & pc_random == 0) { 
      model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for (t in 1:2) {# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
    }
   }
  } else if(length(model_mc_random) != 0 & is_mc_random_c == TRUE & is_int_mc_random_c == FALSE) {
  model_string_jags <- gsub(" + inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i]])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t] <- 1 / ss_g_c_hat[j, t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_g_c_hat[j, t] <- s_g_c_hat[j, t] * s_g_c_hat[j, t] }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c1_random[], mu_g_c_hat[, 1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + inprod(mean_z_c2_random[], mu_g_c_hat[, 2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s] ~ dnorm(mu_g_c_hat[j, 1], tau_g_c_hat[j, 1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g1_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s] ~ dnorm(mu_g_c_hat[j, 2], tau_g_c_hat[j, 2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end g2_c priors costs", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_g_c_hat[j, t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  }
  if(type == "MAR") {
  model_string_jags <- gsub(" + delta_e[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[1] * eff1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * eff2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_e[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[1] * cost1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * cost2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_c[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {# begin mnar priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end mnar priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mnar parameters", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    if(length(model_me_random) != 0) {
    model_string_jags <- gsub("+ d1_e[clus1_me[i]] * eff1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_e[clus2_me[i]] * eff2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_mc_random) != 0) {
    model_string_jags <- gsub("+ d1_c[clus1_mc[i]] * cost1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i]] * cost2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }", "", model_string_jags, fixed = TRUE)
    }
  } else if(type == "MNAR_eff") {
  model_string_jags <- gsub(" + delta_c[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[1] * cost1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * cost2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_c[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_mc_random) != 0) {
    model_string_jags <- gsub("+ d1_c[clus1_mc[i]] * cost1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i]] * cost2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_me_random) != 0 & !("e" %in% model_me_random)) {
    model_string_jags <- gsub("+ d1_e[clus1_me[i]] * eff1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_e[clus2_me[i]] * eff2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_me_random) != 0 & !("e" %in% model_me_random)) {
    model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR_cost") {
  model_string_jags <- gsub(" + delta_e[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[1] * eff1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * eff2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_e[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    if(length(model_me_random) != 0) {
    model_string_jags <- gsub("+ d1_e[clus1_me[i]] * eff1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_e[clus2_me[i]] * eff2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_e_hat[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_e_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_e_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }", "", model_string_jags, fixed = TRUE)
    }
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random)) {
    model_string_jags <- gsub("+ d1_c[clus1_mc[i]] * cost1[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ d2_c[clus2_mc[i]] * cost2[i]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("+ mu_d_c_hat[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_d_c_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_d_c_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_me_random) != 0 & !("e" %in% model_me_random)) {
      model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  } else if(type == "MNAR") {
    if(length(model_me_random) != 0 & !("e" %in% model_me_random)) {
      model_string_jags <- gsub("+ d1_e[clus1_me[i]] * eff1[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_e[clus2_me[i]] * eff2[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_e_hat[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_e_hat[t] <- 1 / ss_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_e_hat[t] <- s_d_e_hat[t] * s_d_e_hat[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_e_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_e_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_me) {d1_e[s] ~ dnorm(mu_d_e_hat[1], tau_d_e_hat[1]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_me) {d2_e[s] ~ dnorm(mu_d_e_hat[2], tau_d_e_hat[2]) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random)) {
      model_string_jags <- gsub("+ d1_c[clus1_mc[i]] * cost1[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ d2_c[clus2_mc[i]] * cost2[i]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ mu_d_c_hat[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_d_c_hat[t] <- 1 / ss_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_d_c_hat[t] <- s_d_c_hat[t] * s_d_c_hat[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_d_c_hat[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_d_c_hat[t] ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {d1_c[s] ~ dnorm(mu_d_c_hat[1], tau_d_c_hat[1]) }", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {d2_c[s] ~ dnorm(mu_d_c_hat[2], tau_d_c_hat[2]) }", "", model_string_jags, fixed = TRUE)
    } 
    if(length(model_mc_random) != 0 & !("c" %in% model_mc_random) & length(model_me_random) != 0 & !("e" %in% model_me_random)) {
      model_string_jags <- gsub("#priors on random effects missing data mechanism", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("#mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("# mean and sd mean mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("for(t in 1:2) {#begin mnar random effects priors", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("} #end mnar random effects priors", "", model_string_jags, fixed = TRUE)
    } 
  }
  if(ind_random == TRUE) {
  model_string_jags <- gsub("+ b1_f[clus1_c[i]] * (eff1[i] - mu_e[1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("+ b2_f[clus2_c[i]] * (eff2[i] - mu_e[2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1_f[s] ~ dnorm(mu_b_f_hat[1], tau_b_f_hat[1]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2_f[s] ~ dnorm(mu_b_f_hat[2], tau_b_f_hat[2]) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {mu_b_f_hat[t] ~ dnorm(0, 0.001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_b_f_hat[t] <- 1 / ss_b_f_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_b_f_hat[t] <- s_b_f_hat[t] * s_b_f_hat[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_b_f_hat[t] ~ dunif(0, 100) }", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("# correlation random effects", "", model_string_jags, fixed = TRUE)
  } 
  if(ind_fixed == TRUE) {
  model_string_jags <- gsub(" + beta_f[1] * (eff1[i] - mu_e[1])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + beta_f[2] * (eff2[i] - mu_e[2])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta_f[t] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
  if(dist_c == "norm") {
  model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
  model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c[1])", "cost1[i] ~ dgamma(mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs1", "tau_c1[i] <- mu_c1[i] / pow(s_c[1], 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i] <- ", "log(mu_c1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c[2])", "cost2[i] ~ dgamma(mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "tau_c2[i] <- mu_c2[i] / pow(s_c[2], 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i] <- ", "log(mu_c2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[t] <- 1 / ss_c[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[t] <- s_c[t] * s_c[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[t] <- exp(ls_c[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_c[t] ~ dunif(-5, 10)", "s_c[t] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c[1])", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c[2])", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
    if(length(model_c_random) == 0) {
    model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1])", "mu_c[1] <- exp(inprod(mean_cov_c1_fixed[], beta[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2])", "mu_c[2] <- exp(inprod(mean_cov_c2_fixed[], beta[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) != 0) {
    model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "mu_c[1] <- exp(inprod(mean_cov_c1_fixed[], beta[, 1]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "mu_c[2] <- exp(inprod(mean_cov_c2_fixed[], beta[, 2]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_c == "lnorm") {
  model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c[1])", "cost1[i] ~ dlnorm(lmu_c1[i], ltau_c[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c1[i] <- ", "lmu_c1[i] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c[2])", "cost2[i] ~ dlnorm(lmu_c2[i], ltau_c[2])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c2[i] <- ", "lmu_c2[i] <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_c[t] <- 1 / ss_c[t]", "ltau_c[t] <- 1 / lss_c[t]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_c[t] <- s_c[t] * s_c[t]", "lss_c[t] <- ls_c[t] * ls_c[t]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_c[t] <- exp(ls_c[t])", "s_c[t] <- sqrt(exp(2 * lmu_c[t] + lss_c[t]) * (exp(lss_c[t]) - 1))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_c[t] ~ dunif(-5, 10)", "ls_c[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#mean for lnorm", "mu_c[t] <- exp(lmu_c[t] + lss_c[t] / 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c[1])", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c[2])", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c[2])", model_string_jags, fixed = TRUE)
    if(length(model_c_random) == 0) {
    model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1])", "lmu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2])", "lmu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2])", model_string_jags, fixed = TRUE)
    }
    if(length(model_c_random) != 0) {
    model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "lmu_c[1] <- inprod(mean_cov_c1_fixed[], beta[, 1]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "lmu_c[2] <- inprod(mean_cov_c2_fixed[], beta[, 2]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2])", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e == "norm") {
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "beta") {
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dbeta(mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- (mu_e1[i] * (1 - mu_e1[i]) / pow(s_e[1], 2) - 1)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dbeta(mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- (mu_e2[i] * (1 - mu_e2[i]) / pow(s_e[2], 2) - 1)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] ~ dunif(0, sqrt(mu_e[t] * (1 - mu_e[t])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- ilogit(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- ilogit(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- ilogit(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- ilogit(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "gamma"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dgamma(mu_e1[i] * tau_e1[i], tau_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- mu_e1[i] / pow(s_e[1], 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dgamma(mu_e2[i] * tau_e2[i], tau_e2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- mu_e2[i] / pow(s_e[2], 2)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.gamma(eff1[i], mu_e1[i] * tau_e1[i], tau_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.gamma(eff2[i], mu_e2[i] * tau_e2[i], tau_e2[i])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "exp"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dexp(1 / mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dexp(1 / mu_e2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] <- mu_e[t]", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.exp(eff1[i], 1/mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.exp(eff2[i], 1/mu_e2[i])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "weibull"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dweib(tau_e1[i], pow(1 / (mu_e1[i] / exp(loggam(1 + 1/tau_e1[i]))), tau_e1[i]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- pow(s_e[1] / mu_e1[i], - 1.086)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dweib(tau_e2[i], pow(1 / (mu_e2[i] / exp(loggam(1 + 1/tau_e2[i]))), tau_e2[i]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- pow(s_e[2] / mu_e2[i], - 1.086)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.weib(eff1[i], tau_e1[i], pow(1 / (mu_e1[i] / exp(loggam(1 + 1/tau_e1[i]))), tau_e1[i]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.weib(eff2[i], tau_e2[i], pow(1 / (mu_e2[i] / exp(loggam(1 + 1/tau_e2[i]))), tau_e2[i]))", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "logis"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dlogis(mu_e1[i], tau_e[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dlogis(mu_e2[i], tau_e[2])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "tau_e[t] <- 1 / sqrt((3 * ss_e[t]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.logis(eff1[i], mu_e1[i], tau_e[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.logis(eff2[i], mu_e2[i], tau_e[2])", model_string_jags, fixed = TRUE)
  } else if(dist_e == "bern"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dbern(mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dbern(mu_e2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] <- sqrt(mu_e[t] * (1 - mu_e[t]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.bern(eff1[i], mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.bern(eff2[i], mu_e2[i])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- ilogit(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- ilogit(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- ilogit(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- ilogit(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "pois"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dpois(mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dpois(mu_e2[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] <- sqrt(mu_e[t])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.pois(eff1[i], mu_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.pois(eff2[i], mu_e2[i])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  } else if(dist_e == "nbinom"){
  model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e[1])", "eff1[i] ~ dnegbin(tau_e[1] / (tau_e[1] + mu_e1[i]), tau_e[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e[2])", "eff2[i] ~ dnegbin(tau_e[2] / (tau_e[2] + mu_e2[i]), tau_e[2])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("tau_e[t] <- 1 / ss_e[t]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ss_e[t] <- s_e[t] * s_e[t]", "s_e[t] <- sqrt((tau_e[t] / (tau_e[t] + mu_e[t])) * tau_e[t]) / (1 - (tau_e[t] / (tau_e[t] + mu_e[t])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("s_e[t] <- exp(ls_e[t])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "tau_e[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.negbin(eff1[i], tau_e[1] / (tau_e[1] + mu_e1[i]), tau_e[1])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.negbin(eff2[i], tau_e[2] / (tau_e[2] + mu_e2[i]), tau_e[2])", model_string_jags, fixed = TRUE)
    if(length(model_e_random) == 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(length(model_e_random) != 0) {
    model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mu_e[1] <- exp(inprod(mean_cov_e1_fixed[], alpha[, 1]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mu_e[2] <- exp(inprod(mean_cov_e2_fixed[], alpha[, 2]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
  }
  if(dist_e %in% c("beta", "gamma", "weibull", "bern", "pois", "exp")) {
     if(dist_c == "gamma"){
     model_string_jags <- gsub("for (t in 1:2) {#begin transformation", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}# end transformation", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
   }
   if(pe_fixed == 1) {
   inprod_e1 <- "X1_e_fixed[i] * alpha[1]"
   inprod_e2 <- "X2_e_fixed[i] * alpha[2]"
   inprod_mean_e1 <- "mean_cov_e1_fixed * alpha[1]"
   inprod_mean_e2 <- "mean_cov_e2_fixed * alpha[2]"
   begin_prior_beta <- "#begin alpha priors effects"
   prior_beta <- "#"
   end_prior_beta <- "#end alpha priors effects"
   prior_beta_e1 <- "alpha[1] ~ dnorm(0, 0.0000001)"
   prior_beta_e2 <- "alpha[2] ~ dnorm(0, 0.0000001)"
   model_string_jags <- gsub("inprod(X1_e_fixed[i, ], alpha[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_e_fixed[i, ], alpha[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_e1_fixed[], alpha[, 1])", inprod_mean_e1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_e2_fixed[], alpha[, 2])", inprod_mean_e2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:pe_fixed) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {alpha[j, t] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("alpha[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_e2, model_string_jags, fixed = TRUE)
   }
   if(length(model_e_random) != 0 & pe_random == 1) {
       model_string_jags <- gsub("inprod(X1_e_random[i, ], a1[, clus1_e[i]])", "X1_e_random[i] * a1[clus1_e[i]]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X2_e_random[i, ], a2[, clus2_e[i]])", "X2_e_random[i] * a2[clus2_e[i]]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "mean_cov_e1_random * mu_a_hat[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "mean_cov_e2_random * mu_a_hat[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pe_random) {#begin a1 priors effects", "#begin a1 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n1_clus_e) {a1[j, s] ~ dnorm(mu_a_hat[j, 1], tau_a_hat[j, 1]) }", "for(s in 1:n1_clus_e) {a1[s] ~ dnorm(mu_a_hat[1], tau_a_hat[1]) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a1 priors effects", "#end a1 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 1:pe_random) {#begin a2 priors effects", "#begin a2 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(s in 1:n2_clus_e) {a2[j, s] ~ dnorm(mu_a_hat[j, 2], tau_a_hat[j, 2]) }", "for(s in 1:n2_clus_e) {a2[s] ~ dnorm(mu_a_hat[2], tau_a_hat[2]) }", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end a2 priors effects", "#end a2 priors effects", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {mu_a_hat[j, t] ~ dnorm(0, 0.001)", "mu_a_hat[t] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_a_hat[j, t] ~ dunif(0, 100) }", "s_a_hat[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe_random) {tau_a_hat[j, t] <- 1 / ss_a_hat[j, t]", "tau_a_hat[t] <- 1 / ss_a_hat[t]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_a_hat[j, t] <- s_a_hat[j, t] * s_a_hat[j, t] }", "ss_a_hat[t] <- s_a_hat[t] * s_a_hat[t]", model_string_jags, fixed = TRUE)
       }
   if(pc_fixed == 1) {
   inprod_c1 <- "X1_c_fixed[i] * beta[1]"
   inprod_c2 <- "X2_c_fixed[i] * beta[2]"
   inprod_mean_c1 <- "mean_cov_c1_fixed * beta[1]"
   inprod_mean_c2 <- "mean_cov_c2_fixed * beta[2]"
   begin_prior_beta <- "#begin beta priors costs"
   prior_beta <- "#"
   end_prior_beta <- "#end beta priors costs"
   prior_beta_c1 <- "beta[1] ~ dnorm(0, 0.0000001)"
   prior_beta_c2 <- "beta[2] ~ dnorm(0, 0.0000001)"
   model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta[, 1])", inprod_mean_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta[, 2])", inprod_mean_c2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:pc_fixed) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {beta[j, t] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end beta priors costs", end_prior_beta,model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("beta[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
   }
   if(length(model_c_random) != 0 & pc_random == 1) {
   model_string_jags <- gsub("inprod(X1_c_random[i, ], b1[, clus1_c[i]])", "X1_c_random[i] * b1[clus1_c[i]]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(X2_c_random[i, ], b2[, clus2_c[i]])", "X2_c_random[i] * b2[clus2_c[i]]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "mean_cov_c1_random * mu_b_hat[1]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "mean_cov_c2_random * mu_b_hat[2]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b1 priors costs", "#begin b1 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n1_clus_c) {b1[j, s] ~ dnorm(mu_b_hat[j, 1], tau_b_hat[j, 1]) }", "for(s in 1:n1_clus_c) {b1[s] ~ dnorm(mu_b_hat[1], tau_b_hat[1]) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b1 priors costs", "#end b1 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 1:pc_random) {#begin b2 priors costs", "#begin b2 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(s in 1:n2_clus_c) {b2[j, s] ~ dnorm(mu_b_hat[j, 2], tau_b_hat[j, 2]) }", "for(s in 1:n2_clus_c) {b2[s] ~ dnorm(mu_b_hat[2], tau_b_hat[2]) }", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end b2 priors costs", "#end b2 priors costs", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j, t] ~ dnorm(0, 0.001)", "mu_b_hat[t] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("s_b_hat[j, t] ~ dunif(0, 100) }", "s_b_hat[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j, t] <- 1 / ss_b_hat[j, t]", "tau_b_hat[t] <- 1 / ss_b_hat[t]", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("ss_b_hat[j, t] <- s_b_hat[j, t] * s_b_hat[j, t] }", "ss_b_hat[t] <- s_b_hat[t] * s_b_hat[t]", model_string_jags, fixed = TRUE)
   }
   if(ze_fixed == 1) {
   inprod_e1 <- "Z1_e_fixed[i] * gamma_e[1]"
   inprod_e2 <- "Z2_e_fixed[i] * gamma_e[2]"
     if(type == "MAR" | type == "MNAR_cost") {
     inprod_mean_e1 <- "mean_z_e1_fixed * gamma_e[1]"
     inprod_mean_e2 <- "mean_z_e2_fixed * gamma_e[2]"
     model_string_jags <- gsub("inprod(mean_z_e1_fixed[], gamma_e[, 1])", inprod_mean_e1, model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("inprod(mean_z_e2_fixed[], gamma_e[, 2])", inprod_mean_e2, model_string_jags, fixed = TRUE)
     }
     if(type == "MNAR_eff" | type == "MNAR") {
     inprod_mean_e1 <- "mean_z_e1_fixed * gamma_e[1] + delta_e[1] * mean(eff1[])"
     inprod_mean_e2 <- "mean_z_e2_fixed * gamma_e[2] + delta_e[2] * mean(eff2[])"
     model_string_jags <- gsub("inprod(mean_z_e1_fixed[], gamma_e[, 1]) + delta_e[1] * mean(eff1[])", inprod_mean_e1, model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("inprod(mean_z_e2_fixed[], gamma_e[, 2]) + delta_e[2] * mean(eff2[])", inprod_mean_e2, model_string_jags, fixed = TRUE)
     }
   begin_prior_gamma <- "#begin gamma priors effects"
   begin_prior_gamma2 <- "#"
   prior_gamma_e1 <- "gamma_e[1] ~ dlogis(0, 1)"
   prior_gamma_e2 <- "gamma_e[2] ~ dlogis(0, 1)"
   end_prior_gamma <- "#end gamma priors effects"
   model_string_jags <- gsub("inprod(Z1_e_fixed[i, ], gamma_e[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("inprod(Z2_e_fixed[i, ], gamma_e[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for (j in 2:ze_fixed) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", prior_gamma_e1, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", prior_gamma_e2, model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
   } 
   if(ze_random == 1 & length(model_me_random) != 0) {
      if(type == "MAR" | type == "MNAR_cost") {
      model_string_jags <- gsub("inprod(mean_z_e1_random[], mu_g_e_hat[, 1])", "mean_z_e1_random * mu_g_e_hat[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e2_random[], mu_g_e_hat[, 2])", "mean_z_e2_random * mu_g_e_hat[2]", model_string_jags, fixed = TRUE)
      }
      if(type == "MNAR_eff" | type == "MNAR") {
      model_string_jags <- gsub("inprod(mean_z_e1_random[], mu_g_e_hat[, 1]) + mu_d_e_hat[1] * mean(eff1[])", "mean_z_e1_random * mu_g_e_hat[1] + mu_d_e_hat[1] * mean(eff1[])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e2_random[], mu_g_e_hat[, 2]) + mu_d_e_hat[2] * mean(eff2[])", "mean_z_e2_random * mu_g_e_hat[2] + mu_d_e_hat[2] * mean(eff2[])", model_string_jags, fixed = TRUE)
      }
      model_string_jags <- gsub("inprod(Z1_e_random[i, ], g1_e[, clus1_me[i]])", "Z1_e_random[i] * g1_e[clus1_me[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_e_random[i, ], g2_e[, clus2_me[i]])", "Z2_e_random[i] * g2_e[clus2_me[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {tau_g_e_hat[j, t] <- 1 / ss_g_e_hat[j, t]", "tau_g_e_hat[t] <- 1 / ss_g_e_hat[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_e_hat[j, t] <- s_g_e_hat[j, t] * s_g_e_hat[j, t] }", "ss_g_e_hat[t] <- s_g_e_hat[t] * s_g_e_hat[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e1_random[], mu_g_e_hat[, 1])", "mean_z_e1_random * mu_g_e_hat[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_e2_random[], mu_g_e_hat[, 2])", "mean_z_e2_random * mu_g_e_hat[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:ze_random) {#begin g1_e priors effects", "#begin g1_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_me) {g1_e[j, s] ~ dnorm(mu_g_e_hat[j, 1], tau_g_e_hat[j, 1]) }", "for(s in 1:n1_clus_me) {g1_e[s] ~ dnorm(mu_g_e_hat[1], tau_g_e_hat[1]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g1_e priors effects", "#end g1_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:ze_random) {#begin g2_e priors effects", "#begin g2_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_me) {g2_e[j, s] ~ dnorm(mu_g_e_hat[j, 2], tau_g_e_hat[j, 2]) }", "for(s in 1:n2_clus_me) {g2_e[s] ~ dnorm(mu_g_e_hat[2], tau_g_e_hat[2]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g2_e priors effects", "#end g2_e priors effects", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:ze_random) {mu_g_e_hat[j, t] ~ dnorm(0, 0.001)", "mu_g_e_hat[t] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_e_hat[j, t] ~ dunif(0, 100) }", "s_g_e_hat[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
   if(zc_fixed == 1) {
      inprod_c1 <- "Z1_c_fixed[i] * gamma_c[1]"
      inprod_c2 <- "Z2_c_fixed[i] * gamma_c[2]"
          if(type == "MAR" | type == "MNAR_eff") {
          inprod_mean_c1 <- "mean_z_c1_fixed * gamma_c[1]"
          inprod_mean_c2 <- "mean_z_c2_fixed * gamma_c[2]"
          model_string_jags <- gsub("inprod(mean_z_c1_fixed[], gamma_c[, 1])", inprod_mean_c1, model_string_jags, fixed = TRUE)
          model_string_jags <- gsub("inprod(mean_z_c2_fixed[], gamma_c[, 2])", inprod_mean_c2, model_string_jags, fixed = TRUE)
          }
          if(type == "MNAR_cost" | type == "MNAR") {
          inprod_mean_c1 <- "mean_z_c1_fixed * gamma_c[1] + delta_c[1] * mean(cost1[])"
          inprod_mean_c2 <- "mean_z_c2_fixed * gamma_c[2] + delta_c[2] * mean(cost2[])"
          model_string_jags <- gsub("inprod(mean_z_c1_fixed[], gamma_c[, 1]) + delta_c[1] * mean(cost1[])", inprod_mean_c1, model_string_jags, fixed = TRUE)
          model_string_jags <- gsub("inprod(mean_z_c2_fixed[], gamma_c[, 2]) + delta_c[2] * mean(cost2[])", inprod_mean_c2, model_string_jags, fixed = TRUE)
          }
      begin_prior_gamma <- "#begin gamma priors costs"
      begin_prior_gamma2 <- "#"
      prior_gamma_c1 <- "gamma_c[1] ~ dlogis(0, 1)"
      prior_gamma_c2 <- "gamma_c[2] ~ dlogis(0, 1)"
      end_prior_gamma <- "#end gamma priors costs"
      model_string_jags <- gsub("inprod(Z1_c_fixed[i, ], gamma_c[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_fixed[i, ], gamma_c[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 2:zc_fixed) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
      }
      if(zc_random == 1 & length(model_mc_random) != 0) {
         if(type == "MAR" | type == "MNAR_eff") {
         model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1])", "mean_z_c1_random * mu_g_c_hat[1]", model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2])", "mean_z_c2_random * mu_g_c_hat[2]", model_string_jags, fixed = TRUE)
         }
         if(type == "MNAR_cost" | type == "MNAR") {
         model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1]) + mu_d_c_hat[1] * mean(cost1[])", "mean_z_c1_random * mu_g_c_hat[1] + mu_d_c_hat[1] * mean(cost1[])", model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2]) + mu_d_c_hat[2] * mean(cost2[])", "mean_z_c2_random * mu_g_c_hat[2] + mu_d_c_hat[2] * mean(cost2[])", model_string_jags, fixed = TRUE)
         }
      model_string_jags <- gsub("inprod(Z1_c_random[i, ], g1_c[, clus1_mc[i]])", "Z1_c_random[i] * g1_c[clus1_mc[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(Z2_c_random[i, ], g2_c[, clus2_mc[i]])", "Z2_c_random[i] * g2_c[clus2_mc[i]]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {tau_g_c_hat[j, t] <- 1 / ss_g_c_hat[j, t]", "tau_g_c_hat[t] <- 1 / ss_g_c_hat[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_g_c_hat[j, t] <- s_g_c_hat[j, t] * s_g_c_hat[j, t] }", "ss_g_c_hat[t] <- s_g_c_hat[t] * s_g_c_hat[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c1_random[], mu_g_c_hat[, 1])", "mean_z_c1_random * mu_g_c_hat[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("inprod(mean_z_c2_random[], mu_g_c_hat[, 2])", "mean_z_c2_random * mu_g_c_hat[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g1_c priors costs", "#begin g1_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n1_clus_mc) {g1_c[j, s] ~ dnorm(mu_g_c_hat[j, 1], tau_g_c_hat[j, 1]) }", "for(s in 1:n1_clus_mc) {g1_c[s] ~ dnorm(mu_g_c_hat[1], tau_g_c_hat[1]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g1_c priors costs", "#end g1_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (j in 1:zc_random) {#begin g2_c priors costs", "#begin g2_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(s in 1:n2_clus_mc) {g2_c[j, s] ~ dnorm(mu_g_c_hat[j, 2], tau_g_c_hat[j, 2]) }", "for(s in 1:n2_clus_mc) {g2_c[s] ~ dnorm(mu_g_c_hat[2], tau_g_c_hat[2]) }", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end g2_c priors costs", "#end g2_c priors costs", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for(j in 1:zc_random) {mu_g_c_hat[j, t] ~ dnorm(0, 0.001)", "mu_g_c_hat[t] ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_g_c_hat[j, t] ~ dunif(0, 100) }", "s_g_c_hat[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      } 
      model_string_jags <- prior_selection(type = type, dist_e = dist_e, dist_c = dist_c, pe_fixed = pe_fixed, pc_fixed = pc_fixed, ze_fixed = ze_fixed, zc_fixed = zc_fixed, 
                                           model_e_random = model_e_random, model_c_random = model_c_random, model_me_random = model_me_random, model_mc_random = model_mc_random,
                                           pe_random = pe_random, pc_random = pc_random, ze_random = ze_random, zc_random = zc_random)
   writeLines(model_string_jags, "selection.txt")
   model_string <- "selection.txt"
   return(model_string)
   }))