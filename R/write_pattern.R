#'An internal function to select which type of pattern mixture model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Pattern mixture models
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param ind_fixed Logical; if TRUE independence between effectiveness and costs is assumed, else correlation is accounted for
#' @param ind_random Logical; if TRUE independence at the level of the random effects between effectiveness and costs is assumed, else correlation is accounted for
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param pe_fixed Number of fixed effects for the effectiveness model
#' @param pc_fixed Number of fixed effects for the cost model
#' @param pe_random Number of random effects for the effectiveness model
#' @param pc_random Number of random effects for the cost model
#' @param model_e_random Random effects formula for the effectiveness model
#' @param model_c_random Random effects formula for the costs model
#' @param d_list Number and type of patterns 
#' @param d1 Pattern indicator in the control 
#' @param d2 Pattern indicator in the intervention
#' @param restriction type of identifying restriction to be imposed
#' @examples
#' # Internal function only
#' # No examples
#' #
#' #


write_pattern <- function(type, dist_e, dist_c, pe_fixed, pc_fixed, ind_fixed, pe_random, pc_random, ind_random, model_e_random, model_c_random, 
                          d_list, d1, d2, restriction) eval.parent(substitute( {
  model_string_jags<-  "
  model {
  
  #control
  for(i in 1:N1) {
  #costs and effects model
  cost1[i] ~ dnorm(mu_c1[i], tau_c_p1[d1[i]])
  eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])
  
  #derive mean and std effects1 
  #derive mean and std costs1
  
  #mean regression in each pattern
  mu_c1[i] <- inprod(X1_c_fixed[i, ], beta_p1[, d1[i]]) + beta_f_p1[d1[i]] * (eff1[i] - meane_p1[d1[i]]) + inprod(X1_c_random[i, ], b1[, clus1_c[i]]) + b1_f[clus1_c[i]] * (eff1[i] - mu_e[1])
  mu_e1[i] <- inprod(X1_e_fixed[i, ], alpha_p1[, d1[i]]) + inprod(X1_e_random[i, ], a1[, clus1_e[i]])

  #patterns model for t=1
  d1[i] ~ dcat(p_prob1[1:n_patterns1])
  
  #loglikelihood
  loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])
  loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c_p1[d1[i]])
  loglik_d1[i] <- logdensity.cat(d1[i], p_prob1[1:n_patterns1])
  }
  
  #intervention
  for(i in 1:N2) {
  #costs and effects model
  cost2[i] ~ dnorm(mu_c2[i], tau_c_p2[d2[i]])
  eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])
  
  #derive mean and std effects2 
  #derive mean and std costs2
  
  #mean regression in each pattern
  mu_c2[i] <- inprod(X2_c_fixed[i, ], beta_p2[, d2[i]]) + beta_f_p2[d2[i]] * (eff2[i] - meane_p2[d2[i]]) + inprod(X2_c_random[i, ], b2[, clus2_c[i]]) + b2_f[clus2_c[i]] * (eff2[i] - mu_e[2])
  mu_e2[i] <- inprod(X2_e_fixed[i, ], alpha_p2[, d2[i]]) + inprod(X2_e_random[i, ], a2[, clus2_e[i]])
  
  #patterns model for t=2
  d2[i] ~ dcat(p_prob2[1:n_patterns2])
  
  #loglikelihood
  loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])
  loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c_p2[d2[i]])
  loglik_d2[i] <- logdensity.cat(d2[i], p_prob2[1:n_patterns2])
  }
  
  #transformation of parameters
  for (d in 1:n_patterns1) {#begin transformation for t=1
  tau_c_p1[d] <- 1 / ss_c_p1[d]
  ss_c_p1[d] <- s_c_p1[d] * s_c_p1[d]
  s_c_p1[d] <- exp(ls_c_p1[d])
  #mean for lnorm for t=1
  tau_e_p1[d] <- 1 / ss_e_p1[d]
  ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]
  s_e_p1[d] <- exp(ls_e_p1[d])
  }# end transformation for t=1

  for (d in 1:n_patterns2) {#begin transformation for t=2
  tau_c_p2[d] <- 1 / ss_c_p2[d]
  ss_c_p2[d] <- s_c_p2[d] * s_c_p2[d]
  s_c_p2[d] <- exp(ls_c_p2[d])
  #mean for lnorm for t=2
  tau_e_p2[d] <- 1 / ss_e_p2[d]
  ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]
  s_e_p2[d] <- exp(ls_e_p2[d])
  }# end transformation for t=2
  
  #transformation of random effects parameters
  for (t in 1:2) {# begin transformation random effects
  for(j in 1:pc_random) {tau_b_hat[j, t] <- 1 / ss_b_hat[j, t]
  ss_b_hat[j, t] <- s_b_hat[j, t] * s_b_hat[j, t] }
  for(j in 1:pe_random) {tau_a_hat[j, t] <- 1 / ss_a_hat[j, t]
  ss_a_hat[j, t] <- s_a_hat[j, t] * s_a_hat[j, t] }
  }#end transformation of random effects 
  
  #calculate means at mean of covariates in each pattern
  for(d in 1:n_patterns1){
  meanc_p1[d] <- inprod(mean_cov_c1_fixed[], beta_p1[, d])
  meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])
  }

  for(d in 1:n_patterns2){
  meanc_p2[d] <- inprod(mean_cov_c2_fixed[], beta_p2[, d])
  meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])
  }

  #add sensitivity parameters to means  
  mu_c_p1[1] <- meanc_p1[1]
  mu_c_p1[2] <- meanc_p1[2]
  mu_c_p1[3] <- meanc_p1[3] + Delta_c[1]
  mu_c_p1[4] <- meanc_p1[4] + Delta_c[1]
  mu_e_p1[1] <- meane_p1[1]
  mu_e_p1[2] <- meane_p1[2] + Delta_e[1]
  mu_e_p1[3] <- meane_p1[3] 
  mu_e_p1[4] <- meane_p1[4] + Delta_e[1]
  
  mu_c_p2[1] <- meanc_p2[1]
  mu_c_p2[2] <- meanc_p2[2]
  mu_c_p2[3] <- meanc_p2[3] + Delta_c[2]
  mu_c_p2[4] <- meanc_p2[4] + Delta_c[2]
  mu_e_p2[1] <- meane_p2[1]
  mu_e_p2[2] <- meane_p2[2] + Delta_e[2]
  mu_e_p2[3] <- meane_p2[3] 
  mu_e_p2[4] <- meane_p2[4] + Delta_e[2]

  #calculate overall means
  mu_c[1] <- sum(mu_c_p1[] * p_prob1[]) + inprod(mean_cov_c1_random[], mu_b_hat[, 1]) 
  mu_c[2] <- sum(mu_c_p2[] * p_prob2[]) + inprod(mean_cov_c2_random[], mu_b_hat[, 2])
  mu_e[1] <- sum(mu_e_p1[] * p_prob1[]) + inprod(mean_cov_e1_random[], mu_a_hat[, 1])
  mu_e[2] <- sum(mu_e_p2[] * p_prob2[]) + inprod(mean_cov_e2_random[], mu_a_hat[, 2])

  #priors
  
  #priors for mean regression coefficients
  for (j in 2:pe_fixed) {#begin alpha priors effects in each pattern
  alpha_p1[j, 1] ~ dnorm(0, 0.0000001)
  alpha_p1[j, 2] <- alpha_p1[j, 1]
  alpha_p1[j, 3] ~ dnorm(0, 0.0000001)
  alpha_p1[j, 4] <- alpha_p1[j, 1]
  alpha_p2[j, 1] ~ dnorm(0, 0.0000001)
  alpha_p2[j, 2] <- alpha_p2[j, 1]
  alpha_p2[j, 3] ~ dnorm(0, 0.0000001)
  alpha_p2[j, 4] <- alpha_p2[j, 1]
  }#end alpha priors effects
  alpha_p1[1, 1] ~ dnorm(0, 0.0000001)
  alpha_p1[1, 2] <- alpha_p1[1, 1]
  alpha_p1[1, 3] ~ dnorm(0, 0.0000001)
  alpha_p1[1, 4] <- alpha_p1[1, 1]
  
  alpha_p2[1, 1] ~ dnorm(0, 0.0000001)
  alpha_p2[1, 2] <- alpha_p2[1, 1]
  alpha_p2[1, 3] ~ dnorm(0, 0.0000001)
  alpha_p2[1, 4] <- alpha_p2[1, 1]

  for (j in 2:pc_fixed) {#begin beta priors costs
  beta_p1[j, 1] ~ dnorm(0, 0.0000001)
  beta_p1[j, 2] ~ dnorm(0, 0.0000001)
  beta_p1[j, 3] <- beta_p1[j, 1]
  beta_p1[j, 4] <- beta_p1[j, 1]
  beta_p2[j, 1] ~ dnorm(0, 0.0000001)
  beta_p2[j, 2] ~ dnorm(0, 0.0000001)
  beta_p2[j, 3] <- beta_p2[j, 1]
  beta_p2[j, 4] <- beta_p2[j, 1]
  }#end beta priors costs
  beta_p1[1, 1] ~ dnorm(0, 0.0000001)
  beta_p1[1, 2] ~ dnorm(0, 0.0000001)
  beta_p1[1, 3] <- beta_p1[1, 1]
  beta_p1[1, 4] <- beta_p1[1, 1]
  
  beta_p2[1, 1] ~ dnorm(0, 0.0000001)
  beta_p2[1, 2] ~ dnorm(0, 0.0000001)
  beta_p2[1, 3] <- beta_p2[1, 1]
  beta_p2[1, 4] <- beta_p2[1, 1]
  
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
  ls_c_p1[1] ~ dunif(-5, 10)
  ls_c_p1[2] ~ dunif(-5, 10)
  ls_c_p1[3] <- ls_c_p1[1]
  ls_c_p1[4] <- ls_c_p1[1]
  
  ls_e_p1[1] ~ dunif(-5, 10)
  ls_e_p1[2] <- ls_e_p1[1]
  ls_e_p1[3] ~ dunif(-5, 10)
  ls_e_p1[4] <- ls_e_p1[1]
  
  ls_c_p2[1] ~ dunif(-5, 10)
  ls_c_p2[2] ~ dunif(-5, 10)
  ls_c_p2[3] <- ls_c_p2[1]
  ls_c_p2[4] <- ls_c_p2[1]
  
  ls_e_p2[1] ~ dunif(-5, 10)
  ls_e_p2[2] <- ls_e_p2[1]
  ls_e_p2[3] ~ dunif(-5, 10)
  ls_e_p2[4] <- ls_e_p2[1]

  #correlation
  beta_f_p1[1] ~ dnorm(0, 0.0000001) 
  beta_f_p1[2] <- beta_f_p1[1]
  beta_f_p1[3] <- beta_f_p1[1] 
  beta_f_p1[4] <- beta_f_p1[1]

  beta_f_p2[1] ~ dnorm(0, 0.0000001)
  beta_f_p2[2] <- beta_f_p2[1]
  beta_f_p2[3] <- beta_f_p2[1]
  beta_f_p2[4] <- beta_f_p2[1]

  for(t in 1:2) {# mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j, t] ~ dnorm(0, 0.001)
  s_b_hat[j, t] ~ dunif(0, 100) }
  for(j in 1:pe_random) {mu_a_hat[j, t] ~ dnorm(0, 0.001)
  s_a_hat[j, t] ~ dunif(0, 100) }
  }# end mean and sd mean regression random coefficients priors

  # correlation random effects
  for(s in 1:n1_clus_c) {b1_f[s] ~ dnorm(mu_b_f_hat[1], tau_b_f_hat[1]) }
  for(s in 1:n2_clus_c) {b2_f[s] ~ dnorm(mu_b_f_hat[2], tau_b_f_hat[2]) }
  for(t in 1:2) {mu_b_f_hat[t] ~ dnorm(0, 0.001)
  tau_b_f_hat[t] <- 1 / ss_b_f_hat[t]
  ss_b_f_hat[t] <- s_b_f_hat[t] * s_b_f_hat[t]
  s_b_f_hat[t] ~ dunif(0, 100) }

  #priors on patterns model
  p_prob1 ~ ddirch(pp1[])
  p_prob2 ~ ddirch(pp2[])

  for (t in 1:2){#begin priors on sensitivity parameters
  Delta_c[t] ~ dunif(range_c[t, 1], range_c[t, 2])  
  Delta_e[t] ~ dunif(range_e[t, 1], range_e[t, 2]) 
  }#end priors on sensitivity parameters

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
  if(length(model_c_random) == 0 & length(model_e_random) == 0 | length(model_e_random) == 0 & pc_random == 0) { 
    model_string_jags <- gsub("for(t in 1:2) {# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("}# end mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("for (t in 1:2) {# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("}#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
  }
  if(type == "MAR") {
    model_string_jags <- gsub(" + Delta_e[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + Delta_e[2]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("Delta_e[t] ~ dunif(range_e[t, 1], range_e[t, 2])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + Delta_c[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + Delta_c[2]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("Delta_c[t] ~ dunif(range_c[t, 1], range_c[t, 2])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (t in 1:2){#begin priors on sensitivity parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end priors on sensitivity parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#add sensitivity parameters to means", "", model_string_jags, fixed = TRUE)
  } else if(type == "MNAR_eff") {
    model_string_jags <- gsub(" + Delta_c[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + Delta_c[2]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("Delta_c[t] ~ dunif(range_c[t, 1], range_c[t, 2])", "", model_string_jags, fixed = TRUE)
  } else if(type == "MNAR_cost") {
    model_string_jags <- gsub(" + Delta_e[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + Delta_e[2]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("Delta_e[t] ~ dunif(range_e[t, 1], range_e[t, 2])", "", model_string_jags, fixed = TRUE)
  }  
  if(restriction == "AC") {
    model_string_jags <- gsub("alpha_p1[j, 2] <- alpha_p1[j, 1]", "alpha_p1[j, 2] <- alpha_p1[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[j, 4] <- alpha_p1[j, 1]", "alpha_p1[j, 4] <- alpha_p1[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 2] <- alpha_p2[j, 1]", "alpha_p2[j, 2] <- alpha_p2[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 4] <- alpha_p2[j, 1]", "alpha_p2[j, 4] <- alpha_p2[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 1]", "alpha_p1[1, 2] <- alpha_p1[1, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 4] <- alpha_p1[1, 1]", "alpha_p1[1, 4] <- alpha_p1[1, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 1]", "alpha_p2[1, 2] <- alpha_p2[1, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 4] <- alpha_p2[1, 1]", "alpha_p2[1, 4] <- alpha_p2[1, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 3] <- beta_p1[j, 1]", "beta_p1[j, 3] <- beta_p1[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 4] <- beta_p1[j, 1]", "beta_p1[j, 4] <- beta_p1[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 3] <- beta_p2[j, 1]", "beta_p2[j, 3] <- beta_p2[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 4] <- beta_p2[j, 1]", "beta_p2[j, 4] <- beta_p2[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 3] <- beta_p1[1, 1]", "beta_p1[1, 3] <- beta_p1[1, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 4] <- beta_p1[1, 1]", "beta_p1[1, 4] <- beta_p1[1, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 3] <- beta_p2[1, 1]", "beta_p2[1, 3] <- beta_p2[1, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 4] <- beta_p2[1, 1]", "beta_p2[1, 4] <- beta_p2[1, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[3] <- ls_c_p1[1]", "ls_c_p1[3] <- ls_c_p1[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[4] <- ls_c_p1[1]", "ls_c_p1[4] <- ls_c_p1[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "ls_e_p1[2] <- ls_e_p1[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "ls_e_p1[4] <- ls_e_p1[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[3] <- ls_c_p2[1]", "ls_c_p2[3] <- ls_c_p2[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[4] <- ls_c_p2[1]", "ls_c_p2[4] <- ls_c_p2[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "ls_e_p2[2] <- ls_e_p2[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "ls_e_p2[4] <- ls_e_p2[3]", model_string_jags, fixed = TRUE)
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
    model_string_jags <- gsub(" + beta_f_p1[d1[i]] * (eff1[i] - meane_p1[d1[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + beta_f_p2[d2[i]] * (eff2[i] - meane_p2[d2[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p1[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p1[2] <- beta_f_p1[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p1[3] <- beta_f_p1[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p1[4] <- beta_f_p1[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p2[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p2[2] <- beta_f_p2[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p2[3] <- beta_f_p2[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p2[4] <- beta_f_p2[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
  if(dist_c == "norm") {
    model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm for t=1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm for t=2", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
    model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c_p1[d1[i]])", "cost1[i] ~ dgamma(mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs1", "tau_c1[i] <- mu_c1[i] / pow(s_c_p1[d1[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c1[i] <- ", "log(mu_c1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c_p2[d2[i]])", "cost2[i] ~ dgamma(mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs2", "tau_c2[i] <- mu_c2[i] / pow(s_c_p2[d2[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c2[i] <- ", "log(mu_c2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p1[d] <- 1 / ss_c_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p1[d] <- s_c_p1[d] * s_c_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p1[d] <- exp(ls_c_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p2[d] <- 1 / ss_c_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p2[d] <- s_c_p2[d] * s_c_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p2[d] <- exp(ls_c_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p1[d] <- inprod(mean_cov_c1_fixed[], beta_p1[, d])", "meanc_p1[d] <- exp(inprod(mean_cov_c1_fixed[], beta_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p2[d] <- inprod(mean_cov_c2_fixed[], beta_p2[, d])", "meanc_p2[d] <- exp(inprod(mean_cov_c2_fixed[], beta_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10)", "s_c_p1[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10)", "s_c_p1[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[3] <- ls_c_p1[1]", "s_c_p1[3] <- s_c_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[4] <- ls_c_p1[1]", "s_c_p1[4] <- s_c_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10)", "s_c_p2[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10)", "s_c_p2[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[3] <- ls_c_p2[1]", "s_c_p2[3] <- s_c_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[4] <- ls_c_p2[1]", "s_c_p2[4] <- s_c_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_c_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "+ exp(inprod(mean_cov_c1_random[], mu_b_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "+ exp(inprod(mean_cov_c2_random[], mu_b_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_c_p1[3] <- ls_c_p1[2]", "s_c_p1[3] <- s_c_p1[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c_p1[4] <- ls_c_p1[2]", "s_c_p1[4] <- s_c_p1[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c_p2[3] <- ls_c_p2[2]", "s_c_p2[3] <- s_c_p2[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c_p2[4] <- ls_c_p2[2]", "s_c_p2[4] <- s_c_p2[2]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("#mean for lnorm for t=1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm for t=2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c_p1[d1[i]])", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c_p2[d2[i]])", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
    
  } else if(dist_c == "lnorm") {
    model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c_p1[d1[i]])", "cost1[i] ~ dlnorm(lmu_c1[i], ltau_c_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c1[i] <- ", "lmu_c1[i] <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c_p2[d2[i]])", "cost2[i] ~ dlnorm(lmu_c2[i], ltau_c_p2[d2[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_c2[i] <- ", "lmu_c2[i] <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p1[d] <- 1 / ss_c_p1[d]", "ltau_c_p1[d] <- 1 / lss_c_p1[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p1[d] <- s_c_p1[d] * s_c_p1[d]", "lss_c_p1[d] <- ls_c_p1[d] * ls_c_p1[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p1[d] <- exp(ls_c_p1[d])", "s_c_p1[d] <- sqrt(exp(2 * lmeanc_p1[d] + lss_c_p1[d]) * (exp(lss_c_p1[d]) - 1))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p2[d] <- 1 / ss_c_p2[d]", "ltau_c_p2[d] <- 1 / lss_c_p2[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p2[d] <- s_c_p2[d] * s_c_p2[d]", "lss_c_p2[d] <- ls_c_p2[d] * ls_c_p2[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p2[d] <- exp(ls_c_p2[d])", "s_c_p2[d] <- sqrt(exp(2 * lmeanc_p2[d] + lss_c_p2[d]) * (exp(lss_c_p2[d]) - 1))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10)", "ls_c_p1[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10)", "ls_c_p1[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10)", "ls_c_p2[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10)", "ls_c_p2[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p1[d] <- inprod(mean_cov_c1_fixed[], beta_p1[, d])", "lmeanc_p1[d] <- inprod(mean_cov_c1_fixed[], beta_p1[, d])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p2[d] <- inprod(mean_cov_c2_fixed[], beta_p2[, d])", "lmeanc_p2[d] <- inprod(mean_cov_c2_fixed[], beta_p2[, d])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm for t=1", "meanc_p1[d] <- exp(lmeanc_p1[d] + lss_c_p1[d] / 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm for t=2", "meanc_p2[d] <- exp(lmeanc_p2[d] + lss_c_p2[d] / 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c_p1[d1[i]])", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c_p2[d2[i]])", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c_p2[d2[i]])", model_string_jags, fixed = TRUE)
     if(length(model_c_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_c1_random[], mu_b_hat[, 1])", "+ exp(inprod(mean_cov_c1_random[], mu_b_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_c2_random[], mu_b_hat[, 2])", "+ exp(inprod(mean_cov_c2_random[], mu_b_hat[, 2]))", model_string_jags, fixed = TRUE)
     }
  }
  if(dist_e == "norm") {
    model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "beta") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dbeta(mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- (mu_e1[i] * (1 - mu_e1[i]) / pow(s_e_p1[d1[i]], 2) - 1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dbeta(mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- (mu_e2[i] * (1 - mu_e2[i]) / pow(s_e_p2[d2[i]], 2) - 1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- ilogit(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- ilogit(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ ilogit(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ ilogit(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "gamma") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dgamma(mu_e1[i] * tau_e1[i], tau_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- mu_e1[i] / pow(s_e_p1[d1[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dgamma(mu_e2[i] * tau_e2[i], tau_e2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- mu_e2[i] / pow(s_e_p2[d2[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- exp(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- exp(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ exp(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ exp(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.gamma(eff1[i], mu_e1[i] * tau_e1[i], tau_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.gamma(eff2[i], mu_e2[i] * tau_e2[i], tau_e2[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "exp") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dexp(1 / mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dexp(1 / mu_e2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- exp(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- exp(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] <- mu_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] <- mu_e_p1[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] <- mu_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] <- mu_e_p2[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ exp(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ exp(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.exp(eff1[i], 1 / mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.exp(eff2[i], 1 / mu_e2[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "weibull") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dweib(tau_e1[i], pow(1 / (mu_e1[i] / exp(loggam(1 + 1/tau_e1[i]))), tau_e1[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- pow(s_e_p1[d1[i]] / mu_e1[i], - 1.086)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dweib(tau_e2[i], pow(1 / (mu_e2[i] / exp(loggam(1 + 1/tau_e2[i]))), tau_e2[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- pow(s_e_p2[d2[i]] / mu_e2[i], - 1.086)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- exp(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- exp(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ exp(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ exp(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.weib(eff1[i], tau_e1[i], pow(1 / (mu_e1[i] / exp(loggam(1 + 1/tau_e1[i]))), tau_e1[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.weib(eff2[i], tau_e2[i], pow(1 / (mu_e2[i] / exp(loggam(1 + 1/tau_e2[i]))), tau_e2[i]))", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "logis") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dlogis(mu_e1[i], tau_e_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dlogis(mu_e2[i], tau_e_p2[d2[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "tau_e_p1[d] <- 1 / sqrt((3 * ss_e_p1[d]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "tau_e_p2[d] <- 1 / sqrt((3 * ss_e_p2[d]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.logis(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.logis(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "bern") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dbern(mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dbern(mu_e2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- ilogit(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- ilogit(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] <- sqrt(mu_e_p1[1] * (1 - mu_e_p1[1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] <- sqrt(mu_e_p1[3] * (1 - mu_e_p1[3]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] <- sqrt(mu_e_p2[1] * (1 - mu_e_p2[1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] <- sqrt(mu_e_p2[3] * (1 - mu_e_p2[3]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ ilogit(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ ilogit(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.bern(eff1[i], mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.bern(eff2[i], mu_e2[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "pois") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dpois(mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dpois(mu_e2[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- exp(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- exp(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "s_e_p1[1] <- sqrt(mu_e_p1[1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "s_e_p1[3] <- sqrt(mu_e_p1[3])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "s_e_p1[4] <- s_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "s_e_p2[1] <- sqrt(mu_e_p2[1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "s_e_p2[3] <- sqrt(mu_e_p2[3])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "s_e_p2[4] <- s_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ exp(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ exp(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "s_e_p1[2] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "s_e_p1[4] <- s_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "s_e_p2[2] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "s_e_p2[4] <- s_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.pois(eff1[i], mu_e1[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.pois(eff2[i], mu_e2[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "nbinom") {
    model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e_p1[d1[i]])", "eff1[i] ~ dnegbin(tau_e_p1[d1[i]] / (tau_e_p1[d1[i]] + mu_e1[i]), tau_e_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e1[i] <- ", "log(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e_p2[d2[i]])", "eff2[i] ~ dnegbin(tau_e_p2[d2[i]] / (tau_e_p2[d2[i]] + mu_e2[i]), tau_e_p2[d2[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("mu_e2[i] <- ", "log(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p1[d] <- 1 / ss_e_p1[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p1[d] <- s_e_p1[d] * s_e_p1[d]", "s_e_p1[d] <- sqrt((tau_e_p1[d] / (tau_e_p1[d] + mu_e_p1[d])) * tau_e_p1[d]) / (1 - (tau_e_p1[d] / (tau_e_p1[d] + mu_e_p1[d])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p1[d] <- exp(ls_e_p1[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p2[d] <- 1 / ss_e_p2[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p2[d] <- s_e_p2[d] * s_e_p2[d]", "s_e_p2[d] <- sqrt((tau_e_p2[d] / (tau_e_p2[d] + mu_e_p2[d])) * tau_e_p2[d]) / (1 - (tau_e_p2[d] / (tau_e_p2[d] + mu_e_p2[d])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p2[d] <- exp(ls_e_p2[d])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p1[d] <- inprod(mean_cov_e1_fixed[], alpha_p1[, d])", "meane_p1[d] <- exp(inprod(mean_cov_e1_fixed[], alpha_p1[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p2[d] <- inprod(mean_cov_e2_fixed[], alpha_p2[, d])", "meane_p2[d] <- exp(inprod(mean_cov_e2_fixed[], alpha_p2[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "tau_e_p1[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "tau_e_p1[2] <- tau_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "tau_e_p1[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "tau_e_p1[4] <- tau_e_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "tau_e_p2[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "tau_e_p2[2] <- tau_e_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "tau_e_p2[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "tau_e_p2[4] <- tau_e_p2[1]", model_string_jags, fixed = TRUE)
    if(length(model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e1_random[], mu_a_hat[, 1])", "+ exp(inprod(mean_cov_e1_random[], mu_a_hat[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("+ inprod(mean_cov_e2_random[], mu_a_hat[, 2])", "+ exp(inprod(mean_cov_e2_random[], mu_a_hat[, 2]))", model_string_jags, fixed = TRUE)
    }
    if(restriction == "AC"){
      model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "tau_e_p1[2] <- tau_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "tau_e_p1[4] <- tau_e_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "tau_e_p2[2] <- tau_e_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "tau_e_p2[4] <- tau_e_p2[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e_p1[d1[i]])", "loglik_e1[i] <- logdensity.negbin(eff1[i], tau_e_p1[d1[i]] / (tau_e_p1[d1[i]] + mu_e1[i]), tau_e_p1[d1[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e_p2[d2[i]])", "loglik_e2[i] <- logdensity.negbin(eff2[i], tau_e_p2[d2[i]] / (tau_e_p2[d2[i]] + mu_e2[i]), tau_e_p2[d2[i]])", model_string_jags, fixed = TRUE)
  }
  if(dist_e %in% c("beta", "gamma", "weibull", "bern", "pois", "exp")) {
    if(dist_c == "gamma") {
    model_string_jags <- gsub("for (d in 1:n_patterns1) {#begin transformation for t=1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}# end transformation for t=1", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (d in 1:n_patterns2) {#begin transformation for t=2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}# end transformation for t=2", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
    }
  }
  if(pe_fixed == 1) {
    inprod_e1 <- "X1_e_fixed[i] * alpha_p1[d1[i]]"
    inprod_e2 <- "X2_e_fixed[i] * alpha_p2[d2[i]]"
    inprod_mean_e1 <- "mean_cov_e1_fixed * alpha_p1[d]"
    inprod_mean_e2 <- "mean_cov_e2_fixed * alpha_p2[d]"
    begin_prior_beta <- "#begin alpha priors effects"
    prior_beta <- "#"
    end_prior_beta <- "#end alpha priors effects"
    prior_beta_e1_1 <- "alpha_p1[1] ~ dnorm(0, 0.0000001)"
    prior_beta_e1_3 <- "alpha_p1[3] ~ dnorm(0, 0.0000001)"
    prior_beta_e2_1 <- "alpha_p2[1] ~ dnorm(0, 0.0000001)"
    prior_beta_e2_3 <- "alpha_p2[3] ~ dnorm(0, 0.0000001)"
    model_string_jags <- gsub("inprod(X1_e_fixed[i, ], alpha_p1[, d1[i]])", inprod_e1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(X2_e_fixed[i, ], alpha_p2[, d2[i]])", inprod_e2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_e1_fixed[], alpha_p1[, d])", inprod_mean_e1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_e2_fixed[], alpha_p2[, d])", inprod_mean_e2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 2:pe_fixed) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[j, 2] <- alpha_p1[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[j, 4] <- alpha_p1[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 2] <- alpha_p2[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[j, 4] <- alpha_p2[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1_1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001)", prior_beta_e1_3, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 1]", "alpha_p1[2] <- alpha_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p1[1, 4] <- alpha_p1[1, 1]", "alpha_p1[4] <- alpha_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e2_1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001)", prior_beta_e2_3, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 1]", "alpha_p2[2] <- alpha_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p2[1, 4] <- alpha_p2[1, 1]", "alpha_p2[4] <- alpha_p2[1]", model_string_jags, fixed = TRUE)
    if(restriction == "AC"){
      model_string_jags <- gsub("alpha_p1[j, 2] <- alpha_p1[j, 3]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p1[j, 4] <- alpha_p1[j, 3]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p2[j, 2] <- alpha_p2[j, 3]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p2[j, 4] <- alpha_p2[j, 3]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 3]", "alpha_p1[2] <- alpha_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p1[1, 4] <- alpha_p1[1, 3]", "alpha_p1[4] <- alpha_p1[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 3]", "alpha_p2[2] <- alpha_p2[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p2[1, 4] <- alpha_p2[1, 3]", "alpha_p2[4] <- alpha_p2[3]", model_string_jags, fixed = TRUE)
    }
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
    inprod_c1 <- "X1_c_fixed[i] * beta_p1[d1[i]]"
    inprod_c2 <- "X2_c_fixed[i] * beta_p2[d2[i]]"
    inprod_mean_c1 <- "mean_cov_c1_fixed * beta_p1[d]"
    inprod_mean_c2 <- "mean_cov_c2_fixed * beta_p2[d]"
    begin_prior_beta <- "#begin beta priors costs"
    prior_beta <- "#"
    end_prior_beta <- "#end beta priors costs"
    prior_beta_c1_1 <- "beta_p1[1] ~ dnorm(0, 0.0000001)"
    prior_beta_c1_2 <- "beta_p1[2] ~ dnorm(0, 0.0000001)"
    prior_beta_c2_1 <- "beta_p2[1] ~ dnorm(0, 0.0000001)"
    prior_beta_c2_2 <- "beta_p2[2] ~ dnorm(0, 0.0000001)"
    model_string_jags <- gsub("inprod(X1_c_fixed[i, ], beta_p1[, d1[i]])", inprod_c1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(X2_c_fixed[i, ], beta_p2[, d2[i]])", inprod_c2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_c1_fixed[], beta_p1[, d])", inprod_mean_c1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_c2_fixed[], beta_p2[, d])", inprod_mean_c2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 2:pc_fixed) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 3] <- beta_p1[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[j, 4] <- beta_p1[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001)", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 3] <- beta_p2[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[j, 4] <- beta_p2[j, 1]", prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end beta priors costs", end_prior_beta, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1_1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_c1_2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 3] <- beta_p1[1, 1]", "beta_p1[3] <- beta_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p1[1, 4] <- beta_p1[1, 1]", "beta_p1[4] <- beta_p1[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c2_1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_c2_2, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 3] <- beta_p2[1, 1]", "beta_p2[3] <- beta_p2[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p2[1, 4] <- beta_p2[1, 1]", "beta_p2[4] <- beta_p2[1]", model_string_jags, fixed = TRUE)
    if(restriction == "AC"){
      model_string_jags <- gsub("beta_p1[j, 3] <- beta_p1[j, 2]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p1[j, 4] <- beta_p1[j, 2]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p2[j, 3] <- beta_p2[j, 2]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p2[j, 4] <- beta_p2[j, 2]", prior_beta, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p1[1, 3] <- beta_p1[1, 2]", "beta_p1[3] <- beta_p1[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p1[1, 4] <- beta_p1[1, 2]", "beta_p1[4] <- beta_p1[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p2[1, 3] <- beta_p2[1, 2]", "beta_p2[3] <- beta_p2[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p2[1, 4] <- beta_p2[1, 2]", "beta_p2[4] <- beta_p2[2]", model_string_jags, fixed = TRUE)
    }
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
  if(restriction == "CC") {
   if(d_list$n_patterns[1] == 3 | d_list$n_patterns[1] == 2) {
    if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p1[4] <- meanc_p1[4]", "", model_string_jags, fixed = TRUE) }
    if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p1[4] <- meane_p1[4]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[4] <- meanc_p1[4] + Delta_c[1]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[4] <- meane_p1[4] + Delta_e[1]", "", model_string_jags, fixed = TRUE) }
    if(pe_fixed > 1) {
      model_string_jags <- gsub("alpha_p1[1, 4] <- alpha_p1[1, 1]", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p1[j, 4] <- alpha_p1[j, 1]", "", model_string_jags, fixed = TRUE)
      }
    if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[4] <- alpha_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(pc_fixed > 1) {
      model_string_jags <- gsub("beta_p1[1, 4] <- beta_p1[1, 1]", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("beta_p1[j, 4] <- beta_p1[j, 1]", "", model_string_jags, fixed = TRUE) 
      }
    if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[4] <- beta_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[4] <- ls_c_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[4] <- s_c_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e %in% c("beta", "gamma", "exp", "weibull", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p1[4] <- s_e_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[4] <- tau_e_p1[1]", "", model_string_jags, fixed = TRUE) }
    if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p1[4] <- beta_f_p1[1]" , "", model_string_jags, fixed = TRUE) }
   }
    if(d_list$n_patterns[1] == 3){
      if(d_list$d1$d1_ec_obs == TRUE & d_list$d1$d1_c_obs == TRUE & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == TRUE) {
        d1 <- ifelse(d1 == 4, 3, d1)
        if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[3] <- meane_p1[3]", "mu_e_p1[3] <- meane_p1[3] + Delta_e[1]", model_string_jags, fixed = TRUE) }
        if(pe_fixed > 1) {
          model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001)", "alpha_p1[1, 3] <- alpha_p1[1, 1]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001)", "alpha_p1[j, 3] <- alpha_p1[j, 1]", model_string_jags, fixed = TRUE) 
        }
        if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001)", "alpha_p1[3] <- alpha_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "ls_e_p1[3] <- ls_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3])))", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000)", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p1[3] <- mu_e_p1[3]", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100)", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3] * (1 - mu_e_p1[3]))", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3])", "s_e_p1[3] <- s_e_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100)", "tau_e_p1[3] <- tau_e_p1[1]", model_string_jags, fixed = TRUE) }
      } else if(d_list$d1$d1_ec_obs == TRUE & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == TRUE & d_list$d1$d1_ec_mis == TRUE) {
        d1 <- ifelse(d1 == 4, 2, d1)
        if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[2] <- meanc_p1[2]", "mu_c_p1[2] <- meanc_p1[2] + Delta_c[1]", model_string_jags, fixed = TRUE) }
        if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001)", "beta_p1[1, 2] <- beta_p1[1, 1]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001)", "beta_p1[j, 2] <- beta_p1[j, 1]", model_string_jags, fixed = TRUE) 
        }
        if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001)", "beta_p1[2] <- beta_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000)", "s_c_p1[2] <- s_c_p1[1]", model_string_jags, fixed = TRUE) }
      } 
    } 
    if(d_list$n_patterns[1] == 2) {
      if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p1[3] <- meanc_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p1[3] <- meane_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[3] <- meanc_p1[3] + Delta_c[1]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[3] <- meane_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      }
      if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p1[1, 3] <- beta_p1[1, 1]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p1[j, 3] <- beta_p1[j, 1]", "", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[3] <- beta_p1[1]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[3] <- ls_c_p1[1]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[3] <- s_c_p1[1]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "beta") {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3])))", "", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "exp") {model_string_jags <- gsub("s_e_p1[3] <- mu_e_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "bern") {model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3] * (1 - mu_e_p1[3]))", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "pois") {model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3])", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
      if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p1[3] <- beta_f_p1[1]" , "", model_string_jags, fixed = TRUE) }
      if(d_list$d1$d1_ec_obs == TRUE & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == TRUE & d_list$d1$d1_ec_mis == FALSE) {
        d1 <- ifelse(d1 == 3, 2, d1)
        if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[2] <- meanc_p1[2]", "mu_c_p1[2] <- meanc_p1[2] + Delta_c[1]", model_string_jags, fixed = TRUE) }
        if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001)", "beta_p1[1, 2] <- beta_p1[1, 1]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001)", "beta_p1[j, 2] <- beta_p1[j, 1]", model_string_jags, fixed = TRUE) 
        }
        if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001)", "beta_p1[2] <- beta_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000)", "s_c_p1[2] <- s_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(pe_fixed > 1) {
          model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 1]", "alpha_p1[1, 2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("alpha_p1[j, 2] <- alpha_p1[j, 1]", "alpha_p1[j, 2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) 
        }
        if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[2] <- alpha_p1[1]", "alpha_p1[2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) }
        if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[1]", "ls_e_p1[2] ~ dunif(-5, 10)", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] ~ dunif(0, sqrt(meane_p1[2] * (1 - meane_p1[2])))", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] <- mu_e_p1[2]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] <- sqrt(mu_e_p1[2] * (1 - mu_e_p1[2]))", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[1]", "s_e_p1[2] <- sqrt(mu_e_p1[2])", model_string_jags, fixed = TRUE) }
        if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[2] <- tau_e_p1[1]", "tau_e_p1[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }
      } else if(d_list$d1$d1_ec_obs == TRUE & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == TRUE) {
        d1 <- ifelse(d1 == 4, 2, d1)
        if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[2] <- meanc_p1[2]", "mu_c_p1[2] <- meanc_p1[2] + Delta_c[1]", model_string_jags, fixed = TRUE) }
        if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001)", "beta_p1[1, 2] <- beta_p1[1, 1]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001)", "beta_p1[j, 2] <- beta_p1[j, 1]", model_string_jags, fixed = TRUE) 
        }
        if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001)", "beta_p1[2] <- beta_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100)", "ls_c_p1[2] <- ls_c_p1[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000)", "s_c_p1[2] <- s_c_p1[1]", model_string_jags, fixed = TRUE) }
      } 
    }
  }
  if(restriction == "AC") {
    if(d_list$n_patterns[1] == 3 | d_list$n_patterns[1] == 2) {
      if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p1[4] <- meanc_p1[4]", "", model_string_jags, fixed = TRUE) }
      if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p1[4] <- meane_p1[4]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[4] <- meanc_p1[4] + Delta_c[1]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[4] <- meane_p1[4] + Delta_e[1]", "", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p1[1, 4] <- alpha_p1[1, 3]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[j, 4] <- alpha_p1[j, 3]", "", model_string_jags, fixed = TRUE)
      }
      if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[4] <- alpha_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p1[1, 4] <- beta_p1[1, 2]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p1[j, 4] <- beta_p1[j, 2]", "", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[4] <- beta_p1[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[4] <- ls_c_p1[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[4] <- s_c_p1[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[4] <- ls_e_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("beta", "gamma", "exp", "weibull", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p1[4] <- s_e_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[4] <- tau_e_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p1[4] <- beta_f_p1[1]" , "", model_string_jags, fixed = TRUE) }
    }
    if(d_list$n_patterns[1] == 3){
      if(d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_c_obs == TRUE & d_list$d1$d1_e_obs == TRUE & d_list$d1$d1_ec_mis == TRUE & ind_fixed == TRUE) {
        d1 <- ifelse(d1 == 4, 1, d1)
        if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[1] <- meane_p1[1]", "mu_e_p1[1] <- meane_p1[1] + Delta_e[1]", model_string_jags, fixed = TRUE) }
        if(pe_fixed > 1) {
          model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001)", "alpha_p1[1, 1] <- alpha_p1[1, 3]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001)", "alpha_p1[j, 1] <- alpha_p1[j, 3]", model_string_jags, fixed = TRUE) 
        }
        if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001)", "alpha_p1[1] <- alpha_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10)", "ls_e_p1[1] <- ls_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1])))", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000)", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p1[1] <- mu_e_p1[1]", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100)", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p1[1] <- sqrt(mu_e_p1[1] * (1 - mu_e_p1[1]))", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p1[1] <- sqrt(mu_e_p1[1])", "s_e_p1[1] <- s_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100)", "tau_e_p1[1] <- tau_e_p1[3]", model_string_jags, fixed = TRUE) }
        if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[1] <- meanc_p1[1]", "mu_c_p1[1] <- meanc_p1[1] + Delta_c[1]", model_string_jags, fixed = TRUE) }
        if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001)", "beta_p1[1, 1] <- beta_p1[1, 2]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001)", "beta_p1[j, 1] <- beta_p1[j, 2]", model_string_jags, fixed = TRUE) 
        }
        if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001)", "beta_p1[1] <- beta_p1[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10)", "ls_c_p1[1] <- ls_c_p1[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100)", "ls_c_p1[1] <- ls_c_p1[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000)", "s_c_p1[1] <- s_c_p1[2]", model_string_jags, fixed = TRUE) }
      }
    } 
    if(d_list$n_patterns[1] == 2) {
      d1 <- ifelse(d1 == 3, 1, d1)
      if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p1[3] <- meanc_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p1[3] <- meane_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[3] <- meanc_p1[3] + Delta_c[1]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p1[3] <- meane_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 3]", "alpha_p1[1, 2] <- alpha_p1[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[j, 2] <- alpha_p1[j, 3]", "alpha_p1[j, 2] <- alpha_p1[j, 1]", model_string_jags, fixed = TRUE) 
      }
      if(pe_fixed == 1) {
        model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p1[2] <- alpha_p1[3]", "alpha_p1[2] <- alpha_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p1[1, 3] <- beta_p1[1, 2]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p1[j, 3] <- beta_p1[j, 2]", "", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {
        model_string_jags <- gsub("beta_p1[3] <- beta_p1[2]", "", model_string_jags, fixed = TRUE) 
      }
      if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[3] <- ls_c_p1[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[3] <- s_c_p1[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {
        model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("ls_e_p1[2] <- ls_e_p1[3]", "ls_e_p1[2] <- ls_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3])))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "exp") {
        model_string_jags <- gsub("s_e_p1[3] <- mu_e_p1[3]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "weibull") {
        model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "bern") {
        model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3] * (1 - mu_e_p1[3]))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "pois") {
        model_string_jags <- gsub("s_e_p1[3] <- sqrt(mu_e_p1[3])", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p1[2] <- s_e_p1[3]", "s_e_p1[2] <- s_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "nbinom") {
        model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("tau_e_p1[2] <- tau_e_p1[3]", "tau_e_p1[2] <- tau_e_p1[1]", model_string_jags, fixed = TRUE) 
      }
      if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p1[3] <- beta_f_p1[1]" , "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p1[1] <- meanc_p1[1]", "mu_c_p1[1] <- meanc_p1[1] + Delta_c[1]", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001)", "beta_p1[1, 1] <- beta_p1[1, 2]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001)", "beta_p1[j, 1] <- beta_p1[j, 2]", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001)", "beta_p1[1] <- beta_p1[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10)", "ls_c_p1[1] <- ls_c_p1[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100)", "ls_c_p1[1] <- ls_c_p1[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000)", "s_c_p1[1] <- s_c_p1[2]", model_string_jags, fixed = TRUE) }
    }
  }
  if(restriction == "CC") {
  if(d_list$n_patterns[2] == 3 | d_list$n_patterns[2] == 2) {
    if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p2[4] <- meanc_p2[4]", "", model_string_jags, fixed = TRUE) }
    if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p2[4] <- meane_p2[4]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[4] <- meanc_p2[4] + Delta_c[2]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[4] <- meane_p2[4] + Delta_e[2]", "", model_string_jags, fixed = TRUE) }
    if(pe_fixed > 1) {
      model_string_jags <- gsub("alpha_p2[1, 4] <- alpha_p2[1, 1]", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p2[j, 4] <- alpha_p2[j, 1]", "", model_string_jags, fixed = TRUE) 
      }
    if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[4] <- alpha_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(pc_fixed > 1) {
      model_string_jags <- gsub("beta_p2[1, 4] <- beta_p2[1, 1]", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("beta_p2[j, 4] <- beta_p2[j, 1]", "", model_string_jags, fixed = TRUE) 
      }
    if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[4] <- beta_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[4] <- ls_c_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[4] <- s_c_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e %in% c("beta", "gamma", "exp", "weibull", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p2[4] <- s_e_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[4] <- tau_e_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p2[4] <- beta_f_p2[1]" , "", model_string_jags, fixed = TRUE) }
  }
  if(d_list$n_patterns[2] == 3){
    if(d_list$d2$d2_ec_obs == TRUE & d_list$d2$d2_c_obs == TRUE & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == TRUE) {
      d2 <- ifelse(d2 == 4, 3, d2)
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[3] <- meane_p2[3]", "mu_e_p2[3] <- meane_p2[3] + Delta_e[2]", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001)", "alpha_p2[1, 3] <- alpha_p2[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001)", "alpha_p2[j, 3] <- alpha_p2[j, 1]", model_string_jags, fixed = TRUE) 
        }
      if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001)", "alpha_p2[3] <- alpha_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "ls_e_p2[3] <- ls_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "beta") {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3])))", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000)", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "exp") {model_string_jags <- gsub("s_e_p2[3] <- mu_e_p2[3]", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100)", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "bern") {model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3] * (1 - mu_e_p2[3]))", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "pois") {model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3])", "s_e_p2[3] <- s_e_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100)", "tau_e_p2[3] <- tau_e_p2[1]", model_string_jags, fixed = TRUE) }
    } else if(d_list$d2$d2_ec_obs == TRUE & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == TRUE & d_list$d2$d2_ec_mis == TRUE) {
      d2 <- ifelse(d2 == 4, 2, d2)
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[2] <- meanc_p2[2]", "mu_c_p2[2] <- meanc_p2[2] + Delta_c[2]", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001)", "beta_p2[1, 2] <- beta_p2[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001)", "beta_p2[j, 2] <- beta_p2[j, 1]", model_string_jags, fixed = TRUE) 
        }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001)", "beta_p2[2] <- beta_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000)", "s_c_p2[2] <- s_c_p2[1]", model_string_jags, fixed = TRUE) }
    } 
  } 
  if(d_list$n_patterns[2] == 2) {
    if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p2[3] <- meanc_p2[3]", "", model_string_jags, fixed = TRUE) }
    if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p2[3] <- meane_p2[3]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[3] <- meanc_p2[3] + Delta_c[2]", "", model_string_jags, fixed = TRUE) }
    if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[3] <- meane_p2[3]", "", model_string_jags, fixed = TRUE) }
    if(pe_fixed > 1) {
      model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      }
    if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) }
    if(pc_fixed > 1) {
      model_string_jags <- gsub("beta_p2[1, 3] <- beta_p2[1, 1]", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("beta_p2[j, 3] <- beta_p2[j, 1]", "", model_string_jags, fixed = TRUE) 
      }
    if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[3] <- beta_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[3] <- ls_c_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[3] <- s_c_p2[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "beta") {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3])))", "", model_string_jags, fixed = TRUE) }
    if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "exp") {model_string_jags <- gsub("s_e_p2[3] <- mu_e_p2[3]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "bern") {model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3] * (1 - mu_e_p2[3]))", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "pois") {model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3])", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
    if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p2[3] <- beta_f_p2[1]" , "", model_string_jags, fixed = TRUE) }
    if(d_list$d2$d2_ec_obs == TRUE & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == TRUE & d_list$d2$d2_ec_mis == FALSE) {
      d2 <- ifelse(d2 == 3, 2, d2)
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[2] <- meanc_p2[2]", "mu_c_p2[2] <- meanc_p2[2] + Delta_c[2]", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001)", "beta_p2[1, 2] <- beta_p2[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001)", "beta_p2[j, 2] <- beta_p2[j, 1]", model_string_jags, fixed = TRUE) 
        }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001)", "beta_p2[2] <- beta_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000)", "s_c_p2[2] <- s_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 1]", "alpha_p2[1, 2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[j, 2] <- alpha_p2[j, 1]", "alpha_p2[j, 2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) 
        }
      if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[2] <- alpha_p2[1]", "alpha_p2[2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[1]", "ls_e_p2[2] ~ dunif(-5, 10)", model_string_jags, fixed = TRUE) }
      if(dist_e == "beta") {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] ~ dunif(0, sqrt(meane_p2[2] * (1 - meane_p2[2])))", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE) }
      if(dist_e == "exp") {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] <- mu_e_p2[2]", model_string_jags, fixed = TRUE) }
      if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }
      if(dist_e == "bern") {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] <- sqrt(mu_e_p2[2] * (1 - mu_e_p2[2]))", model_string_jags, fixed = TRUE) }
      if(dist_e == "pois") {model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[1]", "s_e_p2[2] <- sqrt(mu_e_p2[2])", model_string_jags, fixed = TRUE) }
      if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[2] <- tau_e_p2[1]", "tau_e_p2[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }
    } else if(d_list$d2$d2_ec_obs == TRUE & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == TRUE) {
      d2 <- ifelse(d2 == 4, 2, d2)
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[2] <- meanc_p2[2]", "mu_c_p2[2] <- meanc_p2[2] + Delta_c[2]", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001)", "beta_p2[1, 2] <- beta_p2[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001)", "beta_p2[j, 2] <- beta_p2[j, 1]", model_string_jags, fixed = TRUE) 
        }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001)", "beta_p2[2] <- beta_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100)", "ls_c_p2[2] <- ls_c_p2[1]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000)", "s_c_p2[2] <- s_c_p2[1]", model_string_jags, fixed = TRUE) }
    } 
   }
  }
  if(restriction == "AC") {
    if(d_list$n_patterns[2] == 3 | d_list$n_patterns[2] == 2) {
      if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p2[4] <- meanc_p2[4]", "", model_string_jags, fixed = TRUE) }
      if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p2[4] <- meane_p2[4]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[4] <- meanc_p2[4] + Delta_c[2]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[4] <- meane_p2[4] + Delta_e[2]", "", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p2[1, 4] <- alpha_p2[1, 3]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[j, 4] <- alpha_p2[j, 3]", "", model_string_jags, fixed = TRUE)
      }
      if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[4] <- alpha_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 4] <- beta_p2[1, 2]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 4] <- beta_p2[j, 2]", "", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[4] <- beta_p2[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[4] <- ls_c_p2[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[4] <- s_c_p2[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[4] <- ls_e_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("beta", "gamma", "exp", "weibull", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p2[4] <- s_e_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[4] <- tau_e_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p2[4] <- beta_f_p2[1]" , "", model_string_jags, fixed = TRUE) }
    }
    if(d_list$n_patterns[2] == 3){
      if(d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_c_obs == TRUE & d_list$d2$d2_e_obs == TRUE & d_list$d2$d2_ec_mis == TRUE & ind_fixed == TRUE) {
        d2 <- ifelse(d2 == 4, 1, d2)
        if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[1] <- meane_p2[1]", "mu_e_p2[1] <- meane_p2[1] + Delta_e[2]", model_string_jags, fixed = TRUE) }
        if(pe_fixed > 1) {
          model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001)", "alpha_p2[1, 1] <- alpha_p2[1, 3]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001)", "alpha_p2[j, 1] <- alpha_p2[j, 3]", model_string_jags, fixed = TRUE) 
        }
        if(pe_fixed == 1) {model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001)", "alpha_p2[1] <- alpha_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "norm") {model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10)", "ls_e_p2[1] <- ls_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1])))", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000)", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p2[1] <- mu_e_p2[1]", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weibull") {model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100)", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p2[1] <- sqrt(mu_e_p2[1] * (1 - mu_e_p2[1]))", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p2[1] <- sqrt(mu_e_p2[1])", "s_e_p2[1] <- s_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "nbinom") {model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100)", "tau_e_p2[1] <- tau_e_p2[3]", model_string_jags, fixed = TRUE) }
        if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[1] <- meanc_p2[1]", "mu_c_p2[1] <- meanc_p2[1] + Delta_c[2]", model_string_jags, fixed = TRUE) }
        if(pc_fixed > 1) {
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001)", "beta_p2[1, 1] <- beta_p2[1, 2]", model_string_jags, fixed = TRUE) 
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001)", "beta_p2[j, 1] <- beta_p2[j, 2]", model_string_jags, fixed = TRUE) 
        }
        if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001)", "beta_p2[1] <- beta_p2[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10)", "ls_c_p2[1] <- ls_c_p2[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100)", "ls_c_p2[1] <- ls_c_p2[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000)", "s_c_p2[1] <- s_c_p2[2]", model_string_jags, fixed = TRUE) }
      }
    } 
    if(d_list$n_patterns[2] == 2) {
      d2 <- ifelse(d2 == 3, 1, d2)
      if(type == "MAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_c_p2[3] <- meanc_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_e_p2[3] <- meane_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[3] <- meanc_p2[3] + Delta_c[2]", "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_eff") {model_string_jags <- gsub("mu_e_p2[3] <- meane_p2[3]", "", model_string_jags, fixed = TRUE) }
      if(pe_fixed > 1) {
        model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 3]", "alpha_p2[1, 2] <- alpha_p2[1, 1]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[j, 2] <- alpha_p2[j, 3]", "alpha_p2[j, 2] <- alpha_p2[j, 1]", model_string_jags, fixed = TRUE) 
      }
      if(pe_fixed == 1) {
        model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("alpha_p2[2] <- alpha_p2[3]", "alpha_p2[2] <- alpha_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 3] <- beta_p2[1, 2]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 3] <- beta_p2[j, 2]", "", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {
        model_string_jags <- gsub("beta_p2[3] <- beta_p2[2]", "", model_string_jags, fixed = TRUE) 
      }
      if(dist_c == "norm" | dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[3] <- ls_c_p2[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[3] <- s_c_p2[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {
        model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("ls_e_p2[2] <- ls_e_p2[3]", "ls_e_p2[2] <- ls_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3])))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "exp") {
        model_string_jags <- gsub("s_e_p2[3] <- mu_e_p2[3]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "weibull") {
        model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "bern") {
        model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3] * (1 - mu_e_p2[3]))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "pois") {
        model_string_jags <- gsub("s_e_p2[3] <- sqrt(mu_e_p2[3])", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p2[2] <- s_e_p2[3]", "s_e_p2[2] <- s_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "nbinom") {
        model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("tau_e_p2[2] <- tau_e_p2[3]", "tau_e_p2[2] <- tau_e_p2[1]", model_string_jags, fixed = TRUE) 
      }
      if(ind_fixed == FALSE) {model_string_jags <- gsub("beta_f_p2[3] <- beta_f_p2[1]" , "", model_string_jags, fixed = TRUE) }
      if(type == "MNAR" | type == "MNAR_cost") {model_string_jags <- gsub("mu_c_p2[1] <- meanc_p2[1]", "mu_c_p2[1] <- meanc_p2[1] + Delta_c[2]", model_string_jags, fixed = TRUE) }
      if(pc_fixed > 1) {
        model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001)", "beta_p2[1, 1] <- beta_p2[1, 2]", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001)", "beta_p2[j, 1] <- beta_p2[j, 2]", model_string_jags, fixed = TRUE) 
      }
      if(pc_fixed == 1) {model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001)", "beta_p2[1] <- beta_p2[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "norm") {model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10)", "ls_c_p2[1] <- ls_c_p2[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100)", "ls_c_p2[1] <- ls_c_p2[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000)", "s_c_p2[1] <- s_c_p2[2]", model_string_jags, fixed = TRUE) }
    }
  }
  model_string_jags <- prior_pattern(type = type, dist_e = dist_e, dist_c = dist_c, pe_fixed = pe_fixed, pc_fixed = pc_fixed, 
                                     model_e_random = model_e_random, model_c_random = model_c_random, pe_random = pe_random, pc_random = pc_random, 
                                     d_list = d_list, restriction = restriction)
  writeLines(model_string_jags, "pattern.txt")
  model_string <- "pattern.txt"
  return(model_string)
}))