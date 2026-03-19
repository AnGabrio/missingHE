#'An internal function to select which type of pattern mixture model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Pattern mixture models
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



write_pattern <- function(dist_e , dist_c, type, model_txt_info) {
  model_string_jags<-  "
  model {
  
  #full sample
  for(i in 1:n) {
  #costs and effects model
  cost[i] ~ dnorm(cmu_c[i], tau_c_p[d_mod[i]])
  eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])
  
  #derive mean and std effects 
  #derive mean and std costs
  
  #mean regression in each pattern
  cmu_c[i] <- inprod(X_c_fixed[i, ], beta_p[, d_mod[i]]) + beta_f_p[d_mod[i]] * (eff[i] - meane_p[d_mod[i]]) + inprod(X_c_random[i, ], b[, clus_c[i]]) + b_f[clus_c[i]] * (eff[i] - tmu_e)
  cmu_e[i] <- inprod(X_e_fixed[i, ], alpha_p[, d_mod[i]]) + inprod(X_e_random[i, ], a[, clus_e[i]])

  #patterns model for t=1
  d_mod[i] ~ dcat(p_prob[1:n_patterns])
  
  #loglikelihood
  loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])
  loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c_p[d_mod[i]])
  loglik_d[i] <- logdensity.cat(d_mod[i], p_prob[1:n_patterns])
  }
  
  #transformation of parameters
  for (d in 1:n_patterns) {#begin transformation
  tau_c_p[d] <- 1 / ss_c_p[d]
  ss_c_p[d] <- s_c_p[d] * s_c_p[d]
  #std for lnorm 
  #mean for lnorm
  tau_e_p[d] <- 1 / ss_e_p[d]
  ss_e_p[d] <- s_e_p[d] * s_e_p[d]
  }# end transformation 
  
  #transformation of random effects parameters
  # begin transformation random effects
  for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]
  ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }
  for(j in 1:pe_random) {tau_a_hat[j] <- 1 / ss_a_hat[j]
  ss_a_hat[j] <- s_a_hat[j] * s_a_hat[j] }
  #end transformation of random effects 
  
  #calculate means at mean of covariates in each pattern
  for(d in 1:n_patterns){
  meanc_p[d] <- inprod(mean_cov_c_fixed[], beta_p[, d])
  meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])
  }

  #add sensitivity parameters to means  
  mu_c_p[1] <- meanc_p[1]
  mu_c_p[2] <- meanc_p[2]
  mu_c_p[3] <- meanc_p[3] + delta_c
  mu_c_p[4] <- meanc_p[4] + delta_c
  mu_e_p[1] <- meane_p[1]
  mu_e_p[2] <- meane_p[2] + delta_e
  mu_e_p[3] <- meane_p[3] 
  mu_e_p[4] <- meane_p[4] + delta_e

  #calculate overall means
  tmu_c <- sum(mu_c_p[] * p_prob[]) + inprod(mean_cov_c_random[], mu_b_hat[]) 
  tmu_e <- sum(mu_e_p[] * p_prob[]) + inprod(mean_cov_e_random[], mu_a_hat[])

  #priors
  
  #priors for mean regression coefficients
  for (j in 1:pe_fixed) {#begin alpha priors effects in each pattern
  alpha_p[j, 1] ~ dnorm(0, 0.0000001)
  alpha_p[j, 2] <- alpha_p[j, 1]
  alpha_p[j, 3] ~ dnorm(0, 0.0000001)
  alpha_p[j, 4] <- alpha_p[j, 1]
  }#end alpha priors effects in each pattern

  for (j in 1:pc_fixed) {#begin beta priors costs in each pattern
  beta_p[j, 1] ~ dnorm(0, 0.0000001)
  beta_p[j, 2] ~ dnorm(0, 0.0000001)
  beta_p[j, 3] <- beta_p[j, 1]
  beta_p[j, 4] <- beta_p[j, 1]
  }#end beta priors costs in each pattern
  
  #priors for mean regression random coefficients
  for (j in 1:pe_random) {#begin a priors effects
  for(s in 1:n_clus_e) {a[j, s] ~ dnorm(mu_a_hat[j], tau_a_hat[j]) }
  }#end a priors effects

  for (j in 1:pc_random) {#begin b priors costs
  for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }
  }#end b priors costs
  
  #standard deviation priors
  s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_c_p[3] <- s_c_p[1]
  s_c_p[4] <- s_c_p[1]
  
  s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_e_p[2] <- s_e_p[1]
  s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)
  s_e_p[4] <- s_e_p[1]

  #correlation
  beta_f_p[1] ~ dnorm(0, 0.0000001) 
  beta_f_p[2] <- beta_f_p[1]
  beta_f_p[3] <- beta_f_p[1] 
  beta_f_p[4] <- beta_f_p[1]

  # mean and sd mean regression random coefficients priors
  for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)
  s_b_hat[j] ~ dunif(0, 100) }
  for(j in 1:pe_random) {mu_a_hat[j] ~ dnorm(0, 0.001)
  s_a_hat[j] ~ dunif(0, 100) }  
  # end mean and sd mean regression random coefficients priors

  # correlation random effects
  for(s in 1:n_clus_c) {b_f[s] ~ dnorm(mu_b_f_hat, tau_b_f_hat) }
  mu_b_f_hat ~ dnorm(0, 0.001)
  tau_b_f_hat <- 1 / ss_b_f_hat
  ss_b_f_hat <- s_b_f_hat * s_b_f_hat
  s_b_f_hat ~ dunif(0, 100)  

  #priors on patterns model
  p_prob ~ ddirch(pp[])

  #begin priors on sensitivity parameters
  delta_c ~ dunif(0, 1)  
  delta_e ~ dunif(0, 1) 
  #end priors on sensitivity parameters

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
    if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE) }
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
    if(length(model_txt_info$model_e_random) == 0) { model_string_jags <- gsub("#priors for mean regression random coefficients", "", model_string_jags, fixed = TRUE) }
  }
  if(length(model_txt_info$model_c_random) == 0 & length(model_txt_info$model_e_random) == 0 | length(model_txt_info$model_e_random) == 0 & model_txt_info$pc_random == 0) { 
    model_string_jags <- gsub("# mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# end mean and sd mean regression random coefficients priors", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#transformation of random effects parameters", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("# begin transformation random effects", "", model_string_jags, fixed = TRUE) 
    model_string_jags <- gsub("#end transformation of random effects", "", model_string_jags, fixed = TRUE) 
  }
  if(type == "MAR") {
    model_string_jags <- gsub(" + delta_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_e ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + delta_c", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#begin priors on sensitivity parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#end priors on sensitivity parameters", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#add sensitivity parameters to means", "", model_string_jags, fixed = TRUE)
  } else if(type == "MNAR_eff") {
    model_string_jags <- gsub(" + delta_c", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_c ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  } else if(type == "MNAR_cost") {
    model_string_jags <- gsub(" + delta_e", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("delta_e ~ dunif(0, 1)", "", model_string_jags, fixed = TRUE)
  }  
  if(model_txt_info$restriction == "AC") {
    model_string_jags <- gsub("alpha_p[j, 2] <- alpha_p[j, 1]", "alpha_p[j, 2] <- alpha_p[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p[j, 4] <- alpha_p[j, 1]", "alpha_p[j, 4] <- alpha_p[j, 3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p[j, 3] <- beta_p[j, 1]", "beta_p[j, 3] <- beta_p[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_p[j, 4] <- beta_p[j, 1]", "beta_p[j, 4] <- beta_p[j, 2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[3] <- s_c_p[1]", "s_c_p[3] <- s_c_p[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[4] <- s_c_p[1]", "s_c_p[4] <- s_c_p[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] <- s_e_p[3]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[4] <- s_e_p[1]", "s_e_p[4] <- s_e_p[3]", model_string_jags, fixed = TRUE)
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
    model_string_jags <- gsub(" + beta_f_p[d_mod[i]] * (eff[i] - meane_p[d_mod[i]])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p[2] <- beta_f_p[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p[3] <- beta_f_p[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f_p[4] <- beta_f_p[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
  } 
  if(dist_c == "norm") {
    model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
  }
  if(dist_c == "gamma") {
    model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c_p[d_mod[i]])", "cost[i] ~ dgamma(cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs", "ctau_c[i] <- cmu_c[i] / pow(s_c_p[d_mod[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_c[i] <- ", "log(cmu_c[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p[d] <- 1 / ss_c_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p[d] <- s_c_p[d] * s_c_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p[d] <- inprod(mean_cov_c_fixed[], beta_p[, d])", "meanc_p[d] <- exp(inprod(mean_cov_c_fixed[], beta_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_c_random[], mu_b_hat[])", "+ exp(inprod(mean_cov_c_random[], mu_b_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("#mean for lnorm", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#std for lnorm", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c_p[d_mod[i]])", "loglik_c[i] <- logdensity.gamma(cost[i], cmu_c[i] * ctau_c[i], ctau_c[i])", model_string_jags, fixed = TRUE)
  } else if(dist_c == "lnorm") {
    model_string_jags <- gsub("cost[i] ~ dnorm(cmu_c[i], tau_c_p[d_mod[i]])", "cost[i] ~ dlnorm(clmu_c[i], ltau_c_p[d_mod[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std costs", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_c[i] <- ", "clmu_c[i] <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_c_p[d] <- 1 / ss_c_p[d]", "ltau_c_p[d] <- 1 / lss_c_p[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_c_p[d] <- s_c_p[d] * s_c_p[d]", "lss_c_p[d] <- ls_c_p[d] * ls_c_p[d]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#std for lnorm", "s_c_p[d] <- sqrt(exp(2 * lmeanc_p[d] + lss_c_p[d]) * (exp(lss_c_p[d]) - 1))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c_p[1] ~ dunif(0, 10)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c_p[2] ~ dunif(0, 10)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[3] <- s_c_p[1]", "ls_c_p[3] <- ls_c_p[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[4] <- s_c_p[1]", "ls_c_p[4] <- ls_c_p[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[3] <- s_c_p[2]", "ls_c_p[3] <- ls_c_p[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_c_p[4] <- s_c_p[2]", "ls_c_p[4] <- ls_c_p[2]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meanc_p[d] <- inprod(mean_cov_c_fixed[], beta_p[, d])", "lmeanc_p[d] <- inprod(mean_cov_c_fixed[], beta_p[, d])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#mean for lnorm", "meanc_p[d] <- exp(lmeanc_p[d] + lss_c_p[d] / 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_c[i] <- logdensity.norm(cost[i], cmu_c[i], tau_c_p[d_mod[i]])", "loglik_c[i] <- logdensity.lnorm(cost[i], clmu_c[i], ltau_c_p[d_mod[i]])", model_string_jags, fixed = TRUE)
     if(length(model_txt_info$model_c_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_c_random[], mu_b_hat[])", "+ exp(inprod(mean_cov_c_random[], mu_b_hat[]))", model_string_jags, fixed = TRUE)
     }
  }
  if(dist_e == "norm") {
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "beta") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dbeta(cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- (cmu_e[i] * (1 - cmu_e[i]) / pow(s_e_p[d_mod[i]], 2) - 1)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- ilogit(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ ilogit(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.beta(eff[i], cmu_e[i] * ctau_e[i], (1 - cmu_e[i]) * ctau_e[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "gamma") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dgamma(cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- cmu_e[i] / pow(s_e_p[d_mod[i]], 2)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- exp(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ exp(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.gamma(eff[i], cmu_e[i] * ctau_e[i], ctau_e[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "exp") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dexp(1 / cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- exp(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] <- mu_e_p[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] <- mu_e_p[3]", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ exp(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.exp(eff[i], 1 / cmu_e[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "weib") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dweib(ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "ctau_e[i] <- pow(s_e_p[d_mod[i]] / cmu_e[i], - 1.086)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- exp(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ exp(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.weib(eff[i], ctau_e[i], pow(1 / (cmu_e[i] / exp(loggam(1 + 1/ctau_e[i]))), ctau_e[i]))", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "logis") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dlogis(cmu_e[i], tau_e_p[d_mod[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "tau_e_p[d] <- 1 / sqrt((3 * ss_e_p[d]) / pow(3.14159265 , 2))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.logis(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "bern") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dbern(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "logit(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- ilogit(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] <- sqrt(mu_e_p[1] * (1 - mu_e_p[1]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] <- sqrt(mu_e_p[3] * (1 - mu_e_p[3]))", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ ilogit(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.bern(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "pois") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dpois(cmu_e[i])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- exp(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] <- sqrt(mu_e_p[1])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] <- sqrt(mu_e_p[3])", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ exp(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.pois(eff[i], cmu_e[i])", model_string_jags, fixed = TRUE)
  }
  if(dist_e == "negbin") {
    model_string_jags <- gsub("eff[i] ~ dnorm(cmu_e[i], tau_e_p[d_mod[i]])", "eff[i] ~ dnegbin(tau_e_p[d_mod[i]] / (tau_e_p[d_mod[i]] + cmu_e[i]), tau_e_p[d_mod[i]])", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#derive mean and std effects", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("cmu_e[i] <- ", "log(cmu_e[i]) <- ", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("tau_e_p[d] <- 1 / ss_e_p[d]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_e_p[d] <- s_e_p[d] * s_e_p[d]", "s_e_p[d] <- sqrt((tau_e_p[d] / (tau_e_p[d] + tmu_e_p[d])) * tau_e_p[d]) / (1 - (tau_e_p[d] / (tau_e_p[d] + tmu_e_p[d])))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("meane_p[d] <- inprod(mean_cov_e_fixed[], alpha_p[, d])", "meane_p[d] <- exp(inprod(mean_cov_e_fixed[], alpha_p[, d]))", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e_p[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[2] <- s_e_p[2]", "tau_e_p[2] <- tau_e_p[1]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "tau_e_p[3] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_e_p[4] <- s_e_p[4]", "tau_e_p[4] <- tau_e_p[4]", model_string_jags, fixed = TRUE)
    if(length(model_txt_info$model_e_random) != 0) {
      model_string_jags <- gsub("+ inprod(mean_cov_e_random[], mu_a_hat[])", "+ exp(inprod(mean_cov_e_random[], mu_a_hat[]))", model_string_jags, fixed = TRUE)
    }
    if(model_txt_info$restriction == "AC"){
      model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "tau_e_p[2] <- tau_e_p[3]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[4] <- s_e_p[3]", "tau_e_p[4] <- tau_e_p[3]", model_string_jags, fixed = TRUE)
    }
    model_string_jags <- gsub("loglik_e[i] <- logdensity.norm(eff[i], cmu_e[i], tau_e_p[d_mod[i]])", "loglik_e[i] <- logdensity.negbin(eff[i], tau_e_p[d_mod[i]] / (tau_e_p[d_mod[i]] + cmu_e[i]), tau_e_p[d_mod[i]])", model_string_jags, fixed = TRUE)
  }
  if(dist_e %in% c("beta", "gamma", "weib", "bern", "pois", "exp")) {
    if(dist_c == "gamma") {
    model_string_jags <- gsub("for (d in 1:n_patterns) {#begin transformation", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}# end transformation", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
    }
  }
  if(length(model_txt_info$model_e_random) != 0 & model_txt_info$pe_random == 1) {
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
  if(length(model_txt_info$model_c_random) != 0 & model_txt_info$pc_random == 1) {
    model_string_jags <- gsub("inprod(X_c_random[i, ], b[, clus_c[i]])", "X_c_random[i] * b[clus_c[i]]", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("inprod(mean_cov_c1_random[], mu_b_hat[])", "mean_cov_c_random * mu_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for (j in 1:pc_random) {#begin b priors costs", "#begin b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(s in 1:n_clus_c) {b[j, s] ~ dnorm(mu_b_hat[j], tau_b_hat[j]) }", "for(s in 1:n_clus_c) {b[s] ~ dnorm(mu_b_hat, tau_b_hat) }", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("}#end b priors costs", "#end b priors costs", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {mu_b_hat[j] ~ dnorm(0, 0.001)", "mu_b_hat ~ dnorm(0, 0.001)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("s_b_hat[j] ~ dunif(0, 100) }", "s_b_hat ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("for(j in 1:pc_random) {tau_b_hat[j] <- 1 / ss_b_hat[j]", "tau_b_hat <- 1 / ss_b_hat", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ss_b_hat[j] <- s_b_hat[j] * s_b_hat[j] }", "ss_b_hat <- s_b_hat * s_b_hat", model_string_jags, fixed = TRUE)
  }
  if(model_txt_info$restriction == "CC") {
   if(model_txt_info$n_patterns %in% c(2, 3)) {
    if(type %in% c("MAR", "MNAR_eff")) {model_string_jags <- gsub("mu_c_p[4] <- meanc_p[4]", "", model_string_jags, fixed = TRUE) }
    if(type %in% c("MAR", "MNAR_cost")) {model_string_jags <- gsub("mu_e_p[4] <- meane_p[4]", "", model_string_jags, fixed = TRUE) }
    if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[4] <- meanc_p[4] + delta_c", "", model_string_jags, fixed = TRUE) }
    if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[4] <- meane_p[4] + delta_e", "", model_string_jags, fixed = TRUE) }
      model_string_jags <- gsub("alpha_p[j, 4] <- alpha_p[j, 1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 4] <- beta_p[j, 1]", "", model_string_jags, fixed = TRUE) 
    if(dist_c %in% c("norm", "gamma")) {model_string_jags <- gsub("s_c_p[4] <- s_c_p[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_c == "lnorm") {
      model_string_jags <- gsub("ls_c_p[4] <- ls_c_p[1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c_p[4] <- s_c_p[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "norm") {model_string_jags <- gsub("s_e_p[4] <- s_e_p[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e %in% c("beta", "gamma", "exp", "weib", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p[4] <- s_e_p[1]", "", model_string_jags, fixed = TRUE) }
    if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[4] <- tau_e_p[1]", "", model_string_jags, fixed = TRUE) }
    if(!model_txt_info$ind_fixed) {model_string_jags <- gsub("beta_f_p[4] <- beta_f_p[1]" , "", model_string_jags, fixed = TRUE) }
   }
    if(model_txt_info$n_patterns == 3){
      if(all(model_txt_info$d_or %in% c(1, 2, 4))) {
        if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[3] <- meane_p[3]", "mu_e_p[3] <- meane_p[3] + delta_e", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", "alpha_p[j, 3] <- alpha_p[j, 1]", model_string_jags, fixed = TRUE) 
        if(dist_e == "norm") {model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p[3] <- mu_e_p[3]", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weib") {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3] * (1 - mu_e_p[3]))", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3])", "s_e_p[3] <- s_e_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", "tau_e_p[3] <- tau_e_p[1]", model_string_jags, fixed = TRUE) }
      }
      if(all(model_txt_info$d_or %in% c(1, 3, 4))) {
        if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[2] <- meanc_p[2]", "mu_c_p[2] <- meanc_p[2] + delta_c", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", "beta_p[j, 2] <- beta_p[j, 1]", model_string_jags, fixed = TRUE) 
        if(dist_c == "norm") {model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[2] <- s_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", "ls_c_p[2] <- ls_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", "s_c_p[2] <- s_c_p[1]", model_string_jags, fixed = TRUE) }        
      }
    }
    if(model_txt_info$n_patterns == 2){
      if(type %in% c("MAR", "MNAR_eff")) {model_string_jags <- gsub("mu_c_p[3] <- meanc_p[3]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MAR", "MNAR_cost")) {model_string_jags <- gsub("mu_e_p[3] <- meane_p[3]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[3] <- meanc_p[3] + delta_c", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[3] <- meane_p[3]", "", model_string_jags, fixed = TRUE) }
      model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("beta_p[j, 3] <- beta_p[j, 1]", "", model_string_jags, fixed = TRUE) 
      if(dist_c %in% c("norm", "gamma")) {model_string_jags <- gsub("s_c_p[3] <- s_c_p[1]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {
        model_string_jags <- gsub("ls_c_p[3] <- ls_c_p[1]", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[3] <- s_c_p[1]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "beta") {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", "", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "exp") {model_string_jags <- gsub("s_e_p[3] <- mu_e_p[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "weib") {model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "bern") {model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3] * (1 - mu_e_p[3]))", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "pois") {model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3])", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) }
      if(!model_txt_info$ind_fixed) {model_string_jags <- gsub("beta_f_p[3] <- beta_f_p[1]" , "", model_string_jags, fixed = TRUE) }  
      if(all(model_txt_info$d_or %in% c(1, 3))) {
        if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[2] <- meanc_p[2]", "mu_c_p[2] <- meanc_p[2] + delta_c", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", "beta_p[j, 2] <- beta_p[j, 1]", model_string_jags, fixed = TRUE) 
        if(dist_c == "norm") {model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", "ls_c_p[2] <- ls_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", "ls_c_p[2] <- ls_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", "s_c_p[2] <- s_c_p[1]", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("alpha_p[j, 2] <- alpha_p[j, 1]", "alpha_p[j, 2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE) 
        if(dist_e == "norm") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] ~ dunif(0, sqrt(mu_e_p[2] * (1 - mu_e_p[2])))", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] <- mu_e_p[2]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weib") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] <- sqrt(mu_e_p[2] * (1 - mu_e_p[2]))", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p[2] <- s_e_p[1]", "s_e_p[2] <- sqrt(mu_e_p[2])", model_string_jags, fixed = TRUE) }
        if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[2] <- tau_e_p[1]", "tau_e_p[2] ~ dunif(0, 100)", model_string_jags, fixed = TRUE) }        
      }
      if(all(model_txt_info$d_or %in% c(1, 4))) {
        if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[2] <- meanc_p[2]", "mu_c_p1[2] <- meanc_p[2] + delta_c", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", "beta_p[j, 2] <- beta_p[j, 1]", model_string_jags, fixed = TRUE) 
        if(dist_c == "norm") {model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[2] <- s_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", "ls_c_p[2] <- ls_c_p[1]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", "s_c_p[2] <- s_c_p[1]", model_string_jags, fixed = TRUE) }        
      }
    }
  }
  if(model_txt_info$restriction == "AC") {
    if(model_txt_info$n_patterns %in% c(2, 3)) {
      if(type %in% c("MAR", "MNAR_eff")) {model_string_jags <- gsub("mu_c_p[4] <- meanc_p[4]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MAR", "MNAR_cost")) {model_string_jags <- gsub("mu_e_p[4] <- meane_p[4]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[4] <- meanc_p[4] + delta_c", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[4] <- meane_p[4] + delta_e", "", model_string_jags, fixed = TRUE) }
      model_string_jags <- gsub("alpha_p[j, 4] <- alpha_p[j, 3]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 4] <- beta_p[j, 2]", "", model_string_jags, fixed = TRUE) 
      if(dist_c %in% c("norm", "gamma")) {model_string_jags <- gsub("s_c_p[4] <- s_c_p[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {
        model_string_jags <- gsub("ls_c_p[4] <- ls_c_p[2]", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[4] <- s_c_p[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {model_string_jags <- gsub("s_e_p1[4] <- s_e_p1[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e %in% c("beta", "gamma", "exp", "weib", "logis", "bern", "pois")) {model_string_jags <- gsub("s_e_p[4] <- s_e_p[3]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[4] <- tau_e_p[3]", "", model_string_jags, fixed = TRUE) }
      if(!model_txt_info$ind_fixed) {model_string_jags <- gsub("beta_f_p[4] <- beta_f_p[1]" , "", model_string_jags, fixed = TRUE) }      
    }
    if(model_txt_info$n_patterns == 3){
      if(all(model_txt_info$d_or %in% c(2, 3, 4)) & model_txt_info$ind_fixed) {
        if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[1] <- meane_p[1]", "mu_e_p[1] <- meane_p[1] + delta_e", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", "alpha_p[j, 1] <- alpha_p[j, 3]", model_string_jags, fixed = TRUE) 
        if(dist_e == "norm") {model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "beta") {model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e %in% c("gamma", "logis")) {model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "exp") {model_string_jags <- gsub("s_e_p[1] <- mu_e_p[1]", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "weib") {model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "bern") {model_string_jags <- gsub("s_e_p[1] <- sqrt(mu_e_p[1] * (1 - mu_e_p[1]))", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "pois") {model_string_jags <- gsub("s_e_p[1] <- sqrt(mu_e_p[1])", "s_e_p[1] <- s_e_p[3]", model_string_jags, fixed = TRUE) }
        if(dist_e == "negbin") {model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", "tau_e_p[1] <- tau_e_p[3]", model_string_jags, fixed = TRUE) }        
        if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[1] <- meanc_p[1]", "mu_c_p[1] <- meanc_p[1] + delta_c", model_string_jags, fixed = TRUE) }
        model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", "beta_p[j, 1] <- beta_p[j, 2]", model_string_jags, fixed = TRUE) 
        if(dist_c == "norm") {model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[1] <- s_c_p[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", "ls_c_p[1] <- ls_c_p[2]", model_string_jags, fixed = TRUE) }
        if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", "s_c_p[1] <- s_c_p[2]", model_string_jags, fixed = TRUE) }
      }
    }
    if(model_txt_info$n_patterns == 2){
      if(type %in% c("MAR", "MNAR_eff")) {model_string_jags <- gsub("mu_c_p[3] <- meanc_p[3]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MAR", "MNAR_cost")) {model_string_jags <- gsub("mu_e_p[3] <- meane_p[3]", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[3] <- meanc_p[3] + delta_c", "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_eff")) {model_string_jags <- gsub("mu_e_p[3] <- meane_p[3]", "", model_string_jags, fixed = TRUE) }      
      model_string_jags <- gsub("alpha_p[1, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p[1, 2] <- alpha_p[1, 3]", "alpha_p[1, 2] <- alpha_p[1, 1]", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("alpha_p[j, 2] <- alpha_p[j, 3]", "alpha_p[j, 2] <- alpha_p[j, 1]", model_string_jags, fixed = TRUE) 
      model_string_jags <- gsub("beta_p[j, 3] <- beta_p[j, 2]", "", model_string_jags, fixed = TRUE) 
      if(dist_c %in% c("norm", "gamma")) {model_string_jags <- gsub("s_c_p[3] <- s_c_p[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {
        model_string_jags <- gsub("ls_c_p[3] <- ls_c_p[2]", "", model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[3] <- s_c_p[2]", "", model_string_jags, fixed = TRUE) }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "exp") {
        model_string_jags <- gsub("s_e_p[3] <- mu_e_p[3]", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "weib") {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "bern") {
        model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3] * (1 - mu_e_p[3]))", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "pois") {
        model_string_jags <- gsub("s_e_p[3] <- sqrt(mu_e_p[3])", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("s_e_p[2] <- s_e_p[3]", "s_e_p[2] <- s_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(dist_e == "negbin") {
        model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", "", model_string_jags, fixed = TRUE) 
        model_string_jags <- gsub("tau_e_p[2] <- tau_e_p[3]", "tau_e_p[2] <- tau_e_p[1]", model_string_jags, fixed = TRUE) 
      }
      if(!model_txt_info$ind_fixed) {model_string_jags <- gsub("beta_f_p[3] <- beta_f_p[1]" , "", model_string_jags, fixed = TRUE) }
      if(type %in% c("MNAR", "MNAR_cost")) {model_string_jags <- gsub("mu_c_p[1] <- meanc_p[1]", "mu_c_p[1] <- meanc_p[1] + delta_c", model_string_jags, fixed = TRUE) }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", "beta_p[j, 1] <- beta_p[j, 2]", model_string_jags, fixed = TRUE) 
      if(dist_c == "norm") {model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", "s_c_p[1] <- s_c_p[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "lnorm") {model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", "ls_c_p[1] <- ls_c_p[2]", model_string_jags, fixed = TRUE) }
      if(dist_c == "gamma") {model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", "s_c_p[1] <- s_c_p[2]", model_string_jags, fixed = TRUE) }
    }
  }
  model_string_jags <- prior_pattern(type = type, dist_e = dist_e, dist_c = dist_c, 
                                     model_txt_info = model_txt_info, model_string_jags = model_string_jags)
  writeLines(model_string_jags, "pattern.txt")
  model_string <- "pattern.txt"
  return(model_string)
}