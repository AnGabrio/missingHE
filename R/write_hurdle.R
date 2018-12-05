#'An internal function to select which type of hurdle model to execute for both effectiveness and costs. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of structural value mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Hurdle models
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta')
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param ind Logical; if TRUE independence between effectiveness and costs is assumed, else correlation is accounted for
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR) and Structural At Random (SAR)
#' @param pe Number of covariates for the effectiveness model
#' @param pc Number of cvoariates for the cost model
#' @param ze Number of covariates or the structural indicators model for the effectiveness
#' @param zc Number of covariates or the structural indicators model for the costs
#' @param se Structural value for the effectiveness 
#' @param sc Structural value for the costs
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_hurdle <- function(dist_e , dist_c, ind, type, pe, pc, ze, zc, se, sc) eval.parent(substitute( {
  model_string_jags <- "
  model{

  #control
  for(i in 1:N1) {
  #costs and effects model
  cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])
  eff1[i] ~ dnorm(mu_e1[i], tau_e1[d_eff1[i] + 1])
  
  #derive mean and std effects1
  #derive mean and std costs1 

  #mean regression
  mu_c1[i] <- inprod(X1_c[i, ], beta1[, d_cost1[i] + 1]) + beta_f1[d_cost1[i] + 1] * (eff1[i] - mu_e[1])
  mu_e1[i] <- inprod(X1_e[i, ], alpha1[, d_eff1[i] + 1])

  #structural values mechanism
  d_eff1[i] ~ dbern(pq_1[i])
  logit(pq_1[i]) <- inprod(Z1_e[i, ], gamma_e[, 1])
  d_cost1[i] ~ dbern(pc_1[i])
  logit(pc_1[i]) <- inprod(Z1_c[i, ], gamma_c[, 1])

  #loglikelihood
  loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1[d_eff1[i] + 1])
  loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])
  loglik_de1[i] <- logdensity.bern(d_eff1[i], pq_1[i])
  loglik_dc1[i] <- logdensity.bern(d_cost1[i], pc_1[i])
  }
  

  #intervention
  for(i in 1:N2) {
  #costs and effects model
  cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])
  eff2[i] ~ dnorm(mu_e2[i], tau_e2[d_eff2[i] + 1])

  #derive mean and std effects2 
  #derive mean and std costs2

  #mean regression
  mu_c2[i] <- inprod(X2_c[i, ], beta2[, d_cost2[i] + 1]) + beta_f2[d_cost2[i] + 1] * (eff2[i] - mu_e[2])
  mu_e2[i] <- inprod(X2_e[i, ], alpha2[, d_eff2[i] + 1])

  #structural values mechanism
  d_eff2[i] ~ dbern(pq_2[i])
  logit(pq_2[i]) <- inprod(Z2_e[i, ], gamma_e[, 2])
  d_cost2[i] ~ dbern(pc_2[i])
  logit(pc_2[i]) <- inprod(Z2_c[i, ], gamma_c[, 2])

  #loglikelihood
  loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2[d_eff2[i] + 1])
  loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])
  loglik_de2[i] <- logdensity.bern(d_eff2[i], pq_2[i])
  loglik_dc2[i] <- logdensity.bern(d_cost2[i], pc_2[i])
  }
  
  #transformation of parameters
  for (t in 1:2) {#begin transformation costs
  tau_c1[t] <- 1 / ss_c1[t]
  ss_c1[t] <- s_c1[t] * s_c1[t]
  s_c1[t] <- exp(ls_c1[t])
  tau_c2[t] <- 1 / ss_c2[t]
  ss_c2[t] <- s_c2[t] * s_c2[t]
  s_c2[t] <- exp(ls_c2[t])
  }#end transformation costs
  #mean for lnorm1
  for (t in 1:2) {#begin transformation effects
  tau_e1[t] <- 1 / ss_e1[t]
  ss_e1[t] <- s_e1[t] * s_e1[t]
  s_e1[t] <- exp(ls_e1[t])
  tau_e2[t] <- 1 / ss_e2[t]
  ss_e2[t] <- s_e2[t] * s_e2[t]
  s_e2[t] <- exp(ls_e2[t])
  }#end transformation effects
  #mean for lnorm2
  
  #calculate means at mean of covariates for non-structural values
  nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])
  nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])
  nu_e[1] <- inprod(mean_cov_e1[], alpha1[, 1])
  nu_e[2] <- inprod(mean_cov_e2[], alpha2[, 1])
  
  #obtain coefficients only for non-structural values
  for(j in 1:pe) {#alpha non-structural effects
  alpha[j, 1] <- alpha1[j, 1]
  alpha[j, 2] <- alpha2[j, 1] }
  for(j in 1:pc) {#beta non-structural costs
  beta[j, 1] <- beta1[j, 1]
  beta[j, 2] <- beta2[j, 1] }
  
  #weighted mean parameters
  mu_e[1] <- nu_e[1] * (1 - p_e[1]) + se * p_e[1]
  mu_e[2] <- nu_e[2] * (1 - p_e[2]) + se * p_e[2]
  
  mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]
  mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]
  
  #sd parameters for non-structural values
  s_e[1] <- s_e1[1]
  s_e[2] <- s_e2[1]
  
  s_c[1] <- s_c1[1]
  s_c[2] <- s_c2[1]
  
  #structural values probability
  p_c[1] <- ilogit(inprod(mean_z_c1[], gamma_c[, 1]))
  p_c[2] <- ilogit(inprod(mean_z_c2[], gamma_c[, 2]))
  p_e[1] <- ilogit(inprod(mean_z_e1[], gamma_e[, 1]))
  p_e[2] <- ilogit(inprod(mean_z_e2[], gamma_e[, 2]))
  
  #priors
  
  #priors for mean regression coefficients
  for (j in 2:pe) {#begin alpha priors effects
  alpha1[j, 1] ~ dnorm(0, 0.0000001)
  alpha2[j, 1] ~ dnorm(0, 0.0000001)
  alpha1[j, 2] <- 0
  alpha2[j, 2] <- 0
  }#end alpha priors effects
  alpha1[1, 1] ~ dnorm(0, 0.0000001)
  alpha2[1, 1] ~ dnorm(0, 0.0000001)
  alpha1[1, 2] <- se
  alpha2[1, 2] <- se
  
  for (j in 2:pc) {#begin beta priors costs
  beta1[j, 1] ~ dnorm(0, 0.0000001)
  beta2[j, 1] ~ dnorm(0, 0.0000001)
  beta1[j, 2] <- 0
  beta2[j, 2] <- 0
  }#end beta priors costs
  beta1[1, 1] ~ dnorm(0, 0.0000001)
  beta2[1, 1] ~ dnorm(0, 0.0000001)
  beta1[1, 2] <- sc
  beta2[1, 2] <- sc
  
  #standard deviation priors
  ls_c1[1] ~ dunif(-5, 10)         
  ls_e1[1] ~ dunif(-5, 10)
  ls_c2[1] ~ dunif(-5, 10)         
  ls_e2[1] ~ dunif(-5, 10)
  ls_c1[2] <- sdc         
  ls_e1[2] <- sde
  ls_c2[2] <- sdc         
  ls_e2[2] <- sde

  #correlation
  beta_f1[1] ~ dnorm(0, 0.0000001)
  beta_f2[1] ~ dnorm(0, 0.0000001)
  beta_f1[2] <- 0
  beta_f2[2] <- 0

  beta_f[1] <- beta_f1[1]
  beta_f[2] <- beta_f2[1]
  
  #priors on structural values mechanism
  for (j in 2:ze) {#begin gamma priors effects
  for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }
  }#end gamma priors effects
  gamma_e[1, 1] ~ dlogis(0, 1)
  gamma_e[1, 2] ~ dlogis(0, 1)
  
  for (j in 2:zc) {#begin gamma priors costs
  for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }
  }#end gamma priors costs
  gamma_c[1, 1] ~ dlogis(0, 1)
  gamma_c[1, 2] ~ dlogis(0, 1)
  
}
  "
   if(is.null(se) == FALSE & is.null(sc) == FALSE) {
     if(type == "SCAR") {
       model_string_jags <- gsub("logit(pq_1[i]) <- inprod(Z1_e[i, ], gamma_e[, 1])", "logit(pq_1[i]) <- gamma_e[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("logit(pc_1[i]) <- inprod(Z1_c[i, ], gamma_c[, 1])", "logit(pc_1[i]) <- gamma_c[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("logit(pq_2[i]) <- inprod(Z2_e[i, ], gamma_e[, 2])", "logit(pq_2[i]) <- gamma_e[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("logit(pc_2[i]) <- inprod(Z2_c[i, ], gamma_c[, 2])", "logit(pc_2[i]) <- gamma_c[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_c[1] <- ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", "p_c[1] <- ilogit(gamma_c[1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_c[2] <- ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", "p_c[2] <- ilogit(gamma_c[2])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_e[1] <- ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", "p_e[1] <- ilogit(gamma_e[1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_e[2] <- ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", "p_e[2] <- ilogit(gamma_e[2])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end gamma priors effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", "gamma_e[1] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", "gamma_e[2] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end gamma priors costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", "gamma_c[1] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", "gamma_c[2] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
     }
     if(dist_c == "norm") {
       model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
     } else if(dist_c == "gamma") {
      model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])", "cost1[i] ~ dgamma(mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs1", "tau_c1[i] <- mu_c1[i] / pow(s_c1[d_cost1[i] + 1], 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c1[i] <- ", "log(mu_c1[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])", "cost2[i] ~ dgamma(mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs2", "tau_c2[i] <- mu_c2[i] / pow(s_c2[d_cost2[i] + 1], 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c2[i] <- ", "log(mu_c2[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (t in 1:2) {#begin transformation costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c1[t] <- 1 / ss_c1[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c1[t] <- s_c1[t] * s_c1[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c1[t] <- exp(ls_c1[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c2[t] <- 1 / ss_c2[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c2[t] <- s_c2[t] * s_c2[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c2[t] <- exp(ls_c2[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end transformation costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", "nu_c[1] <- exp(inprod(mean_cov_c1[], beta1[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", "nu_c[2] <- exp(inprod(mean_cov_c2[], beta2[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10)", "s_c1[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10)", "s_c2[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c1[2] <- sdc", "s_c1[2] <- sdc", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c2[2] <- sdc", "s_c2[2] <- sdc", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]", "mu_c[1] <- nu_c[1] * (1 - p_c[1]) + exp(sc) * p_c[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]", "mu_c[2] <- nu_c[2] * (1 - p_c[2]) + exp(sc) * p_c[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
     } else if(dist_c == "lnorm") {
       model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])", "cost1[i] ~ dlnorm(lmu_c1[i], ltau_c1[d_cost1[i] + 1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c1[i] <- ", "lmu_c1[i] <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])", "cost2[i] ~ dlnorm(lmu_c2[i], ltau_c2[d_cost2[i] + 1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c2[i] <- ", "lmu_c2[i] <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c1[t] <- 1 / ss_c1[t]", "ltau_c1[t] <- 1 / lss_c1[t]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c1[t] <- s_c1[t] * s_c1[t]", "lss_c1[t] <- ls_c1[t] * ls_c1[t]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c1[t] <- exp(ls_c1[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c2[t] <- 1 / ss_c2[t]", "ltau_c2[t] <- 1 / lss_c2[t]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c2[t] <- s_c2[t] * s_c2[t]", "lss_c2[t] <- ls_c2[t] * ls_c2[t]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c2[t] <- exp(ls_c2[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", "lnu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", "lnu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]", "mu_c[1] <- nu_c[1] * (1 - p_c[1]) + exp(sc) * p_c[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]", "mu_c[2] <- nu_c[2] * (1 - p_c[2]) + exp(sc) * p_c[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm1", "nu_c[1] <- exp(lnu_c[1] + lss_c1[1] / 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm2", "nu_c[2] <- exp(lnu_c[2] + lss_c2[1] / 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10)", "ls_c1[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10)", "ls_c2[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c[1] <- s_c1[1]", "s_c[1] <- sqrt(exp(2 * lnu_c[1] + lss_c1[1]) * (exp(lss_c1[1]) - 1))",model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c[2] <- s_c2[1]", "s_c[2] <- sqrt(exp(2 * lnu_c[2] + lss_c2[1]) * (exp(lss_c2[1]) - 1))",model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c1[d_cost1[i] + 1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c2[d_cost2[i] + 1])", model_string_jags, fixed = TRUE)
     }
     if(dist_e == "norm") {
       model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
     } else if(dist_e == "beta") {
       model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e1[d_eff1[i] + 1])", "eff1[i] ~ dbeta(mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- (mu_e1[i] * (1 - mu_e1[i]) / pow(s_e1[d_eff1[i] + 1], 2) - 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e1[i] <- ","logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e2[d_eff2[i] + 1])", "eff2[i] ~ dbeta(mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- (mu_e2[i] * (1 - mu_e2[i]) / pow(s_e2[d_eff2[i] + 1], 2) - 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e2[i] <- ","logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (t in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_e1[t] <- 1 / ss_e1[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_e1[t] <- s_e1[t] * s_e1[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_e1[t] <- exp(ls_e1[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_e2[t] <- 1 / ss_e2[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_e2[t] <- s_e2[t] * s_e2[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_e2[t] <- exp(ls_e2[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_e[1] <- inprod(mean_cov_e1[], alpha1[, 1])", "nu_e[1] <- ilogit(inprod(mean_cov_e1[], alpha1[, 1]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_e[2] <- inprod(mean_cov_e2[], alpha2[, 1])", "nu_e[2] <- ilogit(inprod(mean_cov_e2[], alpha2[, 1]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e1[1] ~ dunif(-5, 10)", "s_e1[1] ~ dunif(0, sqrt(nu_e[1] * (1 - nu_e[1])))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e2[1] ~ dunif(-5, 10)", "s_e2[1] ~ dunif(0, sqrt(nu_e[2] * (1 - nu_e[2])))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e1[2] <- sde", "s_e1[2] <- sde", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e2[2] <- sde", "s_e2[2] <- sde", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e[1] <- nu_e[1] * (1 - p_e[1]) + se * p_e[1]", "mu_e[1] <- nu_e[1] * (1 - p_e[1]) + ilogit(se) * p_e[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e[2] <- nu_e[2] * (1 - p_e[2]) + se * p_e[2]", "mu_e[2] <- nu_e[2] * (1 - p_e[2]) + ilogit(se) * p_e[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1[d_eff1[i] + 1])", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2[d_eff2[i] + 1])", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
     }
     if(dist_e == "beta" & dist_c == "gamma") {
       model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
     if(pe == 1) {
         inprod_e1 <- "X1_e[i] * alpha1[d_eff1[i] + 1]"
         inprod_e2 <- "X2_e[i] * alpha2[d_eff2[i] + 1]"
         inprod_mean_e1 <- "mean_cov_e1 * alpha1[1]"
         inprod_mean_e2 <- "mean_cov_e2 * alpha2[1]"
         begin_beta_nons_e <- prior_beta_e1j <- prior_beta_e2j <- prior_beta_e10 <- prior_beta_e20 <- "#"
         beta_nons_e1 <- "alpha[1] <- alpha1[1]"
         beta_nons_e2 <- "alpha[2] <- alpha2[1]"
         begin_prior_beta <- "#begin alpha priors effects"
         prior_beta_e1 <- "alpha1[1] ~ dnorm(0, 0.0000001)"
         prior_beta_e2 <- "alpha2[1] ~ dnorm(0, 0.0000001)"
         prior_beta_e1s <- "alpha1[2] <- se"
         prior_beta_e2s <- "alpha2[2] <- se"
         end_prior_beta <- "#end alpha priors effects"
         model_string_jags <- gsub("inprod(X1_e[i, ], alpha1[, d_eff1[i] + 1])", inprod_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(X2_e[i, ], alpha2[, d_eff2[i] + 1])", inprod_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_e1[], alpha1[, 1])", inprod_mean_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_e2[], alpha2[, 1])", inprod_mean_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(j in 1:pe) {#alpha non-structural effects", begin_beta_nons_e, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha[j, 1] <- alpha1[j, 1]", beta_nons_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha[j, 2] <- alpha2[j, 1] }", beta_nons_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:pe) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_e2j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[j, 2] <- 0", prior_beta_e10, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[j, 2] <- 0", prior_beta_e20, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[1, 2] <- se", prior_beta_e1s, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[1, 2] <- se", prior_beta_e2s, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(pc == 1) {
         inprod_c1 <- "X1_c[i] * beta1[d_cost1[i] + 1]"
         inprod_c2 <- "X2_c[i] * beta2[d_cost2[i] + 1]"
         inprod_mean_c1 <- "mean_cov_c1 * beta1[1]"
         inprod_mean_c2 <- "mean_cov_c2 * beta2[1]"
         begin_beta_nons_c <- prior_beta_c1j <- prior_beta_c2j <- prior_beta_c10 <- prior_beta_c20 <- "#"
         beta_nons_c1 <- "beta[1] <- beta1[1]"
         beta_nons_c2 <- "beta[2] <- beta2[1]"
         begin_prior_beta <- "#begin beta priors costs"
         prior_beta_c1 <- "beta1[1] ~ dnorm(0, 0.0000001)"
         prior_beta_c2 <- "beta2[1] ~ dnorm(0, 0.0000001)"
         prior_beta_c1s <- "beta1[2] <- sc"
         prior_beta_c2s <- "beta2[2] <- sc"
         end_prior_beta <- "#end beta priors costs"
         model_string_jags <- gsub("inprod(X1_c[i, ], beta1[, d_cost1[i] + 1])", inprod_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(X2_c[i, ], beta2[, d_cost2[i] + 1])", inprod_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_c1[], beta1[, 1])", inprod_mean_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_c2[], beta2[, 1])", inprod_mean_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(j in 1:pc) {#beta non-structural costs", begin_beta_nons_c, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta[j, 1] <- beta1[j, 1]", beta_nons_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta[j, 2] <- beta2[j, 1] }", beta_nons_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_c2j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[j, 2] <- 0", prior_beta_c10, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[j, 2] <- 0", prior_beta_c20, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[1, 2] <- sc", prior_beta_c1s, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[1, 2] <- sc", prior_beta_c2s, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end beta priors costs", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(is.null(se) == FALSE) {
     if(ze == 1){
       if(type == "SAR"){
         inprod_e1 <- "Z1_e[i] * gamma_e[1]"
         inprod_e2 <- "Z2_e[i] * gamma_e[2]"
         inprod_mean_e1 <- "ilogit(mean_z_e1 * gamma_e[1])"
         inprod_mean_e2 <- "ilogit(mean_z_e2 * gamma_e[2])"
         begin_prior_gamma <- "#begin gamma priors effects"
         begin_prior_gamma2 <- "#"
         prior_gamma_e1 <- "gamma_e[1] ~ dlogis(0, 1)"
         prior_gamma_e2 <- "gamma_e[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors effects"
         model_string_jags <- gsub("inprod(Z1_e[i, ], gamma_e[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_e[i, ], gamma_e[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", inprod_mean_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", inprod_mean_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", prior_gamma_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", prior_gamma_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
     if(is.null(sc) == FALSE) {
     if(zc == 1) {
       if(type == "SAR") {
         inprod_c1 <- "Z1_c[i] * gamma_c[1]"
         inprod_c2 <- "Z2_c[i] * gamma_c[2]"
         inprod_mean_c1 <- "ilogit(mean_z_c1 * gamma_c[1])"
         inprod_mean_c2 <- "ilogit(mean_z_c2 * gamma_c[2])"
         begin_prior_gamma <- "#begin gamma priors costs"
         begin_prior_gamma2 <- "#"
         prior_gamma_c1 <- "gamma_c[1] ~ dlogis(0, 1)"
         prior_gamma_c2 <- "gamma_c[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors costs"
         model_string_jags <- gsub("inprod(Z1_c[i, ], gamma_c[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_c[i, ], gamma_c[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", inprod_mean_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", inprod_mean_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
   } else if(is.null(se) == TRUE & is.null(sc) == FALSE) {
     model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e1[d_eff1[i] + 1])", "eff1[i] ~ dnorm(mu_e1[i], tau_e1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_e1[i] <- inprod(X1_e[i, ], alpha1[, d_eff1[i] + 1])", "mu_e1[i] <- inprod(X1_e[i, ], alpha1[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("d_eff1[i] ~ dbern(pq_1[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("logit(pq_1[i]) <- inprod(Z1_e[i, ], gamma_e[, 1])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e2[d_eff2[i] + 1])", "eff2[i] ~ dnorm(mu_e2[i], tau_e2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_e2[i] <- inprod(X2_e[i, ], alpha2[, d_eff2[i] + 1])", "mu_e2[i] <- inprod(X2_e[i, ], alpha2[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("d_eff2[i] ~ dbern(pq_2[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("logit(pq_2[i]) <- inprod(Z2_e[i, ], gamma_e[, 2])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("tau_e1[t] <- 1 / ss_e1[t]", "tau_e1 <- 1 / ss_e1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ss_e1[t] <- s_e1[t] * s_e1[t]", "ss_e1 <- s_e1 * s_e1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_e1[t] <- exp(ls_e1[t])", "s_e1 <- exp(ls_e1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("tau_e2[t] <- 1 / ss_e2[t]", "tau_e2 <- 1 / ss_e2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ss_e2[t] <- s_e2[t] * s_e2[t]", "ss_e2 <- s_e2 * s_e2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_e2[t] <- exp(ls_e2[t])", "s_e2 <- exp(ls_e2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for (t in 1:2) {#begin transformation effects", "#begin transformation effects", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end transformation effects", "#end transformation effects", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("nu_e[1] <- inprod(mean_cov_e1[], alpha1[, 1])", "nu_e[1] <- inprod(mean_cov_e1[], alpha1[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("nu_e[2] <- inprod(mean_cov_e2[], alpha2[, 1])", "nu_e[2] <- inprod(mean_cov_e2[], alpha2[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha[j, 1] <- alpha1[j, 1]", "alpha[j, 1] <- alpha1[j]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha[j, 2] <- alpha2[j, 1] }", "alpha[j, 2] <- alpha2[j] }", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_e[1] <- nu_e[1] * (1 - p_e[1]) + se * p_e[1]", "mu_e[1] <- nu_e[1]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_e[2] <- nu_e[2] * (1 - p_e[2]) + se * p_e[2]", "mu_e[2] <- nu_e[2]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_e[1] <- s_e1[1]", "s_e[1] <- s_e1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_e[2] <- s_e2[1]", "s_e[2] <- s_e2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("p_e[1] <- ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("p_e[2] <- ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha1[j, 1] ~ dnorm(0, 0.0000001)", "alpha1[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha2[j, 1] ~ dnorm(0, 0.0000001)", "alpha2[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha1[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha2[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha1[1, 1] ~ dnorm(0, 0.0000001)", "alpha1[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha2[1, 1] ~ dnorm(0, 0.0000001)", "alpha2[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha1[1, 2] <- se", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("alpha2[1, 2] <- se", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_e1[1] ~ dunif(-5, 10)", "ls_e1 ~ dunif(-5, 10)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_e2[1] ~ dunif(-5, 10)", "ls_e2 ~ dunif(-5, 10)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_e1[2] <- sde", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_e2[2] <- sde", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end gamma priors effects", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1[d_eff1[i] + 1])", "loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2[d_eff2[i] + 1])", "loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_de1[i] <- logdensity.bern(d_eff1[i], pq_1[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_de2[i] <- logdensity.bern(d_eff2[i], pq_2[i])", "", model_string_jags, fixed = TRUE)
     if(type == "SCAR") {
       model_string_jags <- gsub("logit(pc_1[i]) <- inprod(Z1_c[i, ], gamma_c[, 1])", "logit(pc_1[i]) <- gamma_c[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("logit(pc_2[i]) <- inprod(Z2_c[i, ], gamma_c[, 2])", "logit(pc_2[i]) <- gamma_c[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_c[1] <- ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", "p_c[1] <- ilogit(gamma_c[1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_c[2] <- ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", "p_c[2] <- ilogit(gamma_c[2])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end gamma priors costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", "gamma_c[1] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", "gamma_c[2] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
     }
    if(dist_c == "norm") {
      model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
    } else if(dist_c == "gamma") {
      model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])", "cost1[i] ~ dgamma(mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs1", "tau_c1[i] <- mu_c1[i] / pow(s_c1[d_cost1[i] + 1], 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c1[i] <- ", "log(mu_c1[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])", "cost2[i] ~ dgamma(mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs2", "tau_c2[i] <- mu_c2[i] / pow(s_c2[d_cost2[i] + 1], 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c2[i] <- ", "log(mu_c2[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("for (t in 1:2) {#begin transformation costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c1[t] <- 1 / ss_c1[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c1[t] <- s_c1[t] * s_c1[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c1[t] <- exp(ls_c1[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c2[t] <- 1 / ss_c2[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c2[t] <- s_c2[t] * s_c2[t]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c2[t] <- exp(ls_c2[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("}#end transformation costs", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", "nu_c[1] <- exp(inprod(mean_cov_c1[], beta1[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", "nu_c[2] <- exp(inprod(mean_cov_c2[], beta2[, 1]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10)", "s_c1[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10)", "s_c2[1] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c1[2] <- sdc", "s_c1[2] <- sdc", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c2[2] <- sdc", "s_c2[2] <- sdc", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]", "mu_c[1] <- nu_c[1] * (1 - p_c[1]) + exp(sc) * p_c[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]", "mu_c[2] <- nu_c[2] * (1 - p_c[2]) + exp(sc) * p_c[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
    } else if(dist_c == "lnorm") {
      model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])", "cost1[i] ~ dlnorm(lmu_c1[i], ltau_c1[d_cost1[i] + 1])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c1[i] <- ", "lmu_c1[i] <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])", "cost2[i] ~ dlnorm(lmu_c2[i], ltau_c2[d_cost2[i] + 1])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c2[i] <- ", "lmu_c2[i] <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c1[t] <- 1 / ss_c1[t]", "ltau_c1[t] <- 1 / lss_c1[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c1[t] <- s_c1[t] * s_c1[t]", "lss_c1[t] <- ls_c1[t] * ls_c1[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c1[t] <- exp(ls_c1[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_c2[t] <- 1 / ss_c2[t]", "ltau_c2[t] <- 1 / lss_c2[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_c2[t] <- s_c2[t] * s_c2[t]", "lss_c2[t] <- ls_c2[t] * ls_c2[t]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c2[t] <- exp(ls_c2[t])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", "lnu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", "lnu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]", "mu_c[1] <- nu_c[1] * (1 - p_c[1]) + exp(sc) * p_c[1]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]", "mu_c[2] <- nu_c[2] * (1 - p_c[2]) + exp(sc) * p_c[2]", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm1", "nu_c[1] <- exp(lnu_c[1] + lss_c1[1] / 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#mean for lnorm2", "nu_c[2] <- exp(lnu_c[2] + lss_c2[1] / 2)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10)", "ls_c1[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10)", "ls_c2[1] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c[1] <- s_c1[1]", "s_c[1] <- sqrt(exp(2 * lnu_c[1] + lss_c1[1]) * (exp(lss_c1[1]) - 1))",model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_c[2] <- s_c2[1]", "s_c[2] <- sqrt(exp(2 * lnu_c[2] + lss_c2[1]) * (exp(lss_c2[1]) - 1))",model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c1[d_cost1[i] + 1])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c2[d_cost2[i] + 1])", model_string_jags, fixed = TRUE)
    }
    if(dist_e == "norm") {
      model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
    } else if(dist_e == "beta") {
      model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e1)", "eff1[i] ~ dbeta(mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- (mu_e1[i] * (1 - mu_e1[i]) / pow(s_e1, 2) - 1)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_e1[i] <- ","logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e2)", "eff2[i] ~ dbeta(mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- (mu_e2[i] * (1 - mu_e2[i]) / pow(s_e2, 2) - 1)", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("mu_e2[i] <- ","logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#begin transformation effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_e1 <- 1 / ss_e1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_e1 <- s_e1 * s_e1", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e1 <- exp(ls_e1)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_e2 <- 1 / ss_e2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ss_e2 <- s_e2 * s_e2", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e2 <- exp(ls_e2)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#end transformation effects", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_e[1] <- inprod(mean_cov_e1[], alpha1[])", "nu_e[1] <- ilogit(inprod(mean_cov_e1[], alpha1[]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("nu_e[2] <- inprod(mean_cov_e2[], alpha2[])", "nu_e[2] <- ilogit(inprod(mean_cov_e2[], alpha2[]))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e1 ~ dunif(-5, 10)", "s_e1 ~ dunif(0, sqrt(nu_e[1] * (1 - nu_e[1])))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("ls_e2 ~ dunif(-5, 10)", "s_e2 ~ dunif(0, sqrt(nu_e[2] * (1 - nu_e[2])))", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1)", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2)", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
    }
     if(dist_e == "beta" & dist_c == "gamma") {
       model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
     if(pe == 1) {
         inprod_e1 <- "X1_e[i] * alpha1[]"
         inprod_e2 <- "X2_e[i] * alpha2[]"
         inprod_mean_e1 <- "mean_cov_e1 * alpha1[1]"
         inprod_mean_e2 <- "mean_cov_e2 * alpha2[1]"
         begin_beta_nons_e <- prior_beta_e1j <- prior_beta_e2j <- prior_beta_e10 <- prior_beta_e20 <- "#"
         beta_nons_e1 <- "alpha[1] <- alpha1[1]"
         beta_nons_e2 <- "alpha[2] <- alpha2[1]"
         begin_prior_beta <- "#begin alpha priors effects"
         prior_beta_e1 <- "alpha1[1] ~ dnorm(0, 0.0000001)"
         prior_beta_e2 <- "alpha2[1] ~ dnorm(0, 0.0000001)"
         end_prior_beta <- "#end beta priors effects"
         model_string_jags <- gsub("inprod(X1_e[i, ], alpha1[])", inprod_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(X2_e[i, ], alpha2[])", inprod_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_e1[], alpha1[])", inprod_mean_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_e2[], alpha2[])", inprod_mean_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(j in 1:pe) {#alpha non-structural effects", begin_beta_nons_e, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha[j, 1] <- alpha1[j]", beta_nons_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha[j, 2] <- alpha2[j] }", beta_nons_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:pe) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[j] ~ dnorm(0, 0.0000001)", prior_beta_e1j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[j] ~ dnorm(0, 0.0000001)", prior_beta_e2j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[j] <- 0", prior_beta_e10, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[j] <- 0", prior_beta_e20, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha1[1] ~ dnorm(0, 0.0000001)", prior_beta_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("alpha2[1] ~ dnorm(0, 0.0000001)", prior_beta_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(pc == 1) {
       inprod_c1 <- "X1_c[i] * beta1[d_cost1[i] + 1]"
       inprod_c2 <- "X2_c[i] * beta2[d_cost2[i] + 1]"
       inprod_mean_c1 <- "mean_cov_c1 * beta1[1]"
       inprod_mean_c2 <- "mean_cov_c2 * beta2[1]"
       begin_beta_nons_c <- prior_beta_c1j <- prior_beta_c2j <- prior_beta_c10 <- prior_beta_c20 <- "#"
       beta_nons_c1 <- "beta[1] <- beta1[1]"
       beta_nons_c2 <- "beta[2] <- beta2[1]"
       begin_prior_beta <- "#begin beta priors costs"
       prior_beta_c1 <- "beta1[1] ~ dnorm(0, 0.0000001)"
       prior_beta_c2 <- "beta2[1] ~ dnorm(0, 0.0000001)"
       prior_beta_c1s <- "beta1[2] <- sc"
       prior_beta_c2s <- "beta2[2] <- sc"
       end_prior_beta <- "#end beta priors costs"
       model_string_jags <- gsub("inprod(X1_c[i, ], beta1[, d_cost1[i] + 1])", inprod_c1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X2_c[i, ], beta2[, d_cost2[i] + 1])", inprod_c2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_c1[], beta1[, 1])", inprod_mean_c1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_c2[], beta2[, 1])", inprod_mean_c2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pc) {#beta non-structural costs", begin_beta_nons_c, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta[j, 1] <- beta1[j, 1]", beta_nons_c1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta[j, 2] <- beta2[j, 1] }", beta_nons_c2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1j, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_c2j, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta1[j, 2] <- 0", prior_beta_c10, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta2[j, 2] <- 0", prior_beta_c20, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta1[1, 2] <- sc", prior_beta_c1s, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("beta2[1, 2] <- sc", prior_beta_c2s, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end beta priors costs", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(is.null(se) == FALSE) {
     if(ze == 1){
       if(type == "SAR"){
         inprod_e1 <- "Z1_e[i] * gamma_e[1]"
         inprod_e2 <- "Z2_e[i] * gamma_e[2]"
         inprod_mean_e1 <- "ilogit(mean_z_e1 * gamma_e[1])"
         inprod_mean_e2 <- "ilogit(mean_z_e2 * gamma_e[2])"
         begin_prior_gamma <- "#begin gamma priors effects"
         begin_prior_gamma2 <- "#"
         prior_gamma_e1 <- "gamma_e[1] ~ dlogis(0, 1)"
         prior_gamma_e2 <- "gamma_e[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors effects"
         model_string_jags <- gsub("inprod(Z1_e[i, ], gamma_e[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_e[i, ], gamma_e[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", inprod_mean_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", inprod_mean_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", prior_gamma_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", prior_gamma_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
     if(is.null(sc) == FALSE) {
     if(zc == 1) {
       if(type == "SAR") {
         inprod_c1 <- "Z1_c[i] * gamma_c[1]"
         inprod_c2 <- "Z2_c[i] * gamma_c[2]"
         inprod_mean_c1 <- "ilogit(mean_z_c1 * gamma_c[1])"
         inprod_mean_c2 <- "ilogit(mean_z_c2 * gamma_c[2])"
         begin_prior_gamma <- "#begin gamma priors costs"
         begin_prior_gamma2 <- "#"
         prior_gamma_c1 <- "gamma_c[1] ~ dlogis(0, 1)"
         prior_gamma_c2 <- "gamma_c[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors costs"
         model_string_jags <- gsub("inprod(Z1_c[i, ], gamma_c[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_c[i, ], gamma_c[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", inprod_mean_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", inprod_mean_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
   } else if(is.null(se) == FALSE & is.null(sc) == TRUE) {
     model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1[d_cost1[i] + 1])", "cost1[i] ~ dnorm(mu_c1[i], tau_c1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_c1[i] <- inprod(X1_c[i, ], beta1[, d_cost1[i] + 1]) + beta_f1[d_cost1[i] + 1] * (eff1[i] - mu_e[1])", "mu_c1[i] <- inprod(X1_c[i, ], beta1[]) + beta_f[1] * (eff1[i] - mu_e[1])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("d_cost1[i] ~ dbern(pc_1[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("logit(pc_1[i]) <- inprod(Z1_c[i, ], gamma_c[, 1])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2[d_cost2[i] + 1])", "cost2[i] ~ dnorm(mu_c2[i], tau_c2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_c2[i] <- inprod(X2_c[i, ], beta2[, d_cost2[i] + 1]) + beta_f2[d_cost2[i] + 1] * (eff2[i] - mu_e[2])", "mu_c2[i] <- inprod(X2_c[i, ], beta2[]) + beta_f[2] * (eff2[i] - mu_e[2])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("d_cost2[i] ~ dbern(pc_2[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("logit(pc_2[i]) <- inprod(Z2_c[i, ], gamma_c[, 2])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("tau_c1[t] <- 1 / ss_c1[t]", "tau_c1 <- 1 / ss_c1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ss_c1[t] <- s_c1[t] * s_c1[t]", "ss_c1 <- s_c1 * s_c1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_c1[t] <- exp(ls_c1[t])", "s_c1 <- exp(ls_c1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("tau_c2[t] <- 1 / ss_c2[t]", "tau_c2 <- 1 / ss_c2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ss_c2[t] <- s_c2[t] * s_c2[t]", "ss_c2 <- s_c2 * s_c2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_c2[t] <- exp(ls_c2[t])", "s_c2 <- exp(ls_c2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for (t in 1:2) {#begin transformation costs", "#begin transformation costs", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end transformation costs", "#end transformation costs", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[, 1])", "nu_c[1] <- inprod(mean_cov_c1[], beta1[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[, 1])", "nu_c[2] <- inprod(mean_cov_c2[], beta2[])", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta[j, 1] <- beta1[j, 1]", "beta[j, 1] <- beta1[j]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta[j, 2] <- beta2[j, 1] }", "beta[j, 2] <- beta2[j] }", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_c[1] <- nu_c[1] * (1 - p_c[1]) + sc * p_c[1]", "mu_c[1] <- nu_c[1]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("mu_c[2] <- nu_c[2] * (1 - p_c[2]) + sc * p_c[2]", "mu_c[2] <- nu_c[2]", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_c[1] <- s_c1[1]", "s_c[1] <- s_c1", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("s_c[2] <- s_c2[1]", "s_c[2] <- s_c2", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("p_c[1] <- ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("p_c[2] <- ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta1[j, 1] ~ dnorm(0, 0.0000001)", "beta1[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta2[j, 1] ~ dnorm(0, 0.0000001)", "beta2[j] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta1[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta2[j, 2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta1[1, 1] ~ dnorm(0, 0.0000001)", "beta1[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta2[1, 1] ~ dnorm(0, 0.0000001)", "beta2[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta1[1, 2] <- sc", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta2[1, 2] <- sc", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10)", "ls_c1 ~ dunif(-5, 10)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10)", "ls_c2 ~ dunif(-5, 10)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_c1[2] <- sdc", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("ls_c2[2] <- sdc", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("}#end gamma priors costs", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f1[1] ~ dnorm(0, 0.0000001)", "beta_f[1] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f2[1] ~ dnorm(0, 0.0000001)", "beta_f[2] ~ dnorm(0, 0.0000001)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f1[2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f2[2] <- 0", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f[1] <- beta_f1[1]", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("beta_f[2] <- beta_f2[1]", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1[d_cost1[i] + 1])", "loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2[d_cost2[i] + 1])", "loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2)", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_dc1[i] <- logdensity.bern(d_cost1[i], pc_1[i])", "", model_string_jags, fixed = TRUE)
     model_string_jags <- gsub("loglik_dc2[i] <- logdensity.bern(d_cost2[i], pc_2[i])", "", model_string_jags, fixed = TRUE)
     if(type == "SCAR") {
       model_string_jags <- gsub("logit(pq_1[i]) <- inprod(Z1_e[i, ], gamma_e[, 1])", "logit(pq_1[i]) <- gamma_e[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("logit(pq_2[i]) <- inprod(Z2_e[i, ], gamma_e[, 2])", "logit(pq_2[i]) <- gamma_e[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_e[1] <- ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", "p_e[1] <- ilogit(gamma_e[1])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("p_e[2] <- ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", "p_e[2] <- ilogit(gamma_e[2])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end gamma priors effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", "gamma_e[1] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", "gamma_e[2] ~ dlogis(0, 1)", model_string_jags, fixed = TRUE)
     }
     if(dist_c == "norm") {
       model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
     } else if(dist_c == "gamma") {
       model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1)", "cost1[i] ~ dgamma(mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs1", "tau_c1[i] <- mu_c1[i] / pow(s_c1, 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c1[i] <- ","log(mu_c1[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2)", "cost2[i] ~ dgamma(mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs2", "tau_c2[i] <- mu_c2[i] / pow(s_c2, 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c2[i] <- ","log(mu_c2[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#begin transformation costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c1 <- 1 / ss_c1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c1 <- s_c1 * s_c1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c1 <- exp(ls_c1)", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c2 <- 1 / ss_c2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c2 <- s_c2 * s_c2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c2 <- exp(ls_c2)", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#end transformation costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[])", "nu_c[1] <- exp(inprod(mean_cov_c1[], beta1[]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[])", "nu_c[2] <- exp(inprod(mean_cov_c2[], beta2[]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c1 ~ dunif(-5, 10)", "s_c1 ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c2 ~ dunif(-5, 10)", "s_c2 ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1)", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2)", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
     } else if(dist_c == "lnorm") {
       model_string_jags <- gsub("cost1[i] ~ dnorm(mu_c1[i], tau_c1)", "cost1[i] ~ dlnorm(lmu_c1[i], ltau_c1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c1[i] <- ","lmu_c1[i] <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("cost2[i] ~ dnorm(mu_c2[i], tau_c2)", "cost2[i] ~ dlnorm(lmu_c2[i], ltau_c2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std costs2", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c2[i] <- ","lmu_c2[i] <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#begin transformation costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c1 <- 1 / ss_c1", "ltau_c1 <- 1 / lss_c1", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c1 <- s_c1 * s_c1", "lss_c1 <- ls_c1 * ls_c1", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c1 <- exp(ls_c1)", "s_c1 <- sqrt(exp(2 * lnu_c[1] + lss_c1) * (exp(lss_c1) - 1))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_c2 <- 1 / ss_c2", "ltau_c2 <- 1 / lss_c2", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_c2 <- s_c2 * s_c2", "lss_c2 <- ls_c2 * ls_c2", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_c2 <- exp(ls_c2)", "s_c2 <- sqrt(exp(2 * lnu_c[2] + lss_c2) * (exp(lss_c2) - 1))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#end transformation costs", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[1] <- inprod(mean_cov_c1[], beta1[])", "lnu_c[1] <- inprod(mean_cov_c1[], beta1[])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_c[2] <- inprod(mean_cov_c2[], beta2[])", "lnu_c[2] <- inprod(mean_cov_c2[], beta2[])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm1", "mu_c[1] <- exp(lnu_c[1] + lss_c1 / 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#mean for lnorm2", "mu_c[2] <- exp(lnu_c[2] + lss_c2 / 2)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c1 ~ dunif(-5, 10)", "ls_c1 ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_c2 ~ dunif(-5, 10)", "ls_c2 ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c[1] <- nu_c[1]","", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_c[2] <- nu_c[2]","", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c1)", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c2)", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c2)", model_string_jags, fixed = TRUE)
     }
     if(dist_e == "norm") {
       model_string_jags <- gsub("#derive mean and std effects1", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects2", "", model_string_jags, fixed = TRUE)
     } else if(dist_e == "beta") {
       model_string_jags <- gsub("eff1[i] ~ dnorm(mu_e1[i], tau_e1[d_eff1[i] + 1])", "eff1[i] ~ dbeta(mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects1", "tau_e1[i] <- (mu_e1[i] * (1 - mu_e1[i]) / pow(s_e1[d_eff1[i] + 1], 2) - 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e1[i] <- ","logit(mu_e1[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("eff2[i] ~ dnorm(mu_e2[i], tau_e2[d_eff2[i] + 1])", "eff2[i] ~ dbeta(mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("#derive mean and std effects2", "tau_e2[i] <- (mu_e2[i] * (1 - mu_e2[i]) / pow(s_e2[d_eff2[i] + 1], 2) - 1)", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e2[i] <- ","logit(mu_e2[i]) <- ", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (t in 1:2) {#begin transformation effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_e1[t] <- 1 / ss_e1[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_e1[t] <- s_e1[t] * s_e1[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_e1[t] <- exp(ls_e1[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("tau_e2[t] <- 1 / ss_e2[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ss_e2[t] <- s_e2[t] * s_e2[t]", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("s_e2[t] <- exp(ls_e2[t])", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end transformation effects", "", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_e[1] <- inprod(mean_cov_e1[], alpha1[, 1])", "nu_e[1] <- ilogit(inprod(mean_cov_e1[], alpha1[, 1]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("nu_e[2] <- inprod(mean_cov_e2[], alpha2[, 1])", "nu_e[2] <- ilogit(inprod(mean_cov_e2[], alpha2[, 1]))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e1[1] ~ dunif(-5, 10)", "s_e1[1] ~ dunif(0, sqrt(nu_e[1] * (1 - nu_e[1])))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e2[1] ~ dunif(-5, 10)", "s_e2[1] ~ dunif(0, sqrt(nu_e[2] * (1 - nu_e[2])))", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e1[2] <- sde", "s_e1[2] <- sde", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("ls_e2[2] <- sde", "s_e2[2] <- sde", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e[1] <- nu_e[1] * (1 - p_e[1]) + se * p_e[1]", "mu_e[1] <- nu_e[1] * (1 - p_e[1]) + ilogit(se) * p_e[1]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("mu_e[2] <- nu_e[2] * (1 - p_e[2]) + se * p_e[2]", "mu_e[2] <- nu_e[2] * (1 - p_e[2]) + ilogit(se) * p_e[2]", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e1[d_eff1[i] + 1])", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e2[d_eff2[i] + 1])", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)       
     }
     if(dist_e == "beta" & dist_c == "gamma") {
       model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
     }
     if(pe == 1) {
       inprod_e1 <- "X1_e[i] * alpha1[d_eff1[i] + 1]"
       inprod_e2 <- "X2_e[i] * alpha2[d_eff2[i] + 1]"
       inprod_mean_e1 <- "mean_cov_e1 * alpha1[1]"
       inprod_mean_e2 <- "mean_cov_e2 * alpha2[1]"
       begin_beta_nons_e <- prior_beta_e1j <- prior_beta_e2j <- prior_beta_e10 <- prior_beta_e20 <- "#"
       beta_nons_e1 <- "alpha[1] <- alpha1[1]"
       beta_nons_e2 <- "alpha[2] <- alpha2[1]"
       begin_prior_beta <- "#begin alpha priors effects"
       prior_beta_e1 <- "alpha1[1] ~ dnorm(0, 0.0000001)"
       prior_beta_e2 <- "alpha2[1] ~ dnorm(0, 0.0000001)"
       prior_beta_e1s <- "alpha1[2] <- se"
       prior_beta_e2s <- "alpha2[2] <- se"
       end_prior_beta <- "#end alpha priors effects"
       model_string_jags <- gsub("inprod(X1_e[i, ], alpha1[, d_eff1[i] + 1])", inprod_e1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(X2_e[i, ], alpha2[, d_eff2[i] + 1])", inprod_e2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e1[], alpha1[, 1])", inprod_mean_e1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("inprod(mean_cov_e2[], alpha2[, 1])", inprod_mean_e2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for(j in 1:pe) {#alpha non-structural effects", begin_beta_nons_e, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha[j, 1] <- alpha1[j, 1]", beta_nons_e1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha[j, 2] <- alpha2[j, 1] }", beta_nons_e2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("for (j in 2:pe) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha1[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1j, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha2[j, 1] ~ dnorm(0, 0.0000001)", prior_beta_e2j, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha1[j, 2] <- 0", prior_beta_e10, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha2[j, 2] <- 0", prior_beta_e20, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha1[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha2[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e2, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha1[1, 2] <- se", prior_beta_e1s, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("alpha2[1, 2] <- se", prior_beta_e2s, model_string_jags, fixed = TRUE)
       model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(pc == 1) {
         inprod_c1 <- "X1_c[i] * beta1[]"
         inprod_c2 <- "X2_c[i] * beta2[]"
         inprod_mean_c1 <- "mean_cov_c1 * beta1[1]"
         inprod_mean_c2 <- "mean_cov_c2 * beta2[1]"
         begin_beta_nons_c <- prior_beta_c1j <- prior_beta_c2j <- prior_beta_c10 <- prior_beta_c20 <- "#"
         beta_nons_c1 <- "beta[1] <- beta1[1]"
         beta_nons_c2 <- "beta[2] <- beta2[1]"
         begin_prior_beta <- "#begin beta priors costs"
         prior_beta_c1 <- "beta1[1] ~ dnorm(0, 0.0000001)"
         prior_beta_c2 <- "beta2[1] ~ dnorm(0, 0.0000001)"
         end_prior_beta <- "#end beta priors costs"
         model_string_jags <- gsub("inprod(X1_c[i, ], beta1[])", inprod_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(X2_c[i, ], beta2[])", inprod_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_c1[], beta1[])", inprod_mean_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(mean_cov_c2[], beta2[])", inprod_mean_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(j in 1:pc) {#beta non-structural costs", begin_beta_nons_c, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta[j, 1] <- beta1[j]", beta_nons_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta[j, 2] <- beta2[j] }", beta_nons_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[j] ~ dnorm(0, 0.0000001)", prior_beta_c1j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[j] ~ dnorm(0, 0.0000001)", prior_beta_c2j, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[j] <- 0", prior_beta_c10, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[j] <- 0", prior_beta_c20, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta1[1] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("beta2[1] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end beta priors costs", end_prior_beta, model_string_jags, fixed = TRUE)
     }
     if(is.null(se) == FALSE) {
     if(ze == 1){
       if(type == "SAR"){
         inprod_e1 <- "Z1_e[i] * gamma_e[1]"
         inprod_e2 <- "Z2_e[i] * gamma_e[2]"
         inprod_mean_e1 <- "ilogit(mean_z_e1 * gamma_e[1])"
         inprod_mean_e2 <- "ilogit(mean_z_e2 * gamma_e[2])"
         begin_prior_gamma <- "#begin gamma priors effects"
         begin_prior_gamma2 <- "#"
         prior_gamma_e1 <- "gamma_e[1] ~ dlogis(0, 1)"
         prior_gamma_e2 <- "gamma_e[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors effects"
         model_string_jags <- gsub("inprod(Z1_e[i, ], gamma_e[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_e[i, ], gamma_e[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", inprod_mean_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", inprod_mean_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", prior_gamma_e1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", prior_gamma_e2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
     if(is.null(sc) == FALSE) {
     if(zc == 1) {
       if(type == "SAR") {
         inprod_c1 <- "Z1_c[i] * gamma_c[1]"
         inprod_c2 <- "Z2_c[i] * gamma_c[2]"
         inprod_mean_c1 <- "ilogit(mean_z_c1 * gamma_c[1])"
         inprod_mean_c2 <- "ilogit(mean_z_c2 * gamma_c[2])"
         begin_prior_gamma <- "#begin gamma priors costs"
         begin_prior_gamma2 <- "#"
         prior_gamma_c1 <- "gamma_c[1] ~ dlogis(0, 1)"
         prior_gamma_c2 <- "gamma_c[2] ~ dlogis(0, 1)"
         end_prior_gamma <- "#end gamma priors costs"
         model_string_jags <- gsub("inprod(Z1_c[i, ], gamma_c[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("inprod(Z2_c[i, ], gamma_c[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", inprod_mean_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", inprod_mean_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
         model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
       }
      }
     }
    }
  if(ind == TRUE) {
    if(is.null(sc) == FALSE) {
    model_string_jags <- gsub(" + beta_f1[d_cost1[i] + 1] * (eff1[i] - mu_e[1])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub(" + beta_f2[d_cost2[i] + 1] * (eff2[i] - mu_e[2])", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f1[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f2[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f1[2] <- 0", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f2[2] <- 0", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[1] <- beta_f1[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("beta_f[2] <- beta_f2[1]", "", model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
    } else if(is.null(sc) == TRUE) {
      model_string_jags <- gsub(" + beta_f[1] * (eff1[i] - mu_e[1])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub(" + beta_f[2] * (eff2[i] - mu_e[2])", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_f[1] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_f[2] ~ dnorm(0, 0.0000001)", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("#correlation", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_f[1] <- beta_f1[1]", "", model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_f[2] <- beta_f2[1]", "", model_string_jags, fixed = TRUE)
    }
  }
  model_string_jags <- prior_hurdle(type = type, dist_e = dist_e, dist_c = dist_c, pe = pe, pc = pc, ze = ze, zc = zc, se = se, sc = sc)
  writeLines(model_string_jags, "hurdle.txt")
  model_string <- "hurdle.txt"
  return(model_string)
}))