#'An internal function to select which type of selection model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Selection models
#' @param dist_e Distribution assumed for the effects. Current available choices are: Normal ('norm') or Beta ('beta')
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param ind Logical; if TRUE independence between effectiveness and costs is assumed, else correlation is accounted for
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param pe Number of covariates for the effectiveness model
#' @param pc Number of covariates for the cost model
#' @param ze Number of covariates or the missingness indicators model for the effectiveness
#' @param zc Number of covariates or the missingness indicators model for the costs
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_selection <- function(dist_e , dist_c, ind, type, pe, pc, ze, zc) eval.parent(substitute( {
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
  mu_c1[i] <- inprod(X1_c[i, ], beta[, 1]) + beta_f[1] * (eff1[i] - mu_e[1])
  mu_e1[i] <- inprod(X1_e[i, ], alpha[, 1])   

  #missing data mechanism
  m_eff1[i] ~ dbern(pq_1[i])
  logit(pq_1[i]) <- inprod(Z1_e[i, ], gamma_e[, 1]) + delta_e[1] * eff1[i]
  m_cost1[i] ~ dbern(pc_1[i])
  logit(pc_1[i]) <- inprod(Z1_c[i, ], gamma_c[, 1]) + delta_c[1] * cost1[i]

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
  mu_c2[i] <- inprod(X2_c[i, ], beta[, 2]) + beta_f[2] * (eff2[i] - mu_e[2])
  mu_e2[i] <- inprod(X2_e[i, ], alpha[, 2])

  #missing data mechanism
  m_eff2[i] ~ dbern(pq_2[i])
  logit(pq_2[i]) <- inprod(Z2_e[i, ], gamma_e[, 2]) + delta_e[2] * eff2[i]
  m_cost2[i] ~ dbern(pc_2[i])
  logit(pc_2[i]) <- inprod(Z2_c[i, ], gamma_c[, 2]) + delta_c[2] * cost2[i]

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
  
  #missingness probability
  p_c[1] <- ilogit(inprod(mean_z_c1[], gamma_c[, 1]) + delta_c[1] * mean(cost1[]))
  p_c[2] <- ilogit(inprod(mean_z_c2[], gamma_c[, 2]) + delta_c[2] * mean(cost2[]))
  p_e[1] <- ilogit(inprod(mean_z_e1[], gamma_e[, 1]) + delta_e[1] * mean(eff1[]))
  p_e[2] <- ilogit(inprod(mean_z_e2[], gamma_e[, 2]) + delta_e[2] * mean(eff2[]))
  
  #calculate means at mean of covariates
  mu_c[1] <- inprod(mean_cov_c1[], beta[, 1])
  mu_c[2] <- inprod(mean_cov_c2[], beta[, 2])
  mu_e[1] <- inprod(mean_cov_e1[], alpha[, 1])
  mu_e[2] <- inprod(mean_cov_e2[], alpha[, 2])
  
  #priors
  
  #priors for mean regression coefficients
  for (j in 2:pe) {#begin alpha priors effects
  for(t in 1:2) {alpha[j, t] ~ dnorm(0, 0.0000001) }
  }#end alpha priors effects
  alpha[1, 1] ~ dnorm(0, 0.0000001)
  alpha[1, 2] ~ dnorm(0, 0.0000001)
  
  for (j in 2:pc) {#begin beta priors costs
  for(t in 1:2) {beta[j, t] ~ dnorm(0, 0.0000001) }
  }#end beta priors costs
  beta[1, 1] ~ dnorm(0, 0.0000001)
  beta[1, 2] ~ dnorm(0, 0.0000001)
  
  #standard deviation priors
  for(t in 1:2) {
  ls_c[t] ~ dunif(-5, 10)
  ls_e[t] ~ dunif(-5, 10)

  #correlation
  beta_f[t] ~ dnorm(0, 0.0000001)
  }
  
  #priors on missing data mechanism
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
  
  #mnar parameters
  for(t in 1:2) {# begin mnar priors
  delta_e[t] ~ dnorm(0, 1)
  delta_c[t] ~ dnorm(0, 1)
  }#end mnar priors
  
}
 "
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
 } else if(type == "MNAR_eff") {
  model_string_jags <- gsub(" + delta_c[1] * mean(cost1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * mean(cost2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[1] * cost1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_c[2] * cost2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_c[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
 } else if(type == "MNAR_cost") {
  model_string_jags <- gsub(" + delta_e[1] * mean(eff1[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * mean(eff2[])", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[1] * eff1[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub(" + delta_e[2] * eff2[i]", "", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("delta_e[t] ~ dnorm(0, 1)", "", model_string_jags, fixed = TRUE)
 }  
 if(ind == TRUE) {
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
  model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1[], beta[, 1])", "mu_c[1] <- exp(inprod(mean_cov_c1[], beta[, 1]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2[], beta[, 2])", "mu_c[2] <- exp(inprod(mean_cov_c2[], beta[, 2]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_c[t] ~ dunif(-5, 10)", "s_c[t] ~ dunif(0, 10000)", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c[1])", "loglik_c1[i] <- logdensity.gamma(cost1[i], mu_c1[i] * tau_c1[i], tau_c1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c[2])", "loglik_c2[i] <- logdensity.gamma(cost2[i], mu_c2[i] * tau_c2[i], tau_c2[i])", model_string_jags, fixed = TRUE)
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
   model_string_jags <- gsub("mu_c[1] <- inprod(mean_cov_c1[], beta[, 1])", "lmu_c[1] <- inprod(mean_cov_c1[], beta[, 1])", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("mu_c[2] <- inprod(mean_cov_c2[], beta[, 2])", "lmu_c[2] <- inprod(mean_cov_c2[], beta[, 2])", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("ls_c[t] ~ dunif(-5, 10)", "ls_c[t] ~ dunif(0, 100)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#mean for lnorm", "mu_c[t] <- exp(lmu_c[t] + lss_c[t] / 2)", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("loglik_c1[i] <- logdensity.norm(cost1[i], mu_c1[i], tau_c[1])", "loglik_c1[i] <- logdensity.lnorm(cost1[i], lmu_c1[i], ltau_c[1])", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("loglik_c2[i] <- logdensity.norm(cost2[i], mu_c2[i], tau_c[2])", "loglik_c2[i] <- logdensity.lnorm(cost2[i], lmu_c2[i], ltau_c[2])", model_string_jags, fixed = TRUE)
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
  model_string_jags <- gsub("mu_e[1] <- inprod(mean_cov_e1[], alpha[, 1])", "mu_e[1] <- ilogit(inprod(mean_cov_e1[], alpha[, 1]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("mu_e[2] <- inprod(mean_cov_e2[], alpha[, 2])", "mu_e[2] <- ilogit(inprod(mean_cov_e2[], alpha[, 2]))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10)", "s_e[t] ~ dunif(0, sqrt(mu_e[t] * (1 - mu_e[t])))", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e1[i] <- logdensity.norm(eff1[i], mu_e1[i], tau_e[1])", "loglik_e1[i] <- logdensity.beta(eff1[i], mu_e1[i] * tau_e1[i], (1 - mu_e1[i]) * tau_e1[i])", model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("loglik_e2[i] <- logdensity.norm(eff2[i], mu_e2[i], tau_e[2])", "loglik_e2[i] <- logdensity.beta(eff2[i], mu_e2[i] * tau_e2[i], (1 - mu_e2[i]) * tau_e2[i])", model_string_jags, fixed = TRUE)
 }
 if(dist_e == "beta" & dist_c == "gamma") {
   model_string_jags <- gsub("for (t in 1:2) {#begin transformation", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("}# end transformation", "", model_string_jags, fixed = TRUE)
   model_string_jags <- gsub("#transformation of parameters", "", model_string_jags, fixed = TRUE)
 }
 if(pe == 1) {
  inprod_e1 <- "X1_e[i] * alpha[1]"
  inprod_e2 <- "X2_e[i] * alpha[2]"
  inprod_mean_e1 <- "mean_cov_e1 * alpha[1]"
  inprod_mean_e2 <- "mean_cov_e2 * alpha[2]"
  begin_prior_beta <- "#begin alpha priors effects"
  prior_beta <- "#"
  end_prior_beta <- "#end alpha priors effects"
  prior_beta_e1 <- "alpha[1] ~ dnorm(0, 0.0000001)"
  prior_beta_e2 <- "alpha[2] ~ dnorm(0, 0.0000001)"
  model_string_jags <- gsub("inprod(X1_e[i, ], alpha[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(X2_e[i, ], alpha[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(mean_cov_e1[], alpha[, 1])", inprod_mean_e1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(mean_cov_e2[], alpha[, 2])", inprod_mean_e2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 2:pe) {#begin alpha priors effects", begin_prior_beta, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {alpha[j, t] ~ dnorm(0, 0.0000001) }",prior_beta, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end alpha priors effects", end_prior_beta, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_e1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("alpha[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_e2, model_string_jags, fixed = TRUE)
 }
 if(pc == 1) {
  inprod_c1 <- "X1_c[i] * beta[1]"
  inprod_c2 <- "X2_c[i] * beta[2]"
  inprod_mean_c1 <- "mean_cov_c1 * beta[1]"
  inprod_mean_c2 <- "mean_cov_c2 * beta[2]"
  begin_prior_beta <- "#begin beta priors costs"
  prior_beta <- "#"
  end_prior_beta <- "#end beta priors costs"
  prior_beta_c1 <- "beta[1] ~ dnorm(0, 0.0000001)"
  prior_beta_c2 <- "beta[2] ~ dnorm(0, 0.0000001)"
  model_string_jags <- gsub("inprod(X1_c[i, ], beta[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(X2_c[i, ], beta[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(mean_cov_c1[], beta[, 1])", inprod_mean_c1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(mean_cov_c2[], beta[, 2])", inprod_mean_c2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {beta[j, t] ~ dnorm(0, 0.0000001) }", prior_beta, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end beta priors costs", end_prior_beta,model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[1, 1] ~ dnorm(0, 0.0000001)", prior_beta_c1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("beta[1, 2] ~ dnorm(0, 0.0000001)", prior_beta_c2, model_string_jags, fixed = TRUE)
 }  
 if(ze == 1) {
  inprod_e1 <- "Z1_e[i] * gamma_e[1]"
  inprod_e2 <- "Z2_e[i] * gamma_e[2]"
  if(type == "MAR" | type == "MNAR_cost") {
    inprod_mean_e1 <- "ilogit(mean_z_e1 * gamma_e[1])"
    inprod_mean_e2 <- "ilogit(mean_z_e2 * gamma_e[2])"
    model_string_jags <- gsub("ilogit(inprod(mean_z_e1[], gamma_e[, 1]))", inprod_mean_e1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ilogit(inprod(mean_z_e2[], gamma_e[, 2]))", inprod_mean_e2, model_string_jags, fixed = TRUE)
  }
  if(type == "MNAR_eff" | type == "MNAR") {
    inprod_mean_e1 <- "ilogit(mean_z_e1 * gamma_e[1] + delta_e[1] * mean(eff1[]))"
    inprod_mean_e2 <- "ilogit(mean_z_e2 * gamma_e[2] + delta_e[2] * mean(eff2[]))"
    model_string_jags <- gsub("ilogit(inprod(mean_z_e1[], gamma_e[, 1]) + delta_e[1] * mean(eff1[]))", inprod_mean_e1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ilogit(inprod(mean_z_e2[], gamma_e[, 2]) + delta_e[2] * mean(eff2[]))", inprod_mean_e2, model_string_jags, fixed = TRUE)
  }
  begin_prior_gamma <- "#begin gamma priors effects"
  begin_prior_gamma2 <- "#"
  prior_gamma_e1 <- "gamma_e[1] ~ dlogis(0, 1)"
  prior_gamma_e2 <- "gamma_e[2] ~ dlogis(0, 1)"
  end_prior_gamma <- "#end gamma priors effects"
  model_string_jags <- gsub("inprod(Z1_e[i, ], gamma_e[, 1])", inprod_e1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(Z2_e[i, ], gamma_e[, 2])", inprod_e2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {gamma_e[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1)", prior_gamma_e1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1)", prior_gamma_e2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags, fixed = TRUE)
 }  
 if(zc == 1) {
  inprod_c1 <- "Z1_c[i] * gamma_c[1]"
  inprod_c2 <- "Z2_c[i] * gamma_c[2]"
  if(type == "MAR" | type == "MNAR_eff") {
    inprod_mean_c1 <- "ilogit(mean_z_c1 * gamma_c[1])"
    inprod_mean_c2 <- "ilogit(mean_z_c2 * gamma_c[2])"
    model_string_jags <- gsub("ilogit(inprod(mean_z_c1[], gamma_c[, 1]))", inprod_mean_c1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ilogit(inprod(mean_z_c2[], gamma_c[, 2]))", inprod_mean_c2, model_string_jags, fixed = TRUE)
  }
  if(type == "MNAR_cost" | type == "MNAR") {
    inprod_mean_c1 <- "ilogit(mean_z_c1 * gamma_c[1] + delta_c[1] * mean(cost1[]))"
    inprod_mean_c2 <- "ilogit(mean_z_c2 * gamma_c[2] + delta_c[2] * mean(cost2[]))"
    model_string_jags <- gsub("ilogit(inprod(mean_z_c1[], gamma_c[, 1]) + delta_c[1] * mean(cost1[]))", inprod_mean_c1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("ilogit(inprod(mean_z_c2[], gamma_c[, 2]) + delta_c[2] * mean(cost2[]))", inprod_mean_c2, model_string_jags, fixed = TRUE)
  }
  begin_prior_gamma <- "#begin gamma priors costs"
  begin_prior_gamma2 <- "#"
  prior_gamma_c1 <- "gamma_c[1] ~ dlogis(0, 1)"
  prior_gamma_c2 <- "gamma_c[2] ~ dlogis(0, 1)"
  end_prior_gamma <- "#end gamma priors costs"
  model_string_jags <- gsub("inprod(Z1_c[i, ], gamma_c[, 1])", inprod_c1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("inprod(Z2_c[i, ], gamma_c[, 2])", inprod_c2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for (j in 2:zc) {#begin gamma priors costs", begin_prior_gamma, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("for(t in 1:2) {gamma_c[j, t] ~ dnorm(0, 0.01) }", begin_prior_gamma2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1)", prior_gamma_c1, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1)", prior_gamma_c2, model_string_jags, fixed = TRUE)
  model_string_jags <- gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags, fixed = TRUE)
 }
 model_string_jags <- prior_selection(type = type, dist_e = dist_e, dist_c = dist_c, pe = pe, pc = pc, ze = ze, zc = zc)
 writeLines(model_string_jags, "selection.txt")
 model_string <- "selection.txt"
 return(model_string)
}))