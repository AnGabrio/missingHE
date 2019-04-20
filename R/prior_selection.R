#' An internal function to change the hyperprior parameters in the selection model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of selection model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm') or Beta ('beta').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param pe Number of covariates for the effectiveness model
#' @param pc Number of covariates for the cost model
#' @param ze Number of covariates or the missingness indicators model for the effectiveness
#' @param zc Number of covariates or the missingness indicators model for the costs
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_selection <- function(type, dist_e, dist_c, pe, pc, ze, zc) eval.parent( substitute( {
  if(type == "MNAR" | type == "MNAR_eff" | type == "MNAR_cost") {
      if(is.null(delta.prior.e) == FALSE & grepl("delta_e[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(delta.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_deltae <- delta.prior.e
        prior_deltae_str <- paste("delta_e[t] ~ dnorm(", prior_deltae[1], ",", prior_deltae[2])
        model_string_jags <- gsub("delta_e[t] ~ dnorm(0, 1", prior_deltae_str, model_string_jags, fixed = TRUE) }
    if(is.null(delta.prior.c) == FALSE & grepl("delta_c[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(delta.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_deltac <- delta.prior.c
      prior_deltac_str <- paste("delta_c[t] ~ dnorm(", prior_deltac[1], ",", prior_deltac[2])
      model_string_jags <- gsub("delta_c[t] ~ dnorm(0, 1", prior_deltac_str, model_string_jags, fixed = TRUE) } 
  }
    if(ze == 1) {
      if(is.null(gamma0.prior.e) == FALSE & grepl("gamma_e[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_e[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(gamma0.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_pe <- gamma0.prior.e
        prior_pe_str <- paste("gamma_e[1] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
        model_string_jags <- gsub("gamma_e[1] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE)
        prior_pe_str <- paste("gamma_e[2] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
        model_string_jags <- gsub("gamma_e[2] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE) }
        } else if(ze > 1) {
          if(is.null(gamma0.prior.e) == FALSE & grepl("gamma_e[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_e[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.e) != 2) {stop("provide correct hyper prior values") }
            prior_pe <- gamma0.prior.e
            prior_pe_str <- paste("gamma_e[1, 1] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
            model_string_jags <- gsub("gamma_e[1, 1] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE)
            prior_pe_str <- paste("gamma_e[1, 2] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
            model_string_jags <- gsub("gamma_e[1, 2] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE) }
          if(is.null(gamma.prior.e) == FALSE & grepl("gamma_e[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma.prior.e) != 2) {stop("provide correct hyper prior values") }
            prior_gammae <- gamma.prior.e
            prior_gammae_str <- paste("gamma_e[j, t] ~ dnorm(", prior_gammae[1], ",", prior_gammae[2])
            model_string_jags <- gsub("gamma_e[j, t] ~ dnorm(0, 0.01", prior_gammae_str, model_string_jags, fixed = TRUE) }
        }
      if(zc == 1) {
          if(is.null(gamma0.prior.c) == FALSE & grepl("gamma_c[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_c[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_pc <- gamma0.prior.c
            prior_pc_str <- paste("gamma_c[1] ~ dlogis(",prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE)
            prior_pc_str <- paste("gamma_c[2] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[2] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE) }
        } else if(zc >1 ) {
          if(is.null(gamma0.prior.c) == FALSE & grepl("gamma_c[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_c[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_pc <- gamma0.prior.c
            prior_pc_str <- paste("gamma_c[1, 1] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1, 1] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE)
            prior_pc_str <- paste("gamma_c[1, 2] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1, 2] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE) }
          if(is.null(gamma.prior.c) == FALSE & grepl("gamma_c[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma.prior.c) != 2){stop("provide correct hyper prior values") }
            prior_gammac <- gamma.prior.c
            prior_gammac_str <- paste("gamma_c[j, t] ~ dnorm(", prior_gammac[1], ",", prior_gammac[2])
            model_string_jags <- gsub("gamma_c[j, t] ~ dnorm(0, 0.01", prior_gammac_str, model_string_jags, fixed = TRUE) }
        }
    if(pe == 1) {
      if(is.null(alpha0.prior) == FALSE & grepl("alpha[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_mue <- alpha0.prior
        prior_mue_str <- paste("alpha[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
        model_string_jags <- gsub("alpha[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE)
        prior_mue_str <- paste("alpha[2] ~ dnorm(", prior_mue[1],",", prior_mue[2])
        model_string_jags <- gsub("alpha[2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE) }
    } else if(pe > 1){
      if(is.null(alpha0.prior) == FALSE & grepl("alpha[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_mue <- alpha0.prior
        prior_mue_str <- paste("alpha[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
        model_string_jags <- gsub("alpha[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE)
        prior_mue_str <- paste("alpha[1, 2] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
        model_string_jags <- gsub("alpha[1, 2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE) }
        if(is.null(alpha.prior) == FALSE & grepl("alpha[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha.prior) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- alpha.prior
        prior_alphae_str <- paste("alpha[j, t] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("alpha[j, t] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
    }
    if(pc == 1) {
      if(is.null(beta0.prior) == FALSE & grepl("beta[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muc <- beta0.prior
        prior_muc_str <- paste("beta[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
        prior_muc_str <- paste("beta[2] ~ dnorm(",prior_muc[1],",", prior_muc[2])
        model_string_jags <- gsub("beta[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
    } else if(pc > 1) {
      if(is.null(beta0.prior) == FALSE & grepl("beta[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muc <- beta0.prior
        prior_muc_str <- paste("beta[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
        prior_muc_str <- paste("beta[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta.prior) != 2){stop("provide correct hyper prior values") }
        prior_betac <- beta.prior
        prior_betac_str <- paste("beta[j, t] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
        model_string_jags <- gsub("beta[j, t] ~ dnorm(0, 0.01", prior_betac_str, model_string_jags, fixed = TRUE) }
    }
  if(dist_e == "norm") {
  if(is.null(sigma.prior.e) == FALSE & grepl("ls_e[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
    if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
    prior_alphae <- sigma.prior.e
    prior_alphae_str <- paste("ls_e[t] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
    model_string_jags <- gsub("ls_e[t] ~ dunif(-5, 10", prior_alphae_str, model_string_jags, fixed = TRUE) }
  } else if(dist_e == "beta") {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e[t] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e[t] ~ dunif(0, sqrt(mu_e[t] * (1 - mu_e[t]))", prior_alphae_str, model_string_jags, fixed = TRUE) }
  }
  if(dist_c == "norm") {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c[t] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c[t] ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed = TRUE) }
  } else if(dist_c == "lnorm") {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c[t] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c[t] ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed = TRUE) }
  } else if(dist_c == "gamma") {
      if(is.null(sigma.prior.c) == FALSE & grepl("s_c[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("s_c[t] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c[t] ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE) }
  }
      if(exists("beta_f.prior") == TRUE) {
      if(is.null(beta_f.prior) == FALSE & grepl("beta_f", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta_f.prior) != 2) {stop("provide correct hyper prior values") }
          prior_beta_f <- beta_f.prior
          prior_beta_f_str <- paste("beta_f[t] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
          model_string_jags <- gsub("beta_f[t] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE)
        }
      }
  model_string_prior <- model_string_jags
    return(model_string_prior)
}))