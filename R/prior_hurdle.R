#' An internal function to change the hyperprior parameters in the hurdle model provided by the user depending on the type of
#' structural value mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of hurdle model selected according 
#' to the type of structural value mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions hurdle models
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#' and Structural At Random (SAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param pe_fixed Number of fixed effects for the effectiveness model
#' @param pc_fixed Number of fixed effects for the cost model
#' @param ze_fixed Number of fixed effects or the structural indicators model for the effectiveness
#' @param zc_fixed Number of fixed effects or the structural indicators model for the costs
#' @param pe_random Number of random effects for the effectiveness model
#' @param pc_random Number of random effects for the cost model
#' @param ze_random Number of random effects or the structural indicators model for the effectiveness
#' @param zc_random Number of random effects or the structural indicators model for the costs
#' @param model_e_random Random effects formula for the effectiveness model
#' @param model_c_random Random effects formula for the costs model
#' @param model_se_random Random effects formula for the structural indicators model for the effectiveness
#' @param model_sc_random Random effects formula for the structural indicators model for the costs
#' @param se Structural value for the effectiveness 
#' @param sc Structural value for the costs
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_hurdle <- function(type, dist_e, dist_c, pe_fixed, pc_fixed , ze_fixed, zc_fixed, model_e_random, model_c_random, 
                         model_se_random, model_sc_random, pe_random, pc_random, ze_random, zc_random, se, sc) eval.parent( substitute( {
  if(is.null(se) == TRUE) {
      if(is.null(gamma0.prior.e) == FALSE | is.null(gamma.prior.e) == FALSE | is.null(mu.g.prior.e) == FALSE | is.null(s.g.prior.e) == FALSE |
         is.null(mu.g0.prior.e) == FALSE | is.null(s.g0.prior.e) == FALSE) {
        stop("priors for structural effects cannot be provided if 'se' is null")
      }
  }
  if(is.null(sc) == TRUE) {
      if(is.null(gamma0.prior.c) == FALSE | is.null(gamma.prior.c) == FALSE | is.null(mu.g.prior.c) == FALSE | is.null(s.g.prior.c) == FALSE |
         is.null(mu.g0.prior.c) == FALSE | is.null(s.g0.prior.c) == FALSE) {
        stop("priors for structural costs cannot be provided if 'sc' is null")
      }
  }                       
  if(ze_fixed == 1) {
      if(is.null(gamma0.prior.e) == FALSE & grepl("gamma_e[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_e[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
         if(length(gamma0.prior.e) != 2) {stop("provide correct hyper prior values") }
         prior_pe <- gamma0.prior.e
         prior_pe_str <- paste("gamma_e[1] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
         model_string_jags <- gsub("gamma_e[1] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE)
         prior_pe_str <- paste("gamma_e[2] ~ dlogis(", prior_pe[1], ",", prior_pe[2])
         model_string_jags <- gsub("gamma_e[2] ~ dlogis(0, 1", prior_pe_str, model_string_jags, fixed = TRUE) }
  } else if(ze_fixed > 1) {
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
  if(length(model_se_random) != 0 & ze_random == 1) {
       if(is.null(mu.g0.prior.e) == FALSE & grepl("mu_g_e_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
         if(length(mu.g0.prior.e) != 2) {stop("provide correct hyper prior values") }
         prior_g0e <- mu.g0.prior.e
         prior_g0e_str <- paste("mu_g_e_hat[t] ~ dnorm(", prior_g0e[1], ",", prior_g0e[2])
         model_string_jags <- gsub("mu_g_e_hat[t] ~ dnorm(0, 0.001", prior_g0e_str, model_string_jags, fixed = TRUE) }
           if(is.null(s.g0.prior.e) == FALSE & grepl("s_g_e_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             if(length(s.g0.prior.e) != 2) {stop("provide correct hyper prior values") }
             prior_g0e <- s.g0.prior.e
             prior_g0e_str <- paste("s_g_e_hat[t] ~ dunif(", prior_g0e[1], ",", prior_g0e[2])
             model_string_jags <- gsub("s_g_e_hat[t] ~ dunif(0, 100", prior_g0e_str, model_string_jags, fixed = TRUE) }
  } else if(length(model_se_random) != 0 & ze_random > 1) {
           if(is.null(mu.g.prior.e) == FALSE & grepl("mu_g_e_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             if(length(mu.g.prior.e) != 2) {stop("provide correct hyper prior values") }
             prior_ge <- mu.g.prior.e
             prior_ge_str <- paste("mu_g_e_hat[j, t] ~ dnorm(", prior_ge[1], ",", prior_ge[2])
             model_string_jags <- gsub("mu_g_e_hat[j, t] ~ dnorm(0, 0.001", prior_ge_str, model_string_jags, fixed = TRUE) }
              if(is.null(s.g.prior.e) == FALSE & grepl("s_g_e_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.g.prior.e) != 2) {stop("provide correct hyper prior values") }
               prior_ge <- s.g.prior.e
               prior_ge_str <- paste("s_g_e_hat[j, t] ~ dunif(", prior_ge[1], ",", prior_ge[2])
               model_string_jags <- gsub("s_g_e_hat[j, t] ~ dunif(0, 100", prior_ge_str, model_string_jags, fixed = TRUE) }
  }
  if(zc_fixed == 1) {
         if(is.null(gamma0.prior.c) == FALSE & grepl("gamma_c[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_c[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
           if(length(gamma0.prior.c) != 2) {stop("provide correct hyper prior values") }
           prior_pc <- gamma0.prior.c
           prior_pc_str <- paste("gamma_c[1] ~ dlogis(",prior_pc[1], ",", prior_pc[2])
           model_string_jags <- gsub("gamma_c[1] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE)
           prior_pc_str <- paste("gamma_c[2] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
           model_string_jags <- gsub("gamma_c[2] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE) }
  } else if(zc_fixed >1 ) {
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
  if(length(model_sc_random) != 0 & zc_random == 1) {
         if(is.null(mu.g0.prior.c) == FALSE & grepl("mu_g_c_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
           if(length(mu.g0.prior.c) != 2) {stop("provide correct hyper prior values") }
           prior_g0c <- mu.g0.prior.c
           prior_g0c_str <- paste("mu_g_c_hat[t] ~ dnorm(", prior_g0c[1], ",", prior_g0c[2])
           model_string_jags <- gsub("mu_g_c_hat[t] ~ dnorm(0, 0.001", prior_g0c_str, model_string_jags, fixed = TRUE) }
            if(is.null(s.g0.prior.c) == FALSE & grepl("s_g_c_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
              if(length(s.g0.prior.c) != 2) {stop("provide correct hyper prior values") }
                 prior_g0c <- s.g0.prior.c
                 prior_g0c_str <- paste("s_g_c_hat[t] ~ dunif(", prior_g0c[1], ",", prior_g0c[2])
                 model_string_jags <- gsub("s_g_c_hat[t] ~ dunif(0, 100", prior_g0c_str, model_string_jags, fixed = TRUE) }
  } else if(length(model_sc_random) != 0 & zc_random > 1) {
         if(is.null(mu.g.prior.c) == FALSE & grepl("mu_g_c_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
           if(length(mu.g.prior.c) != 2) {stop("provide correct hyper prior values") }
             prior_gc <- mu.g.prior.c
             prior_gc_str <- paste("mu_g_c_hat[j, t] ~ dnorm(", prior_gc[1], ",", prior_gc[2])
             model_string_jags <- gsub("mu_g_c_hat[j, t] ~ dnorm(0, 0.001", prior_gc_str, model_string_jags, fixed = TRUE) }
              if(is.null(s.g.prior.c) == FALSE & grepl("s_g_c_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.g.prior.c) != 2) {stop("provide correct hyper prior values") }
               prior_gc <- s.g.prior.c
               prior_gc_str <- paste("s_g_c_hat[j, t] ~ dunif(", prior_gc[1], ",", prior_gc[2])
               model_string_jags <- gsub("s_g_c_hat[j, t] ~ dunif(0, 100", prior_gc_str, model_string_jags, fixed = TRUE) }
  }                            
  if(pe_fixed == 1) {
         if(is.null(alpha0.prior) == FALSE & grepl("alpha1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
           if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
           prior_mue <- alpha0.prior
           prior_mue_str <- paste("alpha1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
           model_string_jags <- gsub("alpha1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE)
           prior_mue_str <- paste("alpha2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
           model_string_jags <- gsub("alpha2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE) }
  } else if(pe_fixed > 1){
         if(is.null(alpha0.prior) == FALSE & grepl("alpha1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
           if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
           prior_mue <- alpha0.prior
           prior_mue_str <- paste("alpha1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
           model_string_jags <- gsub("alpha1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE)
           prior_mue_str <- paste("alpha2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
           model_string_jags <- gsub("alpha2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags, fixed = TRUE) }
            if(is.null(alpha.prior) == FALSE & grepl("alpha1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
              if(length(alpha.prior) != 2) {stop("provide correct hyper prior values") }
              prior_alphae <- alpha.prior
              prior_alphae_str <- paste("alpha1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
              model_string_jags <- gsub("alpha1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE)
              prior_alphae_str <- paste("alpha2[j, 1] ~ dnorm(", prior_alphae[1],",", prior_alphae[2])
              model_string_jags <- gsub("alpha2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
  }
  if(length(model_e_random) != 0 & pe_random == 1) {
          if(is.null(mu.a0.prior) == FALSE & grepl("mu_a_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.a0.prior) != 2) {stop("provide correct hyper prior values") }
            prior_a0 <- mu.a0.prior
            prior_a0_str <- paste("mu_a_hat[t] ~ dnorm(", prior_a0[1], ",", prior_a0[2])
            model_string_jags <- gsub("mu_a_hat[t] ~ dnorm(0, 0.001", prior_a0_str, model_string_jags, fixed = TRUE) }
             if(is.null(s.a0.prior) == FALSE & grepl("s_a_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.a0.prior) != 2) {stop("provide correct hyper prior values") }
               prior_a0 <- s.a0.prior
               prior_a0_str <- paste("s_a_hat[t] ~ dunif(", prior_a0[1], ",", prior_a0[2])
               model_string_jags <- gsub("s_a_hat[t] ~ dunif(0, 100", prior_a0_str, model_string_jags, fixed = TRUE) }
  } else if(length(model_e_random) != 0 & pe_random > 1) {
          if(is.null(mu.a.prior) == FALSE & grepl("mu_a_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.a.prior) != 2) {stop("provide correct hyper prior values") }
            prior_a <- mu.a.prior
            prior_a_str <- paste("mu_a_hat[j, t] ~ dnorm(", prior_a[1], ",", prior_a[2])
            model_string_jags <- gsub("mu_a_hat[j, t] ~ dnorm(0, 0.001", prior_a_str, model_string_jags, fixed = TRUE) }
             if(is.null(s.a.prior) == FALSE & grepl("s_a_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.a.prior) != 2) {stop("provide correct hyper prior values") }
               prior_a <- s.a.prior
               prior_a_str <- paste("s_a_hat[j, t] ~ dunif(", prior_a[1], ",", prior_a[2])
               model_string_jags <- gsub("s_a_hat[j, t] ~ dunif(0, 100", prior_a_str, model_string_jags, fixed = TRUE) }
  }
  if(pc_fixed == 1) {
           if(is.null(beta0.prior) == FALSE & grepl("beta1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
             prior_muc <- beta0.prior
             prior_muc_str <- paste("beta1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
             model_string_jags <- gsub("beta1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
             prior_muc_str <- paste("beta2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
             model_string_jags <- gsub("beta2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
  } else if(pc_fixed > 1) {
          if(is.null(beta.prior) == FALSE & grepl("beta1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(beta.prior) != 2) {stop("provide correct hyper prior values") }
            prior_muc <- beta.prior
            prior_muc_str <- paste("beta1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
            model_string_jags <- gsub("beta1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
            prior_muc_str <- paste("beta2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
            model_string_jags <- gsub("beta2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
          if(is.null(beta0.prior) == FALSE & grepl("beta1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
            prior_muc <- beta0.prior
            prior_muc_str <- paste("beta1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
            model_string_jags <- gsub("beta1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
            prior_muc_str <- paste("beta2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
            model_string_jags <- gsub("beta2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
             if(is.null(beta.prior) == FALSE & grepl("beta1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(beta.prior) != 2){stop("provide correct hyper prior values") }
               prior_betac <- beta.prior
               prior_betac_str <- paste("beta1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
               model_string_jags <- gsub("beta1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE)
               prior_betac_str <- paste("beta2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
               model_string_jags <- gsub("beta2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
  }
  if(length(model_c_random) != 0 & pc_random == 1) {
          if(is.null(mu.b0.prior) == FALSE & grepl("mu_b_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.b0.prior) != 2) {stop("provide correct hyper prior values") }
            prior_b0 <- mu.b0.prior
            prior_b0_str <- paste("mu_b_hat[t] ~ dnorm(", prior_b0[1], ",", prior_b0[2])
            model_string_jags <- gsub("mu_b_hat[t] ~ dnorm(0, 0.001", prior_b0_str, model_string_jags, fixed = TRUE) }
             if(is.null(s.b0.prior) == FALSE & grepl("s_b_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.b0.prior) != 2) {stop("provide correct hyper prior values") }
               prior_b0 <- s.b0.prior
               prior_b0_str <- paste("s_b_hat[t] ~ dunif(", prior_b0[1], ",", prior_b0[2])
               model_string_jags <- gsub("s_b_hat[t] ~ dunif(0, 100", prior_b0_str, model_string_jags, fixed = TRUE) }
  } else if(length(model_c_random) != 0 & pc_random > 1) {
          if(is.null(mu.b.prior) == FALSE & grepl("mu_b_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.b.prior) != 2) {stop("provide correct hyper prior values") }
            prior_b <- mu.b.prior
            prior_b_str <- paste("mu_b_hat[j, t] ~ dnorm(", prior_b[1], ",", prior_b[2])
            model_string_jags <- gsub("mu_b_hat[j, t] ~ dnorm(0, 0.001", prior_b_str, model_string_jags, fixed = TRUE) }
             if(is.null(s.b.prior) == FALSE & grepl("s_b_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               if(length(s.b.prior) != 2) {stop("provide correct hyper prior values") }
               prior_b <- s.b.prior
               prior_b_str <- paste("s_b_hat[j, t] ~ dunif(", prior_b[1], ",", prior_b[2])
               model_string_jags <- gsub("s_b_hat[j, t] ~ dunif(0, 100", prior_b_str, model_string_jags, fixed = TRUE) }
  }  
  if(dist_e == "norm") {
  if(is.null(se) == FALSE) {
  if(is.null(sigma.prior.e) == FALSE & grepl("ls_e1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_e2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
    if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
    prior_alphae <- sigma.prior.e
    prior_alphae_str <- paste("ls_e1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
    model_string_jags <- gsub("ls_e1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags, fixed = TRUE)
    prior_alphae_str <- paste("ls_e2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
    model_string_jags <- gsub("ls_e2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
  } else if(is.null(se) == TRUE) {
    if(is.null(sigma.prior.e) == FALSE & grepl("ls_e1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_e2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e
      prior_alphae_str <- paste("ls_e1 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
      model_string_jags <- gsub("ls_e1 ~ dunif(-5, 10", prior_alphae_str, model_string_jags, fixed = TRUE)
      prior_alphae_str <- paste("ls_e2 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
      model_string_jags <- gsub("ls_e2 ~ dunif(-5, 10", prior_alphae_str, model_string_jags, fixed = TRUE) }
     }
  } else if(dist_e == "beta") {
     if(is.null(se) == FALSE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1[1] ~ dunif(0, sqrt(nu_e[1] * (1 - nu_e[1]))", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2[1] ~ dunif(0, sqrt(nu_e[2] * (1 - nu_e[2]))", prior_alphae_str, model_string_jags, fixed = TRUE) }
    } else if(is.null(se) == TRUE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1 ~ dunif(0, sqrt(nu_e[1] * (1 - nu_e[1]))", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2 ~ dunif(0, sqrt(nu_e[2] * (1 - nu_e[2]))", prior_alphae_str, model_string_jags, fixed = TRUE) }
    }
  } else if(dist_e == "gamma" | dist_e == "logis") {
    if(is.null(se) == FALSE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags, fixed = TRUE) }
    } else if(is.null(se) == TRUE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1 ~ dunif(0, 10000", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2 ~ dunif(0, 10000", prior_alphae_str, model_string_jags, fixed = TRUE) }
    }
  } else if(dist_e == "exp" | dist_e == "bern" | dist_e == "pois") {
        if(is.null(sigma.prior.e) == FALSE) {
          stop("no prior for sigma required for the effects under the 'exp', 'bern', 'pois' distributions") }
  } else if(dist_e == "weibull") {
    if(is.null(se) == FALSE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE) }
    } else if(is.null(se) == TRUE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("s_e1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_e2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("s_e1 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e1 ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("s_e2 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("s_e2 ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE) }
    }
  } else if(dist_e == "nbinom") {
    if(is.null(se) == FALSE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("tau_e1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("tau_e2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("tau_e1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("tau_e1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("tau_e2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("tau_e2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE) }
    } else if(is.null(se) == TRUE) {
      if(is.null(sigma.prior.e) == FALSE & grepl("tau_e1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("tau_e2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
        prior_alphae <- sigma.prior.e
        prior_alphae_str <- paste("tau_e1 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("tau_e1 ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE)
        prior_alphae_str <- paste("tau_e2 ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
        model_string_jags <- gsub("tau_e2 ~ dunif(0, 100", prior_alphae_str, model_string_jags, fixed = TRUE) }
    }
  }
  if(dist_c == "norm") {
  if(is.null(sc) == FALSE) {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_c2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed=TRUE)
      prior_alphac_str <- paste("ls_c2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed=TRUE) }
  } else if(is.null(sc) == TRUE) {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_c2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c1 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c1 ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed = TRUE)
      prior_alphac_str <- paste("ls_c2 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c2 ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed = TRUE) }
    }
  } else if(dist_c == "lnorm") {
    if(is.null(sc) == FALSE) {
      if(is.null(sigma.prior.c) == FALSE & grepl("ls_c1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_c2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("ls_c1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("ls_c1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed=TRUE)
        prior_alphac_str <- paste("ls_c2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("ls_c2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed=TRUE) }
    } else if(is.null(sc) == TRUE) {
      if(is.null(sigma.prior.c) == FALSE & grepl("ls_c1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("ls_c2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("ls_c1 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("ls_c1 ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed = TRUE)
        prior_alphac_str <- paste("ls_c2 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("ls_c2 ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed = TRUE) }
    }
  } else if(dist_c == "gamma") {
    if(is.null(sc) == FALSE) {
      if(is.null(sigma.prior.c) == FALSE & grepl("s_c1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_c2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("s_c1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE)
        prior_alphac_str <- paste("s_c2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE) }
    } else if(is.null(sc) == TRUE) {
      if(is.null(sigma.prior.c) == FALSE & grepl("s_c1 ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("s_c2 ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("s_c1 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c1 ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE)
        prior_alphac_str <- paste("s_c2 ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c2 ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE) }
    }
  }
    if(exists("beta_f.prior") == TRUE) {
      if(is.null(beta_f.prior) == FALSE & grepl("beta_f", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta_f.prior) != 2) {stop("provide correct hyper prior values") }
        if(is.null(se) == FALSE & is.null(sc) == TRUE) {
          prior_beta_f <- beta_f.prior
          prior_beta_f_str <- paste("beta_f[1] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
          model_string_jags <- gsub("beta_f[1] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE)
          prior_beta_f_str <- paste("beta_f[2] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
          model_string_jags <- gsub("beta_f[2] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE)}
        if(is.null(sc) == FALSE) {
          prior_beta_f <- beta_f.prior
          prior_beta_f_str <- paste("beta_f1[1] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
          model_string_jags <- gsub("beta_f1[1] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE)
          prior_beta_f_str <- paste("beta_f2[1] ~ dnorm(", prior_beta_f[1], ",",prior_beta_f[2])
          model_string_jags <- gsub("beta_f2[1] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE) }
        }
    }
    if(exists("mu.b_f.prior") == TRUE) {
     if(is.null(mu.b_f.prior) == FALSE & grepl("mu_b_f_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(mu.b_f.prior) != 2) {stop("provide correct hyper prior values") }
      prior_b_f <- mu.b_f.prior
      prior_b_f_str <- paste("mu_b_f_hat[t] ~ dnorm(", prior_b_f[1], ",", prior_b_f[2])
      model_string_jags <- gsub("mu_b_f_hat[t] ~ dnorm(0, 0.001", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.b_f.prior") == TRUE) {
     if(is.null(s.b_f.prior) == FALSE & grepl("s_b_f_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(s.b_f.prior) != 2) {stop("provide correct hyper prior values") }
      prior_b_f <- s.b_f.prior
      prior_b_f_str <- paste("s_b_f_hat[t] ~ dunif(", prior_b_f[1], ",", prior_b_f[2])
      model_string_jags <- gsub("s_b_f_hat[t] ~ dunif(0, 100", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    }
  model_string_prior <- model_string_jags
    return(model_string_prior)
}))