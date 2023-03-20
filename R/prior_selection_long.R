#' An internal function to change the hyperprior parameters in the selection model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of selection model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_u distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
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
#' #no examples
#' #
#' #

prior_selection_long <- function(type, dist_u, dist_c, pu_fixed, pc_fixed , zu_fixed, zc_fixed, model_u_random, model_c_random, 
                            model_mu_random, model_mc_random, pu_random, pc_random, zu_random, zc_random) eval.parent( substitute( {
  if(type == "MNAR" | type == "MNAR_eff" | type == "MNAR_cost") {
      if(is.null(delta.prior.u) == FALSE & grepl("delta_u[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(delta.prior.u) != 2) {stop("provide correct hyper prior values") }
        prior_deltau <- delta.prior.u
        prior_deltau_str <- paste("delta_u[t, time, mdrop] ~ dnorm(", prior_deltau[1], ",", prior_deltau[2])
        model_string_jags <- gsub("delta_u[t, time, mdrop] ~ dnorm(0, 1", prior_deltau_str, model_string_jags, fixed = TRUE) }
    if(is.null(delta.prior.c) == FALSE & grepl("delta_c[t, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(delta.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_deltac <- delta.prior.c
      prior_deltac_str <- paste("delta_c[t, time, mdrop] ~ dnorm(", prior_deltac[1], ",", prior_deltac[2])
      model_string_jags <- gsub("delta_c[t, time, mdrop] ~ dnorm(0, 1", prior_deltac_str, model_string_jags, fixed = TRUE) } 
    if(is.null(mu.d.prior.u) == FALSE & grepl("mu_d_u_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(mu.d.prior.u) != 2) {stop("provide correct hyper prior values") }
      prior_mudu <- mu.d.prior.u
      prior_mudu_str <- paste("mu_d_u_hat[t, time] ~ dnorm(", prior_mudu[1], ",", prior_mudu[2])
      model_string_jags <- gsub("mu_d_u_hat[t, time] ~ dnorm(0, 1", prior_mudu_str, model_string_jags, fixed = TRUE) }
    if(is.null(s.d.prior.u) == FALSE & grepl("s_d_u_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(s.d.prior.u) != 2) {stop("provide correct hyper prior values") }
      prior_sdu <- s.d.prior.u
      prior_sdu_str <- paste("s_d_u_hat[t, time] ~ dunif(", prior_sdu[1], ",", prior_sdu[2])
      model_string_jags <- gsub("s_d_u_hat[t, time] ~ dunif(0, 1", prior_sdu_str, model_string_jags, fixed = TRUE) }
    if(is.null(mu.d.prior.c) == FALSE & grepl("mu_d_c_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(mu.d.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_mudc <- mu.d.prior.c
      prior_mudc_str <- paste("mu_d_c_hat[t, time] ~ dnorm(", prior_mudc[1], ",", prior_mudc[2])
      model_string_jags <- gsub("mu_d_c_hat[t, time] ~ dnorm(0, 1", prior_mudc_str, model_string_jags, fixed = TRUE) }
    if(is.null(s.d.prior.c) == FALSE & grepl("s_d_c_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(s.d.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_sdc <- s.d.prior.c
      prior_sdc_str <- paste("s_d_c_hat[t, time] ~ dunif(", prior_sdc[1], ",", prior_sdc[2])
      model_string_jags <- gsub("s_d_c_hat[t, time] ~ dunif(0, 1", prior_sdc_str, model_string_jags, fixed = TRUE) }
  }
    if(zu_fixed == 1) {
      if(is.null(gamma0.prior.u) == FALSE & grepl("gamma_u[1, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_u[2, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(gamma0.prior.u) != 2) {stop("provide correct hyper prior values") }
        prior_pu <- gamma0.prior.u
        prior_pu_str <- paste("gamma_u[1, time, mdrop] ~ dlogis(", prior_pu[1], ",", prior_pu[2])
        model_string_jags <- gsub("gamma_u[1, time, mdrop] ~ dlogis(0, 1", prior_pu_str, model_string_jags, fixed = TRUE)
        prior_pu_str <- paste("gamma_u[2, time, mdrop] ~ dlogis(", prior_pu[1], ",", prior_pu[2])
        model_string_jags <- gsub("gamma_u[2, time, mdrop] ~ dlogis(0, 1", prior_pu_str, model_string_jags, fixed = TRUE) }
        } else if(zu_fixed > 1) {
          if(is.null(gamma0.prior.u) == FALSE & grepl("gamma_u[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_u[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_pu <- gamma0.prior.u
            prior_pu_str <- paste("gamma_u[1, 1] ~ dlogis(", prior_pu[1], ",", prior_pu[2])
            model_string_jags <- gsub("gamma_u[1, 1] ~ dlogis(0, 1", prior_pu_str, model_string_jags, fixed = TRUE)
            prior_pu_str <- paste("gamma_u[1, 2] ~ dlogis(", prior_pu[1], ",", prior_pu[2])
            model_string_jags <- gsub("gamma_u[1, 2] ~ dlogis(0, 1", prior_pu_str, model_string_jags, fixed = TRUE) }
          if(is.null(gamma.prior.u) == FALSE & grepl("gamma_u[j, t, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_gammau <- gamma.prior.u
            prior_gammau_str <- paste("gamma_u[j, t, time, mdrop] ~ dnorm(", prior_gammau[1], ",", prior_gammau[2])
            model_string_jags <- gsub("gamma_u[j, t, time, mdrop] ~ dnorm(0, 0.01", prior_gammau_str, model_string_jags, fixed = TRUE) }
       }
        if(length(model_mu_random) != 0 & zu_random == 1) {
          if(is.null(mu.g0.prior.u) == FALSE & grepl("mu_g_u_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.g0.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_g0u <- mu.g0.prior.u
            prior_g0u_str <- paste("mu_g_u_hat[t, time] ~ dnorm(", prior_g0u[1], ",", prior_g0u[2])
            model_string_jags <- gsub("mu_g_u_hat[t, time] ~ dnorm(0, 0.001", prior_g0u_str, model_string_jags, fixed = TRUE) }
          if(is.null(s.g0.prior.u) == FALSE & grepl("s_g_u_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(s.g0.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_g0u <- s.g0.prior.u
            prior_g0u_str <- paste("s_g_u_hat[t, time] ~ dunif(", prior_g0u[1], ",", prior_g0u[2])
            model_string_jags <- gsub("s_g_u_hat[t, time] ~ dunif(0, 100", prior_g0u_str, model_string_jags, fixed = TRUE) }
            } else if(length(model_mu_random) != 0 & zu_random > 1) {
          if(is.null(mu.g.prior.u) == FALSE & grepl("mu_g_u_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.g.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_gu <- mu.g.prior.u
            prior_gu_str <- paste("mu_g_u_hat[j, t, time] ~ dnorm(", prior_gu[1], ",", prior_gu[2])
            model_string_jags <- gsub("mu_g_u_hat[j, t, time] ~ dnorm(0, 0.001", prior_gu_str, model_string_jags, fixed = TRUE) }
          if(is.null(s.g.prior.u) == FALSE & grepl("s_g_u_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(s.g.prior.u) != 2) {stop("provide correct hyper prior values") }
            prior_gu <- s.g.prior.u
            prior_gu_str <- paste("s_g_u_hat[j, t, time] ~ dunif(", prior_gu[1], ",", prior_gu[2])
            model_string_jags <- gsub("s_g_u_hat[j, t, time] ~ dunif(0, 100", prior_gu_str, model_string_jags, fixed = TRUE) }
       }
      if(zc_fixed == 1) {
          if(is.null(gamma0.prior.c) == FALSE & grepl("gamma_c[1, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_c[2, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_pc <- gamma0.prior.c
            prior_pc_str <- paste("gamma_c[1, time, mdrop] ~ dlogis(",prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1, time, mdrop] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE)
            prior_pc_str <- paste("gamma_c[2, time, mdrop] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[2, time, mdrop] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE) }
        } else if(zc_fixed >1 ) {
          if(is.null(gamma0.prior.c) == FALSE & grepl("gamma_c[1, 1, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("gamma_c[1, 2, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_pc <- gamma0.prior.c
            prior_pc_str <- paste("gamma_c[1, 1, time, mdrop] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1, 1, time, mdrop] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE)
            prior_pc_str <- paste("gamma_c[1, 2, time, mdrop] ~ dlogis(", prior_pc[1], ",", prior_pc[2])
            model_string_jags <- gsub("gamma_c[1, 2, time, mdrop] ~ dlogis(0, 1", prior_pc_str, model_string_jags, fixed = TRUE) }
          if(is.null(gamma.prior.c) == FALSE & grepl("gamma_c[j, t, time, mdrop] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(gamma.prior.c) != 2){stop("provide correct hyper prior values") }
            prior_gammac <- gamma.prior.c
            prior_gammac_str <- paste("gamma_c[j, t, time, mdrop] ~ dnorm(", prior_gammac[1], ",", prior_gammac[2])
            model_string_jags <- gsub("gamma_c[j, t, time, mdrop] ~ dnorm(0, 0.01", prior_gammac_str, model_string_jags, fixed = TRUE) }
        }
        if(length(model_mc_random) != 0 & zc_random == 1) {
          if(is.null(mu.g0.prior.c) == FALSE & grepl("mu_g_c_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.g0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_g0c <- mu.g0.prior.c
            prior_g0c_str <- paste("mu_g_c_hat[t, time] ~ dnorm(", prior_g0c[1], ",", prior_g0c[2])
            model_string_jags <- gsub("mu_g_c_hat[t, time] ~ dnorm(0, 0.001", prior_g0c_str, model_string_jags, fixed = TRUE) }
          if(is.null(s.g0.prior.c) == FALSE & grepl("s_g_c_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(s.g0.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_g0c <- s.g0.prior.c
            prior_g0c_str <- paste("s_g_c_hat[t, time] ~ dunif(", prior_g0c[1], ",", prior_g0c[2])
            model_string_jags <- gsub("s_g_c_hat[t, time] ~ dunif(0, 100", prior_g0c_str, model_string_jags, fixed = TRUE) }
            } else if(length(model_mc_random) != 0 & zc_random > 1) {
          if(is.null(mu.g.prior.c) == FALSE & grepl("mu_g_c_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(mu.g.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_gc <- mu.g.prior.c
            prior_gc_str <- paste("mu_g_c_hat[j, t, time] ~ dnorm(", prior_gc[1], ",", prior_gc[2])
            model_string_jags <- gsub("mu_g_c_hat[j, t, time] ~ dnorm(0, 0.001", prior_gc_str, model_string_jags, fixed = TRUE) }
          if(is.null(s.g.prior.c) == FALSE & grepl("s_g_c_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            if(length(s.g.prior.c) != 2) {stop("provide correct hyper prior values") }
            prior_gc <- s.g.prior.c
            prior_gc_str <- paste("s_g_c_hat[j, t, time] ~ dunif(", prior_gc[1], ",", prior_gc[2])
            model_string_jags <- gsub("s_g_c_hat[j, t, time] ~ dunif(0, 100", prior_gc_str, model_string_jags, fixed = TRUE) }
        }                    
    if(pu_fixed == 1) {
      if(is.null(alpha0.prior) == FALSE & grepl("alpha[1, time] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha[2, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muu <- alpha0.prior
        prior_muu_str <- paste("alpha[1, time] ~ dnorm(", prior_muu[1], ",", prior_muu[2])
        model_string_jags <- gsub("alpha[1, time] ~ dnorm(0, 0.0000001", prior_muu_str, model_string_jags,fixed = TRUE)
        prior_muu_str <- paste("alpha[2, time] ~ dnorm(", prior_muu[1],",", prior_muu[2])
        model_string_jags <- gsub("alpha[2, time] ~ dnorm(0, 0.0000001", prior_muu_str, model_string_jags, fixed = TRUE) }
    } else if(pu_fixed > 1){
      if(is.null(alpha0.prior) == FALSE & grepl("alpha[1, 1, time] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("alpha[1, 2, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muu <- alpha0.prior
        prior_muu_str <- paste("alpha[1, 1, time] ~ dnorm(", prior_muu[1], ",", prior_muu[2])
        model_string_jags <- gsub("alpha[1, 1, time] ~ dnorm(0, 0.0000001", prior_muu_str, model_string_jags,fixed = TRUE)
        prior_muu_str <- paste("alpha[1, 2, time] ~ dnorm(", prior_muu[1], ",", prior_muu[2])
        model_string_jags <- gsub("alpha[1, 2, time] ~ dnorm(0, 0.0000001", prior_muu_str, model_string_jags, fixed = TRUE) }
        if(is.null(alpha.prior) == FALSE & grepl("alpha[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha.prior) != 2) {stop("provide correct hyper prior values") }
        prior_alphau <- alpha.prior
        prior_alphau_str <- paste("alpha[j, t, time] ~ dnorm(", prior_alphau[1], ",", prior_alphau[2])
        model_string_jags <- gsub("alpha[j, t, time] ~ dnorm(0, 0.0000001", prior_alphau_str, model_string_jags, fixed = TRUE) }
    }
    if(length(model_u_random) != 0 & pu_random == 1) {
      if(is.null(mu.a0.prior) == FALSE & grepl("mu_a_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.a0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a0 <- mu.a0.prior
        prior_a0_str <- paste("mu_a_hat[t, time] ~ dnorm(", prior_a0[1], ",", prior_a0[2])
        model_string_jags <- gsub("mu_a_hat[t, time] ~ dnorm(0, 0.001", prior_a0_str, model_string_jags, fixed = TRUE) }
      if(is.null(s.a0.prior) == FALSE & grepl("s_a_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.a0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a0 <- s.a0.prior
        prior_a0_str <- paste("s_a_hat[t, time] ~ dunif(", prior_a0[1], ",", prior_a0[2])
        model_string_jags <- gsub("s_a_hat[t, time] ~ dunif(0, 100", prior_a0_str, model_string_jags, fixed = TRUE) }
        } else if(length(model_u_random) != 0 & pu_random > 1) {
      if(is.null(mu.a.prior) == FALSE & grepl("mu_a_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.a.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a <- mu.a.prior
        prior_a_str <- paste("mu_a_hat[j, t, time] ~ dnorm(", prior_a[1], ",", prior_a[2])
        model_string_jags <- gsub("mu_a_hat[j, t, time] ~ dnorm(0, 0.001", prior_a_str, model_string_jags, fixed = TRUE) }
      if(is.null(s.a.prior) == FALSE & grepl("s_a_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.a.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a <- s.a.prior
        prior_a_str <- paste("s_a_hat[j, t, time] ~ dunif(", prior_a[1], ",", prior_a[2])
        model_string_jags <- gsub("s_a_hat[j, t, time] ~ dunif(0, 100", prior_a_str, model_string_jags, fixed = TRUE) }
    }
    if(pc_fixed == 1) {
      if(is.null(beta0.prior) == FALSE & grepl("beta[1, time] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta[2, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muc <- beta0.prior
        prior_muc_str <- paste("beta[1, time] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1, time] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
        prior_muc_str <- paste("beta[2, time] ~ dnorm(",prior_muc[1],",", prior_muc[2])
        model_string_jags <- gsub("beta[2, time] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
    } else if(pc_fixed > 1) {
      if(is.null(beta0.prior) == FALSE & grepl("beta[1, 1, time] ~ ", model_string_jags, fixed = TRUE) == TRUE & grepl("beta[1, 2, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_muc <- beta0.prior
        prior_muc_str <- paste("beta[1, 1, time] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1, 1, time] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE)
        prior_muc_str <- paste("beta[1, 2, time] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
        model_string_jags <- gsub("beta[1, 2, time] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags, fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta.prior) != 2){stop("provide correct hyper prior values") }
        prior_betac <- beta.prior
        prior_betac_str <- paste("beta[j, t, time] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
        model_string_jags <- gsub("beta[j, t, time] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
    }
    if(length(model_c_random) != 0 & pc_random == 1) {
      if(is.null(mu.b0.prior) == FALSE & grepl("mu_b_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b0 <- mu.b0.prior
        prior_b0_str <- paste("mu_b_hat[t, time] ~ dnorm(", prior_b0[1], ",", prior_b0[2])
        model_string_jags <- gsub("mu_b_hat[t, time] ~ dnorm(0, 0.001", prior_b0_str, model_string_jags, fixed = TRUE) }
      if(is.null(s.b0.prior) == FALSE & grepl("s_b_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b0 <- s.b0.prior
        prior_b0_str <- paste("s_b_hat[t, time] ~ dunif(", prior_b0[1], ",", prior_b0[2])
        model_string_jags <- gsub("s_b_hat[t, time] ~ dunif(0, 100", prior_b0_str, model_string_jags, fixed = TRUE) }
        } else if(length(model_c_random) != 0 & pc_random > 1) {
      if(is.null(mu.b.prior) == FALSE & grepl("mu_b_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b <- mu.b.prior
        prior_b_str <- paste("mu_b_hat[j, t, time] ~ dnorm(", prior_b[1], ",", prior_b[2])
        model_string_jags <- gsub("mu_b_hat[j, t, time] ~ dnorm(0, 0.001", prior_b_str, model_string_jags, fixed = TRUE) }
      if(is.null(s.b.prior) == FALSE & grepl("s_b_hat[j, t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b <- s.b.prior
        prior_b_str <- paste("s_b_hat[j, t, time] ~ dunif(", prior_b[1], ",", prior_b[2])
        model_string_jags <- gsub("s_b_hat[j, t, time] ~ dunif(0, 100", prior_b_str, model_string_jags, fixed = TRUE) }
    }                             
  if(dist_u == "norm") {
  if(is.null(sigma.prior.u) == FALSE & grepl("ls_u[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
    if(length(sigma.prior.u) != 2) {stop("provide correct hyper prior values") }
    prior_alphau <- sigma.prior.u
    prior_alphau_str <- paste("ls_u[t, time] ~ dunif(", prior_alphau[1], ",", prior_alphau[2])
    model_string_jags <- gsub("ls_u[t, time] ~ dunif(-5, 10", prior_alphau_str, model_string_jags, fixed = TRUE) }
  } else if(dist_u == "beta") {
      if(is.null(sigma.prior.u) == FALSE & grepl("s_u[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.u) != 2) {stop("provide correct hyper prior values") }
        prior_alphau <- sigma.prior.u
        prior_alphau_str <- paste("s_u[t, time] ~ dunif(", prior_alphau[1], ",", prior_alphau[2])
        model_string_jags <- gsub("s_u[t, time] ~ dunif(0, sqrt(mu_u[t, time] * (1 - mu_u[t, time])))", prior_alphau_str, model_string_jags, fixed = TRUE) }
  } else if(dist_u == "gamma" | dist_u == "logis") {
    if(is.null(sigma.prior.u) == FALSE & grepl("s_u[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.u) != 2) {stop("provide correct hyper prior values") }
      prior_alphau <- sigma.prior.u
      prior_alphau_str <- paste("s_u[t, time] ~ dunif(", prior_alphau[1], ",", prior_alphau[2])
      model_string_jags <- gsub("s_u[t, time] ~ dunif(0, 10000", prior_alphau_str, model_string_jags, fixed = TRUE) }
  } else if(dist_u == "exp" | dist_u == "bern" | dist_u == "pois") {
    if(is.null(sigma.prior.u) == FALSE) {
      stop("no prior for sigma required for the effects under the 'exp', 'bern', 'pois' distributions") }
  } else if(dist_u == "weibull") {
      if(is.null(sigma.prior.u) == FALSE & grepl("s_u[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.u) != 2) {stop("provide correct hyper prior values") }
        prior_alphau <- sigma.prior.u
        prior_alphau_str <- paste("s_u[t, time] ~ dunif(", prior_alphau[1], ",", prior_alphau[2])
        model_string_jags <- gsub("s_u[t, time] ~ dunif(0, 100", prior_alphau_str, model_string_jags, fixed = TRUE) }
  } else if(dist_u == "nbinom") {
    if(is.null(sigma.prior.u) == FALSE & grepl("tau_u[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.u) != 2) {stop("provide correct hyper prior values") }
      prior_alphau <- sigma.prior.u
      prior_alphau_str <- paste("tau_u[t, time] ~ dunif(", prior_alphau[1], ",", prior_alphau[2])
      model_string_jags <- gsub("tau_u[t, time] ~ dunif(0, 100", prior_alphau_str, model_string_jags, fixed = TRUE) }
  }
  if(dist_c == "norm") {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c[t, time] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c[t, time] ~ dunif(-5, 10", prior_alphac_str, model_string_jags, fixed = TRUE) }
  } else if(dist_c == "lnorm") {
    if(is.null(sigma.prior.c) == FALSE & grepl("ls_c[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c
      prior_alphac_str <- paste("ls_c[t, time] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
      model_string_jags <- gsub("ls_c[t, time] ~ dunif(0, 100", prior_alphac_str, model_string_jags, fixed = TRUE) }
  } else if(dist_c == "gamma") {
      if(is.null(sigma.prior.c) == FALSE & grepl("s_c[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
        prior_alphac <- sigma.prior.c
        prior_alphac_str <- paste("s_c[t, time] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
        model_string_jags <- gsub("s_c[t, time] ~ dunif(0, 10000", prior_alphac_str, model_string_jags, fixed = TRUE) }
  }
    if(exists("beta_f.prior") == TRUE) {
      if(is.null(beta_f.prior) == FALSE & grepl("beta_f", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta_f.prior) != 2) {stop("provide correct hyper prior values") }
          prior_beta_f <- beta_f.prior
          prior_beta_f_str <- paste("beta_f[t, time] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
          model_string_jags <- gsub("beta_f[t, time] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags, fixed = TRUE)
       }
     }
    if(exists("mu.b_f.prior") == TRUE) {
      if(is.null(mu.b_f.prior) == FALSE & grepl("mu_b_f_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b_f.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_f <- mu.b_f.prior
        prior_b_f_str <- paste("mu_b_f_hat[t, time] ~ dnorm(", prior_b_f[1], ",", prior_b_f[2])
        model_string_jags <- gsub("mu_b_f_hat[t, time] ~ dnorm(0, 0.001", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.b_f.prior") == TRUE) {
      if(is.null(s.b_f.prior) == FALSE & grepl("s_b_f_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b_f.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_f <- s.b_f.prior
        prior_b_f_str <- paste("s_b_f_hat[t, time] ~ dunif(", prior_b_f[1], ",", prior_b_f[2])
        model_string_jags <- gsub("s_b_f_hat[t, time] ~ dunif(0, 100", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    }
    if(exists("beta_tu.prior") == TRUE) {
      if(is.null(beta_tu.prior) == FALSE & grepl("beta_tu", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_beta_tu <- beta_tu.prior
        prior_beta_tu_str <- paste("beta_tu[t, time] ~ dnorm(", prior_beta_tu[1], ",", prior_beta_tu[2])
        model_string_jags <- gsub("beta_tu[t, time] ~ dnorm(0, 0.0000001", prior_beta_tu_str, model_string_jags, fixed = TRUE)
      }
    }
    if(exists("beta_tc.prior") == TRUE) {
      if(is.null(beta_tc.prior) == FALSE & grepl("beta_tc", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(beta_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_beta_tc <- beta_tc.prior
        prior_beta_tc_str <- paste("beta_tc[t, time] ~ dnorm(", prior_beta_tc[1], ",", prior_beta_tc[2])
        model_string_jags <- gsub("beta_tc[t, time] ~ dnorm(0, 0.0000001", prior_beta_tc_str, model_string_jags, fixed = TRUE)
      }
    }
    if(exists("mu.b_tu.prior") == TRUE) {
      if(is.null(mu.b_tu.prior) == FALSE & grepl("mu_b_tu_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_tu <- mu.b_tu.prior
        prior_b_tu_str <- paste("mu_b_tu_hat[t, time] ~ dnorm(", prior_b_tu[1], ",", prior_b_tu[2])
        model_string_jags <- gsub("mu_b_tu_hat[t, time] ~ dnorm(0, 0.001", prior_b_tu_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.b_tu.prior") == TRUE) {
      if(is.null(s.b_tu.prior) == FALSE & grepl("s_b_tu_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_tu <- s.b_tu.prior
        prior_b_tu_str <- paste("s_b_tu_hat[t, time] ~ dunif(", prior_b_tu[1], ",", prior_b_tu[2])
        model_string_jags <- gsub("s_b_tu_hat[t, time] ~ dunif(0, 100", prior_b_tu_str, model_string_jags, fixed = TRUE)
      }
    } 
    if(exists("mu.b_tc.prior") == TRUE) {
      if(is.null(mu.b_tc.prior) == FALSE & grepl("mu_b_tc_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_tc <- mu.b_tc.prior
        prior_b_tc_str <- paste("mu_b_tc_hat[t, time] ~ dnorm(", prior_b_tc[1], ",", prior_b_tc[2])
        model_string_jags <- gsub("mu_b_tc_hat[t, time] ~ dnorm(0, 0.001", prior_b_tc_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.b_tc.prior") == TRUE) {
      if(is.null(s.b_tc.prior) == FALSE & grepl("s_b_tc_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_tc <- s.b_tc.prior
        prior_b_tc_str <- paste("s_b_tc_hat[t, time] ~ dunif(", prior_b_tc[1], ",", prior_b_tc[2])
        model_string_jags <- gsub("s_b_tc_hat[t, time] ~ dunif(0, 100", prior_b_tc_str, model_string_jags, fixed = TRUE)
      }
    } 
    if(exists("alpha_tu.prior") == TRUE) {
      if(is.null(alpha_tu.prior) == FALSE & grepl("alpha_tu", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_alpha_tu <- alpha_tu.prior
        prior_alpha_tu_str <- paste("alpha_tu[t, time] ~ dnorm(", prior_alpha_tu[1], ",", prior_alpha_tu[2])
        model_string_jags <- gsub("alpha_tu[t, time] ~ dnorm(0, 0.0000001", prior_alpha_tu_str, model_string_jags, fixed = TRUE)
      }
    }
    if(exists("alpha_tc.prior") == TRUE) {
      if(is.null(alpha_tc.prior) == FALSE & grepl("alpha_tc", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(alpha_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_alpha_tc <- alpha_tc.prior
        prior_alpha_tc_str <- paste("alpha_tc[t, time] ~ dnorm(", prior_alpha_tc[1], ",", prior_alpha_tc[2])
        model_string_jags <- gsub("alpha_tc[t, time] ~ dnorm(0, 0.0000001", prior_alpha_tc_str, model_string_jags, fixed = TRUE)
      }
    }
    if(exists("mu.a_tu.prior") == TRUE) {
      if(is.null(mu.a_tu.prior) == FALSE & grepl("mu_a_tu_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.a_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a_tu <- mu.a_tu.prior
        prior_a_tu_str <- paste("mu_a_tu_hat[t, time] ~ dnorm(", prior_a_tu[1], ",", prior_a_tu[2])
        model_string_jags <- gsub("mu_a_tu_hat[t, time] ~ dnorm(0, 0.001", prior_a_tu_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.a_tu.prior") == TRUE) {
      if(is.null(s.a_tu.prior) == FALSE & grepl("s_a_tu_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.a_tu.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a_tu <- s.a_tu.prior
        prior_a_tu_str <- paste("s_a_tu_hat[t, time] ~ dunif(", prior_a_tu[1], ",", prior_a_tu[2])
        model_string_jags <- gsub("s_a_tu_hat[t, time] ~ dunif(0, 100", prior_a_tu_str, model_string_jags, fixed = TRUE)
      }
    } 
    if(exists("mu.a_tc.prior") == TRUE) {
      if(is.null(mu.a_tc.prior) == FALSE & grepl("mu_a_tc_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.a_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a_tc <- mu.a_tc.prior
        prior_a_tc_str <- paste("mu_a_tc_hat[t, time] ~ dnorm(", prior_a_tc[1], ",", prior_a_tc[2])
        model_string_jags <- gsub("mu_a_tc_hat[t, time] ~ dnorm(0, 0.001", prior_a_tc_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.a_tc.prior") == TRUE) {
      if(is.null(s.a_tc.prior) == FALSE & grepl("s_a_tc_hat[t, time] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.a_tc.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a_tc <- s.a_tc.prior
        prior_a_tc_str <- paste("s_a_tc_hat[t, time] ~ dunif(", prior_a_tc[1], ",", prior_a_tc[2])
        model_string_jags <- gsub("s_a_tc_hat[t, time] ~ dunif(0, 100", prior_a_tc_str, model_string_jags, fixed = TRUE)
      }
    } 
  model_string_prior <- model_string_jags
  return(model_string_prior)
}))