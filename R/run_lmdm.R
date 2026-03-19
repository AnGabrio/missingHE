#' An internal function to execute a JAGS longitudinal missing data model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} function and obtain posterior inferences.
#' @param data_model list containing the data for the model to be passed to JAGS.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR).
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param model_info list containing model and MCMC information to be passed to JAGS. 
#' @keywords JAGS Bayesian longitudinal missing data models 
#' @importFrom stats rnorm rbeta rgamma rlnorm rweibull rnbinom rbinom rpois rlogis rexp
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_lmdm <- function(data_model, type, dist_e, dist_c, model_info) {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  n <- data_model$n; n_long <- data_model$n_long; time <- data_model$time; max_time <- max(time, na.rm = TRUE)
  trt_lev <- names(data_model$n_trt); pe_fixed <- data_model$pe_fixed; pc_fixed <- data_model$pc_fixed
  ze_fixed <- data_model$ze_fixed; zc_fixed <- data_model$zc_fixed; mean_cov_e_fixed <- data_model$mean_cov_e_fixed; 
  mean_cov_c_fixed <- data_model$mean_cov_c_fixed; mean_z_e_fixed <- data_model$mean_z_e_fixed; mean_z_c_fixed <- data_model$mean_z_c_fixed
  mean_cov_e_random <- data_model$mean_cov_e_random; mean_cov_c_random <- data_model$mean_cov_c_random; mean_z_e_random <- data_model$mean_z_e_random; mean_z_c_random <- data_model$mean_z_c_random
  pe_random <- data_model$pe_random; pc_random <- data_model$pc_random; ze_random <- data_model$ze_random; zc_random <- data_model$zc_random
  n_clus_e <- data_model$n_clus_e; n_clus_c <- data_model$n_clus_c; n_clus_me <- data_model$n_clus_me; n_clus_mc <- data_model$n_clus_mc
  trt_pos_e <- data_model$trt_pos_e; trt_pos_c <- data_model$trt_pos_c; trt_pos_me <- data_model$trt_pos_me; trt_pos_mc <- data_model$trt_pos_mc
  eff <- data_model$eff; cost <- data_model$cost; eff_long <- data_model$eff_long; cost_long <- data_model$cost_long
  m_eff <- data_model$m_eff; m_cost <- data_model$m_cost; m_eff_long <- data_model$m_eff_long; m_cost_long <- data_model$m_cost_long
  n_patt_e <- max(data_model$m_eff_long, na.rm = TRUE); n_patt_c <- max(data_model$m_cost_long, na.rm = TRUE); 
  X_e_fixed <- as.matrix(data_model$X_e_fixed); X_c_fixed <- as.matrix(data_model$X_c_fixed)
  Z_e_fixed <- as.matrix(data_model$Z_e_fixed); Z_c_fixed <- as.matrix(data_model$Z_c_fixed)
  X_e_random <- data_model$X_e_random; X_c_random <- data_model$X_c_random; Z_e_random <- data_model$Z_e_random; Z_c_random <- data_model$Z_c_random
  if(ze_fixed == 1) { Z_e_fixed <- as.vector(unlist(data_model$Z_e_fixed))} 
  if(zc_fixed == 1) { Z_c_fixed <- as.vector(unlist(data_model$Z_c_fixed))}
  if(ze_fixed > 1) { Z_e_fixed <- as.matrix(data_model$Z_e_fixed)}
  if(zc_fixed > 1) { Z_c_fixed <- as.matrix(data_model$Z_c_fixed)}
  if(pe_random == 1) { X_e_random <- as.vector(unlist(data_model$X_e_random))}
  if(pc_random == 1) { X_c_random <- as.vector(unlist(data_model$X_c_random))}
  if(pe_random > 1) { X_e_random <- as.matrix(data_model$X_e_random)} 
  if(pc_random > 1) { X_c_random <- as.matrix(data_model$X_c_random)}
  if(ze_random == 1) { Z_e_random <- as.vector(unlist(data_model$Z_e_random))}
  if(zc_random == 1) { Z_c_random <- as.vector(unlist(data_model$Z_c_random))}
  if(ze_random > 1) { Z_e_random <- as.matrix(data_model$Z_e_random)} 
  if(zc_random > 1) { Z_c_random <- as.matrix(data_model$Z_c_random)}
  clus_e <- data_model$clus_e; clus_c <- data_model$clus_c; clus_me <- data_model$clus_me; clus_mc <- data_model$clus_mc
  datalist <- list("n"= n, "eff" = eff, "cost" = cost, "m_eff" = m_eff, "m_cost" = m_cost, "max_time" = max_time, 
                   "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "Z_e_fixed" = Z_e_fixed, "Z_c_fixed" = Z_c_fixed, 
                   "mean_cov_e_fixed" = mean_cov_e_fixed, "mean_cov_c_fixed" = mean_cov_c_fixed, "mean_z_e_fixed" = mean_z_e_fixed, "mean_z_c_fixed" = mean_z_c_fixed, 
                   "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "ze_fixed" = ze_fixed, "zc_fixed" = zc_fixed, 
                   "X_e_random" = X_e_random, "mean_cov_e_random" = mean_cov_e_random, "pe_random" = pe_random, 
                   "clus_e" = clus_e, "n_clus_e" = n_clus_e, "X_c_random" = X_c_random, "mean_cov_c_random" = mean_cov_c_random, "pc_random" = pc_random, 
                   "clus_c"= clus_c, "n_clus_c" = n_clus_c, "Z_e_random" = Z_e_random, "mean_z_e_random" = mean_z_e_random, "ze_random" = ze_random, 
                   "clus_me"= clus_me, "n_clus_me" = n_clus_me, "Z_c_random" = Z_c_random,"mean_z_c_random" = mean_z_c_random, "zc_random" = zc_random, 
                   "clus_mc" = clus_mc, "n_clus_mc" = n_clus_mc, "n_patt_e" = n_patt_e, "n_patt_c" = n_patt_c)
  e_random_list <- c("X_e_random","mean_cov_e_random", "pe_random", "clus_e", "n_clus_e")
  c_random_list <- c("X_c_random","mean_cov_c_random", "pc_random", "clus_c", "n_clus_c")
  me_random_list <- c("Z_e_random","mean_z_e_random", "ze_random", "clus_me", "n_clus_me")
  mc_random_list <- c("Z_c_random","mean_z_c_random", "zc_random", "clus_mc", "n_clus_mc")
  inits <- model_info$inits
  n.chains <- model_info$n.chains
  n.iter <- model_info$n.iter
  n.burnin <- model_info$n.burnin
  n.thin <- model_info$n.thin
  DIC <- model_info$dic
  pd <- model_info$pd
  n.iter.pd <- model_info$n.iter.pd
  n.adapt <- model_info$n.adapt
  is_c_random_c <- model_info$is_c_random_c
  is_int_c_random_c <- model_info$is_int_c_random_c
  is_me_random_e <- model_info$is_me_random_e
  is_mc_random_c <- model_info$is_mc_random_c
  is_int_me_random_e <- model_info$is_int_me_random_e
  is_int_mc_random_c <- model_info$is_int_mc_random_c
  ind_random <- model_info$ind_random
  ind_fixed <- model_info$ind_fixed
  ind_time_fixed <- model_info$ind_time_fixed
  if(ze_fixed == 1) { ze_fixed_index <- match("ze_fixed", names(datalist))
  datalist <- datalist[-ze_fixed_index]}
  if(zc_fixed == 1) { zc_fixed_index <- match("zc_fixed", names(datalist))
  datalist <- datalist[-zc_fixed_index]}
  if(length(model_info$model_e_random) != 0) {
    if(pe_random == 1) { 
      pe_random_index <- match("pe_random", names(datalist))
      datalist <- datalist[-pe_random_index]}
  } else if(length(model_info$model_e_random) == 0) { 
    e_random_index <- match(e_random_list, names(datalist))
    datalist <- datalist[-e_random_index]
  }
  if(length(model_info$model_c_random) != 0 & !is_c_random_c) {
    if(pc_random == 1) {
      pc_random_index <- match("pc_random", names(datalist))
      datalist <- datalist[-pc_random_index]}
  } else if(length(model_info$model_c_random) != 0 & is_c_random_c) {
    c_random_index <- match(c_random_list[1:3], names(datalist))
    datalist <- datalist[-c_random_index] 
  } else if(length(model_info$model_c_random) == 0) { 
    c_random_index <- match(c_random_list, names(datalist))
    datalist <- datalist[-c_random_index]
  }
  if(length(model_info$model_me_random) != 0 & !is_me_random_e) {
    if(ze_random == 1 | length(model_info$model_me_random) != 0 & is_me_random_e) {
      ze_random_index <- match("ze_random", names(datalist))
      datalist <- datalist[-ze_random_index]}
  } else if(length(model_info$model_me_random) != 0 & is_me_random_e) {
    me_random_index <- match(me_random_list[1:3], names(datalist))
    datalist <- datalist[-me_random_index] 
  } else if(length(model_info$model_me_random) == 0) { 
    me_random_index <- match(me_random_list, names(datalist))
    datalist <- datalist[-me_random_index] 
  }
  if(length(model_info$model_mc_random) != 0 & !is_mc_random_c) {
    if(zc_random == 1 | length(model_info$model_mc_random) != 0 & is_mc_random_c) {
      zc_random_index <- match("zc_random", names(datalist))
      datalist <- datalist[-zc_random_index]}
  } else if(length(model_info$model_mc_random) != 0 & is_mc_random_c) {
    mc_random_index <- match(mc_random_list[1:3], names(datalist))
    datalist <- datalist[-mc_random_index] 
  } else if(length(model_info$model_mc_random) == 0) { 
    mc_random_index <- match(mc_random_list, names(datalist))
    datalist <- datalist[-mc_random_index] 
  }
  params <- c("eff", "cost", "tmu_e", "tmu_c", "s_e", "s_c", "p_e", "p_c", 
              "beta", "alpha", "gamma_e", "gamma_c", "delta_e", "delta_c")
  params <- c(params, "loglik_e", "loglik_c", "loglik_me", "loglik_mc")
  if(!ind_fixed) { params <- c(params, "beta_f")}
  if(model_info$time_dep == "AR1") { params <- c(params, "beta_te", "beta_tc", "alpha_te", "alpha_tc")}
  if(length(model_info$model_e_random) != 0){ params <- c(params, "a")}
  if(length(model_info$model_e_random) != 0 & !ind_time_fixed & model_info$time_dep == "AR1"){ params <- c(params, "a_te", "a_tc")}
  if(length(model_info$model_c_random) != 0 & !is_c_random_c){ params <- c(params, "b")}
  if(length(model_info$model_c_random) != 0 & !is_c_random_c & !ind_time_fixed & model_info$time_dep == "AR1"){ params <- c(params, "b_te", "b_tc")}
  if(length(model_info$model_c_random) != 0 & is_c_random_c & !ind_time_fixed & model_info$time_dep == "AR1"){
    if(length(model_info$model_c_random) == 1 & "e" %in% model_info$model_c_random) { params <- c(params, "b_te", "b_tc")}
  }
  if(length(model_info$model_me_random) != 0 & !is_me_random_e){ params <- c(params, "g_e")}
  if(length(model_info$model_mc_random) != 0 & !is_mc_random_c){ params <- c(params, "g_c")}
  if(length(model_info$model_c_random) != 0 & !ind_random) { params <- c(params, "b_f")}
  if(ind_random | model_info$time_dep == "none"){
    if(any(c("a_te", "a_tc") %in% params)) {
      a_corr_index <- match(c("a_te", "a_tc"), params)
      params <- params[-a_corr_index]}
    if(any(c("b_te", "b_tc") %in% params)) {
      b_corr_index <- match(c("b_te", "b_tc"), params)
      params <- params[-b_corr_index]}
  }
  if(type %in% c("MNAR_cost", "MAR")) {
    deltae_index <- match("delta_e", params)
    params <- params[-deltae_index] 
  }
  if(type %in% c("MNAR_eff", "MAR")) {
    deltac_index <- match("delta_c", params)
    params <- params[-deltac_index]
  }
  if(type %in% c("MNAR", "MNAR_eff")) {
    if(length(model_info$model_me_random) != 0 & "e" %in% model_info$model_me_random){
      params <- c(params, "d_e")}
  }
  if(type %in% c("MNAR", "MNAR_cost")) {
    if(length(model_info$model_mc_random) != 0 & "c" %in% model_info$model_mc_random){
      params <- c(params, "d_c") 
    }
  }
  if(dist_e %in% c("norm", "negbin", "logis")) {
    mci_e_params <- c("cmu_e", "tau_e") 
  } else if(dist_e %in% c("beta", "gamma", "weib")) {
    mci_e_params <- c("cmu_e", "ctau_e")
  } else if(dist_e %in% c("exp", "bern", "pois")) {
    mci_e_params <- c("cmu_e")
  } 
  if(dist_c == "norm") {
    mci_c_params <- c("cmu_c", "tau_c") 
  } else if(dist_c == "gamma") {
    mci_c_params <- c("cmu_c", "ctau_c")
  } else if(dist_c == "lnorm") {
    mci_c_params <- c("clmu_c", "ltau_c")
  } 
  params <- c(params, mci_e_params, mci_c_params)
  model_txt_info <- list("pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "ze_fixed" = ze_fixed, "zc_fixed" = zc_fixed, 
                         "ind_fixed" = ind_fixed, "pe_random" = pe_random, "pc_random" = pc_random, 
                         "ze_random" = ze_random, "zc_random" = zc_random, "ind_time_fixed" = model_info$ind_time_fixed,
                         "ind_random" = ind_random, "model_e_random" = model_info$model_e_random, "model_c_random" = model_info$model_c_random, 
                         "model_me_random" = model_info$model_me_random, "model_mc_random" = model_info$model_mc_random,
                         "is_c_random_c" = is_c_random_c, "is_int_c_random_c" = is_int_c_random_c, 
                         "is_me_random_e" = is_me_random_e, "is_mc_random_c" = is_mc_random_c,
                         "is_int_me_random_e" = is_int_me_random_e, "is_int_mc_random_c" = is_int_mc_random_c, 
                         "prior" = model_info$prior, "time_dep" = model_info$time_dep)
  filein <- write_lmdm(type = type , dist_e = dist_e, dist_c = dist_c, model_txt_info = model_txt_info)
  model <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, 
                        model.file = filein, n.chains = n.chains, n.iter = n.iter, 
                        n.burnin = n.burnin, DIC = DIC, pD = pd, n.thin = n.thin,
                        n.iter.pd = n.iter.pd, n.adapt = n.adapt)
  tmu_e <- model$BUGSoutput$sims.list$tmu_e
  tmu_c <- model$BUGSoutput$sims.list$tmu_c
  s_e <- model$BUGSoutput$sims.list$s_e
  s_c <- model$BUGSoutput$sims.list$s_c
  alpha <- model$BUGSoutput$sims.list$alpha
  beta <- model$BUGSoutput$sims.list$beta
  p_e <- model$BUGSoutput$sims.list$p_e
  p_c <- model$BUGSoutput$sims.list$p_c
  gamma_e <- model$BUGSoutput$sims.list$gamma_e
  gamma_c <- model$BUGSoutput$sims.list$gamma_c
  colnames(s_e) <- colnames(s_c) <- c(1:max_time) 
  dimnames(alpha) <- list(NULL, colnames(as.data.frame(X_e_fixed)), 1:max_time)
  dimnames(beta) <- list(NULL, colnames(as.data.frame(X_c_fixed)), 1:max_time)
  dimnames(p_e) <- list(NULL, 1:max_time, 1:n_patt_e)
  dimnames(p_c) <- list(NULL, 1:max_time, 1:n_patt_c)
  if(is.vector(Z_e_fixed)) {
    dimnames(gamma_e) <- list(NULL, 1:max_time, 1:n_patt_e)
  } else { dimnames(gamma_e) <- list(NULL, colnames(as.data.frame(Z_e_fixed)), 1:max_time, 1:n_patt_e)}
  if(is.vector(Z_c_fixed)) {
    dimnames(gamma_c) <- list(NULL, 1:max_time, 1:n_patt_c)
  } else { dimnames(gamma_c) <- list(NULL, colnames(as.data.frame(Z_c_fixed)), 1:max_time, 1:n_patt_c)}
  if(dist_e %in% c("norm", "negbin", "logis")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    tau_e <- model$BUGSoutput$sims.list$tau_e
  } else if(dist_e %in% c("beta", "gamma", "weib")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    ctau_e <- model$BUGSoutput$sims.list$ctau_e
  } else if(dist_e %in% c("exp", "bern", "pois")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
  } 
  if(dist_c == "norm") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    tau_c <- model$BUGSoutput$sims.list$tau_c
  } else if(dist_c == "gamma") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    ctau_c <- model$BUGSoutput$sims.list$ctau_c
  } else if(dist_c == "lnorm") {
    clmu_c <- model$BUGSoutput$sims.list$clmu_c
    ltau_c <- model$BUGSoutput$sims.list$ltau_c
  } 
  S <- model_info$n.mci
  trt_index <- c(data_model$trt_index)
  trt_index_long <- c(data_model$trt_index_long)
  mu_e <- mu_c <- array(NA, dim = c(dim(beta)[1], length(data_model$n_trt), max_time), dimnames = list(NULL, trt_lev, 1:max_time))
  for(i in 1:n.iter) {
    for(trt in 1:length(data_model$n_trt)) {
      for(time in 1:max_time) {
      if(dist_e == "norm") {
        mu_e[i, trt, time] <- mean(rnorm(S, mean = cmu_e[i, unlist(data_model$trt_index[trt]), time], sd = s_e[i, time]), na.rm = TRUE)
      }
      if(dist_e == "negbin") {
        mu_e[i, trt, time] <- mean(rnbinom(S, prob = tau_e[i, time]/(tau_e[i, time] + cmu_e[i, unlist(data_model$trt_index[trt]), time]), size = tau_e[i, time]), na.rm = TRUE)
      }  
      if(dist_e == "logis") {
        mu_e[i, trt, time] <- mean(rlogis(S, location = cmu_e[i, unlist(data_model$trt_index[trt]), time], scale = 1 / tau_e[i, time]), na.rm = TRUE)
      }
      if(dist_e == "beta") {
        mu_e[i, trt, time] <- mean(rbeta(S, shape1 = cmu_e[i, unlist(data_model$trt_index[trt]), time] * ctau_e[i, unlist(data_model$trt_index[trt]), time], 
                                   shape2 = (1 - cmu_e[i, unlist(data_model$trt_index[trt]), time]) * ctau_e[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_e == "gamma") {
        mu_e[i, trt, time] <- mean(rgamma(S, shape = cmu_e[i, unlist(data_model$trt_index[trt]), time] * ctau_e[i, unlist(data_model$trt_index[trt]), time], 
                                    rate = ctau_e[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_e == "weib") {
        mu_e[i, trt, time] <- mean(rweibull(S, shape = ctau_e[i, unlist(data_model$trt_index[trt]), time]), scale = cmu_e[i, unlist(data_model$trt_index[trt]), time] / exp(lgamma(1 + 1 / ctau_e[i, unlist(data_model$trt_index[trt]), time])), na.rm = TRUE)
      }
      if(dist_e == "exp") {
        mu_e[i, trt, time] <- mean(rexp(S, rate = 1 / cmu_e[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_e == "bern") {
        mu_e[i, trt, time] <- mean(rbinom(S, size = 1, prob = cmu_e[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_e == "pois") {
        mu_e[i, trt, time] <- mean(rpois(S, lambda = cmu_e[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_c == "norm") {
        mu_c[i, trt, time] <- mean(rnorm(S, mean = cmu_c[i, unlist(data_model$trt_index[trt]), time], sd = s_c[i, time]), na.rm = TRUE)
      }
      if(dist_c == "gamma") {
        mu_c[i, trt, time] <- mean(rgamma(S, shape = cmu_c[i, unlist(data_model$trt_index[trt]), time] * ctau_c[i, unlist(data_model$trt_index[trt]), time], 
                                    rate = ctau_c[i, unlist(data_model$trt_index[trt]), time]), na.rm = TRUE)
      }
      if(dist_c == "lnorm") {
        mu_c[i, trt, time] <- mean(rlnorm(S, meanlog = clmu_c[i, unlist(data_model$trt_index[trt]), time], sdlog = 1 / sqrt(ltau_c[i, time])), na.rm = TRUE)
      }
      }
    }
  }
  dimnames(mu_e) <- dimnames(mu_c) <- list(NULL, trt_lev, 1:max_time)
  a <- b <- g_e <- g_c <- d_e <- d_c <- NULL
  if(length(model_info$model_e_random) != 0) { 
    a <- model$BUGSoutput$sims.list$a
    if(is.vector(X_e_random)) {
      dimnames(a) <- list(NULL, paste("(Intercept).", data_model$clus_e_lev, sep = ""), 1:max_time)
    } else if(length(dim(a)) == 4) {
      dimnames(a)[[2]] <- colnames(as.data.frame(X_e_random))
      dimnames(a)[[3]] <- data_model$clus_e_lev
      dimnames(a)[[4]] <- 1:max_time
    }
  }
  if(length(model_info$model_c_random) != 0 & !is_c_random_c) { 
    b <- model$BUGSoutput$sims.list$b
    if(is.vector(X_c_random)) {
      dimnames(b) <- list(NULL, paste("(Intercept).", data_model$clus_c_lev, sep = ""), 1:max_time)
    } else if(length(dim(b)) == 4) {
      dimnames(b)[[2]] <- colnames(as.data.frame(X_c_random))
      dimnames(b)[[3]] <- data_model$clus_c_lev
      dimnames(b)[[4]] <- 1:max_time
    }
  }
  if(length(model_info$model_me_random) != 0 & !is_me_random_e) { 
    g_e <- model$BUGSoutput$sims.list$g_e
    if(is.vector(Z_e_random)) {
      dimnames(g_e) <- list(NULL, paste("(Intercept).", data_model$clus_me_lev, sep = ""), 1:max_time)
    } else if(length(dim(g_e)) == 4) {
      dimnames(g_e)[[2]] <- colnames(as.data.frame(Z_e_random))
      dimnames(g_e)[[3]] <- data_model$clus_me_lev
      dimnames(g_e)[[4]] <- 1:max_time
    }
  }
  if(length(model_info$model_mc_random) != 0 & !is_mc_random_c) { 
    g_c <- model$BUGSoutput$sims.list$g_c
    if(is.vector(Z_c_random)) {
      dimnames(g_c) <- list(NULL, paste("(Intercept).", data_model$clus_mc_lev, sep = ""), 1:max_time)
    } else if(length(dim(g_c)) == 4) {
      dimnames(g_c)[[2]] <- colnames(as.data.frame(Z_c_random))
      dimnames(g_c)[[3]] <- data_model$clus_mc_lev
      dimnames(g_c)[[4]] <- 1:max_time
    }
  }
  eff_pos <- apply(model$BUGSoutput$sims.list$eff, 1L, c) 
  cost_pos <- apply(model$BUGSoutput$sims.list$cost, 1L, c)
  eff_pos_imp_avg <- apply(eff_pos, 1, mean, na.rm = T)[is.na(data_model$eff_long)] 
  cost_pos_imp_avg <- apply(cost_pos, 1, mean, na.rm = T)[is.na(data_model$cost_long)] 
  eff_pos_imp_ql <- apply(eff_pos, 1, quantile, prob = model_info$prob[1], na.rm = T)[is.na(data_model$eff_long)] 
  eff_pos_imp_qu <- apply(eff_pos, 1, quantile, prob = model_info$prob[2], na.rm = T)[is.na(data_model$eff_long)] 
  cost_pos_imp_ql <- apply(cost_pos, 1, quantile, prob = model_info$prob[1], na.rm = T)[is.na(data_model$cost_long)] 
  cost_pos_imp_qu <- apply(cost_pos, 1, quantile, prob = model_info$prob[2], na.rm = T)[is.na(data_model$cost_long)]   
  efft_pos <- costt_pos <- list()
  efft_pos_imp_avg <- efft_pos_imp_ql <- efft_pos_imp_qu <- list()
  costt_pos_imp_avg <- costt_pos_imp_ql <- costt_pos_imp_qu <- list()
  for(i in data_model$trt_lev) {
    efft_pos[[i]] <- eff_pos[trt_index_long[[i]], ]
    efft_pos_imp_avg[[i]] <- apply(efft_pos[[i]][data_model$m_efft[[i]] == 1, ], 2, mean, na.rm = TRUE)
    efft_pos_imp_ql[[i]] <- apply(efft_pos[[i]][data_model$m_efft[[i]] == 1, ], 2, quantile, prob = model_info$prob[1], na.rm = TRUE)
    efft_pos_imp_qu[[i]] <- apply(efft_pos[[i]][data_model$m_efft[[i]] == 1, ], 2, quantile, prob = model_info$prob[2], na.rm = TRUE)
    costt_pos[[i]] <- cost_pos[trt_index_long[[i]], ]
    costt_pos_imp_avg[[i]] <- apply(costt_pos[[i]][data_model$m_costt[[i]] == 1, ], 2, mean, na.rm = TRUE)
    costt_pos_imp_ql[[i]] <- apply(costt_pos[[i]][data_model$m_costt[[i]] == 1, ], 2, quantile, prob = model_info$prob[1], na.rm = TRUE)
    costt_pos_imp_qu[[i]] <- apply(costt_pos[[i]][data_model$m_costt[[i]] == 1, ], 2, quantile, prob = model_info$prob[2], na.rm = TRUE)
  }
  loglik_e <- model$BUGSoutput$sims.list$loglik_e
  loglik_c <- model$BUGSoutput$sims.list$loglik_c
  loglik_me <- model$BUGSoutput$sims.list$loglik_me
  loglik_mc <- model$BUGSoutput$sims.list$loglik_mc
  if(!ind_fixed & !ind_time_fixed & pc_fixed != 0) {
    beta_f <- model$BUGSoutput$sims.list$beta_f
    colnames(beta_f) <- paste("e.time.", 1:max_time, sep = "")
    if(model_info$time_dep == "AR1") {
      beta_te <- model$BUGSoutput$sims.list$beta_te
      beta_tc <- model$BUGSoutput$sims.list$beta_tc
      alpha_te <- model$BUGSoutput$sims.list$alpha_te
      alpha_tc <- model$BUGSoutput$sims.list$alpha_tc
      colnames(beta_te) <- paste("beta.e.time.", 1:max_time, sep = "")
      colnames(beta_tc) <- paste("beta.c.time.", 1:max_time, sep = "")
      colnames(alpha_te) <- paste("alpha.e.time.", 1:max_time, sep = "")
      colnames(alpha_tc) <- paste("alpha.c.time.", 1:max_time, sep = "")     
      beta <- list("beta" = beta, "beta_f" = beta_f, "beta_te" = beta_te, "beta_tc" = beta_tc)
      alpha <- list("alpha" = alpha, "alpha_te" = alpha_te, "alpha_tc" = alpha_tc)
    } else if(model_info$time_dep %in% c("none", "biv")) {
      beta <- list("beta" = beta, "beta_f" = beta_f)}
   }
   if(!ind_fixed & !ind_time_fixed & pc_fixed == 0) {
      beta_f <- model$BUGSoutput$sims.list$beta_f
      colnames(beta_f) <- paste("e.time.", 1:max_time, sep = "")
      if(model_info$time_dep == "AR1") {
        beta_te <- model$BUGSoutput$sims.list$beta_te
        beta_tc <- model$BUGSoutput$sims.list$beta_tc
        alpha_te <- model$BUGSoutput$sims.list$alpha_te
        alpha_tc <- model$BUGSoutput$sims.list$alpha_tc
        colnames(beta_te) <- paste("beta.e.time.", 1:max_time, sep = "")
        colnames(beta_tc) <- paste("beta.c.time.", 1:max_time, sep = "")
        colnames(alpha_te) <- paste("alpha.e.time.", 1:max_time, sep = "")
        colnames(alpha_tc) <- paste("alpha.c.time.", 1:max_time, sep = "")
        beta <- list("beta_f" = beta_f, "beta_te" = beta_te, "beta_tc" = beta_tc)
        alpha <- list("alpha" = alpha, "alpha_te" = alpha_te, "alpha_tc" = alpha_tc)
      } else if(model_info$time_dep %in% c("none", "biv")) {
        beta <- list("beta_f" = beta_f)}
    }
  if(!ind_fixed & ind_time_fixed & pc_fixed != 0) {
    beta_f <- model$BUGSoutput$sims.list$beta_f
    colnames(beta_f) <- paste("e.time.", 1:max_time, sep = "")
    beta <- list("beta" = beta, "beta_f" = beta_f)
  } else if(!ind_fixed & ind_time_fixed & pc_fixed == 0) {
    beta_f <- model$BUGSoutput$sims.list$beta_f
    colnames(beta_f) <- paste("e.time.", 1:max_time, sep = "")
    beta <- list("beta_f" = beta_f)
  }
  if(length(model_info$model_c_random) != 0 & !ind_random & !ind_time_fixed) {
    b_f <- model$BUGSoutput$sims.list$b_f
    dimnames(b_f) <-  list(NULL, paste("e", data_model$clus_c_lev, sep = "."), 1:max_time)
    if(model_info$time_dep %in% c("biv")) {
    if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & is_c_random_c) {
      b <- list("b_f" = b_f) 
    } else if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & !is_c_random_c) {
      b <- list("b" = b, "b_f" = b_f)
    }
    }
    if(model_info$time_dep %in% c("AR1")) {
      b_te <- model$BUGSoutput$sims.list$b_te
      b_tc <- model$BUGSoutput$sims.list$b_tc
      a_te <- model$BUGSoutput$sims.list$a_te
      a_tc <- model$BUGSoutput$sims.list$a_tc
      dimnames(b_te) <- list(NULL, paste("e", data_model$clus_c_lev, sep = "."), 1:max_time)
      dimnames(b_tc) <- list(NULL, paste("c", data_model$clus_c_lev, sep = "."), 1:max_time)
      dimnames(a_te) <- list(NULL, paste("e", data_model$clus_c_lev, sep = "."), 1:max_time)
      dimnames(a_tc) <- list(NULL, paste("c", data_model$clus_c_lev, sep = "."), 1:max_time)
    if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & is_c_random_c) {
      b <- list("b_f" = b_f, "b_te" = b_te, "b_tc" = b_tc) 
    } else if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & !is_c_random_c) {
      b <- list("b" = b, "b_f" = b_f, "b_te" = b_te, "b_tc" = b_tc)
    }
    if(length(model_info$model_e_random) != 0) {
      a <- list("a" = a, "a_te" = a_te, "a_tc" = a_tc) 
    }
  }
  }
  if(length(model_info$model_c_random) != 0 & !ind_random & model_info$time_dep %in% c("none", "biv")) {
    b_f <- model$BUGSoutput$sims.list$b_f
    dimnames(b_f) <-  list(NULL, paste("e", data_model$clus_c_lev, sep = "."), 1:max_time)
    if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & is_c_random_c) {
      b <- list("b_f" = b_f) 
    } else if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & !is_c_random_c) {
      b <- list("b" = b, "b_f" = b_f)
    }
  }
  if(type == "MNAR") {
    delta_e <- model$BUGSoutput$sims.list$delta_e
    delta_c <- model$BUGSoutput$sims.list$delta_c
    colnames(delta_e) <- paste("delta_e.time.", 1:max_time, sep = "")
    colnames(delta_c) <- paste("delta_c.time.", 1:max_time, sep = "")
    if(length(model_info$model_me_random) != 0 & "e" %in% model_info$model_me_random) {
      d_e <- model$BUGSoutput$sims.list$d_e
      dimnames(d_e) <- list(NULL, paste("d_e_", data_model$clus_me_lev, sep = "."), 1:max_time)
      d_e <- list("d_e" = d_e)}
    if(length(model_info$model_mc_random) != 0 & "c" %in% model_info$model_mc_random) {
      d_c <- model$BUGSoutput$sims.list$d_c
      dimnames(d_c) <- list(NULL, paste("d_c_", data_model$clus_mc_lev, sep = "."), 1:max_time)
      d_c <- list("d_c" = d_c)}
  } else if(type == "MNAR_eff") {
    delta_e <- model$BUGSoutput$sims.list$delta_e
    colnames(delta_e) <- paste("delta_e.time.", 1:max_time, sep = "")
    if(length(model_info$model_me_random) != 0 & "e" %in% model_info$model_me_random) {
      d_e <- model$BUGSoutput$sims.list$d_e
      dimnames(d_e) <- list(NULL, paste("d_e_", data_model$clus_me_lev, sep = "."), 1:max_time)
      d_e <- list("d_e" = d_e)}
  } else if(type == "MNAR_cost") {
    delta_c <- model$BUGSoutput$sims.list$delta_c
    colnames(delta_c) <- paste("delta_c.time.", 1:max_time, sep = "")
    if(length(model_info$model_mc_random) != 0 & "c" %in% model_info$model_mc_random) {
      d_c <- model$BUGSoutput$sims.list$d_c
      dimnames(d_c) <- list(NULL, paste("d_c_", data_model$clus_mc_lev, sep = "."), 1:max_time)
      d_c <- list("d_c" = d_c)}
  }  
  if(n.chains > 1) {
    param_jags_show <- c("eff", "cost", "loglik_e", "loglik_c", "loglik_me", "loglik_mc", 
                         "cmu_c", "cmu_e")
    if(dist_e %in% c("beta", "gamma", "weib")) { param_jags_show <- c(param_jags_show, "ctau_e")}
    if(dist_c %in% c("gamma")) { param_jags_show <- c(param_jags_show, "ctau_c")}
    if(dist_c %in% c("lnorm")) { param_jags_show <- c(param_jags_show, "clmu_c")}
    model_sum <- round(jagsresults(x = model, params = param_jags_show, invert = TRUE), digits = 3)
  } else{ model_sum <- NULL}
  loglik <- list("effects" = loglik_e, "costs" = loglik_c, 
                 "missing indicators effects" = loglik_me, 
                 "missing indicators costs" = loglik_mc)
  imputed_e <- list("avg" = eff_pos_imp_avg, "ql" = eff_pos_imp_ql, "qu" = eff_pos_imp_qu, "index" = which(is.na(data_model$eff)))
  imputed_c <- list("avg" = cost_pos_imp_avg, "ql" = cost_pos_imp_ql, "qu" = cost_pos_imp_qu, "index" = which(is.na(data_model$cost)))
  imputed_et <- list("avg" = efft_pos_imp_avg, "ql" = efft_pos_imp_ql, "qu" = efft_pos_imp_avg, "index" = lapply(data_model$m_efft, function(x) which(x == 1)))
  imputed_ct <- list("avg" = costt_pos_imp_avg, "ql" = costt_pos_imp_ql, "qu" = costt_pos_imp_avg, "index" = lapply(data_model$m_costt, function(x) which(x == 1)))
  imputed <- list("effects" = imputed_e, "costs" = imputed_c)
  imputedt <- list("effects" = imputed_et, "costs" = imputed_ct)
  if(type == "MAR") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "pmis_effects" = p_e, "mis_effects_fixed" = gamma_e, "pmis_costs" = p_c, "mis_costs_fixed" = gamma_c, 
                              "effects_random" = a, "costs_random" = b, "mis_effects_random" = g_e, "mis_costs_random" = g_c, 
                              "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, "type" = "MAR", "method" = "LMDM", 
                              "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein)
  } else if(type == "MNAR") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "pmis_effects" = p_e, "mis_effects_fixed" = gamma_e, "pmis_costs" = p_c, "mis_costs_fixed" = gamma_c, 
                              "mis_effects_mnar_fixed" = delta_e, "mis_costs_mnar_fixed" = delta_c,
                              "effects_random" = a, "costs_random" = b, "mis_effects_random" = g_e, "mis_costs_random" = g_c, 
                              "mis_effects_mnar_random" = d_e, "mis_costs_mnar_random" = d_c, "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR", "method" = "LMDM", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein)
  } else if(type == "MNAR_eff") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "pmis_effects" = p_e, "mis_effects_fixed" = gamma_e, "pmis_costs" = p_c, "mis_costs_fixed" = gamma_c, 
                              "mis_effects_mnar_fixed" = delta_e,
                              "effects_random" = a, "costs_random" = b, "mis_effects_random" = g_e, "mis_costs_random" = g_c, 
                              "mis_effects_mnar_random" = d_e, "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR_eff", "method" = "LMDM", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein)
  } else if(type == "MNAR_cost") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "pmis_effects" = p_e, "mis_effects_fixed" = gamma_e, "pmis_costs" = p_c, "mis_costs_fixed" = gamma_c, 
                              "mis_costs_mnar_fixed" = delta_c,
                              "effects_random" = a, "costs_random" = b, "mis_effects_random" = g_e, "mis_costs_random" = g_c, 
                              "mis_costs_mnar_random" = d_c, "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR_cost", "method" = "LMDM", "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein)
  }
  if(n.chains == 1) { model_output_jags <- model_output_jags[-1]}
  return(model_output_jags = model_output_jags)
}