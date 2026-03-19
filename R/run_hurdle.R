#' An internal function to execute a JAGS hurdle model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} function and obtain posterior inferences.
#' @param data_model list containing the data for the model to be passed to JAGS.
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR), Structural At Random (SAR),
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param model_info list containing model and MCMC information to be passed to JAGS. 
#' @keywords JAGS Bayesian hurdle models 
#' @importFrom stats rnorm rbeta rgamma rlnorm rweibull rnbinom rbinom rpois rlogis qlogis rexp
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_hurdle <- function(data_model, type, dist_e, dist_c, model_info) {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  n <- data_model$n; trt_lev <- names(data_model$n_trt); pe_fixed <- data_model$pe_fixed; pc_fixed <- data_model$pc_fixed
  ze_fixed <- data_model$ze_fixed; zc_fixed <- data_model$zc_fixed; se <- data_model$se; sc <- data_model$sc; sde <- data_model$sde; sdc <- data_model$sdc
  mean_cov_e_fixed <- data_model$mean_cov_e_fixed; mean_cov_c_fixed <- data_model$mean_cov_c_fixed; mean_z_e_fixed <- data_model$mean_z_e_fixed; mean_z_c_fixed <- data_model$mean_z_c_fixed
  mean_cov_e_random <- data_model$mean_cov_e_random; mean_cov_c_random <- data_model$mean_cov_c_random; mean_z_e_random <- data_model$mean_z_e_random; mean_z_c_random <- data_model$mean_z_c_random
  pe_random <- data_model$pe_random; pc_random <- data_model$pc_random; ze_random <- data_model$ze_random; zc_random <- data_model$zc_random
  n_clus_e <- data_model$n_clus_e; n_clus_c <- data_model$n_clus_c; n_clus_se <- data_model$n_clus_se; n_clus_sc <- data_model$n_clus_sc
  trt_pos_e <- data_model$trt_pos_e; trt_pos_c <- data_model$trt_pos_c; trt_pos_se <- data_model$trt_pos_se; trt_pos_sc <- data_model$trt_pos_sc
  eff <- data_model$eff; cost <- data_model$cost; m_eff <- data_model$m_eff; m_cost <- data_model$m_cost
  s_eff <- data_model$s_eff; s_cost <- data_model$s_cost; X_e_fixed <- as.matrix(data_model$X_e_fixed)
  X_c_fixed <- as.matrix(data_model$X_c_fixed); X_e_random <- data_model$X_e_random; X_c_random <- data_model$X_c_random
  Z_e_fixed <- data_model$Z_e_fixed; Z_c_fixed <- data_model$Z_c_fixed
  Z_e_random <- data_model$Z_e_random; Z_c_random <- data_model$Z_c_random  
  if(!is.null(Z_e_fixed)) { Z_e_fixed <- as.matrix(Z_e_fixed)}
  if(!is.null(Z_c_fixed)) { Z_c_fixed <- as.matrix(Z_c_fixed)}
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
  clus_e <- data_model$clus_e; clus_c <- data_model$clus_c
  clus_se <- data_model$clus_se; clus_sc <- data_model$clus_sc
  datalist <- list("n"= n, "eff" = eff, "cost" = cost, "s_eff" = s_eff, "s_cost" = s_cost, "se" = se, "sc" = sc,
                   "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "Z_e_fixed" = Z_e_fixed, "Z_c_fixed" = Z_c_fixed, 
                   "mean_cov_e_fixed" = mean_cov_e_fixed, "mean_cov_c_fixed" = mean_cov_c_fixed, "mean_z_e_fixed" = mean_z_e_fixed, "mean_z_c_fixed" = mean_z_c_fixed, 
                   "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "ze_fixed" = ze_fixed, "zc_fixed" = zc_fixed, 
                   "X_e_random" = X_e_random, "mean_cov_e_random" = mean_cov_e_random, "pe_random" = pe_random, 
                   "clus_e" = clus_e, "n_clus_e" = n_clus_e, "X_c_random" = X_c_random, "mean_cov_c_random" = mean_cov_c_random, "pc_random" = pc_random, 
                   "clus_c"= clus_c, "n_clus_c" = n_clus_c, "Z_e_random" = Z_e_random, "mean_z_e_random" = mean_z_e_random, "ze_random" = ze_random, 
                   "clus_se"= clus_se, "n_clus_se" = n_clus_se, "Z_c_random" = Z_c_random,"mean_z_c_random" = mean_z_c_random, "zc_random" = zc_random, 
                   "clus_sc" = clus_sc, "n_clus_sc" = n_clus_sc, "sde" = sde, "sdc" = sdc)
  se_fixed_list <- c("Z_e_fixed","mean_z_e_fixed", "ze_fixed", "s_eff", "sde", "se")
  sc_fixed_list <- c("Z_c_fixed","mean_z_c_fixed", "zc_fixed", "s_cost", "sdc", "sc")
  e_random_list <- c("X_e_random","mean_cov_e_random", "pe_random", "clus_e", "n_clus_e")
  c_random_list <- c("X_c_random","mean_cov_c_random", "pc_random", "clus_c", "n_clus_c")
  se_random_list <- c("Z_e_random","mean_z_e_random", "ze_random", "clus_se", "n_clus_se")
  sc_random_list <- c("Z_c_random","mean_z_c_random", "zc_random", "clus_sc", "n_clus_sc")  
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
  ind_random <- model_info$ind_random
  ind_fixed <- model_info$ind_fixed
  if(ze_fixed == 0) { sde <- NULL; se <- NULL
  ze_fixed_index <- match(se_fixed_list, names(datalist))
  datalist <- datalist[-ze_fixed_index]}
  if(zc_fixed == 0) { sdc <- NULL; sc <- NULL
  zc_fixed_index <- match(sc_fixed_list, names(datalist))
  datalist <- datalist[-zc_fixed_index]}
  if(ze_fixed == 1) { ze_fixed_index <- match("ze_fixed", names(datalist))
  datalist <- datalist[-ze_fixed_index]}
  if(zc_fixed == 1) { zc_fixed_index <- match("zc_fixed", names(datalist))
  datalist <- datalist[-zc_fixed_index]}
  if(ze_random == 1) { ze_random_index <- match("ze_random", names(datalist))
  datalist <- datalist[-ze_random_index]}
  if(zc_random == 1) { zc_random_index <- match("zc_random", names(datalist))
  datalist <- datalist[-zc_random_index]}
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
  if(length(model_info$model_se_random) == 0) { 
    se_random_index <- match(se_random_list, names(datalist))
    datalist <- datalist[-se_random_index] 
  }
  if(length(model_info$model_sc_random) == 0) { 
    sc_random_index <- match(sc_random_list, names(datalist))
    datalist <- datalist[-sc_random_index] 
  }  
  params <- c("eff", "cost", "s_eff", "s_cost", "tmu_e", "tmu_c",
              "s_e", "s_c", "p_e", "p_c", "beta", "alpha", "gamma_e", "gamma_c")
  params <- c(params, "loglik_e", "loglik_c", "loglik_se", "loglik_sc")
  if(ze_fixed == 0) { params_se_index <- match(c("s_eff", "loglik_se", "gamma_e", "p_e"), params)
    params <- params[-params_se_index]}
  if(zc_fixed == 0) { params_sc_index <- match(c("s_cost", "loglik_sc", "gamma_c", "p_c"), params)
  params <- params[-params_sc_index]}
  if(!ind_fixed) { params <- c(params, "beta_f")}
  if(length(model_info$model_e_random) != 0){ params <- c(params, "a")}
  if(length(model_info$model_c_random) != 0 & !is_c_random_c){ params <- c(params, "b")}
  if(length(model_info$model_se_random) != 0){ params <- c(params, "g_e")}
  if(length(model_info$model_sc_random) != 0){ params <- c(params, "g_c")}
  if(length(model_info$model_c_random) != 0 & !ind_random) { params <- c(params, "b_f")}  
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
                         "ze_random" = ze_random, "zc_random" = zc_random, "sde" = sde, "sdc" = sdc, "se" = se, "sc" = sc,
                         "ind_random" = ind_random, "model_e_random" = model_info$model_e_random, "model_c_random" = model_info$model_c_random, 
                         "model_se_random" = model_info$model_se_random, "model_sc_random" = model_info$model_sc_random,
                         "is_c_random_c" = is_c_random_c, "is_int_c_random_c" = is_int_c_random_c, 
                         "prior" = model_info$prior)
  filein <- write_hurdle(type = type , dist_e = dist_e, dist_c = dist_c, model_txt_info = model_txt_info)
  if(dist_c %in% c("gamma", "lnorm")) {
    if(!is.null(datalist$sc)) {
      if(data_model$sc == 0) { 
        datalist$sc = log(0.0000001)
        index_c_0 <- which(datalist$cost == 0)
        datalist$cost[index_c_0] = 0.0000001
      } else { datalist$sc = log(datalist$sc)}
    }
  }
  if(dist_e %in% c("gamma", "weib", "exp", "pois", "negbin")) {
    if(!is.null(data_model$se)) {
      if(data_model$se == 0) { 
        datalist$se = log(0.0000001)
        index_e_0 <- which(datalist$eff == 0)
        datalist$eff[index_e_0] = 0.0000001
      } else { datalist$se = log(datalist$se)}
    }
  }
  if(dist_e == "beta") {
    if(!is.null(data_model$se)) {
      if(data_model$se == 1) { 
        datalist$se = qlogis(1 - 0.0000001)
        index_e_1 <- which(datalist$eff == 1)
        datalist$eff[index_e_1] = 1 - 0.0000001
      } else if(data_model$se == 0) { 
        datalist$se = qlogis(0 + 0.0000001)
        index_e_0 <- which(datalist$eff == 0)
        datalist$eff[index_e_0] = 0 + 0.0000001
      } else { datalist$se = qlogis(datalist$se)}
    }
  }
  model <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, 
                        model.file = filein, n.chains = n.chains, n.iter = n.iter, 
                        n.burnin = n.burnin, DIC = DIC, pD = pd, n.thin = n.thin,
                        n.iter.pd = n.iter.pd, n.adapt = n.adapt)  
  tmu_e <- model$BUGSoutput$sims.list$tmu_e
  tmu_c <- model$BUGSoutput$sims.list$tmu_c
  s_e <- model$BUGSoutput$sims.list$s_e
  if(dim(s_e)[2] == 1) { colnames(s_e) <- "s_e"}
  if(dim(s_e)[2] == 2) { 
    s_e <- as.matrix(model$BUGSoutput$sims.list$s_e[, 1]) 
    colnames(s_e) <- "s_e"}
  s_c <- model$BUGSoutput$sims.list$s_c
  if(dim(s_c)[2] == 1) { colnames(s_c) <- "s_c"}
  if(dim(s_c)[2] == 2) { 
    s_c <- as.matrix(model$BUGSoutput$sims.list$s_c[, 1])
    colnames(s_c) <- "s_c"}
  alpha <- model$BUGSoutput$sims.list$alpha
  if(length(dim(alpha)) == 2) { colnames(alpha) <- colnames(as.data.frame(X_e_fixed))}
  if(length(dim(alpha)) == 3) { alpha <- as.matrix(alpha[, , 1]); colnames(alpha) <- colnames(as.data.frame(X_e_fixed))}
  beta <- model$BUGSoutput$sims.list$beta
  if(length(dim(beta)) == 2) { colnames(beta) <- colnames(as.data.frame(X_c_fixed))}
  if(length(dim(beta)) == 3) { beta <- as.matrix(beta[, , 1]); colnames(beta) <- colnames(as.data.frame(X_c_fixed))}
  if(ze_fixed != 0) {
    p_e <- model$BUGSoutput$sims.list$p_e
    gamma_e <- model$BUGSoutput$sims.list$gamma_e
    colnames(p_e) <- "p_e"
    if(is.vector(Z_e_fixed)) {
      colnames(gamma_e) <- "(Intercept)"
    } else { colnames(gamma_e) <- colnames(as.data.frame(Z_e_fixed))}
  } else { p_e <- gamma_e <- NULL}
  if(zc_fixed != 0) {
    p_c <- model$BUGSoutput$sims.list$p_c
    gamma_c <- model$BUGSoutput$sims.list$gamma_c
    colnames(p_c) <- "p_c"
    if(is.vector(Z_c_fixed)) {
      colnames(gamma_c) <- "(Intercept)"
    } else { colnames(gamma_c) <- colnames(as.data.frame(Z_c_fixed))}
  } else { p_c <- gamma_c <- NULL}
  if(dist_e %in% c("norm", "negbin", "logis")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    tau_e <- model$BUGSoutput$sims.list$tau_e
    if(dim(tau_e)[2] == 2) { 
      tau_e <- as.matrix(model$BUGSoutput$sims.list$tau_e[, 1])}
  } else if(dist_e %in% c("beta", "gamma", "weib")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    ctau_e <- model$BUGSoutput$sims.list$ctau_e
  } else if(dist_e %in% c("exp", "bern", "pois")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
  } 
  if(dist_c == "norm") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    tau_c <- model$BUGSoutput$sims.list$tau_c
    if(dim(tau_c)[2] == 2) { 
      tau_c <- as.matrix(model$BUGSoutput$sims.list$tau_c[, 1])}
  } else if(dist_c == "gamma") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    ctau_c <- model$BUGSoutput$sims.list$ctau_c
  } else if(dist_c == "lnorm") {
    clmu_c <- model$BUGSoutput$sims.list$clmu_c
    ltau_c <- model$BUGSoutput$sims.list$ltau_c
    if(dim(ltau_c)[2] == 2) { 
      ltau_c <- as.matrix(model$BUGSoutput$sims.list$ltau_c[, 1])}
  } 
  S <- model_info$n.mci
  trt_index <- c(data_model$trt_index)
  mu_e <- mu_c <- matrix(NA, nrow = dim(beta)[1], ncol = length(data_model$n_trt))
  for(i in 1:n.iter) {
    for(trt in 1:length(data_model$n_trt)) {
      if(!is.null(data_model$se)) {
        if(dist_e == "norm") {
          mu_e[i, trt] <- mean(rnorm(S, mean = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se], 
                                     sd = s_e[i]), na.rm = TRUE)
        }
        if(dist_e == "negbin") {
          mu_e[i, trt] <- mean(rnbinom(S, prob = tau_e[i]/(tau_e[i] + cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), 
                                       size = tau_e[i]), na.rm = TRUE)
        }  
        if(dist_e == "logis") {
          mu_e[i, trt] <- mean(rlogis(S, location = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se], 
                                      scale = 1 / tau_e[i]), na.rm = TRUE)
        }
        if(dist_e == "beta") {
          mu_e[i, trt] <- mean(rbeta(S, shape1 = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se] * ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se], 
                                     shape2 = (1 - cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]) * ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), na.rm = TRUE)
        }
        if(dist_e == "gamma") {
          mu_e[i, trt] <- mean(rgamma(S, shape = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se] * ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se], 
                                      rate = ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), na.rm = TRUE)
        }
        if(dist_e == "weib") {
          mu_e[i, trt] <- mean(rweibull(S, shape = ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), 
                               scale = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se] / exp(lgamma(1 + 1 / ctau_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se])), na.rm = TRUE)
        }
        if(dist_e == "exp") {
          mu_e[i, trt] <- mean(rexp(S, rate = 1 / cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), na.rm = TRUE)
        }
        if(dist_e == "bern") {
          mu_e[i, trt] <- mean(rbinom(S, size = 1, prob = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), na.rm = TRUE)
        }
        if(dist_e == "pois") {
          mu_e[i, trt] <- mean(rpois(S, lambda = cmu_e[i, unlist(data_model$trt_index[trt])][cmu_e[i, unlist(data_model$trt_index[trt])] != se]), na.rm = TRUE)
        }
      }
      if(is.null(data_model$se)) {
      if(dist_e == "norm") {
        mu_e[i, trt] <- mean(rnorm(S, mean = cmu_e[i, unlist(data_model$trt_index[trt])], sd = s_e[i]), na.rm = TRUE)
      }
      if(dist_e == "negbin") {
        mu_e[i, trt] <- mean(rnbinom(S, prob = tau_e[i]/(tau_e[i] + cmu_e[i, unlist(data_model$trt_index[trt])]), size = tau_e[i]), na.rm = TRUE)
      }  
      if(dist_e == "logis") {
        mu_e[i, trt] <- mean(rlogis(S, location = cmu_e[i, unlist(data_model$trt_index[trt])], scale = 1 / tau_e[i]), na.rm = TRUE)
      }
      if(dist_e == "beta") {
        mu_e[i, trt] <- mean(rbeta(S, shape1 = cmu_e[i, unlist(data_model$trt_index[trt])] * ctau_e[i, unlist(data_model$trt_index[trt])], 
                                   shape2 = (1 - cmu_e[i, unlist(data_model$trt_index[trt])]) * ctau_e[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
      }
      if(dist_e == "gamma") {
        mu_e[i, trt] <- mean(rgamma(S, shape = cmu_e[i, unlist(data_model$trt_index[trt])] * ctau_e[i, unlist(data_model$trt_index[trt])], 
                                    rate = ctau_e[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
      }
      if(dist_e == "weib") {
        mu_e[i, trt] <- mean(rweibull(S, shape = ctau_e[i, unlist(data_model$trt_index[trt])]), scale = cmu_e[i, unlist(data_model$trt_index[trt])] / exp(lgamma(1 + 1 / ctau_e[i, unlist(data_model$trt_index[trt])])), na.rm = TRUE)
      }
      if(dist_e == "exp") {
        mu_e[i, trt] <- mean(rexp(S, rate = 1 / cmu_e[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
      }
      if(dist_e == "bern") {
        mu_e[i, trt] <- mean(rbinom(S, size = 1, prob = cmu_e[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
      }
      if(dist_e == "pois") {
        mu_e[i, trt] <- mean(rpois(S, lambda = cmu_e[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
      }
      }
      if(!is.null(data_model$sc)) {
        if(dist_c == "norm") {
          mu_c[i, trt] <- mean(rnorm(S, mean = cmu_c[i, unlist(data_model$trt_index[trt])][cmu_c[i, unlist(data_model$trt_index[trt])] != sc], 
                                     sd = s_c[i, 1]), na.rm = TRUE)
        }
        if(dist_c == "gamma") {
          mu_c[i, trt] <- mean(rgamma(S, shape = cmu_c[i, unlist(data_model$trt_index[trt])][cmu_c[i, unlist(data_model$trt_index[trt])] != sc] * ctau_c[i, unlist(data_model$trt_index[trt])][cmu_c[i, unlist(data_model$trt_index[trt])] != sc], 
                                      rate = ctau_c[i, unlist(data_model$trt_index[trt])][cmu_c[i, unlist(data_model$trt_index[trt])] != sc]), na.rm = TRUE)
        }
        if(dist_c == "lnorm") {
          mu_c[i, trt] <- mean(rlnorm(S, meanlog = clmu_c[i, unlist(data_model$trt_index[trt])][clmu_c[i, unlist(data_model$trt_index[trt])] != sc], sdlog = 1 / sqrt(ltau_c[i, 1])), na.rm = TRUE)
        }
      }
      if(is.null(data_model$sc)) {
        if(dist_c == "norm") {
          mu_c[i, trt] <- mean(rnorm(S, mean = cmu_c[i, unlist(data_model$trt_index[trt])], sd = s_c[i]), na.rm = TRUE)
        }
        if(dist_c == "gamma") {
          mu_c[i, trt] <- mean(rgamma(S, shape = cmu_c[i, unlist(data_model$trt_index[trt])] * ctau_c[i, unlist(data_model$trt_index[trt])], 
                                      rate = ctau_c[i, unlist(data_model$trt_index[trt])]), na.rm = TRUE)
        }
        if(dist_c == "lnorm") {
          mu_c[i, trt] <- mean(rlnorm(S, meanlog = clmu_c[i, unlist(data_model$trt_index[trt])], sdlog = 1 / sqrt(ltau_c[i])), na.rm = TRUE)
        }
      }
    }
  }
  if(!is.null(data_model$se)) {
    mu_ns_e <- mu_e
    mu_e <- mu_ns_e * (1 - as.vector(p_e)) + se * as.vector(p_e)
    colnames(mu_ns_e) <- data_model$trt_lev
  } else { mu_ns_e <- NULL}
  if(!is.null(data_model$sc)) {
    mu_ns_c <- mu_c
    mu_c <- mu_ns_c * (1 - as.vector(p_c)) + sc * as.vector(p_c)
    colnames(mu_ns_c) <- data_model$trt_lev
  } else { mu_ns_c <- NULL}
  colnames(mu_e) <- colnames(mu_c) <- data_model$trt_lev
  a <- b <- g_e <- g_c <- NULL 
  if(length(model_info$model_e_random) != 0) { 
    a <- model$BUGSoutput$sims.list$a
    if(is.vector(X_e_random)) {
      colnames(a) <- paste("(Intercept).", data_model$clus_e_lev, sep = "")
    } else if (length(dim(a)) == 2) {
      colnames(a) <- paste(colnames(as.data.frame(X_e_random)), data_model$clus_e_lev, sep = ".")
    } else if(length(dim(a)) == 3) {
      dimnames(a)[[2]] <- colnames(as.data.frame(X_e_random))
      dimnames(a)[[3]] <- data_model$clus_e_lev
    }
  }
  if(length(model_info$model_c_random) != 0 & !is_c_random_c) { 
    b <- model$BUGSoutput$sims.list$b
    if(is.vector(X_c_random)) {
      colnames(b) <- paste("(Intercept).", data_model$clus_c_lev, sep = "")
    } else if (length(dim(b)) == 2) {
      colnames(b) <- paste(colnames(as.data.frame(X_c_random)), data_model$clus_c_lev, sep = ".")
    } else if(length(dim(b)) == 3) {
      dimnames(b)[[2]] <- colnames(as.data.frame(X_c_random))
      dimnames(b)[[3]] <- data_model$clus_c_lev
    }
  }
  if(length(model_info$model_se_random) != 0) { 
    g_e <- model$BUGSoutput$sims.list$g_e
    if(is.vector(Z_e_random)) {
      colnames(g_e) <- paste("(Intercept).", data_model$clus_se_lev, sep = "")
    } else if (length(dim(g_e)) == 2) {
      colnames(g_e) <- paste(colnames(as.data.frame(Z_e_random)), data_model$clus_se_lev, sep = ".")
    } else if(length(dim(g_e)) == 3) {
      dimnames(g_e)[[2]] <- colnames(as.data.frame(Z_e_random))
      dimnames(g_e)[[3]] <- data_model$clus_se_lev
    }
  }
  if(length(model_info$model_mc_random) != 0) { 
    g_c <- model$BUGSoutput$sims.list$g_c
    if(is.vector(Z_c_random)) {
      colnames(g_c) <- paste("(Intercept).", data_model$clus_sc_lev, sep = "")
    } else if (length(dim(g_c)) == 2) {
      colnames(g_c) <- paste(colnames(as.data.frame(Z_c_random)), data_model$clus_sc_lev, sep = ".")
    } else if(length(dim(g_c)) == 3) {
      dimnames(g_c)[[2]] <- colnames(as.data.frame(Z_c_random))
      dimnames(g_c)[[3]] <- data_model$clus_sc_lev
    }
  }  
  eff_pos <- model$BUGSoutput$sims.list$eff 
  cost_pos <- model$BUGSoutput$sims.list$cost 
  eff_pos_imp_avg <- apply(model$BUGSoutput$sims.list$eff, 2, mean, na.rm = T)[is.na(data_model$eff)] 
  cost_pos_imp_avg <- apply(model$BUGSoutput$sims.list$cost, 2, mean, na.rm = T)[is.na(data_model$cost)] 
  eff_pos_imp_ql <- apply(model$BUGSoutput$sims.list$eff, 2, quantile, prob = model_info$prob[1], na.rm = T)[is.na(data_model$eff)] 
  eff_pos_imp_qu <- apply(model$BUGSoutput$sims.list$eff, 2, quantile, prob = model_info$prob[2], na.rm = T)[is.na(data_model$eff)] 
  cost_pos_imp_ql <- apply(model$BUGSoutput$sims.list$cost, 2, quantile, prob = model_info$prob[1], na.rm = T)[is.na(data_model$cost)] 
  cost_pos_imp_qu <- apply(model$BUGSoutput$sims.list$cost, 2, quantile, prob = model_info$prob[2], na.rm = T)[is.na(data_model$cost)]   
  efft_pos <- costt_pos <- list()
  efft_pos_imp_avg <- efft_pos_imp_ql <- efft_pos_imp_qu <- list()
  costt_pos_imp_avg <- costt_pos_imp_ql <- costt_pos_imp_qu <- list()  
  for(i in data_model$trt_lev) {
    efft_pos[[i]] <- eff_pos[, trt_index[[i]]]
    efft_pos_imp_avg[[i]] <- apply(efft_pos[[i]][, data_model$m_efft[[i]] == 1], 2, mean, na.rm = TRUE)
    efft_pos_imp_ql[[i]] <- apply(efft_pos[[i]][, data_model$m_efft[[i]] == 1], 2, quantile, prob = model_info$prob[1], na.rm = TRUE)
    efft_pos_imp_qu[[i]] <- apply(efft_pos[[i]][, data_model$m_efft[[i]] == 1], 2, quantile, prob = model_info$prob[2], na.rm = TRUE)
    costt_pos[[i]] <- cost_pos[, trt_index[[i]]]
    costt_pos_imp_avg[[i]] <- apply(costt_pos[[i]][, data_model$m_costt[[i]] == 1], 2, mean, na.rm = TRUE)
    costt_pos_imp_ql[[i]] <- apply(costt_pos[[i]][, data_model$m_costt[[i]] == 1], 2, quantile, prob = model_info$prob[1], na.rm = TRUE)
    costt_pos_imp_qu[[i]] <- apply(costt_pos[[i]][, data_model$m_costt[[i]] == 1], 2, quantile, prob = model_info$prob[2], na.rm = TRUE)
  }
  loglik_e <- model$BUGSoutput$sims.list$loglik_e
  loglik_c <- model$BUGSoutput$sims.list$loglik_c
  if(ze_fixed != 0) { loglik_se <- model$BUGSoutput$sims.list$loglik_se
  } else { loglik_se <- NULL}
  if(zc_fixed != 0) { loglik_sc <- model$BUGSoutput$sims.list$loglik_sc
  } else { loglik_sc <- NULL}
  if(!ind_fixed) {
    beta_f <- model$BUGSoutput$sims.list$beta_f
    if(dim(beta_f)[2] == 1) { colnames(beta_f) <- "e"}
    if(dim(beta_f)[2] == 2) { 
      beta_f <- as.matrix(model$BUGSoutput$sims.list$beta_f[, 1])
      colnames(beta_f) <- "e"}
    beta <- list("beta" = beta, "beta_f" = beta_f)}
  if(length(model_info$model_c_random) != 0 & !ind_random) {
    b_f <- model$BUGSoutput$sims.list$b_f
    colnames(b_f) <- paste("e", data_model$clus_c_lev, sep = ".")
    if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & is_c_random_c) {
      b <- list("b_f" = b_f) 
    } else if(length(model_info$model_c_random) != 0 & "e" %in% model_info$model_c_random & !is_c_random_c) {
      b <- list("b" = b, "b_f" = b_f)
    }
  }  
  if(n.chains > 1) {
    param_jags_show <- c("eff", "cost", "s_eff", "s_cost", "loglik_e", "loglik_c", "loglik_se", "loglik_sc", 
                         "cmu_c", "cmu_e")
    if(ze_fixed == 0) { params_jags_show_index <- match(c("s_eff", "loglik_se"), param_jags_show)
    param_jags_show <- param_jags_show[-params_jags_show_index]}
    if(zc_fixed == 0) { params_jags_show_index <- match(c("s_cost", "loglik_sc"), param_jags_show)
    param_jags_show <- param_jags_show[-params_jags_show_index]}
    if(dist_e %in% c("beta", "gamma", "weib")) { param_jags_show <- c(param_jags_show, "ctau_e")}
    if(dist_c %in% c("gamma")) { param_jags_show <- c(param_jags_show, "ctau_c")}
    if(dist_c %in% c("lnorm")) { param_jags_show <- c(param_jags_show, "clmu_c")}
    model_sum <- round(jagsresults(x = model, params = param_jags_show, invert = TRUE), digits = 3)
  } else{ model_sum <- NULL} 
  loglik <- list("effects" = loglik_e, "costs" = loglik_c, 
                 "structural indicators effects" = loglik_se, 
                 "structural indicators costs" = loglik_sc)
  imputed_e <- list("avg" = eff_pos_imp_avg, "ql" = eff_pos_imp_ql, "qu" = eff_pos_imp_qu, "index" = which(is.na(data_model$eff)))
  imputed_c <- list("avg" = cost_pos_imp_avg, "ql" = cost_pos_imp_ql, "qu" = cost_pos_imp_qu, "index" = which(is.na(data_model$cost)))
  imputed_et <- list("avg" = efft_pos_imp_avg, "ql" = efft_pos_imp_ql, "qu" = efft_pos_imp_avg, "index" = lapply(data_model$m_efft, function(x) which(x == 1)))
  imputed_ct <- list("avg" = costt_pos_imp_avg, "ql" = costt_pos_imp_ql, "qu" = costt_pos_imp_avg, "index" = lapply(data_model$m_costt, function(x) which(x == 1)))
  imputed <- list("effects" = imputed_e, "costs" = imputed_c)
  imputedt <- list("effects" = imputed_et, "costs" = imputed_ct)
  model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c,
                            "mean_nostr_effects" = mu_ns_e, "mean_nostr_costs" = mu_ns_c,
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "pstr_effects" = p_e, "str_effects_fixed" = gamma_e, "pstr_costs" = p_c, "str_costs_fixed" = gamma_c, 
                              "effects_random" = a, "costs_random" = b, "str_effects_random" = g_e, "str_costs_random" = g_c, 
                              "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, "type" = type, "method" = "HURDLE", 
                              "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, 
                              "se" = datalist$se, "sc" = datalist$sc, "filein" = filein)
  if(n.chains == 1) { model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}