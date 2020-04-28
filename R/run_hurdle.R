#' An internal function to execute a JAGS hurdle model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#' and Structural At Random (SAR).
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param se Structural value to be found in the effect data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters
#' @param sde hyper-prior value for the standard deviation of the distribution of the structural effects. The default value is
#' \code{1.0E-6} to approximate a point mass at the structural value provided by the user.
#' @param sdc hyper-prior value for the standard deviation of the distribution of the structural costs. The default value is
#' \code{1.0E-6} to approximate a point mass at the structural value provided by the user.
#' @param ppc Logical. If \code{ppc} is \code{TRUE}, the estimates of the parameters that can be used to generate replications from the model are saved.
#' @keywords JAGS Bayesian hurdle models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_hurdle <- function(type, dist_e, dist_c, inits, se, sc, sde, sdc, ppc) eval.parent(substitute( {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  if(!dist_e %in% c("norm", "beta", "exp", "bern", "nbinom", "weibull", "gamma", "logis", "pois") | !dist_c %in% c("norm", "gamma", "lnorm")) {
    stop("Distributions available for use are 'norm', 'beta', 'gamma', 'weibull', 'logis', 'exp', 'bern', 'pois', 'nbinom' for the effects and 'norm', 'gamma', 'lnorm' for the costs")
  }
  if(!type %in% c("SCAR", "SAR")) {
    stop("Types available for use are 'SCAR', 'SAR'")
  }
  if(is.null(inits) == FALSE) {inits = inits }
  model <- write_hurdle(type = type , dist_e = dist_e, dist_c = dist_c, pe_fixed = pe_fixed, pc_fixed = pc_fixed, ze_fixed = ze_fixed, zc_fixed = zc_fixed, ind_fixed = ind_fixed,
                        pe_random = pe_random, pc_random = pc_random, ze_random = ze_random, zc_random = zc_random, ind_random = ind_random, 
                        model_e_random = model_e_random, model_c_random = model_c_random, model_se_random = model_se_random, model_sc_random = model_sc_random, se = se, sc = sc)
  filein <- model
  if(dist_e %in% c("norm")) {sde <- log(sde) }
  if(dist_c %in% c("norm")) {sdc <- log(sdc) }
  if(dist_c %in% c("gamma", "lnorm")) {
  if(is.null(sc) == FALSE) {
    if(sc == 0) {
      sc = log(0.0000001)
      if(any(which(cost1 == 0)) == TRUE) {
        index_c1_0 <- which(cost1 == 0)
        cost1[index_c1_0] = 0.0000001
      }
      if(any(which(cost2 == 0)) == TRUE) {
        index_c2_0 <- which(cost2 == 0)
        cost2[index_c2_0] = 0.0000001
      }
    } else {
      sc = log(sc)
    }
   }
  }
  if(dist_e %in% c("gamma", "weibull", "exp", "pois", "nbinom")) {
    if(is.null(se) == FALSE) {
      if(se == 0) {
        se = log(0.0000001)
        if(any(which(eff1 == 0)) == TRUE) {
          index_e1_0 <- which(eff1 == 0)
          if(dist_e %in% c("gamma", "weibull", "exp")){
          eff1[index_e1_0] = 0.0000001
          }
        }
        if(any(which(eff2 == 0)) == TRUE) {
          index_e2_0 <- which(eff2 == 0)
          if(dist_e %in% c("gamma", "weibull", "exp")){
          eff2[index_e2_0] = 0.0000001
          }
        }
      } else {
        se = log(se)
      }
    }
  }
  if(dist_e == "beta") {
    if(is.null(se) == FALSE) {
      if(se == 1) {
        se = qlogis(1 - 0.0000001)
        if(any(which(eff1 == 1)) == TRUE) {
          index_e1_1 <- which(eff1 == 1)
          eff1[index_e1_1] = 1 - 0.0000001
        }
        if(any(which(eff2 == 1)) == TRUE) {
          index_e2_1 <- which(eff2 == 1)
          eff2[index_e2_1] = 1 - 0.0000001
        }
      } else if(se == 0){
        se = qlogis(0 + 0.0000001)
        if(any(which(eff1 == 0)) == TRUE) {
          index_e1_0 <- which(eff1 == 0)
          eff1[index_e1_0] = 0 + 0.0000001
        }
        if(any(which(eff2 == 0)) == TRUE) {
          index_e2_0 <- which(eff2 == 0)
          eff2[index_e2_0] = 0 + 0.0000001
        }
      } else {
        se = qlogis(se)
      }
    }
  }
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "d_eff1", "d_eff2", "d_cost1", "d_cost2", 
                   "X1_e_fixed", "X2_e_fixed", "X1_c_fixed", "X2_c_fixed", "Z1_e_fixed", "Z2_e_fixed", "Z1_c_fixed", "Z2_c_fixed", 
                   "mean_cov_e1_fixed", "mean_cov_e2_fixed", "mean_cov_c1_fixed", "mean_cov_c2_fixed", "mean_z_e1_fixed", "mean_z_e2_fixed", "mean_z_c1_fixed", "mean_z_c2_fixed", 
                   "pe_fixed", "pc_fixed", "ze_fixed", "zc_fixed", "X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", 
                   "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e", "X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", 
                   "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c", "Z1_e_random", "Z2_e_random","mean_z_e1_random", "mean_z_e2_random", "ze_random", 
                   "clus1_se", "clus2_se", "n1_clus_se", "n2_clus_se", "Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", 
                   "clus1_sc", "clus2_sc", "n1_clus_sc", "n2_clus_sc", "se", "sc", "sde", "sdc")
  e_random_list <- c("X1_e_random", "X2_e_random","mean_cov_e1_random", "mean_cov_e2_random", "pe_random", "clus1_e", "clus2_e", "n1_clus_e", "n2_clus_e")
  c_random_list <- c("X1_c_random", "X2_c_random","mean_cov_c1_random", "mean_cov_c2_random", "pc_random", "clus1_c", "clus2_c", "n1_clus_c", "n2_clus_c")
  se_fixed_list <- c("Z1_e_fixed", "Z2_e_fixed","mean_z_e1_fixed", "mean_z_e2_fixed", "ze_fixed")
  sc_fixed_list <- c("Z1_c_fixed", "Z2_c_fixed","mean_z_c1_fixed", "mean_z_c2_fixed", "zc_fixed")
  se_random_list <- c("Z1_e_random", "Z2_e_random","mean_z_e1_random", "mean_z_e2_random", "ze_random", "clus1_se", "clus2_se", "n1_clus_se", "n2_clus_se")
  sc_random_list <- c("Z1_c_random", "Z2_c_random","mean_z_c1_random", "mean_z_c2_random", "zc_random", "clus1_sc", "clus2_sc", "n1_clus_sc", "n2_clus_sc")
  if(pe_fixed == 1) {pe_fixed_index <- match("pe_fixed", datalist)
  datalist <- datalist[-pe_fixed_index] }
  if(pc_fixed == 1) {pc_fixed_index <- match("pc_fixed", datalist)
  datalist <- datalist[-pc_fixed_index] }
  if(ze_fixed == 0) {ze_fixed_index <- match(se_fixed_list, datalist)
  datalist <- datalist[-ze_fixed_index] }
  if(zc_fixed == 0) {zc_fixed_index <- match(sc_fixed_list, datalist)
  datalist <- datalist[-zc_fixed_index] }
  if(ze_fixed == 1) {ze_fixed_index <- match("ze_fixed", datalist)
  datalist <- datalist[-ze_fixed_index] }
  if(zc_fixed == 1) {zc_fixed_index <- match("zc_fixed", datalist)
  datalist <- datalist[-zc_fixed_index] }
  if(length(model_e_random) != 0) {
    if(pe_random == 1) {pe_random_index <- match("pe_random", datalist)
    datalist <- datalist[-pe_random_index] }
  } else if(length(model_e_random) == 0) { e_random_index <- match(e_random_list, datalist)
  datalist <- datalist[-e_random_index] }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE) {
    if(pc_random == 1 | length(model_c_random) != 0 & is_c_random_c == TRUE) {
      pc_random_index <- match("pc_random", datalist)
      datalist <- datalist[-pc_random_index] }
  } else if(length(model_c_random) != 0 & is_c_random_c == TRUE) {
    c_random_index <- match(c_random_list[1:5], datalist)
    datalist <- datalist[-c_random_index] 
  } else if(length(model_c_random) == 0) { 
    c_random_index <- match(c_random_list, datalist)
    datalist <- datalist[-c_random_index] 
  }
  if(length(model_se_random) != 0) {
    if(ze_random == 1) {ze_random_index <- match("ze_random", datalist)
    datalist <- datalist[-ze_random_index] }
  } else if(length(model_se_random) == 0) { se_random_index <- match(se_random_list, datalist)
  datalist <- datalist[-se_random_index] }
  if(length(model_sc_random) != 0) {
    if(zc_random == 1) {zc_random_index <- match("zc_random", datalist)
    datalist <- datalist[-zc_random_index] }
  } else if(length(model_sc_random) == 0) { sc_random_index <- match(sc_random_list, datalist)
  datalist <- datalist[-sc_random_index] }
  if(is.null(se) == TRUE) {
    d_eff1_index <- match("d_eff1", datalist)
    d_eff2_index <- match("d_eff2", datalist)
    se_index <- match("se", datalist)
    sde_index <- match("sde", datalist)
    datalist <- datalist[-c(d_eff1_index, d_eff2_index, se_index, sde_index)]
  }
  if(is.null(sc) == TRUE) {
    d_cost1_index <- match("d_cost1", datalist)
    d_cost2_index <- match("d_cost2", datalist)
    sc_index <- match("sc", datalist)
    sdc_index <- match("sdc", datalist)
    datalist <- datalist[-c(d_cost1_index, d_cost2_index, sc_index, sdc_index)]
  }
  DIC <- TRUE
  if(is.null(se) == TRUE) {params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "s_e", "s_c", "p_c", "beta", "alpha", "gamma_c", 
                                       "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_dc1", "loglik_dc2") }
  if(is.null(sc) == TRUE) {params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "s_e", "s_c", "p_e", "beta", "alpha", "gamma_e", 
                                       "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_de1", "loglik_de2") }
  if(is.null(se) == FALSE & is.null(sc) == FALSE) {
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "s_e", "s_c", "p_e", "p_c", "beta", "alpha", "gamma_e", "gamma_c", 
              "loglik_e1", "loglik_e2", "loglik_c1", "loglik_c2", "loglik_de1", "loglik_de2", "loglik_dc1", "loglik_dc2") 
  }
  if(ind_fixed == FALSE) {params <- c(params, "beta_f") }
  if(length(model_e_random) != 0){params <- c(params, "a1", "a2") }
  if(length(model_c_random) != 0 & is_c_random_c == FALSE){params <- c(params, "b1", "b2") }
  if(length(model_se_random) != 0){params <- c(params, "g1_e", "g2_e") }
  if(length(model_sc_random) != 0){params <- c(params, "g1_c", "g2_c") }
  if(length(model_c_random) != 0 & ind_random == FALSE) {params <- c(params, "b1_f", "b2_f") }
  if(ppc == TRUE) { 
    if(dist_e %in% c("norm", "nbinom", "logis")) {
      ppc_e_params <- c("mu_e1", "mu_e2", "tau_e1", "tau_e2") 
    } else if(dist_e %in% c("beta", "gamma", "weibull")) {
      ppc_e_params <- c("mu_e1", "tau_e1", "mu_e2", "tau_e2")
    } else if(dist_e %in% c("exp", "bern", "pois")) {
      ppc_e_params <- c("mu_e1", "mu_e2")
    }
    if(dist_c %in% c("norm")) {
      ppc_c_params <- c("mu_c1", "mu_c2", "tau_c1", "tau_c2") 
    } else if(dist_c %in% c("gamma")) {
      ppc_c_params <- c("mu_c1", "tau_c1", "mu_c2", "tau_c2")
    } else if(dist_c %in% c("lnorm")) {
      ppc_c_params <- c("lmu_c1", "lmu_c2", "ltau_c1", "ltau_c2")
    } 
    params <- c(params, ppc_e_params, ppc_c_params)
  }
  modelN1 <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, model.file = filein, n.chains = n.chains, 
                          n.iter = n.iter, n.burnin = n.burnin, DIC = DIC, n.thin = n.thin)
    mu_e <- modelN1$BUGSoutput$sims.list$mu_e
    mu_c <- modelN1$BUGSoutput$sims.list$mu_c
    s_e <- modelN1$BUGSoutput$sims.list$s_e
    s_c <- modelN1$BUGSoutput$sims.list$s_c
    alpha <- modelN1$BUGSoutput$sims.list$alpha
    beta <- modelN1$BUGSoutput$sims.list$beta
    loglik_e1 <- modelN1$BUGSoutput$sims.list$loglik_e1
    loglik_e2 <- modelN1$BUGSoutput$sims.list$loglik_e2
    loglik_c1 <- modelN1$BUGSoutput$sims.list$loglik_c1
    loglik_c2 <- modelN1$BUGSoutput$sims.list$loglik_c2
    a <- b <- g_e <- g_c <- NULL
    if(length(model_e_random) != 0) { 
      a1 <- modelN1$BUGSoutput$sims.list$a1
      a2 <- modelN1$BUGSoutput$sims.list$a2 
      a <- list("a1" = a1, "a2" = a2) 
    }
    if(length(model_c_random) != 0 & is_c_random_c == FALSE) { 
      b1 <- modelN1$BUGSoutput$sims.list$b1
      b2 <- modelN1$BUGSoutput$sims.list$b2 
      b <- list("b1" = b1, "b2" = b2) 
    }
    if(is.null(se) == TRUE & is.null(sc) == FALSE) {
      p_c <- modelN1$BUGSoutput$sims.list$p_c
      gamma_c <- modelN1$BUGSoutput$sims.list$gamma_c
      loglik_dc1 <- modelN1$BUGSoutput$sims.list$loglik_dc1
      loglik_dc2 <- modelN1$BUGSoutput$sims.list$loglik_dc2
      if(length(model_sc_random) != 0) { 
        g1_c <- modelN1$BUGSoutput$sims.list$g1_c
        g2_c <- modelN1$BUGSoutput$sims.list$g2_c
        g_c <- list("g1_c" = g1_c, "g2_c" = g2_c) 
      }
    } else if(is.null(sc) == TRUE & is.null(se) == FALSE) {
      p_e <- modelN1$BUGSoutput$sims.list$p_e
      gamma_e <- modelN1$BUGSoutput$sims.list$gamma_e
      loglik_de1 <- modelN1$BUGSoutput$sims.list$loglik_de1
      loglik_de2 <- modelN1$BUGSoutput$sims.list$loglik_de2
      if(length(model_se_random) != 0) { 
        g1_e <- modelN1$BUGSoutput$sims.list$g1_e
        g2_e <- modelN1$BUGSoutput$sims.list$g2_e
        g_e <- list("g1_e" = g1_e, "g2_e" = g2_e) 
      }
    } else if(is.null(se) == FALSE & is.null(sc) == FALSE) {
      p_c <- modelN1$BUGSoutput$sims.list$p_c
      gamma_c <- modelN1$BUGSoutput$sims.list$gamma_c
      p_e <- modelN1$BUGSoutput$sims.list$p_e
      gamma_e <- modelN1$BUGSoutput$sims.list$gamma_e
      loglik_dc1 <- modelN1$BUGSoutput$sims.list$loglik_dc1
      loglik_dc2 <- modelN1$BUGSoutput$sims.list$loglik_dc2
      loglik_de1 <- modelN1$BUGSoutput$sims.list$loglik_de1
      loglik_de2 <- modelN1$BUGSoutput$sims.list$loglik_de2
      if(length(model_se_random) != 0) { 
        g1_e <- modelN1$BUGSoutput$sims.list$g1_e
        g2_e <- modelN1$BUGSoutput$sims.list$g2_e
        g_e <- list("g1_e" = g1_e, "g2_e" = g2_e) 
      }
      if(length(model_sc_random) != 0) { 
        g1_c <- modelN1$BUGSoutput$sims.list$g1_c
        g2_c <- modelN1$BUGSoutput$sims.list$g2_c
        g_c <- list("g1_c" = g1_c, "g2_c" = g2_c) 
      }
    }
    eff1_pos <- matrix(eff1, N1, 3)
    cost1_pos <- matrix(cost1, N1, 3)
    eff2_pos <- matrix(eff2, N2, 3)
    cost2_pos <- matrix(cost2, N2, 3)
    eff1_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, mean)
    eff1_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, quantile, probs = prob[1])
    eff1_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$eff1, 2, quantile, probs = prob[2])
    eff2_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, mean)
    eff2_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, quantile, probs = prob[1])
    eff2_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$eff2, 2, quantile, probs = prob[2])
    cost1_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, mean)
    cost1_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, quantile, probs = prob[1])
    cost1_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$cost1, 2, quantile, probs = prob[2])
    cost2_pos[, 1] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, mean)
    cost2_pos[, 2] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, quantile, probs = prob[1])
    cost2_pos[, 3] <- apply(modelN1$BUGSoutput$sims.list$cost2, 2, quantile, probs = prob[2])
    if(ind_fixed == FALSE) {
      beta_f <- modelN1$BUGSoutput$sims.list$beta_f
      beta <- list("beta" = beta, "beta_f" = beta_f)
    }
    if(length(model_c_random) != 0 & ind_random == FALSE) {
      b1_f <- modelN1$BUGSoutput$sims.list$b1_f
      b2_f <- modelN1$BUGSoutput$sims.list$b2_f
      if(length(model_c_random) != 0 & "e" %in% model_c_random & is_c_random_c == TRUE) {
        b <- list("b1_f" = b1_f, "b2_f" = b2_f) 
      } else if(length(model_c_random) != 0 & "e" %in% model_c_random & is_c_random_c == FALSE) {
        b <- list("b1" = b1, "b1_f" = b1_f, "b2" = b2, "b2_f" = b2_f)
      }
    }
  if(is.null(se) == FALSE & is.null(sc) == FALSE){
   if(n.chains >1 ) {
     model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_e1', 'loglik_e2',
                                                            'loglik_c1', 'loglik_c2', 'loglik_de1', 'loglik_de2',
                                                            'loglik_dc1', 'loglik_dc2'), invert = TRUE), digits = 3)
   } else{model_sum <- NULL }
  } else if(is.null(se) == TRUE & is.null(sc) == FALSE) {
    if(n.chains >1 ) {
      model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_e1', 'loglik_e2',
                                                             'loglik_c1', 'loglik_c2','loglik_dc1', 'loglik_dc2'), invert = TRUE), digits = 3)
    } else{model_sum <- NULL }
  } else if(is.null(se) == FALSE & is.null(sc) == TRUE) {
    if(n.chains >1 ) {
      model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2', 'loglik_e1', 'loglik_e2',
                                                             'loglik_c1', 'loglik_c2','loglik_de1', 'loglik_de2'), invert = TRUE), digits = 3)
    } else{model_sum <- NULL }
  }
    loglik_e <- list("control" = loglik_e1, "intervention" = loglik_e2)
    loglik_c <- list("control" = loglik_c1, "intervention" = loglik_c2)
    loglik_de <- NULL
    loglik_dc <- NULL
    if(is.null(se) == FALSE) {
      loglik_de <- list("control" = loglik_de1, "intervention" = loglik_de2)
    }
    if(is.null(sc) == FALSE) {
      loglik_dc <- list("control" = loglik_dc1, "intervention" = loglik_dc2)
    }
    loglik <- list("effects" = loglik_e, "costs" = loglik_c, "structural indicators effects" = loglik_de, "structural indicators costs" = loglik_dc)
    colnames(eff1_pos) <- c("mean", "LB", "UB")
    colnames(eff2_pos) <- c("mean", "LB", "UB")
    colnames(cost1_pos) <- c("mean", "LB", "UB")
    colnames(cost2_pos) <- c("mean", "LB", "UB")
    imputed <- list("effects1" = eff1_pos, "effects2" = eff2_pos, "costs1" = cost1_pos, "costs2" = cost2_pos)
    if(is.null(se) == TRUE) {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "structural_probability_costs" = p_c, 
                                "structural_parameter_costs_fixed" = gamma_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b,
                                "structural_parameter_costs_random" = g_c, "loglik" = loglik, "imputed" = imputed, "type" = "HURDLE_c", "ind_fixed" = ind_fixed, "ind_random" = ind_random,
                                "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    } else if(is.null(sc) == TRUE) {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "structural_probability_effects" = p_e, 
                                "structural_parameter_effects_fixed" = gamma_e, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b,
                                "structural_parameter_effects_random" = g_e, "loglik" = loglik, "imputed" = imputed, "type" = "HURDLE_e", "ind_fixed" = ind_fixed, "ind_random" = ind_random,
                                "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    } else if(is.null(se) == FALSE & is.null(sc) == FALSE) {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects_fixed" = alpha, "covariate_parameter_costs_fixed" = beta, "structural_probability_effects" = p_e, "structural_probability_costs" = p_c,
                                "structural_parameter_effects_fixed" = gamma_e, "structural_parameter_costs_fixed" = gamma_c, "covariate_parameter_effects_random" = a, "covariate_parameter_costs_random" = b,
                                "structural_parameter_effects_random" = g_e, "structural_parameter_costs_random" = g_c, "loglik" = loglik, "imputed" = imputed, "type" = "HURDLE_ec", "ind_fixed" = ind_fixed, "ind_random" = ind_random,
                                "ppc" = ppc, "dist_e" = dist_e, "dist_c" = dist_c)
    }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))