#' An internal function to execute a JAGS pattern mixture model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} function and obtain posterior inferences.
#' @param data_model list containing the data for the model to be passed to JAGS.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR).
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm').
#' @param model_info list containing model and MCMC information to be passed to JAGS. 
#' @keywords JAGS Bayesian pattern mixture models 
#' @importFrom stats rnorm rbeta rgamma rlnorm rweibull rnbinom rbinom rpois rlogis rexp
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #

run_pattern <- function(data_model, type, dist_e, dist_c, model_info) {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  n <- data_model$n; trt_lev <- names(data_model$n_trt); pe_fixed <- data_model$pe_fixed; pc_fixed <- data_model$pc_fixed
  mean_cov_e_fixed <- data_model$mean_cov_e_fixed; mean_cov_c_fixed <- data_model$mean_cov_c_fixed
  mean_cov_e_random <- data_model$mean_cov_e_random; mean_cov_c_random <- data_model$mean_cov_c_random
  pe_random <- data_model$pe_random; pc_random <- data_model$pc_random
  n_clus_e <- data_model$n_clus_e; n_clus_c <- data_model$n_clus_c; n_patterns <- data_model$n_patterns
  trt_pos_e <- data_model$trt_pos_e; trt_pos_c <- data_model$trt_pos_c; n_patternst <- data_model$n_patternst
  eff <- data_model$eff; cost <- data_model$cost; m_eff <- data_model$m_eff; m_cost <- data_model$m_cost
  d_or <- data_model$d; dt <- data_model$dt; X_e_fixed <- as.matrix(data_model$X_e_fixed); X_c_fixed <- as.matrix(data_model$X_c_fixed)
  X_e_random <- data_model$X_e_random; X_c_random <- data_model$X_c_random
  if(pe_random == 1) { X_e_random <- as.vector(unlist(data_model$X_e_random))}
  if(pc_random == 1) { X_c_random <- as.vector(unlist(data_model$X_c_random))}
  if(pe_random > 1) { X_e_random <- as.matrix(data_model$X_e_random)} 
  if(pc_random > 1) { X_c_random <- as.matrix(data_model$X_c_random)}
  clus_e <- data_model$clus_e; clus_c <- data_model$clus_c
  pp <- rep(1, n_patterns); d_mod <- d_or
  if(all(unique(d_or) %in% c(1, 4)) & n_patterns == 2) { d_mod[d_mod == 4] <- 2}
  if(all(unique(d_or) %in% c(1, 3)) & n_patterns == 2) { d_mod[d_mod == 3] <- 2}
  if(all(unique(d_or) %in% c(2, 4)) & n_patterns == 2) { d_mod[d_mod == 2] <- 1; d_mod[d_mod == 4] <- 2}
  if(all(unique(d_or) %in% c(2, 3)) & n_patterns == 2) { d_mod[d_mod == 2] <- 1; d_mod[d_mod == 3] <- 2}
  if(all(unique(d_or) %in% c(3, 4)) & n_patterns == 2) { d_mod[d_mod == 3] <- 1; d_mod[d_mod == 4] <- 2}
  if(all(unique(d_or) %in% c(1, 2, 4)) & n_patterns == 3) { d_mod[d_mod == 4] <- 3}
  if(all(unique(d_or) %in% c(1, 3, 4)) & n_patterns == 3) { d_mod[d_mod == 3] <- 2; d_mod[d_mod == 4] <- 3}
  if(all(unique(d_or) %in% c(2, 3, 4)) & n_patterns == 4) { d_mod[d_mod == 2] <- 1; d_mod[d_mod == 3] <- 2; d_mod[d_mod == 4] <- 3}
  datalist <- list("n"= n, "eff" = eff, "cost" = cost, 
                   "X_e_fixed" = X_e_fixed, "X_c_fixed" = X_c_fixed, "n_patterns" = n_patterns,
                   "mean_cov_e_fixed" = mean_cov_e_fixed, "mean_cov_c_fixed" = mean_cov_c_fixed, 
                   "pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "pp" = pp, "d_mod" = d_mod,
                   "X_e_random" = X_e_random, "mean_cov_e_random" = mean_cov_e_random, "pe_random" = pe_random, 
                   "clus_e" = clus_e, "n_clus_e" = n_clus_e, "X_c_random" = X_c_random, "mean_cov_c_random" = mean_cov_c_random, "pc_random" = pc_random, 
                   "clus_c"= clus_c, "n_clus_c" = n_clus_c)
  e_random_list <- c("X_e_random","mean_cov_e_random", "pe_random", "clus_e", "n_clus_e")
  c_random_list <- c("X_c_random","mean_cov_c_random", "pc_random", "clus_c", "n_clus_c")  
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
  if(!is.null(model_info$patterns.prior)) {
    pp_index <- match("pp", names(datalist))
    datalist <- datalist[-pp_index]
  }
  params <- c("eff", "cost", "tmu_e", "tmu_c", "s_e_p", "s_c_p", "p_prob", 
              "beta_p", "alpha_p", "delta_e", "delta_c")
  params <- c(params, "loglik_e", "loglik_c", "loglik_d")
  if(!ind_fixed) { params <- c(params, "beta_f_p")}
  if(length(model_info$model_e_random) != 0){ params <- c(params, "a")}
  if(length(model_info$model_c_random) != 0 & !is_c_random_c){ params <- c(params, "b")}
  if(length(model_info$model_c_random) != 0 & !ind_random) { params <- c(params, "b_f")} 
  if(type %in% c("MNAR_cost", "MAR")) {
    deltae_index <- match("delta_e", params)
    params <- params[-deltae_index] 
  }
  if(type %in% c("MNAR_eff", "MAR")) {
    deltac_index <- match("delta_c", params)
    params <- params[-deltac_index]
  }
  if(dist_e %in% c("norm", "negbin", "logis")) {
    mci_e_params <- c("cmu_e", "tau_e_p") 
  } else if(dist_e %in% c("beta", "gamma", "weib")) {
    mci_e_params <- c("cmu_e", "ctau_e")
  } else if(dist_e %in% c("exp", "bern", "pois")) {
    mci_e_params <- c("cmu_e")
  } 
  if(dist_c == "norm") {
    mci_c_params <- c("cmu_c", "tau_c_p") 
  } else if(dist_c == "gamma") {
    mci_c_params <- c("cmu_c", "ctau_c")
  } else if(dist_c == "lnorm") {
    mci_c_params <- c("clmu_c", "ltau_c_p")
  }
  params <- c(params, mci_e_params, mci_c_params)
  model_txt_info <- list("pe_fixed" = pe_fixed, "pc_fixed" = pc_fixed, "n_patterns" = n_patterns,
                         "ind_fixed" = ind_fixed, "pe_random" = pe_random, "pc_random" = pc_random,
                         "ind_random" = ind_random, "model_e_random" = model_info$model_e_random, 
                         "model_c_random" = model_info$model_c_random, "is_c_random_c" = is_c_random_c, 
                         "is_int_c_random_c" = is_int_c_random_c, "prior" = model_info$prior, 
                         "restriction" = model_info$restriction, "d_or" = unique(d_or))
  filein <- write_pattern(type = type , dist_e = dist_e, dist_c = dist_c, model_txt_info = model_txt_info)
  model <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, 
                        model.file = filein, n.chains = n.chains, n.iter = n.iter, 
                        n.burnin = n.burnin, DIC = DIC, pD = pd, n.thin = n.thin,
                        n.iter.pd = n.iter.pd, n.adapt = n.adapt)
  tmu_e <- model$BUGSoutput$sims.list$tmu_e
  tmu_c <- model$BUGSoutput$sims.list$tmu_c
  s_c <- model$BUGSoutput$sims.list$s_c_p
  s_e <- model$BUGSoutput$sims.list$s_e_p
  alpha <- model$BUGSoutput$sims.list$alpha_p
  beta <- model$BUGSoutput$sims.list$beta_p
  p_prob <- model$BUGSoutput$sims.list$p_prob
  colnames(s_c) <- paste0("s_c_p", 1:model_txt_info$n_patterns)
  colnames(s_e) <- paste0("s_e_p", 1:model_txt_info$n_patterns)    
  dimnames(alpha) <- list(NULL, colnames(as.data.frame(X_e_fixed)), 1:model_txt_info$n_patterns)
  dimnames(beta) <- list(NULL, colnames(as.data.frame(X_c_fixed)), 1:model_txt_info$n_patterns)
  colnames(p_prob) <- paste0("p", 1:model_txt_info$n_patterns)
  if(dist_e %in% c("norm", "negbin", "logis")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    tau_e <- model$BUGSoutput$sims.list$tau_e_p
  } else if(dist_e %in% c("beta", "gamma", "weib")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
    ctau_e <- model$BUGSoutput$sims.list$ctau_e
  } else if(dist_e %in% c("exp", "bern", "pois")) {
    cmu_e <- model$BUGSoutput$sims.list$cmu_e
  } 
  if(dist_c == "norm") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    tau_c <- model$BUGSoutput$sims.list$tau_c_p
  } else if(dist_c == "gamma") {
    cmu_c <- model$BUGSoutput$sims.list$cmu_c
    ctau_c <- model$BUGSoutput$sims.list$ctau_c
  } else if(dist_c == "lnorm") {
    clmu_c <- model$BUGSoutput$sims.list$clmu_c      
    ltau_c <- model$BUGSoutput$sims.list$ltau_c_p
  } 
  S <- model_info$n.mci
  if(S != n.iter) { 
  S <- n.iter
  warning("Value for 'n.mci' set equal to the number of MCMC iterations to compute average estimates across patterns.")}
  trt_index <- c(data_model$trt_index)
  p_index <- lapply(1:data_model$n_patterns, function(target) which(d_mod == target))
  mu_ep <- mu_cp <- array(NA, dim = c(S, length(data_model$n_trt), data_model$n_patterns), 
                          dimnames = list(NULL, data_model$trt_lev, 1:data_model$n_patterns))
  for(i in 1:n.iter) {
    for(trt in 1:length(data_model$n_trt)) {
      for(p in 1:data_model$n_patterns){
        if(dist_e == "norm") {
          mu_ep[i, trt, p] <- mean(rnorm(S, mean = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], sd = s_e[i, p]), na.rm = TRUE)
        }
        if(dist_e == "negbin") {
          mu_ep[i, trt, p] <- mean(rnbinom(S, prob = tau_e[i]/(tau_e[i] + cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), size = tau_e[i, p]), na.rm = TRUE)
        }  
        if(dist_e == "logis") {
          mu_ep[i, trt, p] <- mean(rlogis(S, location = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], scale = 1 / tau_e[i, p]), na.rm = TRUE)
        }
        if(dist_e == "beta") {
          mu_ep[i, trt, p] <- mean(rbeta(S, shape1 = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))] * ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], 
                                     shape2 = (1 - cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]) * ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_e == "gamma") {
          mu_ep[i, trt, p] <- mean(rgamma(S, shape = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))] * ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], 
                                      rate = ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_e == "weib") {
          mu_ep[i, trt, p] <- mean(rweibull(S, shape = ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), scale = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))] / exp(lgamma(1 + 1 / ctau_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))])), na.rm = TRUE)
        }
        if(dist_e == "exp") {
          mu_ep[i, trt, p] <- mean(rexp(S, rate = 1 / cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_e == "bern") {
          mu_ep[i, trt, p] <- mean(rbinom(S, size = 1, prob = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_e == "pois") {
          mu_ep[i, trt, p] <- mean(rpois(S, lambda = cmu_e[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_c == "norm") {
          mu_cp[i, trt, p] <- mean(rnorm(S, mean = cmu_c[i, unlist(data_model$trt_index[trt])], sd = s_c[i, p]), na.rm = TRUE)
        }
        if(dist_c == "gamma") {
          mu_cp[i, trt, p] <- mean(rgamma(S, shape = cmu_c[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))] * ctau_c[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], 
                                      rate = ctau_c[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))]), na.rm = TRUE)
        }
        if(dist_c == "lnorm") {
          mu_cp[i, trt, p] <- mean(rlnorm(S, meanlog = clmu_c[i, intersect(unlist(data_model$trt_index[trt]), unlist(p_index[p]))], sdlog = 1 / sqrt(ltau_c[i, p])), na.rm = TRUE)
        }
      }
    }
  } 
  mu_e <- mu_c <- matrix(NA, nrow = S, ncol = length(data_model$n_trt), dimnames = list(NULL, data_model$trt_lev))
  delta_ep <- matrix(0, nrow = S, ncol = dim(mu_ep)[3], dimnames = list(NULL, 1:dim(mu_ep)[3]))
  delta_cp <- matrix(0, nrow = S, ncol = dim(mu_cp)[3], dimnames = list(NULL, 1:dim(mu_cp)[3]))
  if(type == "MNAR") {
    delta_e <- model$BUGSoutput$sims.list$delta_e
    delta_c <- model$BUGSoutput$sims.list$delta_c
    delta_ep[, 1:dim(mu_ep)[3]] <- delta_e
    delta_cp[, 1:dim(mu_cp)[3]] <- delta_c
    if(all(unique(d_or) %in% c(1, 2)) & n_patterns == 2) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 3)) & n_patterns == 2) { delta_ep[, 1:2] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 4)) & n_patterns == 2) { delta_ep[, 1] <- delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3)) & n_patterns == 2) { delta_ep[, 2] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(3, 4)) & n_patterns == 2) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 4)) & n_patterns == 2) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3)) & n_patterns == 3) { delta_ep[, c(1, 3)] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 2, 4)) & n_patterns == 3) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 3, 4)) & n_patterns == 3) { delta_ep[, 1:2] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3, 4)) & n_patterns == 3) { delta_ep[, 2] <- rep(0, dim(delta_ep)[1]); delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3, 4)) & n_patterns == 4) { delta_ep[, c(1, 3)] <- rep(0, dim(delta_ep)[1]); delta_cp[, c(1, 2)] <- rep(0, dim(delta_ep)[1])}
    colnames(delta_e) <- "delta_e"
    colnames(delta_c) <- "delta_c"
  } else if(type == "MNAR_eff") {
    delta_e <- model$BUGSoutput$sims.list$delta_e
    delta_ep[, 1:dim(mu_ep)[3]] <- delta_e
    if(all(unique(d_or) %in% c(1, 2)) & n_patterns == 2) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 3)) & n_patterns == 2) { delta_ep[, 1:2] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 4)) & n_patterns == 2) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3)) & n_patterns == 2) { delta_ep[, 2] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(3, 4)) & n_patterns == 2) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3)) & n_patterns == 3) { delta_ep[, c(1, 3)] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 4)) & n_patterns == 3) { delta_ep[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 3, 4)) & n_patterns == 3) { delta_ep[, 1:2] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3, 4)) & n_patterns == 3) { delta_ep[, 2] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3, 4)) & n_patterns == 4) { delta_ep[, c(1, 3)] <- rep(0, dim(delta_ep)[1])}
    colnames(delta_e) <- "delta_e"
  } else if(type == "MNAR_cost") {
    delta_c <- model$BUGSoutput$sims.list$delta_c
    delta_cp[, 1:dim(mu_cp)[3]] <- delta_c
    if(all(unique(d_or) %in% c(1, 2)) & n_patterns == 2) { delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 3)) & n_patterns == 2) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 4)) & n_patterns == 2) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3)) & n_patterns == 2) { delta_cp[, 1] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(2, 4)) & n_patterns == 2) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3)) & n_patterns == 3) { delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 2, 4)) & n_patterns == 3) { delta_cp[, 1:2] <- rep(0, dim(delta_cp)[1])}
    if(all(unique(d_or) %in% c(1, 3, 4)) & n_patterns == 3) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(2, 3, 4)) & n_patterns == 3) { delta_cp[, 1] <- rep(0, dim(delta_ep)[1])}
    if(all(unique(d_or) %in% c(1, 2, 3, 4)) & n_patterns == 4) { delta_cp[, c(1, 2)] <- rep(0, dim(delta_ep)[1])}
    colnames(delta_c) <- "delta_c"
  } 
  for(trt in 1:length(data_model$n_trt)) {
    mu_e[, trt] <-  rowSums(mu_ep[, trt, ] * p_prob[, ] + delta_ep)
    mu_c[, trt] <-  rowSums(mu_cp[, trt, ] * p_prob[, ] + delta_cp)
  }
  a <- b <- NULL
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
  loglik_d <- model$BUGSoutput$sims.list$loglik_d
  if(!ind_fixed) {
    beta_f <- model$BUGSoutput$sims.list$beta_f_p
    colnames(beta_f) <- paste0("e", 1:model_txt_info$n_patterns)
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
    param_jags_show <- c("eff", "cost", "loglik_e", "loglik_c", "loglik_d", "cmu_c", "cmu_e")
    if(dist_e %in% c("beta", "gamma", "weib")) { param_jags_show <- c(param_jags_show, "ctau_e")}
    if(dist_c %in% c("gamma")) { param_jags_show <- c(param_jags_show, "ctau_c")}
    if(dist_c %in% c("lnorm")) { param_jags_show <- c(param_jags_show, "clmu_c")}
    model_sum <- round(jagsresults(x = model, params = param_jags_show, invert = TRUE), digits = 3)
  } else{ model_sum <- NULL}  
  loglik <- list("effects" = loglik_e, "costs" = loglik_c, "pattern indicators" = loglik_d)
  imputed_e <- list("avg" = eff_pos_imp_avg, "ql" = eff_pos_imp_ql, "qu" = eff_pos_imp_qu, "index" = which(is.na(data_model$eff)))
  imputed_c <- list("avg" = cost_pos_imp_avg, "ql" = cost_pos_imp_ql, "qu" = cost_pos_imp_qu, "index" = which(is.na(data_model$cost)))
  imputed_et <- list("avg" = efft_pos_imp_avg, "ql" = efft_pos_imp_ql, "qu" = efft_pos_imp_avg, "index" = lapply(data_model$m_efft, function(x) which(x == 1)))
  imputed_ct <- list("avg" = costt_pos_imp_avg, "ql" = costt_pos_imp_ql, "qu" = costt_pos_imp_avg, "index" = lapply(data_model$m_costt, function(x) which(x == 1)))
  imputed <- list("effects" = imputed_e, "costs" = imputed_c)
  imputedt <- list("effects" = imputed_et, "costs" = imputed_ct)  
  if(type == "MAR") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "p_pattern" = p_prob, "effects_random" = a, "costs_random" = b, 
                              "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, "type" = "MAR", "method" = "PATTERN", 
                              "ind_fixed" = ind_fixed, "ind_random" = ind_random, "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein, "d_mod" = d_mod)
  } else if(type == "MNAR") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "p_pattern" = p_prob, "effects_sens_par" = delta_e, "costs_sens_par" = delta_c,
                              "effects_random" = a, "costs_random" = b, "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR", "method" = "PATTERN", "ind_fixed" = ind_fixed, "ind_random" = ind_random, 
                              "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein, "d_mod" = d_mod)
  } else if(type == "MNAR_eff") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "p_pattern" = p_prob, "effects_sens_par" = delta_e, "effects_random" = a, "costs_random" = b,
                              "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR_eff", "method" = "PATTERN", "ind_fixed" = ind_fixed, "ind_random" = ind_random, 
                              "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein, "d_mod" = d_mod)
  } else if(type == "MNAR_cost") {
    model_output_jags <- list("summary" = model_sum, "model" = model, "mean_effects" = mu_e, "mean_costs" = mu_c, 
                              "sd_effects" = s_e, "sd_costs" = s_c, "effects_fixed" = alpha, "costs_fixed" = beta, 
                              "p_pattern" = p_prob, "costs_sens_par" = delta_c, "effects_random" = a, "costs_random" = b, 
                              "imputed" = imputed, "imputed_trt" = imputedt, "loglik" = loglik, 
                              "type" = "MNAR_cost", "method" = "PATTERN", "ind_fixed" = ind_fixed, "ind_random" = ind_random, 
                              "dist_e" = dist_e, "dist_c" = dist_c, "filein" = filein, "d_mod" = d_mod)
  }  
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}