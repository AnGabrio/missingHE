#' An internal function to execute a JAGS selection model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters
#' @keywords JAGS Bayesian hurdle models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_selection <- function(type, dist_e, dist_c, inits) eval.parent(substitute( {
  if(!isTRUE(requireNamespace("R2jags", quietly = TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  if(!dist_e %in% c("norm","beta") | !dist_c %in% c("norm", "gamma")) {
    stop("Distributions available for use are 'norm','beta' for the effects and 'norm','gamma' for the costs")
  }
  if(!type %in% c("MAR", "MNAR", "MNAR_eff", "MNAR_cost")) {
    stop("Types available for use are 'MAR','MNAR_eff','MNAR_cost'and 'MNAR'")
  }
  if(is.null(inits) == FALSE) {inits = inits }
  model <- write_selection(type = type , dist_e = dist_e, dist_c = dist_c, pe = pe, pc = pc, ze = ze, zc = zc, ind = ind)
  filein <- model
  datalist <- list("N1", "N2", "eff1", "eff2", "cost1", "cost2", "m_eff1", "m_eff2", "m_cost1", "m_cost2", 
                   "X1_e", "X2_e", "X1_c", "X2_c", "Z1_e", "Z2_e", "Z1_c", "Z2_c", "mean_cov_e1", "mean_cov_e2", "mean_cov_c1", 
                   "mean_cov_c2", "mean_z_e1", "mean_z_e2", "mean_z_c1", "mean_z_c2", "pe", "pc", "ze", "zc")
  if(pe == 1) {pe_index <- match("pe", datalist)
  datalist <- datalist[-pe_index] }
  if(pc == 1) {pc_index <- match("pc", datalist)
  datalist <- datalist[-pc_index] }
  if(ze == 1) {ze_index <- match("ze", datalist)
  datalist <- datalist[-ze_index] }
  if(zc == 1) {zc_index <- match("zc", datalist)
  datalist <- datalist[-zc_index] }
  DIC <- TRUE
  params <- c("eff1", "eff2", "cost1", "cost2", "mu_e", "mu_c", "s_e", "s_c", "p_e", "p_c", "beta", "alpha", "gamma_e", "gamma_c", "delta_e", "delta_c")
  if(ind == FALSE) {params <- c(params, "beta_f") }
  if(type == "MNAR_cost" | type == "MAR") {deltae_index <- match("delta_e", params)
  params <- params[-deltae_index] }
  if(type == "MNAR_eff" | type == "MAR") {deltac_index <- match("delta_c", params)
  params <- params[-deltac_index] }
  modelN1 <- R2jags::jags(data = datalist, inits = inits, parameters.to.save = params, model.file = filein, n.chains = n.chains, 
                          n.iter = n.iter, n.burnin = n.burnin, DIC = DIC, n.thin = n.thin)
    mu_e <- modelN1$BUGSoutput$sims.list$mu_e
    mu_c <- modelN1$BUGSoutput$sims.list$mu_c
    s_e <- modelN1$BUGSoutput$sims.list$s_e
    s_c <- modelN1$BUGSoutput$sims.list$s_c
    alpha <- modelN1$BUGSoutput$sims.list$alpha
    beta <- modelN1$BUGSoutput$sims.list$beta
    p_e <- modelN1$BUGSoutput$sims.list$p_e
    p_c <- modelN1$BUGSoutput$sims.list$p_c
    gamma_e <- modelN1$BUGSoutput$sims.list$gamma_e
    gamma_c <- modelN1$BUGSoutput$sims.list$gamma_c
    eff1_pos <- matrix(eff1, N1 ,3)
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
  if(ind == FALSE & dist_e == "norm" & dist_c == "norm") {
    beta_f <- modelN1$BUGSoutput$sims.list$beta_f
    beta <- cbind(beta, beta_f)
  }
  if(type == "MNAR") {
    delta_e <- modelN1$BUGSoutput$sims.list$delta_e
    delta_c <- modelN1$BUGSoutput$sims.list$delta_c
  } else if(type == "MNAR_eff") {
    delta_e <- modelN1$BUGSoutput$sims.list$delta_e
  } else if(type == "MNAR_cost") {
    delta_c <- modelN1$BUGSoutput$sims.list$delta_c
  }
   if(n.chains > 1) {model_sum <- round(jagsresults(x = modelN1, params = c('eff1', 'eff2', 'cost1', 'cost2'), invert = TRUE), digits = 3)
   } else{model_sum <- NULL }
    colnames(eff1_pos) <- c("mean", "LB", "UB")
    colnames(eff2_pos) <- c("mean", "LB", "UB")
    colnames(cost1_pos) <- c("mean", "LB", "UB")
    colnames(cost2_pos) <- c("mean", "LB", "UB")
    imputed <- list("effects1" = eff1_pos, "effects2" = eff2_pos, "costs1" = cost1_pos, "costs2" = cost2_pos)
    if(type == "MAR") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects" = alpha, "covariate_parameter_costs" = beta, "missingness_probability_effects" = p_e, 
                                "missingness_parameter_effects" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs" = gamma_c, 
                                "imputed" = imputed, "type" = "SELECTION", "ind" = ind)
    } else if(type == "MNAR") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects" = alpha, "covariate_parameter_costs" = beta, "missingness_probability_effects" = p_e, 
                                "missingness_parameter_effects" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs" = gamma_c, 
                                "mnar_parameter_effects" = delta_e, "mnar_parameter_costs" = delta_c, "imputed" = imputed, "type" = "SELECTION_ec", "ind" = ind)
    } else if(type == "MNAR_eff") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects" = alpha, "covariate_parameter_costs" = beta, "missingness_probability_effects" = p_e, 
                                "missingness_parameter_effects" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs" = gamma_c, 
                                "mnar_parameter_effects" = delta_e, "imputed" = imputed, "type" = "SELECTION_e", "ind" = ind)
    } else if(type == "MNAR_cost") {
      model_output_jags <- list("summary" = model_sum, "model summary" = modelN1, "mean_effects" = mu_e, "mean_costs" = mu_c, "sd_effects" = s_e, "sd_costs" = s_c, 
                                "covariate_parameter_effects" = alpha, "covariate_parameter_costs" = beta, "missingness_probability_effects" = p_e, 
                                "missingness_parameter_effects" = gamma_e, "missingness_probability_costs" = p_c, "missingness_parameter_costs" = gamma_c, 
                                "mnar_parameter_costs" = delta_c, "imputed" = imputed, "type" = "SELECTION_c", "ind" = ind)
    }
  if(n.chains == 1) {model_output_jags <- model_output_jags[-1] }
  return(model_output_jags = model_output_jags)
}))