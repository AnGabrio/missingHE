#' A function to read and re-arrange the data in different ways for the hurdle model
#'
#' This internal function imports the data and outputs only those variables that are needed to run the hurdle model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.se}, \code{model.sc} (model formulas for the structural effect and cost models) . Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model.
#' @param model.se A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural effects are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the strcutural effects through a logistic-linear model.
#' @param model.sc A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural costs are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the strcutural costs through a logistic-linear model.
#' @param se Structural value to be found in the effect data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param type Type of structural value mechanism assumed, either 'SCAR' (Structural Completely At Random) or 'SAR' (Strcutural At Random).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @keywords read data hurdle models
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response
#' @export
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

data_read_hurdle <- function(data, model.eff, model.cost, model.se, model.sc, se, sc, type = type, center) {
  if(any(names(data) == "e") == TRUE & any(names(data) == "c") == TRUE) {
    e <- as.name("e")
    c <- as.name("c")
  }
  cov_matrix <- subset(data, select = -c(e, c) )
  if(any(is.na(cov_matrix)) == TRUE) {
    stop("no missing covariate or treatment indicator is allowed")
  }
  is.formula<-function (x) { inherits(x, "formula") }
  if(is.formula(model.eff) == FALSE | is.formula(model.cost) == FALSE) {
    stop("model.eff and/or model.cost must be formula objects")
  }
  if(all(names(model.frame(model.eff)) %in% names(data)) == FALSE | 
     all(names(model.frame(model.cost)) %in% names(data)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if("e" %in% labels(terms(model.eff)) | "c" %in% labels(terms(model.cost))) {
    stop("please remove 'e' from the right hand side of model.eff and/or 'c' from the right hand side of model.cost")
  }
  if(names(model.frame(model.eff)[1]) != "e") {
    stop("you must set 'e' as the response in the formula model.eff")
  }
  if("c" %in% names(model.frame(model.eff))) {
    stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
  }
  if(names(model.frame(model.cost)[1]) != "c") {
    stop("you must set 'c' as the response in the formula model.cost")
  }
  if("t" %in% names(model.frame(model.cost)) | "t" %in% names(model.frame(model.eff))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  index_mis_e <- index_mis2_e <- which(is.na(data$e))
  index_mis_c <- index_mis2_c <- which(is.na(data$c))
  data2 <- data
  data$e[is.na(data$e) == TRUE] <- -999999
  data$c[is.na(data$c) == TRUE] <- -999999
  mf_e <- model.frame(formula = model.eff, data = data)
  mf_c <- model.frame(formula = model.cost, data = data)
  terms <- NULL
  x_e <- model.matrix(attr(mf_e, "terms"), data = mf_e)
  x_c <- model.matrix(attr(mf_c, "terms"), data = mf_c)
  if("e" %in% names(mf_c)) {
    mf_c$e[index_mis_e] <- NA
  }
  y_e <- model.response(mf_e)
  y_c <- model.response(mf_c)
  y_e[index_mis_e] <- NA
  y_c[index_mis_c] <- NA
  data$e[index_mis_e] <- NA
  data$c[index_mis_c] <- NA
  N1 <- N2 <- c() 
  N1 <- sum(data$t == 1) 
  N2 <- length(data$t) - N1 
  N <- c(N1, N2)
  m_eff <- rep(0, length(data$e))
  m_eff[index_mis_e] <- 1
  m_cost <- rep(0, length(data$c))
  m_cost[index_mis_c] <- 1
  m_eff1 <- m_eff2 <- m_cost1 <- m_cost2 <- c() 
  t1_index <- which(data$t == 1)
  t2_index <- which(data$t == 2)
  eff1 <- y_e[t1_index] 
  eff2 <- y_e[t2_index] 
  eff <- list(eff1, eff2)
  cost1 <- y_c[t1_index] 
  cost2 <- y_c[t2_index] 
  cost <- list(cost1, cost2)
  m_eff1 <- m_eff[t1_index]
  m_eff2 <- m_eff[t2_index]
  m_eff <- list(m_eff1, m_eff2) 
  m_cost1 <- m_cost[t1_index]
  m_cost2 <- m_cost[t2_index]
  m_cost <- list(m_cost1, m_cost2)  
  N1_cc <- N2_cc <- N1_mis <- N2_mis <- c() 
  N1_cc[1] <- length(na.omit(eff1)) 
  N1_cc[2] <- length(na.omit(cost1)) 
  N2_cc[1] <- length(na.omit(eff2)) 
  N2_cc[2] <- length(na.omit(cost2)) 
  N_cc <- cbind(N1_cc, N2_cc) 
  N1_mis <- N1 - N1_cc 
  N2_mis <- N2 - N2_cc 
  N_mis <- cbind(N1_mis, N2_mis) 
  effects <- list(eff1, eff2) 
  costs <- list(cost1, cost2) 
  eff1_cc <- eff2_cc <- cost1_cc <- cost2_cc <- c() 
  eff1_cc <- na.omit(eff1) 
  eff2_cc <- na.omit(eff2) 
  eff_cc <- list(eff1_cc, eff2_cc) 
  cost1_cc <- na.omit(cost1) 
  cost2_cc <- na.omit(cost2) 
  cost_cc <- list(cost1_cc, cost2_cc)
  cov1_e <- as.data.frame(x_e[t1_index, ])
  names(cov1_e) <- colnames(x_e)
  cov2_e <- as.data.frame(x_e[t2_index, ])
  names(cov2_e) <- colnames(x_e)
  cov_e <- list(cov1_e, cov2_e)
  x_c_hold <- x_c
  if("e" %in% colnames(x_c_hold)) {
    x_c <- subset(x_c_hold, select = -c(e))
  }
  cov1_c <- as.data.frame(x_c[t1_index, ])
  names(cov1_c) <- colnames(x_c)
  cov2_c <- as.data.frame(x_c[t2_index, ])
  names(cov2_c) <- colnames(x_c)
  cov_c <- list(cov1_c, cov2_c)
  cove <- list(cov1_e, cov2_e) 
  mean_cov_e <- list(apply(as.matrix(cov1_e), 2, mean), apply(as.matrix(cov2_e), 2, mean))
  covc <- list(cov1_c, cov2_c) 
  mean_cov_c <- list(apply(as.matrix(cov1_c), 2, mean), apply(as.matrix(cov2_c), 2, mean))
  cov1_e_center <- as.data.frame(scale(cov1_e, scale = FALSE))
  cov2_e_center <- as.data.frame(scale(cov2_e, scale = FALSE))
  cov1_e_center[, 1] <- rep(1, nrow(cov1_e))
  cov2_e_center[, 1] <- rep(1, nrow(cov2_e))
  cov_e_center <- list(cov1_e_center, cov2_e_center)
  mean_cov_e_center <- list(apply(as.matrix(cov1_e_center), 2, mean), apply(as.matrix(cov2_e_center), 2, mean))
  cov1_c_center <- as.data.frame(scale(cov1_c, scale = FALSE))
  cov2_c_center <- as.data.frame(scale(cov2_c, scale = FALSE))
  cov1_c_center[, 1] <- rep(1, nrow(cov1_c))
  cov2_c_center[, 1] <- rep(1, nrow(cov2_c))
  cov_c_center <- list(cov1_c_center, cov2_c_center)
  mean_cov_c_center <- list(apply(as.matrix(cov1_c_center), 2, mean), apply(as.matrix(cov2_c_center), 2, mean))
  if(center == TRUE) {
    cov_e <- cov_e_center
    cov_c <- cov_c_center
    mean_cov_e <- mean_cov_e_center
    mean_cov_c <- mean_cov_c_center
  }
  data2$e[is.na(data2$e) == TRUE] <- -999999
  data2$c[is.na(data2$c) == TRUE] <- -999999
  if(type == "SCAR" | type == "SAR") {
    data2$se <- c(m_eff1, m_eff2)
    data2$sc <- c(m_cost1, m_cost2)
  }
  if(is.formula(model.se) == FALSE | is.formula(model.sc) == FALSE) {
    stop("model.se and/or model.sc must be formula objects")
  }
  if(all(names(model.frame(model.se, data = data2)) %in% names(data2)) == FALSE | 
     all(names(model.frame(model.sc, data = data2)) %in% names(data2)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if(names(model.frame(model.se, data = data2)[1]) != "se") {
    stop("you must set 'se' as the response in the formula model.se")
  }
  if(names(model.frame(model.sc, data = data2)[1]) != "sc") {
    stop("you must set 'sc' as the response in the formula model.sc")
  }
  if("t" %in% names(model.frame(model.sc, data = data2)) | "t" %in% names(model.frame(model.se, data = data2))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.se' and/or 'model.sc'")
  }
  if(any(c("c", "e") %in% names(model.frame(model.se, data = data2))) == TRUE | 
     any(c("e", "c") %in% names(model.frame(model.sc, data = data2))) == TRUE) {
    stop("please remove 'e' and/or 'c' from model.se and/or model.sc")
  }
  zf_e <- model.frame(formula = model.se, data = data2)
  zf_c <- model.frame(formula = model.sc, data = data2)
  z_e <- model.matrix(attr(zf_e, "terms"), data = zf_e)
  z_c <- model.matrix(attr(zf_c, "terms"), data = zf_c)
  y_e <- model.response(mf_e)
  y_c <- model.response(mf_c)
  y_e[index_mis2_e] <- NA
  y_c[index_mis2_c] <- NA
  data2$e[index_mis2_e] <- NA
  data2$c[index_mis2_c] <- NA
  if(is.null(se) == TRUE & is.null(sc) == TRUE) { stop("Structural values in at least one outcome variable are required, please provide the structural value") }
  if(is.null(se) == FALSE) {
  index_str_e <- ifelse(round(y_e, digits = 4) == se, 1, 0)
  if(any(na.omit(index_str_e) == 1) == FALSE) { stop("Provide structural values that are present in the data") }
  d_eff1 <- d_eff2 <- c()
  d_eff1 <- index_str_e[t1_index]
  d_eff2 <- index_str_e[t2_index]
  d_eff <- list(d_eff1, d_eff2) 
  covz1_e <- as.data.frame(z_e[t1_index, ])
  names(covz1_e) <- colnames(z_e)
  covz2_e <- as.data.frame(z_e[t2_index, ])
  names(covz2_e) <- colnames(z_e)
  covz_e <- list(covz1_e, covz2_e)
  covze <- list(covz1_e, covz2_e) 
  mean_covz_e <- list(apply(as.matrix(covz1_e), 2, mean), apply(as.matrix(covz2_e), 2, mean))
  covz1_e_center <- as.data.frame(scale(covz1_e, scale = FALSE))
  covz2_e_center <- as.data.frame(scale(covz2_e, scale = FALSE))
  covz1_e_center[, 1] <- rep(1, nrow(covz1_e))
  covz2_e_center[, 1] <- rep(1, nrow(covz2_e))
  covz_e_center <- list(covz1_e_center, covz2_e_center)
  mean_covz_e_center <- list(apply(as.matrix(covz1_e_center), 2, mean), apply(as.matrix(covz2_e_center), 2, mean))
  if(center == TRUE) {
    covz_e <- covz_e_center
    mean_covz_e <- mean_covz_e_center
  }
  names(covz_e) <- names(mean_covz_e) <- c("Control", "Intervention")
  names(d_eff) <- c("Control", "Intervention")
  }
  if(is.null(sc) == FALSE) {
  index_str_c <- ifelse(round(y_c, digits = 0) == sc, 1, 0)
  if(any(na.omit(index_str_c) == 1) == FALSE){ stop("Provide structural values that are present in the data") }
  d_cost1 <- d_cost2 <- c() 
  d_cost1 <- index_str_c[t1_index]
  d_cost2 <- index_str_c[t2_index]
  d_cost <- list(d_cost1, d_cost2)  
  covz1_c <- as.data.frame(z_c[t1_index, ])
  names(covz1_c) <- colnames(z_c)
  covz2_c <- as.data.frame(z_c[t2_index, ])
  names(covz2_c) <- colnames(z_c)
  covz_c <- list(covz1_c, covz2_c)
  covzc <- list(covz1_c, covz2_c) 
  mean_covz_c <- list(apply(as.matrix(covz1_c), 2, mean), apply(as.matrix(covz2_c), 2, mean))
  mean_covz_c <- list(apply(as.matrix(covz1_c), 2, mean), apply(as.matrix(covz2_c), 2, mean))
  covz1_c_center <- as.data.frame(scale(covz1_c, scale = FALSE))
  covz2_c_center <- as.data.frame(scale(covz2_c, scale = FALSE))
  covz1_c_center[, 1] <- rep(1, nrow(covz1_c))
  covz2_c_center[, 1] <- rep(1, nrow(covz2_c))
  covz_c_center <- list(covz1_c_center, covz2_c_center)
  mean_covz_c_center <- list(apply(as.matrix(covz1_c_center), 2, mean), apply(as.matrix(covz2_c_center), 2, mean))
  if(center == TRUE) {
    covz_c <- covz_c_center
    mean_covz_c <- mean_covz_c_center
  }
  names(covz_c) <- names(mean_covz_c) <- c("Control", "Intervention")
  names(d_cost) <- c("Control", "Intervention")
  }
  names(cov_e) <- names(cov_c) <- names(mean_cov_e) <- names(mean_cov_c) <- c("Control", "Intervention")
  names(m_eff) <- names(m_cost) <- c("Control", "Intervention")
  names(effects) <- names(costs) <- names(eff_cc) <- names(cost_cc) <- c("Control", "Intervention")
  if(is.null(se) == FALSE & is.null(sc) == TRUE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects" = cov_e, "covariates_costs" = cov_c, "mean_cov_effects" = mean_cov_e, "mean_cov_costs" = mean_cov_c, 
                     "covariates_structural_effects" = covz_e, "mean_cov_structural_effects" = mean_covz_e, "structural_effects" = d_eff, "data_ind" = data2) 
  } else if(is.null(se) == TRUE & is.null(sc) == FALSE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N,
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects" = cov_e, "covariates_costs" = cov_c, "mean_cov_effects" = mean_cov_e, "mean_cov_costs" = mean_cov_c, 
                     "covariates_structural_costs" = covz_c, "mean_cov_structural_costs" = mean_covz_c, "structural_costs" = d_cost, "data_ind" = data2) 
  }else if(is.null(se) == FALSE & is.null(sc) == FALSE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects" = cov_e,"covariates_costs" = cov_c, "mean_cov_effects" = mean_cov_e, "mean_cov_costs" = mean_cov_c, 
                     "covariates_structural_effects" = covz_e, "mean_cov_structural_effects" = mean_covz_e, "covariates_structural_costs" = covz_c, 
                     "mean_cov_structural_costs" = mean_covz_c, "structural_effects" = d_eff, "structural_costs" = d_cost, "data_ind" = data2) 
  }
  return(data_raw) 
}


