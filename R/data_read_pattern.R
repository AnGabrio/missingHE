#' A function to read and re-arrange the data in different ways
#'
#' This internal function imports the data and outputs only those variables that are needed to run the model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @keywords read data
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response
#' @export
#' @examples
#' # Internal function only
#' # no examples
#' #
#' #


data_read_pattern <- function(data, model.eff, model.cost, type = type, center) {
  if(is.data.frame(data) == FALSE) {
    stop("object data must be provided as data frame")
  }
  if(any(names(data) == "e") == TRUE & any(names(data) == "c") == TRUE) {
    e <- as.name("e")
    c <- as.name("c")
  }
  cov_matrix <- subset(data, select = -c(e, c))
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  is.formula<-function (x) { inherits(x, "formula") }
  if(is.formula(model.eff) == FALSE | is.formula(model.cost) == FALSE) {
    stop("model.eff and/or model.cost must be formula objects")
  }
  if(all(names(model.frame(model.eff, data = data)) %in% c("e", names(cov_matrix))) == FALSE | 
     all(names(model.frame(model.cost, data = data)) %in% c("c", "e", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the model")
  }
  if(all(names(model.frame(model.eff, data = data)) %in% names(data)) == FALSE | 
     all(names(model.frame(model.cost, data = data)) %in% names(data)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if("e" %in% labels(terms(model.eff)) | "c" %in% labels(terms(model.cost))) {
    stop("please remove 'e' from the right hand side of model.eff and/or 'c' from the right hand side of model.cost")
  }
  if(names(model.frame(model.eff, data = data)[1]) != "e") {
    stop("you must set 'e' as the response in the formula model.eff")
  }
  if("c" %in% names(model.frame(model.eff, data = data))) {
    stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
  }
  if(names(model.frame(model.cost, data = data)[1]) != "c") {
    stop("you must set 'c' as the response in the formula model.cost")
  }
  if("t" %in% names(model.frame(model.cost, data = data)) | "t" %in% names(model.frame(model.eff, data = data))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  if(is.logical(center) == FALSE) { stop("center must be either TRUE or FALSE") }
  index_mis_e <- which(is.na(data$e))
  index_mis_c <- which(is.na(data$c))
  data2 <- data
  data$e[is.na(data$e) == TRUE] <- -999999
  data$c[is.na(data$c) == TRUE] <- -999999
  mf_e <- model.frame(formula = model.eff, data = data)
  mf_c <- model.frame(formula = model.cost, data = data)
  terms <- NULL
  x_e <- model.matrix(attr(mf_e, "terms"), data = mf_e)
  x_c <- model.matrix(attr(mf_c, "terms"), data = mf_c)
  if("e" %in% names(mf_c)){
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
  d1 <- d2 <- c()
  d1[is.na(eff1) == FALSE & is.na(cost1) == FALSE] <- 1
  d1[is.na(eff1) == TRUE & is.na(cost1) == FALSE] <- 2
  d1[is.na(eff1) == FALSE & is.na(cost1) == TRUE] <- 3
  d1[is.na(eff1) == TRUE & is.na(cost1) == TRUE] <- 4
  d2[is.na(eff2) == FALSE & is.na(cost2) == FALSE] <- 1
  d2[is.na(eff2) == TRUE & is.na(cost2) == FALSE] <- 2
  d2[is.na(eff2) == FALSE & is.na(cost2) == TRUE] <- 3
  d2[is.na(eff2) == TRUE & is.na(cost2) == TRUE] <- 4
  n_patterns1 <- length(unique(d1))
  n_patterns2 <- length(unique(d2))
  d1_list <- list()
  d1_list[[1]] <- any(d1 == 1)
  d1_list[[2]] <- any(d1 == 2)
  d1_list[[3]] <- any(d1 == 3)
  d1_list[[4]] <- any(d1 == 4)
  names(d1_list) <- c("d1_ec_obs", "d1_c_obs", "d1_e_obs", "d1_ec_mis")
  d2_list<-list()
  d2_list[[1]] <- any(d2 == 1)
  d2_list[[2]] <- any(d2 == 2)
  d2_list[[3]] <- any(d2 == 3)
  d2_list[[4]] <- any(d2 == 4)
  names(d2_list) <- c("d2_ec_obs", "d2_c_obs", "d2_e_obs", "d2_ec_mis")
  d_list <- list()
  d_list[[1]] <- c(n_patterns1, n_patterns2)
  d_list[[2]] <- d1_list
  d_list[[3]] <- d2_list
  names(d_list) <- c("n_patterns", "d1", "d2")
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
  data2$me <- c(m_eff1, m_eff2)
  data2$mc <- c(m_cost1, m_cost2)
  d <- list(d1, d2)
  names(cov_e) <- names(cov_c) <- names(mean_cov_e) <- names(mean_cov_c) <- c("Control", "Intervention")
  names(m_eff) <- names(m_cost) <- c("Control", "Intervention")
  names(d) <- c("Control", "Intervention")
  names(effects) <- names(costs) <- names(eff_cc) <- names(cost_cc) <- c("Control", "Intervention")
  data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                   "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                   "covariates_effects" = cov_e, "covariates_costs" = cov_c, "mean_cov_effects" = mean_cov_e, "mean_cov_costs" = mean_cov_c, 
                   "data_ind" = data2, "patterns_list" = d_list, "patterns" = d) 
  return(data_raw) 
}


