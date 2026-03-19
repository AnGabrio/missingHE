#' Extract regression coefficient estimates from objects in the class \code{missingHE}
#'
#' Produces a table printout with summary statistics for the regression coefficients of the health economic evaluation probabilistic model
#' run using the function \code{\link{selection}}, \code{\link{pattern}},  \code{\link{hurdle}} or \code{\link{lmdm}}.
#' @param object A \code{missingHE} object containing the results of the Bayesian modelling and the economic evaluation
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#' CI sample quantiles to be calculated and returned for the estimates.
#' @param random Logical. If \code{random} is \code{TRUE}, the estimates of the random effects parameters are printed, when available.
#' @param digits Integer indicating the number of decimal places to be used for rounding (default = 3).
#' @param ... Additional arguments affecting the summary produced.
#' @return Prints a table with some summary statistics, including posterior mean, standard deviation and lower and upper quantiles based on the
#' values specified in \code{prob}, for the posterior distributions of the regression coefficients of the effects and costs models run using the 
#' function \code{selection}, \code{pattern}, \code{hurdle} or \code{\link{lmdm}}. 
#' @seealso \code{\link{selection}} \code{\link{pattern}} \code{\link{hurdle}} \code{\link{lmdm}} \code{\link{diagnostic}} \code{\link{plot.missingHE}}
#' @author Andrea Gabrio
#' @importFrom stats quantile
#' @export
#' @examples  
#' # For examples see the function \code{\link{selection}}, \code{\link{pattern}}, 
#' # \code{\link{hurdle}} or \code{\link{lmdm}}
#' # 
#' # 

coef.missingHE <- function(object, prob = c(0.025, 0.975), 
                           random = FALSE, digits = 3, ...) {
   exArgs <- list(...)
   if(!inherits(object,"missingHE")) {
    stop("Only objects of class 'missingHE' can be used")}
   if(!is.vector(prob) | length(prob) != 2 | !is.numeric(prob) | any(prob <= 0 | prob >= 1)) {
     stop("You must provide valid lower/upper quantiles for the coefficients distributions")}
   if(object$data_format == "wide") {
     if(length(dim(object$model_output$effects_fixed)) == 2) {
       cov_e_fixed_name <- colnames(object$model_output$effects_fixed)
       p_e_fixed <- length(cov_e_fixed_name)
       mean_cov_e_fixed <- apply(object$model_output$effects_fixed, 2 , mean)
       sd_cov_e_fixed <- apply(object$model_output$effects_fixed, 2 , sd)
       ql_cov_e_fixed <- apply(object$model_output$effects_fixed, 2 , quantile, prob = prob[1])
       qu_cov_e_fixed <- apply(object$model_output$effects_fixed, 2 , quantile, prob = prob[2])
     }
     if(length(dim(object$model_output$effects_fixed)) == 3) {
       cov_e_fixed_name <- colnames(cbind.data.frame(object$model_output$effects_fixed))
       p_e_fixed <- length(cov_e_fixed_name)
       mean_cov_e_fixed <- c(t(apply(object$model_output$effects_fixed, 2:3 , mean)))
       sd_cov_e_fixed <- c(t(apply(object$model_output$effects_fixed, 2:3 , sd)))
       ql_cov_e_fixed <- c(t(apply(object$model_output$effects_fixed, 2:3 , quantile, prob = prob[1])))
       qu_cov_e_fixed <- c(t(apply(object$model_output$effects_fixed, 2:3 , quantile, prob = prob[2])))
     }     
    if(object$model_output$ind_fixed) {
      if(length(dim(object$model_output$costs_fixed)) == 2) {
        cov_c_fixed_name <- colnames(object$model_output$costs_fixed)
        p_c_fixed <- length(cov_c_fixed_name)
        mean_cov_c_fixed <- apply(object$model_output$costs_fixed, 2 , mean)
        sd_cov_c_fixed <- apply(object$model_output$costs_fixed, 2 , sd)
        ql_cov_c_fixed <- apply(object$model_output$costs_fixed, 2 , quantile, prob = prob[1])
        qu_cov_c_fixed <- apply(object$model_output$costs_fixed, 2 , quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$costs_fixed)) == 3) {
        cov_c_fixed_name <- colnames(cbind.data.frame(object$model_output$costs_fixed))
        p_c_fixed <- length(cov_c_fixed_name)
        mean_cov_c_fixed <- c(t(apply(object$model_output$costs_fixed, 2:3 , mean)))
        sd_cov_c_fixed <- c(t(apply(object$model_output$costs_fixed, 2:3 , sd)))
        ql_cov_c_fixed <- c(t(apply(object$model_output$costs_fixed, 2:3 , quantile, prob = prob[1])))
        qu_cov_c_fixed <- c(t(apply(object$model_output$costs_fixed, 2:3 , quantile, prob = prob[2])))
      }       
    }
    if(!object$model_output$ind_fixed) {
      costs_fixed <- cbind.data.frame(object$model_output$costs_fixed$beta, object$model_output$costs_fixed$beta_f)
      cov_c_fixed_name <- colnames(costs_fixed)
      p_c_fixed <- length(cov_c_fixed_name)
      mean_cov_c_fixed <- apply(costs_fixed, 2 , mean)
      sd_cov_c_fixed <- apply(costs_fixed, 2 , sd)
      ql_cov_c_fixed <- apply(costs_fixed, 2 , quantile, prob = prob[1])
      qu_cov_c_fixed <- apply(costs_fixed, 2 , quantile, prob = prob[2])
    }
   cov_e_random_name <- NULL
   p_e_random <- 0
   if(!is.null(object$model_output$effects_random)) {
     if(length(dim(object$model_output$effects_random)) == 2) {
     cov_e_random_name <- colnames(object$model_output$effects_random)
     p_e_random <- length(cov_e_random_name)
     mean_cov_e_random <- apply(object$model_output$effects_random, 2 , mean)
     sd_cov_e_random <- apply(object$model_output$effects_random, 2 , sd)
     ql_cov_e_random <- apply(object$model_output$effects_random, 2 , quantile, prob = prob[1])
     qu_cov_e_random <- apply(object$model_output$effects_random, 2 , quantile, prob = prob[2])
     }
     if(length(dim(object$model_output$effects_random)) == 3) {
       cov_e_random_name <- colnames(cbind.data.frame(object$model_output$effects_random))
       p_e_random <- length(cov_e_random_name)
       mean_cov_e_random <- c(t(apply(object$model_output$effects_random, 2:3 , mean)))
       sd_cov_e_random <- c(t(apply(object$model_output$effects_random, 2:3 , sd)))
       ql_cov_e_random <- c(t(apply(object$model_output$effects_random, 2:3 , quantile, prob = prob[1])))
       qu_cov_e_random <- c(t(apply(object$model_output$effects_random, 2:3 , quantile, prob = prob[2])))
     } 
   }
   cov_c_random_name <- NULL
   p_c_random <- 0
   if(!is.null(object$model_output$costs_random)) {
    if(object$model_output$ind_random) {
      if(length(dim(object$model_output$costs_random)) == 2) {
     cov_c_random_name <- colnames(object$model_output$costs_random)
     p_c_random <- length(cov_c_random_name)
     mean_cov_c_random <- apply(object$model_output$costs_random, 2 , mean)
     sd_cov_c_random <- apply(object$model_output$costs_random, 2 , sd)
     ql_cov_c_random <- apply(object$model_output$costs_random, 2 , quantile, prob = prob[1])
     qu_cov_c_random <- apply(object$model_output$costs_random, 2 , quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$costs_random)) == 3) {
        cov_c_random_name <- colnames(cbind.data.frame(object$model_output$costs_random))
        p_c_random <- length(cov_c_random_name)
        mean_cov_c_random <- c(t(apply(object$model_output$costs_random, 2:3 , mean)))
        sd_cov_c_random <- c(t(apply(object$model_output$costs_random, 2:3 , sd)))
        ql_cov_c_random <- c(t(apply(object$model_output$costs_random, 2:3 , quantile, prob = prob[1])))
        qu_cov_c_random <- c(t(apply(object$model_output$costs_random, 2:3 , quantile, prob = prob[2])))
      }      
    }
    if(!object$model_output$ind_random) {
      if(length(object$model_output$costs_random) == 2) {
        costs_random <- cbind.data.frame(object$model_output$costs_random$b, object$model_output$costs_random$b_f)
        cov_c_random_name <- colnames(costs_random)
        p_c_random <- length(cov_c_random_name)
        mean_cov_c_random <- apply(costs_random, 2 , mean)
        sd_cov_c_random <- apply(costs_random, 2 , sd)
        ql_cov_c_random <- apply(costs_random, 2 , quantile, prob = prob[1])
        qu_cov_c_random <- apply(costs_random, 2 , quantile, prob = prob[2])
      }
      if(length(object$model_output$costs_random) == 1) {
        cov_c_random_name <- colnames(object$model_output$costs_random$b_f)
        p_c_random <- length(cov_c_random_name)
        mean_cov_c_random <- apply(object$model_output$costs_random$b_f, 2 , mean)
        sd_cov_c_random <- apply(object$model_output$costs_random$b_f, 2 , sd)
        ql_cov_c_random <- apply(object$model_output$costs_random$b_f, 2 , quantile, prob = prob[1])
        qu_cov_c_random <- apply(object$model_output$costs_random$b_f, 2 , quantile, prob = prob[2])
      }
     }
   }
 }
 if(object$data_format == "long") {
  max_time <- max(object$data_set$data_raw$time_long, na.rm = TRUE)
  if(is.null(names(object$model_output$effects_fixed))) {
    dimnames(object$model_output$effects_fixed)[[3]] <- paste0("time.", dimnames(object$model_output$effects_fixed)[[3]])
    effects_fixed <- as.data.frame(object$model_output$effects_fixed)
  }
  if(!is.null(names(object$model_output$effects_fixed))) {
    dimnames(object$model_output$effects_fixed$alpha)[[3]] <- paste0("time.", dimnames(object$model_output$effects_fixed$alpha)[[3]])
    if(!all(c("alpha_te", "alpha_tc") %in% names(object$model_output$effects_fixed))) {
      effects_fixed <- as.data.frame(object$model_output$effects_fixed$alpha)
    }
    if(all(c("alpha_te", "alpha_tc") %in% names(object$model_output$effects_fixed))) { 
      effects_fixed <- cbind(as.data.frame(object$model_output$effects_fixed$alpha), as.data.frame(object$model_output$effects_fixed$alpha_te), 
                             as.data.frame(object$model_output$effects_fixed$alpha_tc))
    }
  }
  cov_e_fixed_name <- colnames(effects_fixed)
  p_e_fixed <- length(cov_e_fixed_name)
  mean_cov_e_fixed <- apply(effects_fixed, 2 , mean)
  sd_cov_e_fixed <- apply(effects_fixed, 2 , sd)
  ql_cov_e_fixed <- apply(effects_fixed, 2 , quantile, prob = prob[1])
  qu_cov_e_fixed <- apply(effects_fixed, 2 , quantile, prob = prob[2])
  if(is.null(names(object$model_output$costs_fixed))) {
    dimnames(object$model_output$costs_fixed)[[3]] <- paste0("time.", dimnames(object$model_output$costs_fixed)[[3]])
    costs_fixed <- as.data.frame(object$model_output$costs_fixed)
  }
  if(!is.null(names(object$model_output$costs_fixed))) {
    dimnames(object$model_output$costs_fixed$beta)[[3]] <- paste0("time.", dimnames(object$model_output$costs_fixed$beta)[[3]])
    if(!all(c("beta_f", "beta_te", "beta_tc") %in% names(object$model_output$costs_fixed))) {
      costs_fixed <- as.data.frame(object$model_output$costs_fixed$beta)
    }
    if(!all(c("beta_te", "beta_tc") %in% names(object$model_output$costs_fixed)) & "beta_f" %in% names(object$model_output$costs_fixed)) {
      costs_fixed <- cbind(as.data.frame(object$model_output$costs_fixed$beta), as.data.frame(object$model_output$costs_fixed$beta_f))
    }
    if(all(c("beta_te", "beta_tc") %in% names(object$model_output$costs_fixed)) & !"beta_f" %in% names(object$model_output$costs_fixed)) {
      costs_fixed <- cbind(as.data.frame(object$model_output$costs_fixed$beta), 
                           as.data.frame(object$model_output$costs_fixed$beta_te),
                           as.data.frame(object$model_output$costs_fixed$beta_tc))
    }
    if(all(c("beta_f", "beta_te", "beta_tc") %in% names(object$model_output$costs_fixed))) {
      costs_fixed <- cbind(as.data.frame(object$model_output$costs_fixed$beta), 
                           as.data.frame(object$model_output$costs_fixed$beta_f),
                           as.data.frame(object$model_output$costs_fixed$beta_te),
                           as.data.frame(object$model_output$costs_fixed$beta_tc))
    }
  }
  cov_c_fixed_name <- colnames(costs_fixed)
  p_c_fixed <- length(cov_c_fixed_name)
  mean_cov_c_fixed <- apply(costs_fixed, 2 , mean)
  sd_cov_c_fixed <- apply(costs_fixed, 2 , sd)
  ql_cov_c_fixed <- apply(costs_fixed, 2 , quantile, prob = prob[1])
  qu_cov_c_fixed <- apply(costs_fixed, 2 , quantile, prob = prob[2])
  cov_e_random_name <- NULL
  p_e_random <- 0
  if(!is.null(object$model_output$effects_random)) {
  if(is.null(names(object$model_output$effects_random))) {
    dimnames(object$model_output$effects_random)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random)[[3]])
    effects_random <- as.data.frame(object$model_output$effects_random)
  }
  if(!is.null(names(object$model_output$effects_random))) {
    if(!all(c("a_te", "a_tc") %in% names(object$model_output$effects_random)) & "a" %in% names(object$model_output$effects_random)) {
      dimnames(object$model_output$effects_random$a)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a)[[3]])
      effects_random <- as.data.frame(object$model_output$effects_random$a)
    }
    if(all(c("a_te", "a_tc") %in% names(object$model_output$effects_random)) & !"a" %in% names(object$model_output$effects_random)) {
      dimnames(object$model_output$effects_random$a_te)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a_te)[[3]])
      dimnames(object$model_output$effects_random$a_tc)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a_tc)[[3]])
      effects_random <- cbind(as.data.frame(object$model_output$effects_random$a_te), 
                              as.data.frame(object$model_output$effects_random$a_tc))
    }
    if(all(c("a", "a_te", "a_tc") %in% names(object$model_output$effects_random))) { 
      dimnames(object$model_output$effects_random$a)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a)[[3]])
      dimnames(object$model_output$effects_random$a_te)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a_te)[[3]])
      dimnames(object$model_output$effects_random$a_tc)[[3]] <- paste0("time.", dimnames(object$model_output$effects_random$a_tc)[[3]])
      effects_random <- cbind(as.data.frame(object$model_output$effects_random$a), as.data.frame(object$model_output$effects_random$a_te), 
                             as.data.frame(object$model_output$effects_random$a_tc))
    }
  }    
   cov_e_random_name <- colnames(effects_random)
   p_e_random <- length(cov_e_random_name)
   mean_cov_e_random <- apply(effects_random, 2 , mean)
   sd_cov_e_random <- apply(effects_random, 2 , sd)
   ql_cov_e_random <- apply(effects_random, 2 , quantile, prob = prob[1])
   qu_cov_e_random <- apply(effects_random, 2 , quantile, prob = prob[2])    
  }
  cov_c_random_name <- NULL
  p_c_random <- 0
  if(!is.null(object$model_output$costs_random)) {
    if(is.null(names(object$model_output$costs_random))) {
      dimnames(object$model_output$costs_random)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random)[[3]])
      costs_random <- as.data.frame(object$model_output$costs_random)
    }
    if(!is.null(names(object$model_output$costs_random))) {
      if(!all(c("b_te", "b_tc", "b_f") %in% names(object$model_output$costs_random)) & "b" %in% names(object$model_output$costs_random)) {
        dimnames(object$model_output$costs_random$b)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b)[[3]])
        costs_random <- as.data.frame(object$model_output$costs_random$b)
      }
      if(!all(c("b_te", "b_tc") %in% names(object$model_output$costs_random)) & all(c("b", "b_f") %in% names(object$model_output$costs_random))) {
        dimnames(object$model_output$costs_random$b)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b)[[3]])
        dimnames(object$model_output$costs_random$b_f)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_f)[[3]])
        costs_random <- cbind(as.data.frame(object$model_output$costs_random$b), as.data.frame(object$model_output$costs_random$b_f))
      }
      if(all(c("b_te", "b_tc") %in% names(object$model_output$costs_random)) & !all(c("b", "b_f") %in% names(object$model_output$costs_random))) {
        dimnames(object$model_output$costs_random$b_te)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_te)[[3]])
        dimnames(object$model_output$costs_random$b_tc)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_tc)[[3]])
        costs_random <- cbind(as.data.frame(object$model_output$costs_random$b_te), 
                                as.data.frame(object$model_output$costs_random$b_tc))
      }
      if(all(c("b_te", "b_tc", "b_f") %in% names(object$model_output$costs_random)) & !"b" %in% names(object$model_output$costs_random)) {
        dimnames(object$model_output$costs_random$b_f)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_f)[[3]])
        dimnames(object$model_output$costs_random$b_te)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_te)[[3]])
        dimnames(object$model_output$costs_random$b_tc)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_tc)[[3]])
        costs_random <- cbind(as.data.frame(object$model_output$costs_random$b_te), 
                              as.data.frame(object$model_output$costs_random$b_tc))
      }
      if(all(c("b", "b_te", "b_tc") %in% names(object$model_output$costs_random)) & !"b_f" %in% names(object$model_output$costs_random)) { 
        dimnames(object$model_output$costs_random$b)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b)[[3]])
        dimnames(object$model_output$costs_random$b_te)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_te)[[3]])
        dimnames(object$model_output$costs_random$b_tc)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_tc)[[3]])
        costs_random <- cbind(as.data.frame(object$model_output$costs_random$b), as.data.frame(object$model_output$costs_random$b_te), 
                              as.data.frame(object$model_output$costs_random$b_tc))
      }
      if(all(c("b", "b_te", "b_tc", "b_f") %in% names(object$model_output$costs_random))) { 
        dimnames(object$model_output$costs_random$b)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b)[[3]])
        dimnames(object$model_output$costs_random$b_f)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_f)[[3]])
        dimnames(object$model_output$costs_random$b_te)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_te)[[3]])
        dimnames(object$model_output$costs_random$b_tc)[[3]] <- paste0("time.", dimnames(object$model_output$costs_random$b_tc)[[3]])
        costs_random <- cbind(as.data.frame(object$model_output$costs_random$b), as.data.frame(object$model_output$costs_random$b_f), 
                              as.data.frame(object$model_output$costs_random$b_te), as.data.frame(object$model_output$costs_random$b_tc))
      }
    }    
    cov_c_random_name <- colnames(costs_random)
    p_c_random <- length(cov_c_random_name)
    mean_cov_c_random <- apply(costs_random, 2 , mean)
    sd_cov_c_random <- apply(costs_random, 2 , sd)
    ql_cov_c_random <- apply(costs_random, 2 , quantile, prob = prob[1])
    qu_cov_c_random <- apply(costs_random, 2 , quantile, prob = prob[2])    
  }
 } 
   tbl_colname <- c("Mean", "SD", "QL", "QU")
   table_e_fixed <- matrix(NA, nrow = p_e_fixed, ncol = 4)
   rownames(table_e_fixed) <- cov_e_fixed_name
   colnames(table_e_fixed) <- tbl_colname
   table_e_fixed[, 1] <- mean_cov_e_fixed
   table_e_fixed[, 2] <- sd_cov_e_fixed
   table_e_fixed[, 3] <- ql_cov_e_fixed
   table_e_fixed[, 4] <- qu_cov_e_fixed
   table_e_fixed <- round(table_e_fixed, digits = digits)
   table_c_fixed <- matrix(NA, nrow = p_c_fixed, ncol = 4)
   rownames(table_c_fixed) <- cov_c_fixed_name
   colnames(table_c_fixed) <- tbl_colname
   table_c_fixed[, 1] <- mean_cov_c_fixed
   table_c_fixed[, 2] <- sd_cov_c_fixed
   table_c_fixed[, 3] <- ql_cov_c_fixed
   table_c_fixed[, 4] <- qu_cov_c_fixed
   table_c_fixed <- round(table_c_fixed, digits = digits)
   table_list_random <- NULL
   table_e_random <- table_c_random <- NULL
   clus_e <- clus_c <- NULL
   if(p_e_random != 0) {
     mean_cov_e_random <- mean_cov_e_random
     sd_cov_e_random <- sd_cov_e_random
     ql_cov_e_random <- ql_cov_e_random
     qu_cov_e_random <- qu_cov_e_random
     table_e_random <- matrix(NA, nrow = p_e_random, ncol = 4)
     rownames(table_e_random) <- cov_e_random_name
     colnames(table_e_random) <- tbl_colname
     table_e_random[, 1] <- mean_cov_e_random
     table_e_random[, 2] <- sd_cov_e_random
     table_e_random[, 3] <- ql_cov_e_random
     table_e_random[, 4] <- qu_cov_e_random
     table_e_random <- round(table_e_random, digits = digits)
   }
   if(p_c_random != 0) {
     mean_cov_c_random <- mean_cov_c_random
     sd_cov_c_random <- sd_cov_c_random
     ql_cov_c_random <- ql_cov_c_random
     qu_cov_c_random <- qu_cov_c_random
     table_c_random <- matrix(NA, nrow = p_c_random, ncol = 4)
     rownames(table_c_random) <- cov_c_random_name
     colnames(table_c_random) <- tbl_colname
     table_c_random[, 1] <- mean_cov_c_random
     table_c_random[, 2] <- sd_cov_c_random
     table_c_random[, 3] <- ql_cov_c_random
     table_c_random[, 4] <- qu_cov_c_random
     table_c_random <- round(table_c_random, digits = digits)
   }
  table_list_fixed <- list("Effects" = table_e_fixed, "Costs" = table_c_fixed)
  table_list_random <- list("Effects" = table_e_random, "Costs" = table_c_random)
  if(random & p_c_random == 0 & p_e_random == 0) {
     stop("No random effects estimates found")}   
 if(!random) { print(table_list_fixed)}
 if(random) { print(table_list_random)}
}