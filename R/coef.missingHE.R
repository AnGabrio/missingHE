#' Extract regression coefficient estimates from objects in the class \code{missingHE}
#'
#' Produces a table printout with summary statistics for the regression coefficients of the health economic evaluation probabilistic model
#' run using the function \code{\link{selection}}, \code{\link{selection_long}}, \code{\link{pattern}} or \code{\link{hurdle}}.
#' @param object A \code{missingHE} object containing the results of the Bayesian modelling and the economic evaluation
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the estimates.
#' @param random Logical. If \code{random} is \code{TRUE}, the estimates of the random effects parameters are printed, when available.
#' @param time A number indicating the time point at which posterior results for the model coefficients should be reported (only for longitudinal models).
#' @param digits Number of digits to be displayed for each estimate.
#' @param ... Additional arguments affecting the summary produced.
#' @return Prints a table with some summary statistics, including posterior mean, standard deviation and lower and upper quantiles based on the
#' values specified in \code{prob}, for the posterior distributions of the regression coefficients of the effects and costs models run using the 
#' function \code{selection}, \code{\link{selection_long}}, \code{pattern} or \code{hurdle}. 
#' @seealso \code{\link{selection}} \code{\link{selection_long}} \code{\link{pattern}} \code{\link{hurdle}} \code{\link{diagnostic}} \code{\link{plot.missingHE}}
#' @author Andrea Gabrio
#' @importFrom stats quantile
#' @export
#' @examples  
#' # For examples see the function \code{\link{selection}}, \code{\link{selection_long}},
#' # \code{\link{pattern}} or \code{\link{hurdle}}
#' # 
#' # 

coef.missingHE <- function(object, prob = c(0.025, 0.975), random = FALSE, time = 1, digits = 3, ...) {
   exArgs <- list(...)
   if(!inherits(object, what = "missingHE")) {
    stop("Only objects of class 'missingHE' can be used")
   }
   if(length(prob) != 2 | is.numeric(prob) == FALSE | any(prob < 0) != FALSE | any(prob > 1) != FALSE) {
    stop("You must provide valid lower/upper quantiles for the coefficients distributions")
   }
   if(object$data_format == "wide") {
   if(length(grep("^SELECTION", object$model_output$type)) == 1 | length(grep("^HURDLE", object$model_output$type)) == 1) {
   cov_e_fixed <- names(object$data_set$covariates_effects_fixed$Control)
   p_e_fixed <- length(cov_e_fixed)
   if(length(dim(object$model_output$covariate_parameter_effects_fixed)) == 2) {
    mean_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, 2 , mean)
    sd_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, 2 , sd)
    ql_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, 2 , quantile, prob = prob[1])
    qu_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, 2 , quantile, prob = prob[2])
   } else if(length(dim(object$model_output$covariate_parameter_effects_fixed)) == 3) {
    mean_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, c(2, 3) , mean)
    sd_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, c(2, 3) , sd)
    ql_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, c(2, 3) , quantile, prob = prob[1])
    qu_cov_e_fixed <- apply(object$model_output$covariate_parameter_effects_fixed, c(2, 3) , quantile, prob = prob[2])
   }
   cov_c_fixed <- names(object$data_set$covariates_costs_fixed$Control)
   p_c_fixed <- length(cov_c_fixed)
   if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 2 & object$model_output$ind_fixed == TRUE) {
    mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, 2 , mean)
    sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, 2 , sd)
    ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, 2 , quantile, prob = prob[1])
    qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, 2 , quantile, prob = prob[2])
   } else if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 3 & object$model_output$ind_fixed == TRUE) {
    mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, c(2, 3) , mean)
    sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, c(2, 3) , sd)
    ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, c(2, 3) , quantile, prob = prob[1])
    qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed, c(2, 3) , quantile, prob = prob[2])
   }
   if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 0) {
    if(length(dim(object$model_output$covariate_parameter_costs_fixed$beta)) == 2 & object$model_output$ind_fixed == FALSE) {
    mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, 2 , mean)
    sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, 2 , sd)
    ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, 2 , quantile, prob = prob[1])
    qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_costs_fixed$beta)) == 3 & object$model_output$ind_fixed == FALSE) {
    mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, c(2, 3) , mean)
    sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, c(2, 3) , sd)
    ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, c(2, 3) , quantile, prob = prob[1])
    qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta, c(2, 3) , quantile, prob = prob[2])
    }
   }
   if(object$model_output$ind_fixed == FALSE) {
    dep_c_fixed <- "e"
    mean_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f, 2, mean)
    sd_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f, 2 , sd)
    ql_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f, 2 , quantile, prob = prob[1])
    qu_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f, 2 , quantile, prob = prob[2])
   } else {dep_c_fixed <- NULL }
   cov_e_random <- names(object$data_set$covariates_effects_random$Control)
   p_e_random <- length(cov_e_random)
   if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 2) {
    mean_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , mean)
    sd_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , sd)
    ql_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , quantile, prob = prob[1])
    qu_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , quantile, prob = prob[2])
    mean_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , mean)
    sd_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , sd)
    ql_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , quantile, prob = prob[1])
    qu_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , quantile, prob = prob[2])
   } else if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 3) {
    mean_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , mean)
    sd_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , sd)
    ql_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , quantile, prob = prob[1])
    qu_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , quantile, prob = prob[2])
    mean_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , mean)
    sd_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , sd)
    ql_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , quantile, prob = prob[1])
    qu_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , quantile, prob = prob[2])
   }
   cov_c_random <- names(object$data_set$covariates_costs_random$Control)
   p_c_random <- length(cov_c_random)
   if(length(dim(object$model_output$covariate_parameter_costs_random)) == 2) {
    mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , mean)
    sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , sd)
    ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , quantile, prob = prob[1])
    qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , quantile, prob = prob[2])
   } else if(length(dim(object$model_output$covariate_parameter_costs_random)) == 3) {
    mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , mean)
    sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , sd)
    ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , quantile, prob = prob[1])
    qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , quantile, prob = prob[2])
   }
   if(length(dim(object$model_output$covariate_parameter_costs_random)) == 0 & is.null(object$model_output$covariate_parameter_costs_random) == FALSE) {
    if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 2) {
    mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , mean)
    sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , sd)
    ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , quantile, prob = prob[1])
    qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , quantile, prob = prob[2])
    mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , mean)
    sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , sd)
    ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , quantile, prob = prob[1])
    qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 3) {
    mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , mean)
    sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , sd)
    ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , quantile, prob = prob[1])
    qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , quantile, prob = prob[2])
    mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , mean)
    sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , sd)
    ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , quantile, prob = prob[1])
    qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , quantile, prob = prob[2])
    }
   }
   if(object$model_output$ind_random == FALSE) {
    dep_c_random <- "e"
    mean_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2, mean)
    sd_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , sd)
    ql_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , quantile, prob = prob[1])
    qu_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , quantile, prob = prob[2])
    mean_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2, mean)
    sd_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , sd)
    ql_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , quantile, prob = prob[1])
    qu_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , quantile, prob = prob[2])
   } else {dep_c_random <- NULL }
   if(object$model_output$ind_fixed == FALSE) {
     cov_c_fixed <- c(cov_c_fixed, dep_c_fixed)
     p_c_fixed <- length(cov_c_fixed)
     mean_cov_c_fixed <- rbind(mean_cov_c_fixed, mean_dep_c_fixed)
     sd_cov_c_fixed <- rbind(sd_cov_c_fixed, sd_dep_c_fixed)
     ql_cov_c_fixed <- rbind(ql_cov_c_fixed, ql_dep_c_fixed)
     qu_cov_c_fixed <- rbind(qu_cov_c_fixed, qu_dep_c_fixed)
   }
   if(is.matrix(mean_cov_e_fixed) == FALSE) {
     mean_cov_e_fixed_arm1 <- mean_cov_e_fixed[1]
     sd_cov_e_fixed_arm1 <- sd_cov_e_fixed[1]
     ql_cov_e_fixed_arm1 <- ql_cov_e_fixed[1]
     qu_cov_e_fixed_arm1 <- qu_cov_e_fixed[1]
     mean_cov_e_fixed_arm2 <- mean_cov_e_fixed[2]
     sd_cov_e_fixed_arm2 <- sd_cov_e_fixed[2]
     ql_cov_e_fixed_arm2 <- ql_cov_e_fixed[2]
     qu_cov_e_fixed_arm2 <- qu_cov_e_fixed[2]
   } else if(is.matrix(mean_cov_e_fixed) == TRUE) {
     mean_cov_e_fixed_arm1 <- mean_cov_e_fixed[, 1]
     sd_cov_e_fixed_arm1 <- sd_cov_e_fixed[, 1]
     ql_cov_e_fixed_arm1 <- ql_cov_e_fixed[, 1]
     qu_cov_e_fixed_arm1 <- qu_cov_e_fixed[, 1]
     mean_cov_e_fixed_arm2 <- mean_cov_e_fixed[, 2]
     sd_cov_e_fixed_arm2 <- sd_cov_e_fixed[, 2]
     ql_cov_e_fixed_arm2 <- ql_cov_e_fixed[, 2]
     qu_cov_e_fixed_arm2 <- qu_cov_e_fixed[, 2]
   } 
   if(is.matrix(mean_cov_c_fixed) == FALSE) {
     mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[1]
     sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[1]
     ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[1]
     qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[1]
     mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[2]
     sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[2]
     ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[2]
     qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[2]
   } else if(is.matrix(mean_cov_c_fixed) == TRUE) {
     mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[, 1]
     sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[, 1]
     ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[, 1]
     qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[, 1]
     mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[, 2]
     sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[, 2]
     ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[, 2]
     qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[, 2]
   }
  }
  if(length(grep("^PATTERN", object$model_output$type)) == 1) {
    cov_e_fixed <- names(object$data_set$covariates_effects_fixed$Control)
    p_e_fixed <- length(cov_e_fixed)
    if(length(dim(object$model_output$covariate_parameter_effects_fixed_pattern$control)) == 2) {
      mean_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, 2 , mean)
      sd_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, 2 , sd)
      ql_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, 2 , quantile, prob = prob[1])
      qu_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_effects_fixed_pattern$control)) == 3) {
      mean_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, c(2, 3) , mean)
      sd_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, c(2, 3) , sd)
      ql_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, c(2, 3) , quantile, prob = prob[1])
      qu_cov_e1_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$control, c(2, 3) , quantile, prob = prob[2])
    }
    if(length(dim(object$model_output$covariate_parameter_effects_fixed_pattern$intervention)) == 2) {
      mean_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, 2 , mean)
      sd_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, 2 , sd)
      ql_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, 2 , quantile, prob = prob[1])
      qu_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_effects_fixed_pattern$intervention)) == 3) {
      mean_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, c(2, 3) , mean)
      sd_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, c(2, 3) , sd)
      ql_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, c(2, 3) , quantile, prob = prob[1])
      qu_cov_e2_fixed <- apply(object$model_output$covariate_parameter_effects_fixed_pattern$intervention, c(2, 3) , quantile, prob = prob[2])
    }
    cov_c_fixed <- names(object$data_set$covariates_costs_fixed$Control)
    p_c_fixed <- length(cov_c_fixed)
    if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$control)) == 2) {
      mean_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, 2 , mean)
      sd_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, 2 , sd)
      ql_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, 2 , quantile, prob = prob[1])
      qu_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$control)) == 3) {
      mean_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, c(2, 3) , mean)
      sd_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, c(2, 3) , sd)
      ql_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control, c(2, 3) , quantile, prob = prob[2])
    }
    if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$intervention)) == 2) {
      mean_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, 2 , mean)
      sd_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, 2 , sd)
      ql_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, 2 , quantile, prob = prob[1])
      qu_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$intervention)) == 3) {
      mean_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, c(2, 3) , mean)
      sd_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, c(2, 3) , sd)
      ql_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention, c(2, 3) , quantile, prob = prob[2])
    }
    if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$control)) == 0) {
     if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1)) == 2) {
      mean_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, 2 , mean)
      sd_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, 2 , sd)
      ql_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, 2 , quantile, prob = prob[1])
      qu_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, 2 , quantile, prob = prob[2])
     } else if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1)) == 3) {
      mean_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, c(2, 3) , mean)
      sd_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, c(2, 3) , sd)
      ql_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_p1, c(2, 3) , quantile, prob = prob[2])
     }
    }
    if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$intervention)) == 0) {
     if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2)) == 2) {
      mean_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, 2 , mean)
      sd_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, 2 , sd)
      ql_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, 2 , quantile, prob = prob[1])
      qu_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, 2 , quantile, prob = prob[2])
     } else if(length(dim(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2)) == 3) {
      mean_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, c(2, 3) , mean)
      sd_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, c(2, 3) , sd)
      ql_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_p2, c(2, 3) , quantile, prob = prob[2])
     }
    }
    if(object$model_output$ind_fixed == FALSE) {
      dep_c_fixed <- "e"
      mean_dep_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_f_p1, 2, mean)
      sd_dep_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_f_p1, 2 , sd)
      ql_dep_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_f_p1, 2 , quantile, prob = prob[1])
      qu_dep_c1_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$control$beta_f_p1, 2 , quantile, prob = prob[2])
      mean_dep_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_f_p2, 2, mean)
      sd_dep_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_f_p2, 2 , sd)
      ql_dep_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_f_p2, 2 , quantile, prob = prob[1])
      qu_dep_c2_fixed <- apply(object$model_output$covariate_parameter_costs_fixed_pattern$intervention$beta_f_p2, 2 , quantile, prob = prob[2])
    } else {dep_c_fixed <- NULL }
    cov_e_random <- names(object$data_set$covariates_effects_random$Control)
    p_e_random <- length(cov_e_random)
    if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 2) {
      mean_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , mean)
      sd_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , sd)
      ql_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , quantile, prob = prob[1])
      qu_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, 2 , quantile, prob = prob[2])
      mean_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , mean)
      sd_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , sd)
      ql_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , quantile, prob = prob[1])
      qu_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 3) {
      mean_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , mean)
      sd_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , sd)
      ql_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , quantile, prob = prob[1])
      qu_cov_e1_random <- apply(object$model_output$covariate_parameter_effects_random$a1, c(2, 3) , quantile, prob = prob[2])
      mean_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , mean)
      sd_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , sd)
      ql_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , quantile, prob = prob[1])
      qu_cov_e2_random <- apply(object$model_output$covariate_parameter_effects_random$a2, c(2, 3) , quantile, prob = prob[2])
    }
    cov_c_random <- names(object$data_set$covariates_costs_random$Control)
    p_c_random <- length(cov_c_random)
    if(length(dim(object$model_output$covariate_parameter_costs_random)) == 2) {
      mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , mean)
      sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , sd)
      ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , quantile, prob = prob[1])
      qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, 2 , quantile, prob = prob[2])
    } else if(length(dim(object$model_output$covariate_parameter_costs_random)) == 3) {
      mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , mean)
      sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , sd)
      ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random, c(2, 3) , quantile, prob = prob[2])
    }
    if(length(dim(object$model_output$covariate_parameter_costs_random)) == 0 & is.null(object$model_output$covariate_parameter_costs_random) == FALSE) {
     if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 2) {
      mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , mean)
      sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , sd)
      ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , quantile, prob = prob[1])
      qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, 2 , quantile, prob = prob[2])
      mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , mean)
      sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , sd)
      ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , quantile, prob = prob[1])
      qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, 2 , quantile, prob = prob[2])
     } else if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 3) {
      mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , mean)
      sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , sd)
      ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1, c(2, 3) , quantile, prob = prob[2])
      mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , mean)
      sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , sd)
      ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , quantile, prob = prob[1])
      qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2, c(2, 3) , quantile, prob = prob[2])
     }
    }
    if(object$model_output$ind_random == FALSE) {
      dep_c_random <- "e"
      mean_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2, mean)
      sd_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , sd)
      ql_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , quantile, prob = prob[1])
      qu_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f, 2 , quantile, prob = prob[2])
      mean_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2, mean)
      sd_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , sd)
      ql_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , quantile, prob = prob[1])
      qu_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f, 2 , quantile, prob = prob[2])
    } else {dep_c_random <- NULL }
    if(object$model_output$ind_fixed == FALSE) {
      cov_c_fixed <- c(cov_c_fixed, dep_c_fixed)
      p_c_fixed <- length(cov_c_fixed)
      mean_cov_c1_fixed <- rbind(mean_cov_c1_fixed, mean_dep_c1_fixed)
      sd_cov_c1_fixed <- rbind(sd_cov_c1_fixed, sd_dep_c1_fixed)
      ql_cov_c1_fixed <- rbind(ql_cov_c1_fixed, ql_dep_c1_fixed)
      qu_cov_c1_fixed <- rbind(qu_cov_c1_fixed, qu_dep_c1_fixed)
      mean_cov_c2_fixed <- rbind(mean_cov_c2_fixed, mean_dep_c2_fixed)
      sd_cov_c2_fixed <- rbind(sd_cov_c2_fixed, sd_dep_c2_fixed)
      ql_cov_c2_fixed <- rbind(ql_cov_c2_fixed, ql_dep_c2_fixed)
      qu_cov_c2_fixed <- rbind(qu_cov_c2_fixed, qu_dep_c2_fixed)
    }
    mean_cov_e_fixed_arm1 <- mean_cov_e1_fixed
    sd_cov_e_fixed_arm1 <- sd_cov_e1_fixed
    ql_cov_e_fixed_arm1 <- ql_cov_e1_fixed
    qu_cov_e_fixed_arm1 <- qu_cov_e1_fixed
    mean_cov_e_fixed_arm2 <- mean_cov_e2_fixed
    sd_cov_e_fixed_arm2 <- sd_cov_e2_fixed
    ql_cov_e_fixed_arm2 <- ql_cov_e2_fixed
    qu_cov_e_fixed_arm2 <- qu_cov_e2_fixed
    mean_cov_c_fixed_arm1 <- mean_cov_c1_fixed
    sd_cov_c_fixed_arm1 <- sd_cov_c1_fixed
    ql_cov_c_fixed_arm1 <- ql_cov_c1_fixed
    qu_cov_c_fixed_arm1 <- qu_cov_c1_fixed
    mean_cov_c_fixed_arm2 <- mean_cov_c2_fixed
    sd_cov_c_fixed_arm2 <- sd_cov_c2_fixed
    ql_cov_c_fixed_arm2 <- ql_cov_c2_fixed
    qu_cov_c_fixed_arm2 <- qu_cov_c2_fixed
  }
   if(length(grep("^PATTERN", object$model_output$type)) == 1) {
     pat_arm1 <- length(unique(as.numeric(object$data_set$`patterns in comparator arm`)))
     pat_arm2 <- length(unique(as.numeric(object$data_set$`patterns in reference arm`)))
     pat_e1_index <- sort(rep(seq(1:pat_arm1), p_e_fixed))
     pat_e2_index <- sort(rep(seq(1:pat_arm2), p_e_fixed))
     table_e1_fixed <- matrix(NA, nrow = length(pat_e1_index), ncol = 4)
     rownames(table_e1_fixed) <- paste(rep(cov_e_fixed, pat_arm1), pat_e1_index, sep = " pattern")
     colnames(table_e1_fixed) <- c("mean", "sd", "lower", "upper")
     table_e1_fixed[, 1] <- c(mean_cov_e_fixed_arm1)
     table_e1_fixed[, 2] <- c(sd_cov_e_fixed_arm1)
     table_e1_fixed[, 3] <- c(ql_cov_e_fixed_arm1)
     table_e1_fixed[, 4] <- c(qu_cov_e_fixed_arm1)
     table_e1_fixed <- round(table_e1_fixed, digits = digits)
     table_e2_fixed <- matrix(NA, nrow = length(pat_e2_index), ncol = 4)
     rownames(table_e2_fixed) <- paste(rep(cov_e_fixed, pat_arm2), pat_e2_index, sep = " pattern")
     colnames(table_e2_fixed) <- c("mean", "sd", "lower", "upper")
     table_e2_fixed[, 1] <- c(mean_cov_e_fixed_arm2)
     table_e2_fixed[, 2] <- c(sd_cov_e_fixed_arm2)
     table_e2_fixed[, 3] <- c(ql_cov_e_fixed_arm2)
     table_e2_fixed[, 4] <- c(qu_cov_e_fixed_arm2)
     table_e2_fixed <- round(table_e2_fixed, digits = digits)
     pat_c1_index <- sort(rep(seq(1:pat_arm1), p_c_fixed))
     pat_c2_index <- sort(rep(seq(1:pat_arm2), p_c_fixed))
     table_c1_fixed <- matrix(NA, nrow = length(pat_c1_index), ncol = 4)
     rownames(table_c1_fixed) <- paste(rep(cov_c_fixed, pat_arm1), pat_c1_index, sep = " pattern")
     colnames(table_c1_fixed) <- c("mean", "sd", "lower", "upper")
     table_c1_fixed[, 1] <- c(mean_cov_c_fixed_arm1)
     table_c1_fixed[, 2] <- c(sd_cov_c_fixed_arm1)
     table_c1_fixed[, 3] <- c(ql_cov_c_fixed_arm1)
     table_c1_fixed[, 4] <- c(qu_cov_c_fixed_arm1)
     table_c1_fixed <- round(table_c1_fixed, digits = digits)
     table_c2_fixed <- matrix(NA, nrow = length(pat_c2_index), ncol = 4)
     rownames(table_c2_fixed) <- paste(rep(cov_c_fixed, pat_arm2), pat_c2_index, sep = " pattern")
     colnames(table_c2_fixed) <- c("mean", "sd", "lower", "upper")
     table_c2_fixed[, 1] <- c(mean_cov_c_fixed_arm2)
     table_c2_fixed[, 2] <- c(sd_cov_c_fixed_arm2)
     table_c2_fixed[, 3] <- c(ql_cov_c_fixed_arm2)
     table_c2_fixed[, 4] <- c(qu_cov_c_fixed_arm2)
     table_c2_fixed <- round(table_c2_fixed, digits = digits)
  } 
  if(length(grep("^SELECTION", object$model_output$type)) == 1 | length(grep("^HURDLE", object$model_output$type)) == 1) {
  table_e1_fixed <- matrix(NA, nrow = p_e_fixed, ncol = 4)
  rownames(table_e1_fixed) <- cov_e_fixed
  colnames(table_e1_fixed) <- c("mean", "sd", "lower", "upper")
  table_e1_fixed[, 1] <- c(mean_cov_e_fixed_arm1)
  table_e1_fixed[, 2] <- c(sd_cov_e_fixed_arm1)
  table_e1_fixed[, 3] <- c(ql_cov_e_fixed_arm1)
  table_e1_fixed[, 4] <- c(qu_cov_e_fixed_arm1)
  table_e1_fixed <- round(table_e1_fixed, digits = digits)
  table_e2_fixed <- matrix(NA, nrow = p_e_fixed, ncol = 4)
  rownames(table_e2_fixed) <- cov_e_fixed
  colnames(table_e2_fixed) <- c("mean", "sd", "lower", "upper")
  table_e2_fixed[, 1] <- c(mean_cov_e_fixed_arm2)
  table_e2_fixed[, 2] <- c(sd_cov_e_fixed_arm2)
  table_e2_fixed[, 3] <- c(ql_cov_e_fixed_arm2)
  table_e2_fixed[, 4] <- c(qu_cov_e_fixed_arm2)
  table_e2_fixed <- round(table_e2_fixed, digits = digits) 
  table_c1_fixed <- matrix(NA, nrow = p_c_fixed, ncol = 4)
  rownames(table_c1_fixed) <- cov_c_fixed
  colnames(table_c1_fixed) <- c("mean", "sd", "lower", "upper")
  table_c1_fixed[, 1] <- c(mean_cov_c_fixed_arm1)
  table_c1_fixed[, 2] <- c(sd_cov_c_fixed_arm1)
  table_c1_fixed[, 3] <- c(ql_cov_c_fixed_arm1)
  table_c1_fixed[, 4] <- c(qu_cov_c_fixed_arm1)
  table_c1_fixed <- round(table_c1_fixed, digits = digits)
  table_c2_fixed <- matrix(NA, nrow = p_c_fixed, ncol = 4)
  rownames(table_c2_fixed) <- cov_c_fixed
  colnames(table_c2_fixed) <- c("mean", "sd", "lower", "upper")
  table_c2_fixed[, 1] <- c(mean_cov_c_fixed_arm2)
  table_c2_fixed[, 2] <- c(sd_cov_c_fixed_arm2)
  table_c2_fixed[, 3] <- c(ql_cov_c_fixed_arm2)
  table_c2_fixed[, 4] <- c(qu_cov_c_fixed_arm2)
  table_c2_fixed <- round(table_c2_fixed, digits = digits)
  }
  table_list_fixed <- list("Comparator" = list("Effects" = table_e1_fixed, "Costs" = table_c1_fixed),
                           "Reference" = list("Effects" = table_e2_fixed, "Costs" = table_c2_fixed))
  table_list_random <- NULL
  table_e1_random <- table_c1_random <- NULL
  table_e2_random <- table_c2_random <- NULL
  table_list_random <- list("Comparator" = list("Effects" = table_e1_random, "Costs" = table_c1_random),
                           "Reference" = list("Effects" = table_e2_random, "Costs" = table_c2_random))
  clus_e_arm1 <- clus_e_arm2 <- NULL
  clus_c_arm1 <- clus_c_arm2 <- NULL
  if(p_c_random != 0 & object$model_output$ind_random == FALSE) {
    cov_c_random <- c(cov_c_random, dep_c_random)
    p_c_random <- length(cov_c_random)
    mean_cov_c1_random <- rbind(mean_cov_c1_random, mean_dep_c1_random)
    sd_cov_c1_random <- rbind(sd_cov_c1_random, sd_dep_c1_random)
    ql_cov_c1_random <- rbind(ql_cov_c1_random, ql_dep_c1_random)
    qu_cov_c1_random <- rbind(qu_cov_c1_random, qu_dep_c1_random)
    mean_cov_c2_random <- rbind(mean_cov_c2_random, mean_dep_c2_random)
    sd_cov_c2_random <- rbind(sd_cov_c2_random, sd_dep_c2_random)
    ql_cov_c2_random <- rbind(ql_cov_c2_random, ql_dep_c2_random)
    qu_cov_c2_random <- rbind(qu_cov_c2_random, qu_dep_c2_random)
  } else if(p_c_random == 0 & object$model_output$ind_random == FALSE) {
    cov_c_random <- dep_c_random
    p_c_random <- length(cov_c_random)
    mean_cov_c1_random <- mean_dep_c1_random
    sd_cov_c1_random <- sd_dep_c1_random
    ql_cov_c1_random <- ql_dep_c1_random
    qu_cov_c1_random <- qu_dep_c1_random
    mean_cov_c2_random <- mean_dep_c2_random
    sd_cov_c2_random <- sd_dep_c2_random
    ql_cov_c2_random <- ql_dep_c2_random
    qu_cov_c2_random <- qu_dep_c2_random
    pc_random <- 1
  }
  if(p_e_random != 0) {
    mean_cov_e_random_arm1 <- mean_cov_e1_random
    sd_cov_e_random_arm1 <- sd_cov_e1_random
    ql_cov_e_random_arm1 <- ql_cov_e1_random
    qu_cov_e_random_arm1 <- qu_cov_e1_random
    mean_cov_e_random_arm2 <- mean_cov_e2_random
    sd_cov_e_random_arm2 <- sd_cov_e2_random
    ql_cov_e_random_arm2 <- ql_cov_e2_random
    qu_cov_e_random_arm2 <- qu_cov_e2_random
    clus_e_arm1 <- max(as.numeric(object$data_set$clus_effects$Control))
    clus_e_arm2 <- max(as.numeric(object$data_set$clus_effects$Intervention))
    clus_e1_index <- sort(rep(seq(1:clus_e_arm1), p_e_random))
    clus_e2_index <- sort(rep(seq(1:clus_e_arm2), p_e_random))
    table_e1_random <- matrix(NA, nrow = length(clus_e1_index), ncol = 4)
    rownames(table_e1_random) <- paste(rep(cov_e_random, clus_e_arm1), clus_e1_index, sep = " ")
    colnames(table_e1_random) <- c("mean", "sd", "lower", "upper")
    table_e1_random[, 1] <- c(mean_cov_e_random_arm1)
    table_e1_random[, 2] <- c(sd_cov_e_random_arm1)
    table_e1_random[, 3] <- c(ql_cov_e_random_arm1)
    table_e1_random[, 4] <- c(qu_cov_e_random_arm1)
    table_e1_random <- round(table_e1_random, digits = digits)
    table_e2_random <- matrix(NA, nrow = length(clus_e2_index), ncol = 4)
    rownames(table_e2_random) <- paste(rep(cov_e_random, clus_e_arm2), clus_e2_index, sep = " ")
    colnames(table_e2_random) <- c("mean", "sd", "lower", "upper")
    table_e2_random[, 1] <- c(mean_cov_e_random_arm2)
    table_e2_random[, 2] <- c(sd_cov_e_random_arm2)
    table_e2_random[, 3] <- c(ql_cov_e_random_arm2)
    table_e2_random[, 4] <- c(qu_cov_e_random_arm2)
    table_e2_random <- round(table_e2_random, digits = digits)
    table_list_random$Comparator$Effects <- table_e1_random
    table_list_random$Reference$Effects <- table_e2_random
  }
  if(p_c_random != 0 | object$model_output$ind_random == FALSE) {
    mean_cov_c_random_arm1 <- mean_cov_c1_random
    sd_cov_c_random_arm1 <- sd_cov_c1_random
    ql_cov_c_random_arm1 <- ql_cov_c1_random
    qu_cov_c_random_arm1 <- qu_cov_c1_random
    mean_cov_c_random_arm2 <- mean_cov_c2_random
    sd_cov_c_random_arm2 <- sd_cov_c2_random
    ql_cov_c_random_arm2 <- ql_cov_c2_random
    qu_cov_c_random_arm2 <- qu_cov_c2_random
    clus_c_arm1 <- max(as.numeric(object$data_set$clus_costs$Control))
    clus_c_arm2 <- max(as.numeric(object$data_set$clus_costs$Intervention))
    clus_c1_index <- sort(rep(seq(1:clus_c_arm1), p_c_random))
    clus_c2_index <- sort(rep(seq(1:clus_c_arm2), p_c_random))
    table_c1_random <- matrix(NA, nrow = length(clus_c1_index), ncol = 4)
    rownames(table_c1_random) <- paste(rep(cov_c_random, clus_c_arm1), clus_c1_index, sep = " ")
    colnames(table_c1_random) <- c("mean", "sd", "lower", "upper")
    table_c1_random[, 1] <- c(mean_cov_c_random_arm1)
    table_c1_random[, 2] <- c(sd_cov_c_random_arm1)
    table_c1_random[, 3] <- c(ql_cov_c_random_arm1)
    table_c1_random[, 4] <- c(qu_cov_c_random_arm1)
    table_c1_random <- round(table_c1_random, digits = digits)
    table_c2_random <- matrix(NA, nrow = length(clus_c2_index), ncol = 4)
    rownames(table_c2_random) <- paste(rep(cov_c_random, clus_c_arm2), clus_c2_index, sep = " ")
    colnames(table_c2_random) <- c("mean", "sd", "lower", "upper")
    table_c2_random[, 1] <- c(mean_cov_c_random_arm2)
    table_c2_random[, 2] <- c(sd_cov_c_random_arm2)
    table_c2_random[, 3] <- c(ql_cov_c_random_arm2)
    table_c2_random[, 4] <- c(qu_cov_c_random_arm2)
    table_c2_random <- round(table_c2_random, digits = digits)
    table_list_random$Comparator$Costs <- table_c1_random
    table_list_random$Reference$Costs <- table_c2_random
  }
   if(random == TRUE & p_c_random == 0 & p_e_random == 0) {
     stop("No random effects estimates found")
   }
  }
  if(object$data_format == "long") {
    max_time <- dim(object$data_set$effects$Control)[2]
    time_dep <- object$time_dep
    if(time < 1 | time > max_time | is.numeric(time) == FALSE | isTRUE(all.equal(time, as.integer(time))) == FALSE) {
      stop("Time must be a numeric integer. Minimum time is 1 (baseline) and maximum time is the latest follow-up time in the data")
    }
    if(length(grep("^SELECTION", object$model_output$type)) == 1) {
      cov_u_fixed <- names(object$data_set$covariates_effects_fixed$Control)
      p_u_fixed <- length(cov_u_fixed)
      if(length(dim(object$model_output$covariate_parameter_effects_fixed)) == 3) {
        mean_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , time], 2, mean)
        sd_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , time], 2, sd)
        ql_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , time], 2, quantile, prob = prob[1])
        qu_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , time], 2, quantile, prob = prob[2])
      } else if(length(dim(object$model_output$covariate_parameter_effects_fixed)) == 4) {
        mean_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , , time], c(2, 3), mean)
        sd_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , , time], c(2, 3), sd)
        ql_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , , time], c(2, 3), quantile, prob = prob[1])
        qu_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed[, , , time], c(2, 3), quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$covariate_parameter_effects_fixed)) == 0) {
        if(length(dim(object$model_output$covariate_parameter_effects_fixed$alpha)) == 3 & object$model_output$ind_fixed == FALSE) {
          mean_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , time], 2, mean)
          sd_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , time], 2, sd)
          ql_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , time], 2, quantile, prob = prob[1])
          qu_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , time], 2, quantile, prob = prob[2])
        } else if(length(dim(object$model_output$covariate_parameter_effects_fixed$alpha)) == 4 & object$model_output$ind_fixed == FALSE) {
          mean_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , , time], c(2, 3), mean)
          sd_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , , time], c(2, 3), sd)
          ql_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , , time], c(2, 3), quantile, prob = prob[1])
          qu_cov_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha[, , , time], c(2, 3), quantile, prob = prob[2])
        }
      }
      if(object$model_output$ind_fixed == FALSE) {
        dep_dep_utime_u_fixed <- NULL 
        dep_dep_ctime_u_fixed <- NULL
        if(object$model_output$ind_time_fixed == FALSE) {
          if(time >= 2 & time_dep == "AR1") {
          dep_dep_utime_u_fixed <- "u(j-1)"
          mean_dep_utime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tu[, , time], 2, mean)
          sd_dep_utime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tu[, , time], 2, sd)
          ql_dep_utime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tu[, , time], 2, quantile, prob = prob[2])
          dep_dep_ctime_u_fixed <- "c(j-1)"
          mean_dep_ctime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tc[, , time], 2, mean)
          sd_dep_ctime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tc[, , time], 2, sd)
          ql_dep_ctime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_u_fixed <- apply(object$model_output$covariate_parameter_effects_fixed$alpha_tc[, , time], 2, quantile, prob = prob[2])
          }
        }
      } else {
        dep_dep_utime_u_fixed <- NULL 
        dep_dep_ctime_u_fixed <- NULL 
      }
      cov_c_fixed <- names(object$data_set$covariates_costs_fixed$Control)
      p_c_fixed <- length(cov_c_fixed)
      if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 3 & object$model_output$ind_fixed == TRUE) {
        mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , time], 2, mean)
        sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , time], 2, sd)
        ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , time], 2, quantile, prob = prob[1])
        qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , time], 2, quantile, prob = prob[2])
      } else if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 4 & object$model_output$ind_fixed == TRUE) {
        mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , , time], c(2, 3), mean)
        sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , , time], c(2, 3), sd)
        ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , , time], c(2, 3), quantile, prob = prob[1])
        qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed[, , , time], c(2, 3), quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$covariate_parameter_costs_fixed)) == 0) {
        if(length(dim(object$model_output$covariate_parameter_costs_fixed$beta)) == 3 & object$model_output$ind_fixed == FALSE) {
          mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , time], 2, mean)
          sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , time], 2, sd)
          ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , time], 2, quantile, prob = prob[1])
          qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , time], 2, quantile, prob = prob[2])
        } else if(length(dim(object$model_output$covariate_parameter_costs_fixed$beta)) == 4 & object$model_output$ind_fixed == FALSE) {
          mean_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , , time], c(2, 3), mean)
          sd_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , , time], c(2, 3), sd)
          ql_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , , time], c(2, 3), quantile, prob = prob[1])
          qu_cov_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta[, , , time], c(2, 3), quantile, prob = prob[2])
        }
      }
      if(object$model_output$ind_fixed == FALSE) {
        dep_c_fixed <- "u"
        dep_dep_utime_c_fixed <- NULL
        dep_dep_ctime_c_fixed <- NULL
        mean_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f[, , time], 2, mean)
        sd_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f[, , time], 2, sd)
        ql_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f[, , time], 2, quantile, prob = prob[1])
        qu_dep_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_f[, , time], 2, quantile, prob = prob[2])
        if(object$model_output$ind_time_fixed == FALSE) {
          if(time >= 2 & time_dep == "AR1") {
          dep_dep_utime_c_fixed <- "u(j-1)"
          mean_dep_utime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tu[, , time], 2, mean)
          sd_dep_utime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tu[, , time], 2, sd)
          ql_dep_utime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tu[, , time], 2, quantile, prob = prob[2])
          dep_dep_ctime_c_fixed <- "c(j-1)"
          mean_dep_ctime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tc[, , time], 2, mean)
          sd_dep_ctime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tc[, , time], 2, sd)
          ql_dep_ctime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_c_fixed <- apply(object$model_output$covariate_parameter_costs_fixed$beta_tc[, , time], 2, quantile, prob = prob[2])
          }
        }
      } else {
        dep_c_fixed <- NULL 
        dep_dep_utime_c_fixed <- NULL
        dep_dep_ctime_c_fixed <- NULL
      }
      cov_u_random <- names(object$data_set$covariates_effects_random$Control)
      p_u_random <- length(cov_u_random)
      if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 3) {
        mean_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, mean)
        sd_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, sd)
        ql_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, quantile, prob = prob[1])
        qu_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, quantile, prob = prob[2])
        mean_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, mean)
        sd_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, sd)
        ql_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, quantile, prob = prob[1])
        qu_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, quantile, prob = prob[2])
      } else if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 4) {
        mean_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3), mean)
        sd_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3), sd)
        ql_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3), quantile, prob = prob[1])
        qu_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3), quantile, prob = prob[2])
        mean_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3), mean)
        sd_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3), sd)
        ql_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3), quantile, prob = prob[1])
        qu_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3), quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$covariate_parameter_effects_random)) == 0 & is.null(object$model_output$covariate_parameter_effects_random) == FALSE) {
        if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 3) {
          mean_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, mean)
          sd_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, sd)
          ql_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, quantile, prob = prob[1])
          qu_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , time], 2, quantile, prob = prob[2])
          mean_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, mean)
          sd_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, sd)
          ql_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, quantile, prob = prob[1])
          qu_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , time], 2, quantile, prob = prob[2])
        } else if(length(dim(object$model_output$covariate_parameter_effects_random$a1)) == 4) {
          mean_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3) , mean)
          sd_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3) , sd)
          ql_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3) , quantile, prob = prob[1])
          qu_cov_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1[, , , time], c(2, 3) , quantile, prob = prob[2])
          mean_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3) , mean)
          sd_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3) , sd)
          ql_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3) , quantile, prob = prob[1])
          qu_cov_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2[, , , time], c(2, 3) , quantile, prob = prob[2])
        }
      }
      if(object$model_output$ind_random == FALSE) {
        dep_dep_utime_u_random <- NULL 
        dep_dep_ctime_u_random <- NULL
        if(is.null(object$model_output$covariate_parameter_effects_random$a1_tu) == FALSE) {
          if(time >= 2 & time_dep == "AR1") {
          dep_dep_utime_u_random <- "u(j-1)"
          mean_dep_utime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tu[, , time], 2, mean)
          sd_dep_utime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tu[, , time], 2, sd)
          ql_dep_utime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tu[, , time], 2, quantile, prob = prob[2])
          mean_dep_utime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tu[, , time], 2, mean)
          sd_dep_utime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tu[, , time], 2, sd)
          ql_dep_utime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tu[, , time], 2, quantile, prob = prob[2])
          dep_dep_ctime_u_random <- "c(j-1)"
          mean_dep_ctime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tc[, , time], 2, mean)
          sd_dep_ctime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tc[, , time], 2, sd)
          ql_dep_ctime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_u1_random <- apply(object$model_output$covariate_parameter_effects_random$a1_tc[, , time], 2, quantile, prob = prob[2])
          mean_dep_ctime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tc[, , time], 2, mean)
          sd_dep_ctime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tc[, , time], 2, sd)
          ql_dep_ctime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_u2_random <- apply(object$model_output$covariate_parameter_effects_random$a2_tc[, , time], 2, quantile, prob = prob[2])
          }
        } else if(is.null(object$model_output$covariate_parameter_costs_random$a1_tu) == TRUE) {
          dep_dep_utime_u_random <- NULL
          dep_dep_ctime_u_random <- NULL
        }
      } else {
        dep_dep_utime_u_random <- NULL 
        dep_dep_ctime_u_random <- NULL 
      }
      cov_c_random <- names(object$data_set$covariates_costs_random$Control)
      p_c_random <- length(cov_c_random)
      if(length(dim(object$model_output$covariate_parameter_costs_random)) == 3) {
        mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , time], 2, mean)
        sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , time], 2, sd)
        ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , time], 2, quantile, prob = prob[1])
        qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , time], 2, quantile, prob = prob[2])
      } else if(length(dim(object$model_output$covariate_parameter_costs_random)) == 4) {
        mean_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , , time], c(2, 3), mean)
        sd_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , , time], c(2, 3), sd)
        ql_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , , time], c(2, 3), quantile, prob = prob[1])
        qu_cov_c_random <- apply(object$model_output$covariate_parameter_costs_random[, , , time], c(2, 3), quantile, prob = prob[2])
      }
      if(length(dim(object$model_output$covariate_parameter_costs_random)) == 0 & is.null(object$model_output$covariate_parameter_costs_random) == FALSE) {
        if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 3) {
          mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , time], 2, mean)
          sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , time], 2, sd)
          ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , time], 2, quantile, prob = prob[1])
          qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , time], 2, quantile, prob = prob[2])
          mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , time], 2, mean)
          sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , time], 2, sd)
          ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , time], 2, quantile, prob = prob[1])
          qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , time], 2, quantile, prob = prob[2])
        } else if(length(dim(object$model_output$covariate_parameter_costs_random$b1)) == 4) {
          mean_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , , time], c(2, 3), mean)
          sd_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , , time], c(2, 3), sd)
          ql_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , , time], c(2, 3), quantile, prob = prob[1])
          qu_cov_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1[, , , time], c(2, 3), quantile, prob = prob[2])
          mean_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , , time], c(2, 3), mean)
          sd_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , , time], c(2, 3), sd)
          ql_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , , time], c(2, 3), quantile, prob = prob[1])
          qu_cov_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2[, , , time], c(2, 3), quantile, prob = prob[2])
        }
      }
      if(object$model_output$ind_random == FALSE) {
        dep_c_random <- "u"
        mean_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f[, , time], 2, mean)
        sd_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f[, , time], 2, sd)
        ql_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f[, , time], 2, quantile, prob = prob[1])
        qu_dep_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_f[, , time], 2, quantile, prob = prob[2])
        mean_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f[, , time], 2, mean)
        sd_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f[, , time], 2, sd)
        ql_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f[, , time], 2, quantile, prob = prob[1])
        qu_dep_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_f[, , time], 2, quantile, prob = prob[2])
        dep_dep_utime_c_random <- NULL
        dep_dep_ctime_c_random <- NULL
        if(object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == FALSE) {
          if(time >= 2) {
          dep_dep_utime_c_random <- "u(j-1)"
          mean_dep_utime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tu[, , time], 2, mean)
          sd_dep_utime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tu[, , time], 2, sd)
          ql_dep_utime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tu[, , time], 2, quantile, prob = prob[2])
          mean_dep_utime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tu[, , time], 2, mean)
          sd_dep_utime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tu[, , time], 2, sd)
          ql_dep_utime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tu[, , time], 2, quantile, prob = prob[1])
          qu_dep_utime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tu[, , time], 2, quantile, prob = prob[2])
          dep_dep_ctime_c_random <- "c(j-1)"
          mean_dep_ctime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tc[, , time], 2, mean)
          sd_dep_ctime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tc[, , time], 2, sd)
          ql_dep_ctime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_c1_random <- apply(object$model_output$covariate_parameter_costs_random$b1_tc[, , time], 2, quantile, prob = prob[2])
          mean_dep_ctime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tc[, , time], 2, mean)
          sd_dep_ctime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tc[, , time], 2, sd)
          ql_dep_ctime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tc[, , time], 2, quantile, prob = prob[1])
          qu_dep_ctime_c2_random <- apply(object$model_output$covariate_parameter_costs_random$b2_tc[, , time], 2, quantile, prob = prob[2])
          }
        } else if(object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == TRUE) {
          dep_dep_utime_c_random <- NULL
          dep_dep_ctime_c_random <- NULL
        }
      } else {
        dep_c_random <- NULL
        dep_dep_utime_c_random <- NULL
        dep_dep_ctime_c_random <- NULL
      }
      if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == TRUE & p_c_fixed != 0) {
        cov_c_fixed <- c(cov_c_fixed, dep_c_fixed)
        p_c_fixed <- length(cov_c_fixed)
        mean_cov_c_fixed <- rbind(mean_cov_c_fixed, mean_dep_c_fixed)
        sd_cov_c_fixed <- rbind(sd_cov_c_fixed, sd_dep_c_fixed)
        ql_cov_c_fixed <- rbind(ql_cov_c_fixed, ql_dep_c_fixed)
        qu_cov_c_fixed <- rbind(qu_cov_c_fixed, qu_dep_c_fixed)
      } else if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == TRUE & p_c_fixed == 0) {
        cov_c_fixed <- c(dep_c_fixed)
        p_c_fixed <- length(cov_c_fixed)
        mean_cov_c_fixed <- rbind(mean_dep_c_fixed)
        sd_cov_c_fixed <- rbind(sd_dep_c_fixed)
        ql_cov_c_fixed <- rbind(ql_dep_c_fixed)
        qu_cov_c_fixed <- rbind(qu_dep_c_fixed)
      } else if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == FALSE & p_c_fixed != 0) {
        if(time == 1) {
          cov_c_fixed <- c(cov_c_fixed, dep_c_fixed)
          p_c_fixed <- length(cov_c_fixed)
          mean_cov_c_fixed <- rbind(mean_cov_c_fixed, mean_dep_c_fixed)
          sd_cov_c_fixed <- rbind(sd_cov_c_fixed, sd_dep_c_fixed)
          ql_cov_c_fixed <- rbind(ql_cov_c_fixed, ql_dep_c_fixed)
          qu_cov_c_fixed <- rbind(qu_cov_c_fixed, qu_dep_c_fixed)
        } else if(time > 1) {
          if(time_dep == "AR1") {
          cov_c_fixed <- c(cov_c_fixed, dep_c_fixed, dep_dep_utime_c_fixed, dep_dep_ctime_c_fixed)
          p_c_fixed <- length(cov_c_fixed)
          mean_cov_c_fixed <- rbind(mean_cov_c_fixed, mean_dep_c_fixed, mean_dep_utime_c_fixed, mean_dep_ctime_c_fixed)
          sd_cov_c_fixed <- rbind(sd_cov_c_fixed, sd_dep_c_fixed, sd_dep_utime_c_fixed, sd_dep_ctime_c_fixed)
          ql_cov_c_fixed <- rbind(ql_cov_c_fixed, ql_dep_c_fixed, ql_dep_utime_c_fixed, ql_dep_ctime_c_fixed)
          qu_cov_c_fixed <- rbind(qu_cov_c_fixed, qu_dep_c_fixed, qu_dep_utime_c_fixed, qu_dep_ctime_c_fixed)
          } else {
            cov_c_fixed <- c(cov_c_fixed, dep_c_fixed)
            p_c_fixed <- length(cov_c_fixed)
            mean_cov_c_fixed <- rbind(mean_cov_c_fixed, mean_dep_c_fixed)
            sd_cov_c_fixed <- rbind(sd_cov_c_fixed, sd_dep_c_fixed)
            ql_cov_c_fixed <- rbind(ql_cov_c_fixed, ql_dep_c_fixed)
            qu_cov_c_fixed <- rbind(qu_cov_c_fixed, qu_dep_c_fixed)
          }
        }
      } else if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == FALSE & p_c_fixed == 0) {
        if(time == 1) {
          cov_c_fixed <- c(dep_c_fixed)
          p_c_fixed <- length(cov_c_fixed)
          mean_cov_c_fixed <- rbind(mean_dep_c_fixed)
          sd_cov_c_fixed <- rbind(sd_dep_c_fixed)
          ql_cov_c_fixed <- rbind(ql_dep_c_fixed)
          qu_cov_c_fixed <- rbind(qu_dep_c_fixed)
        } else if(time > 1) {
          if(time_dep == "AR1") {
            cov_c_fixed <- c(dep_c_fixed, dep_dep_utime_c_fixed, dep_dep_ctime_c_fixed)
            p_c_fixed <- length(cov_c_fixed)
            mean_cov_c_fixed <- rbind(mean_dep_c_fixed, mean_dep_utime_c_fixed, mean_dep_ctime_c_fixed)
            sd_cov_c_fixed <- rbind(sd_dep_c_fixed, sd_dep_utime_c_fixed, sd_dep_ctime_c_fixed)
            ql_cov_c_fixed <- rbind(ql_dep_c_fixed, ql_dep_utime_c_fixed, ql_dep_ctime_c_fixed)
            qu_cov_c_fixed <- rbind(qu_dep_c_fixed, qu_dep_utime_c_fixed, qu_dep_ctime_c_fixed)
          } else {
            cov_c_fixed <- c(dep_c_fixed)
            p_c_fixed <- length(cov_c_fixed)
            mean_cov_c_fixed <- rbind(mean_dep_c_fixed)
            sd_cov_c_fixed <- rbind(sd_dep_c_fixed)
            ql_cov_c_fixed <- rbind(ql_dep_c_fixed)
            qu_cov_c_fixed <- rbind(qu_dep_c_fixed)
          }
        }
      } 
      if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == FALSE) {
        if(time == 1) {
          cov_u_fixed <- cov_u_fixed
          p_u_fixed <- length(cov_u_fixed)
          mean_cov_u_fixed <- mean_cov_u_fixed
          sd_cov_u_fixed <- sd_cov_u_fixed
          ql_cov_u_fixed <- ql_cov_u_fixed
          qu_cov_u_fixed <- qu_cov_u_fixed
        } else if(time > 1) {
          if(time_dep == "AR1") {
            cov_u_fixed <- c(cov_u_fixed, dep_dep_utime_u_fixed, dep_dep_ctime_u_fixed)
            p_u_fixed <- length(cov_u_fixed)
            mean_cov_u_fixed <- rbind(mean_cov_u_fixed, mean_dep_utime_u_fixed, mean_dep_ctime_u_fixed)
            sd_cov_u_fixed <- rbind(sd_cov_u_fixed, sd_dep_utime_u_fixed, sd_dep_ctime_u_fixed)
            ql_cov_u_fixed <- rbind(ql_cov_u_fixed, ql_dep_utime_u_fixed, ql_dep_ctime_u_fixed)
            qu_cov_u_fixed <- rbind(qu_cov_u_fixed, qu_dep_utime_u_fixed, qu_dep_ctime_u_fixed)
          } else {
            cov_u_fixed <- cov_u_fixed
            p_u_fixed <- length(cov_u_fixed)
            mean_cov_u_fixed <- mean_cov_u_fixed
            sd_cov_u_fixed <- sd_cov_u_fixed
            ql_cov_u_fixed <- ql_cov_u_fixed
            qu_cov_u_fixed <- qu_cov_u_fixed
          }
        }
      }
      if(object$model_output$ind_fixed == TRUE & p_u_fixed == 1) {
        mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[1]
        sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[1]
        ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[1]
        qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[1]
        mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[2]
        sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[2]
        ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[2]
        qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[2]
      } else if(object$model_output$ind_fixed == TRUE & p_u_fixed > 1) {
        mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[, 1]
        sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[, 1]
        ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[, 1]
        qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[, 1]
        mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[, 2]
        sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[, 2]
        ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[, 2]
        qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[, 2]
      } 
      if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == TRUE & p_u_fixed == 1) {
        mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[1]
        sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[1]
        ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[1]
        qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[1]
        mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[2]
        sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[2]
        ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[2]
        qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[2]
      } else if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == TRUE & p_u_fixed > 1) {
        mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[, 1]
        sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[, 1]
        ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[, 1]
        qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[, 1]
        mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[, 2]
        sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[, 2]
        ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[, 2]
        qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[, 2]
      }
      if(object$model_output$ind_fixed == FALSE & object$model_output$ind_time_fixed == FALSE) {
        if(time == 1) {
        mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[1]
        sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[1]
        ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[1]
        qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[1]
        mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[2]
        sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[2]
        ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[2]
        qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[2]
        } else if(time > 1) {
          if(time_dep == "AR1") {
            mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[, 1]
            sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[, 1]
            ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[, 1]
            qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[, 1]
            mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[, 2]
            sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[, 2]
            ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[, 2]
            qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[, 2]
          } else {
            mean_cov_u_fixed_arm1 <- mean_cov_u_fixed[1]
            sd_cov_u_fixed_arm1 <- sd_cov_u_fixed[1]
            ql_cov_u_fixed_arm1 <- ql_cov_u_fixed[1]
            qu_cov_u_fixed_arm1 <- qu_cov_u_fixed[1]
            mean_cov_u_fixed_arm2 <- mean_cov_u_fixed[2]
            sd_cov_u_fixed_arm2 <- sd_cov_u_fixed[2]
            ql_cov_u_fixed_arm2 <- ql_cov_u_fixed[2]
            qu_cov_u_fixed_arm2 <- qu_cov_u_fixed[2]
          }
        }
      }
      if(object$model_output$ind_fixed == TRUE & p_c_fixed == 1) {
        mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[1]
        sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[1]
        ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[1]
        qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[1]
        mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[2]
        sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[2]
        ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[2]
        qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[2]
      } else if(object$model_output$ind_fixed == TRUE & p_c_fixed > 1) {
        mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[, 1]
        sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[, 1]
        ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[, 1]
        qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[, 1]
        mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[, 2]
        sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[, 2]
        ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[, 2]
        qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[, 2]
      } 
      if(object$model_output$ind_fixed == FALSE & p_c_fixed != 0) {
        mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[, 1]
        sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[, 1]
        ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[, 1]
        qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[, 1]
        mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[, 2]
        sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[, 2]
        ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[, 2]
        qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[, 2]
      } else if(object$model_output$ind_fixed == FALSE & p_c_fixed == 0) {
        mean_cov_c_fixed_arm1 <- mean_cov_c_fixed[1]
        sd_cov_c_fixed_arm1 <- sd_cov_c_fixed[1]
        ql_cov_c_fixed_arm1 <- ql_cov_c_fixed[1]
        qu_cov_c_fixed_arm1 <- qu_cov_c_fixed[1]
        mean_cov_c_fixed_arm2 <- mean_cov_c_fixed[2]
        sd_cov_c_fixed_arm2 <- sd_cov_c_fixed[2]
        ql_cov_c_fixed_arm2 <- ql_cov_c_fixed[2]
        qu_cov_c_fixed_arm2 <- qu_cov_c_fixed[2]
      }
    }
    if(length(grep("^SELECTION", object$model_output$type)) == 1) {
      table_u1_fixed <- matrix(NA, nrow = p_u_fixed, ncol = 4)
      rownames(table_u1_fixed) <- cov_u_fixed
      colnames(table_u1_fixed) <- c("mean", "sd", "lower", "upper")
      table_u1_fixed[, 1] <- c(mean_cov_u_fixed_arm1)
      table_u1_fixed[, 2] <- c(sd_cov_u_fixed_arm1)
      table_u1_fixed[, 3] <- c(ql_cov_u_fixed_arm1)
      table_u1_fixed[, 4] <- c(qu_cov_u_fixed_arm1)
      table_u1_fixed <- round(table_u1_fixed, digits = digits)
      table_u2_fixed <- matrix(NA, nrow = p_u_fixed, ncol = 4)
      rownames(table_u2_fixed) <- cov_u_fixed
      colnames(table_u2_fixed) <- c("mean", "sd", "lower", "upper")
      table_u2_fixed[, 1] <- c(mean_cov_u_fixed_arm2)
      table_u2_fixed[, 2] <- c(sd_cov_u_fixed_arm2)
      table_u2_fixed[, 3] <- c(ql_cov_u_fixed_arm2)
      table_u2_fixed[, 4] <- c(qu_cov_u_fixed_arm2)
      table_u2_fixed <- round(table_u2_fixed, digits = digits) 
      table_c1_fixed <- matrix(NA, nrow = p_c_fixed, ncol = 4)
      rownames(table_c1_fixed) <- cov_c_fixed
      colnames(table_c1_fixed) <- c("mean", "sd", "lower", "upper")
      table_c1_fixed[, 1] <- c(mean_cov_c_fixed_arm1)
      table_c1_fixed[, 2] <- c(sd_cov_c_fixed_arm1)
      table_c1_fixed[, 3] <- c(ql_cov_c_fixed_arm1)
      table_c1_fixed[, 4] <- c(qu_cov_c_fixed_arm1)
      table_c1_fixed <- round(table_c1_fixed, digits = digits)
      table_c2_fixed <- matrix(NA, nrow = p_c_fixed, ncol = 4)
      rownames(table_c2_fixed) <- cov_c_fixed
      colnames(table_c2_fixed) <- c("mean", "sd", "lower", "upper")
      table_c2_fixed[, 1] <- c(mean_cov_c_fixed_arm2)
      table_c2_fixed[, 2] <- c(sd_cov_c_fixed_arm2)
      table_c2_fixed[, 3] <- c(ql_cov_c_fixed_arm2)
      table_c2_fixed[, 4] <- c(qu_cov_c_fixed_arm2)
      table_c2_fixed <- round(table_c2_fixed, digits = digits)
    }
    table_list_fixed <- list("Comparator" = list("Effects" = table_u1_fixed, "Costs" = table_c1_fixed),
                             "Reference" = list("Effects" = table_u2_fixed, "Costs" = table_c2_fixed))
    table_list_random <- NULL
    table_u1_random <- table_c1_random <- NULL
    table_u2_random <- table_c2_random <- NULL
    table_list_random <- list("Comparator" = list("Effects" = table_u1_random, "Costs" = table_c1_random),
                              "Reference" = list("Effects" = table_u2_random, "Costs" = table_c2_random))
    clus_u_arm1 <- clus_u_arm2 <- NULL
    clus_c_arm1 <- clus_c_arm2 <- NULL
    if(p_c_random != 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == TRUE) {
      cov_c_random <- c(cov_c_random, dep_c_random)
      p_c_random <- length(cov_c_random)
      mean_cov_c1_random <- rbind(mean_cov_c1_random, mean_dep_c1_random)
      sd_cov_c1_random <- rbind(sd_cov_c1_random, sd_dep_c1_random)
      ql_cov_c1_random <- rbind(ql_cov_c1_random, ql_dep_c1_random)
      qu_cov_c1_random <- rbind(qu_cov_c1_random, qu_dep_c1_random)
      mean_cov_c2_random <- rbind(mean_cov_c2_random, mean_dep_c2_random)
      sd_cov_c2_random <- rbind(sd_cov_c2_random, sd_dep_c2_random)
      ql_cov_c2_random <- rbind(ql_cov_c2_random, ql_dep_c2_random)
      qu_cov_c2_random <- rbind(qu_cov_c2_random, qu_dep_c2_random)
    } else if(p_c_random != 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == FALSE) {
      if(time == 1) {
        cov_c_random <- c(cov_c_random, dep_c_random)
        p_c_random <- length(cov_c_random)
        mean_cov_c1_random <- rbind(mean_cov_c1_random, mean_dep_c1_random)
        sd_cov_c1_random <- rbind(sd_cov_c1_random, sd_dep_c1_random)
        ql_cov_c1_random <- rbind(ql_cov_c1_random, ql_dep_c1_random)
        qu_cov_c1_random <- rbind(qu_cov_c1_random, qu_dep_c1_random)
        mean_cov_c2_random <- rbind(mean_cov_c2_random, mean_dep_c2_random)
        sd_cov_c2_random <- rbind(sd_cov_c2_random, sd_dep_c2_random)
        ql_cov_c2_random <- rbind(ql_cov_c2_random, ql_dep_c2_random)
        qu_cov_c2_random <- rbind(qu_cov_c2_random, qu_dep_c2_random)
      } else if(time > 1) {
       cov_c_random <- c(cov_c_random, dep_c_random, dep_dep_utime_c_random, dep_dep_ctime_c_random)
       p_c_random <- length(cov_c_random)
       mean_cov_c1_random <- rbind(mean_cov_c1_random, mean_dep_c1_random, mean_dep_utime_c1_random, mean_dep_ctime_c1_random)
       sd_cov_c1_random <- rbind(sd_cov_c1_random, sd_dep_c1_random, sd_dep_utime_c1_random, sd_dep_ctime_c1_random)
       ql_cov_c1_random <- rbind(ql_cov_c1_random, ql_dep_c1_random, ql_dep_utime_c1_random, ql_dep_ctime_c1_random)
       qu_cov_c1_random <- rbind(qu_cov_c1_random, qu_dep_c1_random, qu_dep_utime_c1_random, qu_dep_ctime_c1_random)
       mean_cov_c2_random <- rbind(mean_cov_c2_random, mean_dep_c2_random, mean_dep_utime_c2_random, mean_dep_ctime_c2_random)
       sd_cov_c2_random <- rbind(sd_cov_c2_random, sd_dep_c2_random, sd_dep_utime_c2_random, sd_dep_ctime_c2_random)
       ql_cov_c2_random <- rbind(ql_cov_c2_random, ql_dep_c2_random, ql_dep_utime_c2_random, ql_dep_ctime_c2_random)
       qu_cov_c2_random <- rbind(qu_cov_c2_random, qu_dep_c2_random, qu_dep_utime_c2_random, qu_dep_ctime_c2_random)
       }
     } else if(p_c_random == 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == TRUE) {
      cov_c_random <- dep_c_random
      p_c_random <- length(cov_c_random)
      mean_cov_c1_random <- mean_dep_c1_random
      sd_cov_c1_random <- sd_dep_c1_random
      ql_cov_c1_random <- ql_dep_c1_random
      qu_cov_c1_random <- qu_dep_c1_random
      mean_cov_c2_random <- mean_dep_c2_random
      sd_cov_c2_random <- sd_dep_c2_random
      ql_cov_c2_random <- ql_dep_c2_random
      qu_cov_c2_random <- qu_dep_c2_random
      pc_random <- 1
     } else if(p_c_random == 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$b1_tu) == FALSE) {
       if(is.null(dep_dep_utime_c_random) == FALSE) {
       if(time == 1) {
         cov_c_random <- dep_c_random
         p_c_random <- length(cov_c_random)
         mean_cov_c1_random <- mean_dep_c1_random
         sd_cov_c1_random <- sd_dep_c1_random
         ql_cov_c1_random <- ql_dep_c1_random
         qu_cov_c1_random <- qu_dep_c1_random
         mean_cov_c2_random <- mean_dep_c2_random
         sd_cov_c2_random <- sd_dep_c2_random
         ql_cov_c2_random <- ql_dep_c2_random
         qu_cov_c2_random <- qu_dep_c2_random
       } else if(time > 1) {
         cov_c_random <- c(dep_c_random, dep_dep_utime_c_random, dep_dep_ctime_c_random)
         p_c_random <- length(cov_c_random)
         mean_cov_c1_random <- rbind(mean_dep_c1_random, mean_dep_utime_c1_random, mean_dep_ctime_c1_random)
         sd_cov_c1_random <- rbind(sd_dep_c1_random, sd_dep_utime_c1_random, sd_dep_ctime_c1_random)
         ql_cov_c1_random <- rbind(ql_dep_c1_random, ql_dep_utime_c1_random, ql_dep_ctime_c1_random)
         qu_cov_c1_random <- rbind(qu_dep_c1_random, qu_dep_utime_c1_random, qu_dep_ctime_c1_random)
         mean_cov_c2_random <- rbind(mean_dep_c2_random, mean_dep_utime_c2_random, mean_dep_ctime_c2_random)
         sd_cov_c2_random <- rbind(sd_dep_c2_random, sd_dep_utime_c2_random, sd_dep_ctime_c2_random)
         ql_cov_c2_random <- rbind(ql_dep_c2_random, ql_dep_utime_c2_random, ql_dep_ctime_c2_random)
         qu_cov_c2_random <- rbind(qu_dep_c2_random, qu_dep_utime_c2_random, qu_dep_ctime_c2_random)
       }
       } else if(is.null(dep_dep_utime_c_random) == TRUE) {
         if(time == 1) {
           cov_c_random <- dep_c_random
           p_c_random <- length(cov_c_random)
           mean_cov_c1_random <- mean_dep_c1_random
           sd_cov_c1_random <- sd_dep_c1_random
           ql_cov_c1_random <- ql_dep_c1_random
           qu_cov_c1_random <- qu_dep_c1_random
           mean_cov_c2_random <- mean_dep_c2_random
           sd_cov_c2_random <- sd_dep_c2_random
           ql_cov_c2_random <- ql_dep_c2_random
           qu_cov_c2_random <- qu_dep_c2_random
         } else if(time > 1) {
           cov_c_random <- dep_c_random
           p_c_random <- length(cov_c_random)
           mean_cov_c1_random <- mean_dep_c1_random
           sd_cov_c1_random <- sd_dep_c1_random
           ql_cov_c1_random <- ql_dep_c1_random
           qu_cov_c1_random <- qu_dep_c1_random
           mean_cov_c2_random <- mean_dep_c2_random
           sd_cov_c2_random <- sd_dep_c2_random
           ql_cov_c2_random <- ql_dep_c2_random
           qu_cov_c2_random <- qu_dep_c2_random
         }
      }
     }
     if(p_u_random != 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$a1_tu) == FALSE) {
      if(time > 1 & time_dep == "AR1") {
        cov_u_random <- c(cov_u_random, dep_dep_utime_u_random, dep_dep_ctime_u_random)
        p_u_random <- length(cov_u_random)
        mean_cov_u1_random <- rbind(mean_cov_u1_random, mean_dep_utime_u1_random, mean_dep_ctime_u1_random)
        sd_cov_u1_random <- rbind(sd_cov_u1_random, sd_dep_utime_u1_random, sd_dep_ctime_u1_random)
        ql_cov_u1_random <- rbind(ql_cov_u1_random, ql_dep_utime_u1_random, ql_dep_ctime_u1_random)
        qu_cov_u1_random <- rbind(qu_cov_u1_random, qu_dep_utime_u1_random, qu_dep_ctime_u1_random)
        mean_cov_u2_random <- rbind(mean_cov_u2_random, mean_dep_utime_u2_random, mean_dep_ctime_u2_random)
        sd_cov_u2_random <- rbind(sd_cov_u2_random, sd_dep_utime_u2_random, sd_dep_ctime_u2_random)
        ql_cov_u2_random <- rbind(ql_cov_u2_random, ql_dep_utime_u2_random, ql_dep_ctime_u2_random)
        qu_cov_u2_random <- rbind(qu_cov_u2_random, qu_dep_utime_u2_random, qu_dep_ctime_u2_random)
      }
     } else if (p_u_random != 0 & object$model_output$ind_random == FALSE & is.null(object$model_output$covariate_parameter_costs_random$a1_tu) == TRUE) {
       if(time > 1 & time_dep == "AR1") {
         cov_u_random <- c(cov_u_random, dep_dep_utime_u_random, dep_dep_ctime_u_random)
         p_u_random <- length(cov_u_random)
         mean_cov_u1_random <- rbind(mean_cov_u1_random, mean_dep_utime_u1_random, mean_dep_ctime_u1_random)
         sd_cov_u1_random <- rbind(sd_cov_u1_random, sd_dep_utime_u1_random, sd_dep_ctime_u1_random)
         ql_cov_u1_random <- rbind(ql_cov_u1_random, ql_dep_utime_u1_random, ql_dep_ctime_u1_random)
         qu_cov_u1_random <- rbind(qu_cov_u1_random, qu_dep_utime_u1_random, qu_dep_ctime_u1_random)
         mean_cov_u2_random <- rbind(mean_cov_u2_random, mean_dep_utime_u2_random, mean_dep_ctime_u2_random)
         sd_cov_u2_random <- rbind(sd_cov_u2_random, sd_dep_utime_u2_random, sd_dep_ctime_u2_random)
         ql_cov_u2_random <- rbind(ql_cov_u2_random, ql_dep_utime_u2_random, ql_dep_ctime_u2_random)
         qu_cov_u2_random <- rbind(qu_cov_u2_random, qu_dep_utime_u2_random, qu_dep_ctime_u2_random)
       } 
    }
    if(p_u_random != 0) {
      mean_cov_u_random_arm1 <- mean_cov_u1_random
      sd_cov_u_random_arm1 <- sd_cov_u1_random
      ql_cov_u_random_arm1 <- ql_cov_u1_random
      qu_cov_u_random_arm1 <- qu_cov_u1_random
      mean_cov_u_random_arm2 <- mean_cov_u2_random
      sd_cov_u_random_arm2 <- sd_cov_u2_random
      ql_cov_u_random_arm2 <- ql_cov_u2_random
      qu_cov_u_random_arm2 <- qu_cov_u2_random
      clus_u_arm1 <- max(as.numeric(object$data_set$clus_effects$Control))
      clus_u_arm2 <- max(as.numeric(object$data_set$clus_effects$Intervention))
      clus_u1_index <- sort(rep(seq(1:clus_u_arm1), p_u_random))
      clus_u2_index <- sort(rep(seq(1:clus_u_arm2), p_u_random))
      table_u1_random <- matrix(NA, nrow = length(clus_u1_index), ncol = 4)
      rownames(table_u1_random) <- paste(rep(cov_u_random, clus_u_arm1), clus_u1_index, sep = " ")
      colnames(table_u1_random) <- c("mean", "sd", "lower", "upper")
      table_u1_random[, 1] <- c(mean_cov_u_random_arm1)
      table_u1_random[, 2] <- c(sd_cov_u_random_arm1)
      table_u1_random[, 3] <- c(ql_cov_u_random_arm1)
      table_u1_random[, 4] <- c(qu_cov_u_random_arm1)
      table_u1_random <- round(table_u1_random, digits = digits)
      table_u2_random <- matrix(NA, nrow = length(clus_u2_index), ncol = 4)
      rownames(table_u2_random) <- paste(rep(cov_u_random, clus_u_arm2), clus_u2_index, sep = " ")
      colnames(table_u2_random) <- c("mean", "sd", "lower", "upper")
      table_u2_random[, 1] <- c(mean_cov_u_random_arm2)
      table_u2_random[, 2] <- c(sd_cov_u_random_arm2)
      table_u2_random[, 3] <- c(ql_cov_u_random_arm2)
      table_u2_random[, 4] <- c(qu_cov_u_random_arm2)
      table_u2_random <- round(table_u2_random, digits = digits)
      table_list_random$Comparator$Effects <- table_u1_random
      table_list_random$Reference$Effects <- table_u2_random
    }
    if(p_c_random != 0 | object$model_output$ind_random == FALSE) {
      mean_cov_c_random_arm1 <- mean_cov_c1_random
      sd_cov_c_random_arm1 <- sd_cov_c1_random
      ql_cov_c_random_arm1 <- ql_cov_c1_random
      qu_cov_c_random_arm1 <- qu_cov_c1_random
      mean_cov_c_random_arm2 <- mean_cov_c2_random
      sd_cov_c_random_arm2 <- sd_cov_c2_random
      ql_cov_c_random_arm2 <- ql_cov_c2_random
      qu_cov_c_random_arm2 <- qu_cov_c2_random
      clus_c_arm1 <- max(as.numeric(object$data_set$clus_costs$Control))
      clus_c_arm2 <- max(as.numeric(object$data_set$clus_costs$Intervention))
      clus_c1_index <- sort(rep(seq(1:clus_c_arm1), p_c_random))
      clus_c2_index <- sort(rep(seq(1:clus_c_arm2), p_c_random))
      table_c1_random <- matrix(NA, nrow = length(clus_c1_index), ncol = 4)
      rownames(table_c1_random) <- paste(rep(cov_c_random, clus_c_arm1), clus_c1_index, sep = " ")
      colnames(table_c1_random) <- c("mean", "sd", "lower", "upper")
      table_c1_random[, 1] <- c(mean_cov_c_random_arm1)
      table_c1_random[, 2] <- c(sd_cov_c_random_arm1)
      table_c1_random[, 3] <- c(ql_cov_c_random_arm1)
      table_c1_random[, 4] <- c(qu_cov_c_random_arm1)
      table_c1_random <- round(table_c1_random, digits = digits)
      table_c2_random <- matrix(NA, nrow = length(clus_c2_index), ncol = 4)
      rownames(table_c2_random) <- paste(rep(cov_c_random, clus_c_arm2), clus_c2_index, sep = " ")
      colnames(table_c2_random) <- c("mean", "sd", "lower", "upper")
      table_c2_random[, 1] <- c(mean_cov_c_random_arm2)
      table_c2_random[, 2] <- c(sd_cov_c_random_arm2)
      table_c2_random[, 3] <- c(ql_cov_c_random_arm2)
      table_c2_random[, 4] <- c(qu_cov_c_random_arm2)
      table_c2_random <- round(table_c2_random, digits = digits)
      table_list_random$Comparator$Costs <- table_c1_random
      table_list_random$Reference$Costs <- table_c2_random
    }
    if(random == TRUE & p_c_random == 0 & p_u_random == 0) {
      stop("No random effects estimates found for effect and cost models")
    }
  }   
  if(random == FALSE) {
  print(table_list_fixed)  
  } else if(random == TRUE) {
  print(table_list_random)  
  }
}