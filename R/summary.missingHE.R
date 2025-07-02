#' Summary method for objects in the class \code{missingHE}
#'
#' Produces a table printout with some summary results of the health economic evaluation probabilistic model
#' run using the function \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}.
#' @param object A \code{missingHE} object containing the results of the Bayesian modelling and the economic evaluation
#' @param ... Additional arguments affecting the summary produced.
#' @return Prints a table with some information on the health economic model based on the assumption
#' selected for the missingness using the function \code{selection}, \code{pattern} or \code{hurdle}. 
#' Summary information on the main parameters of interests is provided.
#' @seealso \code{\link{selection}} \code{\link{pattern}} \code{\link{hurdle}} \code{\link{diagnostic}} \code{\link{plot.missingHE}}
#' @author Andrea Gabrio
#' @references 
#' Baio, G.(2012). \emph{Bayesian Methods in Health Economcis}. CRC/Chapman Hall, London.
#' @importFrom stats quantile
#' @export
#' @examples 
#' #For examples see the function selection, pattern or hurdle
#' #
#' #

summary.missingHE <- function(object, ...) {
  exArgs <- list(...)
  if(!inherits(object, "missingHE")) {
    stop("Only objects of class 'missingHE' can be used")
  }
  mu_eff1 <- round(mean(object$cea$e[,1]), digits=3)
  mu_eff2 <- round(mean(object$cea$e[,2]), digits = 3)
  mu_cost1 <- round(mean(object$cea$c[,1]), digits = 3)
  mu_cost2 <- round(mean(object$cea$c[,2]), digits = 3)
  mu_delta_e <- round(mean(unlist(object$cea$delta.e)), digits = 3)
  mu_delta_c <- round(mean(unlist(object$cea$delta.c)), digits = 3)
  icer <- round(object$cea$ICER, digits = 3)
  sd_mu_eff1 <- round(sd(object$cea$e[,1]), digits = 3)
  sd_mu_eff2 <- round(sd(object$cea$e[,2]), digits = 3)
  sd_mu_cost1 <- round(sd(object$cea$c[,1]), digits = 3)
  sd_mu_cost2 <- round(sd(object$cea$c[,2]), digits = 3)
  sd_delta_e <- round(sd(unlist(object$cea$delta.e)), digits = 3)
  sd_delta_c <- round(sd(unlist(object$cea$delta.c)), digits = 3)
  lower_mu_eff1 <- round(quantile(object$model_output$mean_effects[,1], probs = 0.05), digits = 3)
  lower_mu_eff2 <- round(quantile(object$model_output$mean_effects[,2], probs = 0.05), digits = 3)
  lower_mu_cost1 <- round(quantile(object$model_output$mean_costs[,1], probs = 0.05), digits = 3)
  lower_mu_cost2 <- round(quantile(object$model_output$mean_costs[,2], probs = 0.05), digits = 3)
  lower_mu_delta_e <- round(quantile(object$cea$delta.e, probs = 0.05), digits = 3)
  lower_mu_delta_c <- round(quantile(object$cea$delta.c, probs = 0.05), digits = 3)
  upper_mu_eff1 <- round(quantile(object$model_output$mean_effects[,1], probs = 0.95), digits = 3)
  upper_mu_eff2 <- round(quantile(object$model_output$mean_effects[,2], probs = 0.95), digits = 3)
  upper_mu_cost1 <- round(quantile(object$model_output$mean_costs[,1], probs = 0.95), digits = 3)
  upper_mu_cost2 <- round(quantile(object$model_output$mean_costs[,2], probs = 0.95), digits = 3)
  upper_mu_delta_e <- round(quantile(unlist(object$cea$delta.e), probs = 0.95), digits = 3)
  upper_mu_delta_c <- round(quantile(unlist(object$cea$delta.c), probs = 0.95), digits = 3)
  rownames_v <- c("mean effects (t = 1)", "mean costs (t = 1)", "mean effects (t = 2)", "mean costs (t = 2)",
                "delta effects", "delta costs", "ICER")
  colnames_v <- c("mean", "sd", "LB", "UB")
  mean_v <- c(mu_eff1, mu_cost1, mu_eff2, mu_cost2, mu_delta_e, mu_delta_c, icer)
  sd_v <- c(sd_mu_eff1, sd_mu_cost1, sd_mu_eff2, sd_mu_cost2, sd_delta_e, sd_delta_c, "")
  lower_v <- c(lower_mu_eff1, lower_mu_cost1, lower_mu_eff2, lower_mu_cost2, lower_mu_delta_e, lower_mu_delta_c, "")
  upper_v <- c(upper_mu_eff1, upper_mu_cost1, upper_mu_eff2, upper_mu_cost2, upper_mu_delta_e, upper_mu_delta_c, "")
  cea_table <- cbind(mean_v, sd_v, lower_v, upper_v)
  colnames(cea_table) <- colnames_v
  rownames(cea_table) <- rownames_v
  cea_table_df <- as.data.frame(cea_table)
  table1 <- cea_table_df[1:2, ]
  table2 <- cea_table_df[3:4, ]
  table3 <- cea_table_df[5:7, ]
  cea_table_list <- list(table1, table2, table3)
  cat(paste("\n Cost-effectiveness analysis summary \n \n", "Comparator intervention:", object$cea$interventions[1], "\n","Reference intervention:", 
            object$cea$interventions[2], "\n \n", "Parameter estimates under", object$type, "assumption"))
  cat("\n \n Comparator intervention \n")
  print(table1, quote = FALSE, digits = 3, justify = "center")
  cat("\n Reference intervention \n")
  print(table2, quote = FALSE, digits = 3, justify = "center")
  cat("\n Incremental results \n")
  print(table3, quote = FALSE, digits = 3, justify = "center")
 }