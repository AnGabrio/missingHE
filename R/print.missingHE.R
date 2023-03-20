#' Print method for the posterior results contained in the objects of class \code{missingHE}
#' 
#' Prints the summary table for the model fitted, with the estimate of the parameters and/or missing values.
#' @keywords print JAGS missing data  
#' @param x A \code{missingHE} object containing the results of the Bayesian model run using the function \code{\link{selection}}, \code{\link{selection_long}}, \code{\link{pattern}} or \code{\link{hurdle}}.
#' @param value.mis Logical. If \code{value.mis} is \code{TRUE}, the model results displayed contain also the imputed values,
#' else if \code{value.mis} is \code{FALSE} the missing values are hidden.
#' @param only.means Logical. If \code{only.means} is \code{TRUE}, then the \code{print} function only shows the summary
#' statistics for the mean effectiveness and costs. Defaults at \code{FALSE} (in which case, shows the summary statistics
#' for all parameters in the model).
#' @param ... additional arguments affecting the printed output produced. For example: \code{digits=} number of 
#' significant digits to be shown in the printed table (default=3). Not available if \code{value.mis=}TRUE.
#' @author Andrea Gabrio
#' @export
#' @examples  
#' # For examples see the function \code{\link{selection}}, \code{\link{selection_long}}, 
#' # \code{\link{pattern}} or \code{\link{hurdle}}
#' # 
#' # 


print.missingHE <- function(x, value.mis = FALSE, only.means = TRUE, ...) {
  exArgs <- list(...)
  if(exists("digits", where = exArgs)) {digits = exArgs$digits} else {digits = 3}
  if(!inherits(x, what = "missingHE")) { 
    stop("Only objects of class 'missingHE' can be used") 
  }
  if(x$model_output$`model summary`$BUGSoutput$n.chains == 1){
    stop("no output is available if n.chain = 1")
  }
  if(x$data_format == "wide") {
  x_print_sum <- x$model_output$summary[, c(1:3, 7:9)]
  x_print_sum2 <- jagsresults(x = x$model_output$`model summary`, params = c('loglik_e1', 'loglik_e2',
                  'loglik_c1', 'loglik_c2', 'loglik_me1', 'loglik_me2','loglik_mc1', 'loglik_mc2',
                  'loglik_de1', 'loglik_de2', 'loglik_dc1', 'loglik_dc2', 'loglik_d1', 'loglik_d2'), invert = TRUE)
  x_print_sum2 <- x_print_sum2[, c(1:3, 7:9)]
    if(x$model_output$`model summary`$BUGSoutput$n.chains > 1) {
        if(value.mis == FALSE & only.means == FALSE) {
          print(x_print_sum, digits = digits)
        } else if(value.mis == TRUE & only.means == FALSE) {
          print(x_print_sum2, digits = digits)
        }
      if(only.means == TRUE) {
        print(x_print_sum[grep("mu_[ec]\\[", rownames(x_print_sum)),])
      }
    }
  }
  if(x$data_format == "long") {
    x_print_sum <- x$model_output$summary[, c(1:3, 7:9)]
    x_print_sum2 <- jagsresults(x = x$model_output$`model summary`, params = c('loglik_u1', 'loglik_u2',
                                                                               'loglik_c1', 'loglik_c2', 'loglik_mu1', 'loglik_mu2','loglik_mc1', 'loglik_mc2',
                                                                               'loglik_du1', 'loglik_du2', 'loglik_dc1', 'loglik_dc2', 'loglik_d1', 'loglik_d2'), invert = TRUE)
    x_print_sum2 <- x_print_sum2[, c(1:3, 7:9)]
    if(x$model_output$`model summary`$BUGSoutput$n.chains > 1) {
      if(value.mis == FALSE & only.means == FALSE) {
        print(x_print_sum, digits = digits)
      } else if(value.mis == TRUE & only.means == FALSE) {
        print(x_print_sum2, digits = digits)
      }
      if(only.means == TRUE) {
        print(x_print_sum[grep("mu_[uc]\\[", rownames(x_print_sum)),])
      }
    }
  }
}
