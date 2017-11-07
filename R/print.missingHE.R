#' Print a summary of the \code{JAGS} model fitted by \code{\link{selection}} or \code{\link{hurdle}}
#' 
#' Prints the summary table for the model fitted, with the estimate of the parameters and/or missing values.
#' @keywords print JAGS missing data  
#' @param x A \code{missingHE} object containing the results of the Bayesian model run using the function \code{\link{selection}} or \code{\link{hurdle}}.
#' @param value.mis Logical. If \code{value.mis} is \code{TRUE}, the model results displayed contain also the imputed values,
#' else if \code{value.mis} is \code{FALSE} the missing values are hidden.
#' @param ... additional arguments affecting the printed output produced. For example: \code{digits=} number of 
#' significant digits to be shown in the printed table (default=3). Not available if \code{value.mis=}TRUE.
#' @author Andrea Gabrio
#' @export
#' @examples  
#' #For examples see the function selection or hurdle 
#' # 
#' # 


print.missingHE<-function(x,value.mis=FALSE,...){
  #define additional inputs as a list
  exArgs <- list(...)
  #number of digits as additional option
  if(exists("digits",where=exArgs)) {digits=exArgs$digits} else {digits=3}
  #x can only be object of class missingHE 
  if(class(x)!="missingHE"){ 
    stop("Only objects of class 'missingHE' can be used") 
  }
  if(x$model_output$`model summary`$BUGSoutput$n.chains==1){
    stop("no output is available if n.chain=1")
  }
  #redefine output
  x_print_sum<-x$model_output$summary[,c(1:3,7:9)]
  x_print_sum2<-x$model_output$`model summary`$BUGSoutput$summary[,c(1:3,7:9)]
    if(x$model_output$`model summary`$BUGSoutput$n.chains>1){
        if(value.mis==FALSE){
          print(x_print_sum,digits = digits)
        }else if(value.mis==TRUE){
          print(x_print_sum2,digits = digits)
        }
      }
}
