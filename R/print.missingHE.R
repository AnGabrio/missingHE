#' Print a summary of the \code{BUGS}/\code{JAGS} model fitted by \code{\link{run_model}}
#' 
#' Prints the summary table for the model fitted, with the estimate of the parameters and/or missing values.
#' @keywords print BUGS JAGS missing data  
#' @param x A \code{missingHE} object containing the results of the Bayesian modelling.
#' @param value.mis Logical. If \code{value.mis} is \code{TRUE}, the model results displayed contain also the imputed values,
#' else if \code{value.mis} is \code{FALSE} the missing values are hidden.
#' @param ... additional arguments affecting the printed output produced. For example: \code{digits=} number of 
#' significant digits to be shown in the printed table (default=3). Not available if \code{value.mis=}TRUE.
#' @details If the argument \code{value.mis} is set to \code{TRUE}, the resulting table will differ based on whether 
#' \code{JAGS} or \code{BUGS} were used to fit the model. Under \code{JAGS} the tables report the posterior results for both observed
#' and imputed data, while under \code{BUGS} only the imputed data posterior results are shown.
#' @author Andrea Gabrio
#' @export
#' @examples  
#' #For examples see the function run_model 
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
  #if forward sampling selected no missing values 
  if(x$model_class=="forward"){
      stop("No missing values if model is run in forward sampling")
  }
  if(x$model_output$type=="JAGS"){
      if(x$model_output$`model summary`$BUGSoutput$n.chains>1){
        if(value.mis==FALSE){
          print(x$model_output$summary,digits = digits)
        }else if(value.mis==TRUE){
          print(x$model_output$`model summary`,digits = digits)
        }
      }else if(x$model_output$`model summary`$BUGSoutput$n.chains==1){
        stop("no output is available if n.chain=1")
      }
    }else if(x$model_output$type=="BUGS"){
      if(x$model_output$`model summary`$n.chains>1){
        if(value.mis==FALSE){
          print(x$model_output$summary,digits = digits)
        }else if(value.mis==TRUE){
          print(x$model_output$`model summary`,digits = digits)
        }
      }else if(x$model_output$`model summary`$n.chains==1){
          stop("no output is available if n.chain=1")
    }
  }
}
