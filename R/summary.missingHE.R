#' Summary method for objects in the class \code{missingHE}
#'
#' Produces a table printout with some summary results of the health economic evaluation simulation model
#'  run using the function \code{\link{run_model}}.
#' @param object A \code{missingHE} object containing the results of the Bayesian modelling and the economic evaluation
#' @param ... additional arguments affecting the summary produced.
#' @return Prints a summary table with some information on the health economic simulation model run with the 
#' function \code{run_model}, and synthetic information on the main parameters of interests.
#' @seealso \code{\link{run_model}} \code{\link{diagnostic_checks}} \code{\link{plot.missingHE}}
#' @author Andrea Gabrio
#' @references 
#' Baio, G.(2012). \emph{Bayesian Methods in Health Economcis}. CRC/Chapman Hall, London.
#' @importFrom stats quantile
#' @export
#' @examples 
#' #For examples see the function run_model
#' #
#' #

summary.missingHE<-function(object,...){
  #define additional inputs as a list
  exArgs <- list(...)
  #x can only be object of class missing
  if(class(object)!="missingHE"){
    stop("Only objects of class 'missing' can be used")
  }
  #not possible if forward sampling selected
  if(object$model_class=="forward"){
    stop("summary results not available for forward sampling output")
  }
  #create table with summary CEA results
  #obtain results to plug in
  #mean column
  mu_eff1<-round(mean(object$cea$e[,1]),digits=3)
  mu_eff2<-round(mean(object$cea$e[,2]),digits = 3)
  mu_cost1<-round(mean(object$cea$c[,1]),digits = 3)
  mu_cost2<-round(mean(object$cea$c[,2]),digits = 3)
  sigma_eff1<-round(mean(object$model_output$sd_effects[,1]),digits = 3)
  sigma_eff2<-round(mean(object$model_output$sd_effects[,2]),digits = 3)
  sigma_cost1<-round(mean(object$model_output$sd_costs[,1]),digits = 3)
  sigma_cost2<-round(mean(object$model_output$sd_costs[,2]),digits = 3)
  mu_delta_e<-round(mean(object$cea$delta.e),digits = 3)
  mu_delta_c<-round(mean(object$cea$delta.c),digits = 3)
  icer<-round(object$cea$ICER,digits = 3)
  #sd column
  sd_mu_eff1<-round(sd(object$cea$e[,1]),digits = 3)
  sd_mu_eff2<-round(sd(object$cea$e[,2]),digits = 3)
  sd_mu_cost1<-round(sd(object$cea$c[,1]),digits = 3)
  sd_mu_cost2<-round(sd(object$cea$c[,2]),digits = 3)
  sd_sigma_eff1<-round(sd(object$model_output$sd_effects[,1]),digits = 3)
  sd_sigma_eff2<-round(sd(object$model_output$sd_effects[,2]),digits = 3)
  sd_sigma_cost1<-round(sd(object$model_output$sd_costs[,1]),digits = 3)
  sd_sigma_cost2<-round(sd(object$model_output$sd_costs[,2]),digits = 3)
  sd_delta_e<-round(sd(object$cea$delta.e),digits = 3)
  sd_delta_c<-round(sd(object$cea$delta.c),digits = 3)
  #lower bound CI column
  lower_mu_eff1<-round(quantile(object$model_output$mean_effects[,1],probs = 0.05),digits = 3)
  lower_mu_eff2<-round(quantile(object$model_output$mean_effects[,2],probs = 0.05),digits = 3)
  lower_mu_cost1<-round(quantile(object$model_output$mean_costs[,1],probs = 0.05),digits = 3)
  lower_mu_cost2<-round(quantile(object$model_output$mean_costs[,2],probs = 0.05),digits = 3)
  lower_sigma_eff1<-round(quantile(object$model_output$sd_effects[,1],probs = 0.05),digits = 3)
  lower_sigma_eff2<-round(quantile(object$model_output$sd_effects[,2],probs = 0.05),digits = 3)
  lower_sigma_cost1<-round(quantile(object$model_output$sd_costs[,1],probs = 0.05),digits = 3)
  lower_sigma_cost2<-round(quantile(object$model_output$sd_costs[,2],probs = 0.05),digits = 3)
  lower_mu_delta_e<-round(quantile(object$cea$delta.e,probs = 0.05),digits = 3)
  lower_mu_delta_c<-round(quantile(object$cea$delta.c,probs = 0.05),digits = 3)
  #upper bound CI column
  upper_mu_eff1<-round(quantile(object$model_output$mean_effects[,1],probs = 0.95),digits = 3)
  upper_mu_eff2<-round(quantile(object$model_output$mean_effects[,2],probs = 0.95),digits = 3)
  upper_mu_cost1<-round(quantile(object$model_output$mean_costs[,1],probs = 0.95),digits = 3)
  upper_mu_cost2<-round(quantile(object$model_output$mean_costs[,2],probs = 0.95),digits = 3)
  upper_sigma_eff1<-round(quantile(object$model_output$sd_effects[,1],probs = 0.95),digits = 3)
  upper_sigma_eff2<-round(quantile(object$model_output$sd_effects[,2],probs = 0.95),digits = 3)
  upper_sigma_cost1<-round(quantile(object$model_output$sd_costs[,1],probs = 0.95),digits = 3)
  upper_sigma_cost2<-round(quantile(object$model_output$sd_costs[,2],probs = 0.95),digits = 3)
  upper_mu_delta_e<-round(quantile(object$cea$delta.e,probs = 0.95),digits = 3)
  upper_mu_delta_c<-round(quantile(object$cea$delta.c,probs = 0.95),digits = 3)
  #names of rows
  rownames_v<-c("mean effects","mean costs","sd effects","sd costs","mean effects","mean costs","sd effects","sd costs",
                "delta effects","delta costs","ICER")
  #names of columns
  colnames_v<-c("mean","sd","LB","UB")
  mean_v<-c(mu_eff1,mu_cost1,sigma_eff1,sigma_cost1,mu_eff2,mu_cost2,sigma_eff2,sigma_cost2,mu_delta_e,mu_delta_c,icer)
  sd_v<-c(sd_mu_eff1,sd_mu_cost1,sd_sigma_eff1,sd_sigma_cost1,sd_mu_eff2,sd_mu_cost2,sd_sigma_eff2,sd_sigma_cost2,sd_delta_e,sd_delta_c,"")
  lower_v<-c(lower_mu_eff1,lower_mu_cost1,lower_sigma_eff1,lower_sigma_cost1,lower_mu_eff2,lower_mu_cost2,lower_sigma_eff2,lower_sigma_cost2,lower_mu_delta_e,lower_mu_delta_c,"")
  upper_v<-c(upper_mu_eff1,upper_mu_cost1,upper_sigma_eff1,upper_sigma_cost1,upper_mu_eff2,upper_mu_cost2,upper_sigma_eff2,upper_sigma_cost2,upper_mu_delta_e,upper_mu_delta_c,"")
  #create output table
  cea_table<-cbind(mean_v,sd_v,lower_v,upper_v)
  colnames(cea_table)<-colnames_v
  rownames(cea_table)<-rownames_v
  #convert into data frame
  cea_table_df<-as.data.frame(cea_table)
  #split into multiple tables for better display in R console
  table1<-cea_table_df[1:4,]
  table2<-cea_table_df[5:8,]
  table3<-cea_table_df[9:11,]
  #create list
  cea_table_list<-list(table1,table2,table3)
  #write sentences to be placed before the table 
  cat(paste("\n Cost-effectiveness analysis summary \n \n","Comparator intervention:",object$cea$interventions[1],"\n","Reference intervention:", object$cea$interventions[2],"\n \n",
            "Model of Analysis (MoA) parameter estimates under",object$type,"assumption"))
  cat("\n \n Comparator intervention \n")
  print(table1,quote = F, digits = 3, justify = "center")
  cat("\n Reference intervention \n")
  print(table2,quote = F, digits = 3, justify = "center")
  cat("\n Incremental results \n")
  print(table3,quote = F, digits = 3, justify = "center")
  return("cea_table"=cea_table_list)
  }