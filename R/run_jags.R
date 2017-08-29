#' An internal function to execute a JAGS model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of missingness mechanism assumed. Choices are Missing Completely At Random (MCAR),
#'  Missing At Random (MAR), Missing Not At Random (MNAR). For 'MNAR' alternative versions are available 
#'  depending on whether the mechanism for only one variable is condiered, that is for the effects (MNAR_eff)
#'  or the costs (MNAR_cost), or also covariates are included either for the effects (MNAR_eff_cov),
#'  the costs (MNAR_cost), or both (MNAR_cov).
#' @param dist_e effect data matrix where rows represent individuals and columns represent the arms (only 2 supported)
#' @param dist_c cost data where rows represent individuals and columns represent the arms (only 2 supprted)
#' @param forward Logical. If \code{forward} is \code{TRUE}, the model is run in forward sampling mode
#'  without providing any data, else if \code{forward} is \code{FALSE} standard sampling mode is selected.
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, BUGS will generate initial values for parameters
#' @keywords JAGS Bayesian model 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_jags<-function(type,dist_e,dist_c,forward,inits)eval.parent(substitute({
  #load namespace of R2jags else ask to install the package
  if(!isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  #check distributions for effects and costs are among those available
  if(!dist_e %in% c("norm","beta")|!dist_c %in% c("norm","gamma")) {
    stop("Distributions available for use are 'norm','beta' for the effects and 'norm','gamma' for the costs")
  }
  #check type of missingness mechanism is among those available
  if(!type %in% c("MCAR","MAR","MNAR","MNAR_eff","MNAR_cost","MNAR_cov","MNAR_eff_cov","MNAR_cost_cov")) {
    stop("Types available for use are 'MCAR','MAR','MNAR_eff','MNAR_cost','MNAR','MNAR_eff_cov','MNAR_cost_cov','MNAR_cov'")
  }
  #define objects that are common to all types of distributions and mechanisms
  #null initialised parameter values by default or user-provided
  if(is.null(inits)==FALSE){inits=inits}
  #write model and create the txt file in current WD
  model<-write_model(type = type ,dist_e = dist_e,dist_c = dist_c,program = "JAGS")
  #define model file from output of write_model
  filein<-model$model_string
  #define data list based on type of MoA and MoM
  #assuming forward sampling (no data) or standard sampling (with data)
  #if no covariate included
  if(type=="MCAR"|type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
    if(stand==TRUE){
      if(forward==TRUE){
        datalist<-list("N1","N2","mean_eff","mean_cost","sd_eff","sd_cost")
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1_s","eff2_s","cost1_s","cost2_s","mean_eff",
                       "sd_eff","mean_cost","sd_cost","m_eff1","m_eff2","m_cost1","m_cost2")
      }
    }else if(stand==FALSE){
      if(forward==TRUE){
        datalist<-list("N1","N2")
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1","eff2","cost1","cost2","m_eff1","m_eff2","m_cost1","m_cost2")
      }
    }
  }
  #if covariate inlcuded
  if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
    if(stand==TRUE){
      if(forward==TRUE){
        datalist<-list("N1","N2","mean_eff","mean_cost","sd_eff","sd_cost",
                       "pe","pc","X1_es","X2_es","X1_cs","X2_cs","mean_cov_e1_t","mean_cov_e2_t","mean_cov_c1_t","mean_cov_c2_t")
        if(pe==1){
          pe_index<-match("pe",datalist)
          datalist<-datalist[-pe_index]
        }
        if(pc==1){
          pc_index<-match("pc",datalist)
          datalist<-datalist[-pc_index]
        }
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1_s","eff2_s","cost1_s","cost2_s","X1_es","X2_es","X1_cs","X2_cs","pe","pc",
                       "mean_eff","sd_eff","mean_cost","sd_cost","m_eff1","m_eff2","m_cost1","m_cost2","mean_cov_e1_t","mean_cov_e2_t","mean_cov_c1_t","mean_cov_c2_t")
        if(pe==1){
          pe_index<-match("pe",datalist)
          datalist<-datalist[-pe_index]
        }
        if(pc==1){
          pc_index<-match("pc",datalist)
          datalist<-datalist[-pc_index]
        }
      }
    }else if(stand==FALSE){
      if(forward==TRUE){
        datalist<-list("N1","N2","pe","pc","X1_e","X2_e","X1_c","X2_c","mean_cov_e1","mean_cov_e2","mean_cov_c1","mean_cov_c2")
        if(pe==1){
          pe_index<-match("pe",datalist)
          datalist<-datalist[-pe_index]
        }
        if(pc==1){
          pc_index<-match("pc",datalist)
          datalist<-datalist[-pc_index]
        }
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1","eff2","cost1","cost2","m_eff1","m_eff2","m_cost1","m_cost2",
                       "pe","pc","X1_e","X2_e","X1_c","X2_c","mean_cov_e1","mean_cov_e2","mean_cov_c1","mean_cov_c2")
        if(pe==1){
          pe_index<-match("pe",datalist)
          datalist<-datalist[-pe_index]
        }
        if(pc==1){
          pc_index<-match("pc",datalist)
          datalist<-datalist[-pc_index]
        }
      }
    }
  }
  #DIC is set to FALSE as no data provided
  DIC<-TRUE
  if(forward==FALSE){DIC=TRUE}else{DIC=FALSE}
  #define all parameters to monitor depending on MoM type selected
  if(type=="MCAR"){params<-c("mu_e","mu_c","s_e","s_c")}
  if(type=="MNAR_eff"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_e","delta_e")}
  if(type=="MNAR_cost"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_c","delta_c")}
  if(type=="MNAR"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_e","gamma0_c","delta_e","delta_c")}
  if(type=="MAR"){params<-c("mu_e","mu_c","s_e","s_c","beta_e","beta_c")}
  if(type=="MNAR_eff_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta_e","gamma_e","delta_e")}
  if(type=="MNAR_cost_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta_c","gamma_c","delta_c")}
  if(type=="MNAR_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta_e","beta_c","gamma_e","gamma_c","delta_e","delta_c")}
  if(ind==FALSE & dist_e=="norm" & dist_c=="norm"){params<-c(params,"theta")}
  if(forward==FALSE){params<-c(params,"eff1","cost1","eff2","cost2")}
  #run model
  modelN1<-R2jags::jags(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=n.chains,
                        n.iter=n.iter,n.burnin = n.burnin,DIC=DIC,n.thin=n.thin)
  #call jags function to perform sampling using all inputs previously defined
  #save parameters simulations into list object to be returned
  #parameters common to all models 
  if(forward==TRUE){
  mu_e<-modelN1$BUGSoutput$sims.list$mu_e
  mu_c<-modelN1$BUGSoutput$sims.list$mu_c
  s_e<-modelN1$BUGSoutput$sims.list$s_e
  s_c<-modelN1$BUGSoutput$sims.list$s_c
  }else if(forward==FALSE){
    mu_e<-modelN1$BUGSoutput$sims.list$mu_e
    mu_c<-modelN1$BUGSoutput$sims.list$mu_c
    s_e<-modelN1$BUGSoutput$sims.list$s_e
    s_c<-modelN1$BUGSoutput$sims.list$s_c
    eff1_pos<-matrix(eff1,N1,3)
    cost1_pos<-matrix(cost1,N1,3)
    eff2_pos<-matrix(eff2,N2,3)
    cost2_pos<-matrix(cost2,N2,3)
    eff1_pos[,1]<-apply(modelN1$BUGSoutput$sims.list$eff1,2,mean)
    eff1_pos[,2]<-apply(modelN1$BUGSoutput$sims.list$eff1,2,quantile,probs=prob[1])
    eff1_pos[,3]<-apply(modelN1$BUGSoutput$sims.list$eff1,2,quantile,probs=prob[2])
    eff2_pos[,1]<-apply(modelN1$BUGSoutput$sims.list$eff2,2,mean)
    eff2_pos[,2]<-apply(modelN1$BUGSoutput$sims.list$eff2,2,quantile,probs=prob[1])
    eff2_pos[,3]<-apply(modelN1$BUGSoutput$sims.list$eff2,2,quantile,probs=prob[2])
    cost1_pos[,1]<-apply(modelN1$BUGSoutput$sims.list$cost1,2,mean)
    cost1_pos[,2]<-apply(modelN1$BUGSoutput$sims.list$cost1,2,quantile,probs=prob[1])
    cost1_pos[,3]<-apply(modelN1$BUGSoutput$sims.list$cost1,2,quantile,probs=prob[2])
    cost2_pos[,1]<-apply(modelN1$BUGSoutput$sims.list$cost2,2,mean)
    cost2_pos[,2]<-apply(modelN1$BUGSoutput$sims.list$cost2,2,quantile,probs=prob[1])
    cost2_pos[,3]<-apply(modelN1$BUGSoutput$sims.list$cost2,2,quantile,probs=prob[2])
  }
  if(ind==FALSE & dist_e=="norm" & dist_c=="norm"){theta<-modelN1$BUGSoutput$sims.list$theta}
  if(type=="MAR"){
    #beta0_e<-modelN1$BUGSoutput$sims.list$beta0_e
    #beta0_c<-modelN1$BUGSoutput$sims.list$beta0_c
    beta_e<-modelN1$BUGSoutput$sims.list$beta_e
    beta_c<-modelN1$BUGSoutput$sims.list$beta_c
  }
  if(type=="MNAR"){
    gamma0_e<-modelN1$BUGSoutput$sims.list$gamma0_e
    delta_e<-modelN1$BUGSoutput$sims.list$delta_e
    gamma0_c<-modelN1$BUGSoutput$sims.list$gamma0_c
    delta_c<-modelN1$BUGSoutput$sims.list$delta_c
  }
  if(type=="MNAR_eff"){
    gamma0_e<-modelN1$BUGSoutput$sims.list$gamma0_e
    delta_e<-modelN1$BUGSoutput$sims.list$delta_e
  }
  if(type=="MNAR_cost"){
    gamma0_c<-modelN1$BUGSoutput$sims.list$gamma0_c
    delta_c<-modelN1$BUGSoutput$sims.list$delta_c
  }
  if(type=="MNAR_cov"){
    #beta0_e<-modelN1$BUGSoutput$sims.list$beta0_e
    #beta0_c<-modelN1$BUGSoutput$sims.list$beta0_c
    #gamma0_e<-modelN1$BUGSoutput$sims.list$gamma0_e
    beta_e<-modelN1$BUGSoutput$sims.list$beta_e
    gamma_e<-modelN1$BUGSoutput$sims.list$gamma_e
    delta_e<-modelN1$BUGSoutput$sims.list$delta_e
    #gamma0_c<-modelN1$BUGSoutput$sims.list$gamma0_c
    beta_c<-modelN1$BUGSoutput$sims.list$beta_c
    gamma_c<-modelN1$BUGSoutput$sims.list$gamma_c
    delta_c<-modelN1$BUGSoutput$sims.list$delta_c
  }
  if(type=="MNAR_eff_cov"){
    #beta0_e<-modelN1$BUGSoutput$sims.list$beta0_e
    #gamma0_e<-modelN1$BUGSoutput$sims.list$gamma0_e
    beta_e<-modelN1$BUGSoutput$sims.list$beta_e
    gamma_e<-modelN1$BUGSoutput$sims.list$gamma_e
    delta_e<-modelN1$BUGSoutput$sims.list$delta_e
  }
  if(type=="MNAR_cost_cov"){
    #beta0_c<-modelN1$BUGSoutput$sims.list$beta0_c
    #gamma0_c<-modelN1$BUGSoutput$sims.list$gamma0_c
    beta_c<-modelN1$BUGSoutput$sims.list$beta_c
    gamma_c<-modelN1$BUGSoutput$sims.list$gamma_c
    delta_c<-modelN1$BUGSoutput$sims.list$delta_c
  }
  #hide constant variables data if standard sampling
  if(forward==FALSE){
   if(n.chains>1){model_sum<-round(jagsresults(x=modelN1, params=c('eff1','eff2','cost1','cost2'), invert=TRUE),digits = 3)
   }else{model_sum<-NULL}
  }else if(forward==TRUE){model_sum<-"default"}
  #save imputed outcome data if standard sampling else set it as null
  if(forward==FALSE){
    #set colnames
    colnames(eff1_pos)<-c("mean","LB","UB")
    colnames(eff2_pos)<-c("mean","LB","UB")
    colnames(cost1_pos)<-c("mean","LB","UB")
    colnames(cost2_pos)<-c("mean","LB","UB")
    imputed<-list("effects1"=eff1_pos,"effects2"=eff2_pos,"costs1"=cost1_pos,"costs2"=cost2_pos)
  }else{imputed<-NULL}
  #define model output list
  if(type=="MCAR"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,"imputed"=imputed,"type"="JAGS")
  }
  if(type=="MAR"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,"imputed"=imputed,"type"="JAGS")
  }
  if(type=="MNAR"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                                "baseline_parameter_effects"=gamma0_e,"baseline_parameter_costs"=gamma0_c,
                                "MNAR_parameter_effects"=delta_e,"MNAR_parameter_costs"=delta_c,"type"="JAGS")
  }
  if(type=="MNAR_eff"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                                "baseline_parameter_effects"=gamma0_e,"MNAR_parameter_effects"=delta_e,"type"="JAGS")
  }
  if(type=="MNAR_cost"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                                "baseline_parameter_costs"=gamma0_c,"MNAR_parameter_costs"=delta_c,"type"="JAGS")
  }
  if(type=="MNAR_cov"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                                "covariate_parameter_miss_effects"=gamma_e,"covariate_parameter_miss_costs"=gamma_c,
                                "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,
                                "MNAR_parameter_effects"=delta_e,"MNAR_parameter_costs"=delta_c,"type"="JAGS")
  }
  if(type=="MNAR_eff_cov"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "covariate_parameter_miss_effects"=gamma_e,"covariate_parameter_effects"=beta_e,"MNAR_parameter_effects"=delta_e,"type"="JAGS")
  }
  if(type=="MNAR_cost_cov"){
    model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
    "covariate_parameter_miss_costs"=gamma_c,"covariate_parameter_costs"=beta_c,"MNAR_parameter_costs"=delta_c,"type"="JAGS")
  }
  if(forward==TRUE){model_output_jags$summary<-NULL}
  if(n.chains==1){model_output_jags<-model_output_jags[-1]}
  return(model_output_jags=model_output_jags)
  #model output list
}))