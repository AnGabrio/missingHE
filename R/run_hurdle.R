#' An internal function to execute a JAGS hurdle model and get posterior results
#'
#' This function fits a JAGS using the \code{\link[R2jags]{jags}} funciton and obtain posterior inferences.
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#' and Structural At Random (MNAR).
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @param se Structural value to be found in the effect data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param inits a list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the BUGS model, 
#' or a function creating (possibly random) initial values. If inits is NULL, JAGS will generate initial values for parameters
#' @param sde hyper-prior value for the standard deviation of the distribution of the structural effects. The default value is
#' \code{1.0E-6} to approximate a point mass at the structural value provided by the user.
#' @param sdc hyper-prior value for the standard deviation of the distribution of the structural costs. The default value is
#' \code{1.0E-6} to approximate a point mass at the structural value provided by the user.
#' @keywords JAGS Bayesian hurdle models 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_hurdle<-function(type,dist_e,dist_c,inits,se=se,sc=sc,sde=sde,sdc=sdc)eval.parent(substitute({
  #load namespace of R2jags else ask to install the package
  if(!isTRUE(requireNamespace("R2jags",quietly=TRUE))) {
    stop("You need to install the R package 'R2jags'. Please run in your R terminal:\n install.packages('R2jags')")
  }
  #check distributions for effects and costs are among those available
  if(!dist_e %in% c("norm","beta")|!dist_c %in% c("norm","gamma")) {
    stop("Distributions available for use are 'norm' or 'beta' for the effects and 'norm'or'gamma' for the costs")
  }
  #check type of structural values mechanism is among those available
  if(!type %in% c("SCAR","SAR")) {
    stop("Types available for use are 'SCAR','SAR'")
  }
  #define objects that are common to all types of distributions and mechanisms
  #null initialised parameter values by default or user-provided
  if(is.null(inits)==FALSE){inits=inits}
  #write model and create the txt file in current WD
  model<-write_hurdle(type = type ,dist_e = dist_e,dist_c = dist_c,se=se,sc=sc)
  #define model file from output of write_model
  filein<-model$model_string
  #sd for structural values for normal on log scale
  if(dist_e=="norm"){sde<-log(sde)}
  if(dist_c=="norm"){sdc<-log(sdc)}
  #for gamma avoid zero if no hurdle selected for costs
  if(dist_c=="gamma"){
  if(is.null(sc)==TRUE){
    cost1<-cost1+0.01
    cost2<-cost2+0.01
  }
  if(is.null(sc)==FALSE){
    if(sc==0){
      if(any(which(cost1==0))==TRUE){
        index_c1_0<-which(cost1==0)
        cost1[index_c1_0]=0.000001
      }
      if(any(which(cost2==0))==TRUE){
        index_c2_0<-which(cost2==0)
        cost2[index_c2_0]=0.000001
      }
    }
   }
  }
  #for beta avoid ones if no hurdle selected for effects
  if(dist_e=="beta"){
    if(is.null(se)==TRUE){
      eff1<-eff1-0.01
      eff2<-eff2-0.01
    }
    if(is.null(se)==FALSE){
      if(se==1){
        if(any(which(eff1==1))==TRUE){
          index_e1_1<-which(eff1==1)
          eff1[index_e1_1]=1-0.0000001
        }
        if(any(which(eff2==1))==TRUE){
          index_e2_1<-which(eff2==1)
          eff2[index_e2_1]=1-0.0000001
        }
      }
    }
  }
  #define data list based on covariates included in models and mechanisms
  if(type=="SCAR"){
        datalist<-list("N1","N2","eff1","eff2","cost1","cost2","d_eff1","d_eff2","d_cost1","d_cost2",
                       "X1_e","X2_e","X1_c","X2_c","mean_cov_e1","mean_cov_e2","mean_cov_c1","se","sc",
                       "mean_cov_c2","pe","pc","sde","sdc")
        if(pe==1){pe_index<-match("pe",datalist)
          datalist<-datalist[-pe_index]}
        if(pc==1){pc_index<-match("pc",datalist)
        datalist<-datalist[-pc_index]}
        if(is.null(se)==TRUE){
          d_eff1_index<-match("d_eff1",datalist)
          d_eff2_index<-match("d_eff2",datalist)
          se_index<-match("se",datalist)
          sde_index<-match("sde",datalist)
          datalist<-datalist[-c(d_eff1_index,d_eff2_index,se_index,sde_index)]
        }
        if(is.null(sc)==TRUE){
          d_cost1_index<-match("d_cost1",datalist)
          d_cost2_index<-match("d_cost2",datalist)
          sc_index<-match("sc",datalist)
          sdc_index<-match("sdc",datalist)
          datalist<-datalist[-c(d_cost1_index,d_cost2_index,sc_index,sdc_index)]
        }
  }
  if(type=="SAR"){
    datalist<-list("N1","N2","eff1","eff2","cost1","cost2","d_eff1","d_eff2","d_cost1","d_cost2",
                   "X1_e","X2_e","X1_c","X2_c","Z1_e","Z2_e","Z1_c","Z2_c","mean_cov_e1","mean_cov_e2","mean_cov_c1","sc","se",
                   "mean_cov_c2","mean_z_e1","mean_z_e2","mean_z_c1","mean_z_c2","pe","pc","ze","zc","sde","sdc")
    if(pe==1){pe_index<-match("pe",datalist)
    datalist<-datalist[-pe_index]}
    if(pc==1){pc_index<-match("pc",datalist)
    datalist<-datalist[-pc_index]}
    if(is.null(ze)==FALSE){
    if(ze==1){ze_index<-match("ze",datalist)
    datalist<-datalist[-ze_index]}
    }else if(is.null(ze)==TRUE){
      ze_index<-match("ze",datalist)
      datalist<-datalist[-ze_index]
      Z1_e_index<-match("Z1_e",datalist)
      datalist<-datalist[-Z1_e_index]
      Z2_e_index<-match("Z2_e",datalist)
      datalist<-datalist[-Z2_e_index]
      mean_z_e1_index<-match("mean_z_e1",datalist)
      datalist<-datalist[-mean_z_e1_index]
      mean_z_e2_index<-match("mean_z_e2",datalist)
      datalist<-datalist[-mean_z_e2_index]
    }
    if(is.null(zc)==FALSE){
      if(zc==1){zc_index<-match("zc",datalist)
      datalist<-datalist[-zc_index]}
    }else if(is.null(zc)==TRUE){
      zc_index<-match("zc",datalist)
      datalist<-datalist[-zc_index]
      Z1_c_index<-match("Z1_c",datalist)
      datalist<-datalist[-Z1_c_index]
      Z2_c_index<-match("Z2_c",datalist)
      datalist<-datalist[-Z2_c_index]
      mean_z_c1_index<-match("mean_z_c1",datalist)
      datalist<-datalist[-mean_z_c1_index]
      mean_z_c2_index<-match("mean_z_c2",datalist)
      datalist<-datalist[-mean_z_c2_index]
    }
    if(is.null(se)==TRUE){
      d_eff1_index<-match("d_eff1",datalist)
      d_eff2_index<-match("d_eff2",datalist)
      se_index<-match("se",datalist)
      sde_index<-match("sde",datalist)
      datalist<-datalist[-c(d_eff1_index,d_eff2_index,se_index,sde_index)]
      if(any(datalist=="ze")==TRUE){
        ze_index2<-match("ze",datalist)
        datalist<-datalist[-ze_index2]
      }
    }
    if(is.null(sc)==TRUE){
      d_cost1_index<-match("d_cost1",datalist)
      d_cost2_index<-match("d_cost2",datalist)
      sc_index<-match("sc",datalist)
      sdc_index<-match("sdc",datalist)
      datalist<-datalist[-c(d_cost1_index,d_cost2_index,sc_index,sdc_index)]
      if(any(datalist=="zc")==TRUE){
        zc_index2<-match("zc",datalist)
        datalist<-datalist[-zc_index2]
      }
    }
  }
  #DIC
  DIC<-TRUE
  #define all parameters to monitor
  if(is.null(se)==TRUE){params<-c("eff1","eff2","cost1","cost2","mu_e","mu_c","s_e","s_c","p_c","beta_c","beta_e","gamma_c")}
  if(is.null(sc)==TRUE){params<-c("eff1","eff2","cost1","cost2","mu_e","mu_c","s_e","s_c","p_e","beta_c","beta_e","gamma_e")}
  if(is.null(se)==FALSE & is.null(sc)==FALSE){
  params<-c("eff1","eff2","cost1","cost2","mu_e","mu_c","s_e","s_c","p_e","p_c","beta_c","beta_e","gamma_e","gamma_c")}
  if(ind==FALSE){params<-c(params,"rho")}
  #run model
  modelN1<-R2jags::jags(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=n.chains,
                        n.iter=n.iter,n.burnin = n.burnin,DIC=DIC,n.thin=n.thin)
  #call jags function to perform sampling using all inputs previously defined
  #save parameters simulations into list object to be returned
  #parameters common to all models 
    mu_e<-modelN1$BUGSoutput$sims.list$mu_e
    mu_c<-modelN1$BUGSoutput$sims.list$mu_c
    s_e<-modelN1$BUGSoutput$sims.list$s_e
    s_c<-modelN1$BUGSoutput$sims.list$s_c
    beta_e<-modelN1$BUGSoutput$sims.list$beta_e
    beta_c<-modelN1$BUGSoutput$sims.list$beta_c
    if(is.null(se)==TRUE & is.null(sc)==FALSE){
      p_c<-modelN1$BUGSoutput$sims.list$p_c
      gamma_c<-modelN1$BUGSoutput$sims.list$gamma_c
    }else if(is.null(sc)==TRUE & is.null(se)==FALSE){
      p_e<-modelN1$BUGSoutput$sims.list$p_e
      gamma_e<-modelN1$BUGSoutput$sims.list$gamma_e
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      p_c<-modelN1$BUGSoutput$sims.list$p_c
      gamma_c<-modelN1$BUGSoutput$sims.list$gamma_c
      p_e<-modelN1$BUGSoutput$sims.list$p_e
      gamma_e<-modelN1$BUGSoutput$sims.list$gamma_e
    }
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
  if(ind==FALSE){rho<-modelN1$BUGSoutput$sims.list$rho}
  #hide constant variables
   if(n.chains>1){model_sum<-round(jagsresults(x=modelN1, params=c('eff1','eff2','cost1','cost2'), invert=TRUE),digits = 3)
   }else{model_sum<-NULL}
  #save imputed outcome data
    #set colnames
    colnames(eff1_pos)<-c("mean","LB","UB")
    colnames(eff2_pos)<-c("mean","LB","UB")
    colnames(cost1_pos)<-c("mean","LB","UB")
    colnames(cost2_pos)<-c("mean","LB","UB")
    imputed<-list("effects1"=eff1_pos,"effects2"=eff2_pos,"costs1"=cost1_pos,"costs2"=cost2_pos)
  #define model output list
    if(is.null(se)==TRUE){
      model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                              "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,"structural_probability_costs"=p_c,
                              "structural_parameter_costs"=gamma_c,"imputed"=imputed,"type"="HURDLE_c","ind"=ind)
    }else if(is.null(sc)==TRUE){
      model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                              "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,"structural_probability_effects"=p_e,
                              "structural_parameter_effects"=gamma_e,"imputed"=imputed,"type"="HURDLE_e","ind"=ind)
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      model_output_jags<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                              "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,"structural_probability_effects"=p_e,
                              "structural_parameter_effects"=gamma_e,"structural_probability_costs"=p_c,"structural_parameter_costs"=gamma_c,
                              "imputed"=imputed,"type"="HURDLE_ec","ind"=ind)
    }
  if(n.chains==1){model_output_jags<-model_output_jags[-1]}
  return(model_output_jags=model_output_jags)
  #model output list
}))