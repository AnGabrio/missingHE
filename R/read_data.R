#' A function to read and re-arrange the data in different ways
#'
#' This internal function imports the data and outputs only those variables that are needed to run the model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model.
#' @keywords read data
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response
#' @export
#' @examples
#' \dontrun{
#' #create a data set which respects the requirements specified in "data" (see Arguments)
#' N1<-150
#' N2<-100
#' m_eff1<-m_cost1<-rbinom(N1,1,0.25)
#' m_eff2<-m_cost2<-rbinom(N2,1,0.25)
#' m_cost1<-m_cost1<-rbinom(N1,1,0.25)
#' m_cost2<-m_cost2<-rbinom(N2,1,0.25)
#' eff1<-rnorm(N1,0.5,0.5)
#' eff2<-rnorm(N2,0.5,0.5)
#' cost1<-rnorm(N1,90,20)
#' cost2<-rnorm(N2,90,20)
#' #introduce missingness
#' eff1[m_eff1==1]<-NA
#' eff2[m_eff2==1]<-NA
#' cost1[m_cost1==1]<-NA
#' cost2[m_cost2==1]<-NA
#' #arrange data frame
#' e<-c(eff1,eff2)
#' c<-c(cost1,cost2)
#' m_eff<-c(m_eff1,m_eff2)
#' m_cost<-c(m_cost1,m_cost2)
#' t<-c(t1,t2)
#' data<-data.frame(e,c,t)
#' #run the function
#' date_rearranged<-read_data(data=data)
#' }
#' #
#' #


read_data<-function(data,model.eff=model.eff,model.cost=model.cost){
  if(!any(c("e","c","t")%in% names(data))==TRUE){
    stop("Please rename or provide variables in the data frame as 'e', 'c' and 't' for the effectiveness, cost and treatment indicator")
  }
  if(any(names(data)=="e")==TRUE & any(names(data)=="c")==TRUE){
    e<-as.name("e")
    c<-as.name("c")
  }
  cov_matrix<-subset(data, select = -c(e,c) )
  if(any(is.na(cov_matrix))==TRUE){
    stop("no missing covariate or treatment indicator is allowed")
  }
  if(any(levels(as.factor(cov_matrix$t))!=c("1","2"))==TRUE){
    stop("A two arm indicator variable must be provided")
  }
  index_mis_e<-which(is.na(data$e))
  index_mis_c<-which(is.na(data$c))
  data$e[is.na(data$e)==TRUE]<--999999
  data$c[is.na(data$c)==TRUE]<--999999
  mf_e <- model.frame(formula=model.eff, data=data)
  if("c"%in%names(mf_e)){
    stop("only dependence of costs on effectiveness is allowed. Please remove 'c' from 'model.eff'")
  }
  if(!any(names(mf_e)%in% names(data))==TRUE){
    stop("you must provide names in the formula that correspond to those in the data")
  }
  mf_c <- model.frame(formula=model.cost, data=data)
  if(!any(names(mf_c)%in% names(data))==TRUE){
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if("t"%in%names(mf_e)|"t"%in%names(mf_c)){
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  x_e <- model.matrix(attr(mf_e, "terms"), data=mf_e)
  x_c <- model.matrix(attr(mf_c, "terms"), data=mf_c)
  if("e"%in%names(mf_c)){
    mf_c$e[index_mis_e]<-NA
  }
  y_e <- model.response(mf_e)
  y_c <- model.response(mf_c)
  y_e[index_mis_e]<-NA
  y_c[index_mis_c]<-NA
  data$e[index_mis_e]<-NA
  data$c[index_mis_c]<-NA
  #divide by treatment arm
  N1<-N2<-c() 
  N1<-sum(data$t==1) 
  N2<-length(data$t)-N1 
  N<-c(N1,N2)
  m_eff<-rep(0,length(data$e))
  m_eff[index_mis_e]<-1
  m_cost<-rep(0,length(data$c))
  m_cost[index_mis_c]<-1
  m_eff1<-m_eff2<-m_cost1<-m_cost2<-c() 
  t1_index<-which(data$t==1)
  t2_index<-which(data$t==2)
  eff1<-y_e[t1_index] 
  eff2<-y_e[t2_index] 
  eff<-list(eff1,eff2)
  cost1<-y_c[t1_index] 
  cost2<-y_c[t2_index] 
  cost<-list(cost1,cost2)
  m_eff1<-m_eff[t1_index]
  m_eff2<-m_eff[t2_index]
  m_eff<-list(m_eff1,m_eff2) 
  m_cost1<-m_cost[t1_index]
  m_cost2<-m_cost[t2_index]
  m_cost<-list(m_cost1,m_cost2) 
  #create and define standardised variables (mean 0 and sd 0.5) for each arm 
  stand_05<-function(x){ 
    w=NULL
    w=(x-mean(x,na.rm = TRUE))/sd(x,na.rm = TRUE)/2 
    return(w) 
  } 
  eff1_stand<-lapply(eff,stand_05)[[1]] 
  eff2_stand<-lapply(eff,stand_05)[[2]] 
  cost1_stand<-lapply(cost,stand_05)[[1]] 
  cost2_stand<-lapply(cost,stand_05)[[2]] 
  #define mean and sd for each variable and arm 
  #used to backtransform on natural scale after modelling 
  mean_eff<-mean_cost<-sd_eff<-sd_cost<-c() 
  mean_eff = c(mean(eff1,na.rm=TRUE),mean(eff2,na.rm=TRUE)) 
  mean_cost = c(mean(cost1,na.rm=TRUE),mean(cost2,na.rm=TRUE)) 
  sd_eff<-c(sd(eff1,na.rm=TRUE),sd(eff2,na.rm=TRUE)) 
  sd_cost<-c(sd(cost1,na.rm=TRUE),sd(cost2,na.rm=TRUE)) 
  #define variables on standardised scale for each arm 
  eff1_s<-eff2_s<-cost1_s<-cost2_s<-c() 
  eff1_s<-eff1_stand 
  eff2_s<-eff2_stand 
  eff_s<-list(eff1_s,eff2_s) 
  cost1_s<-cost1_stand 
  cost2_s<-cost2_stand 
  cost_s<-list(cost1_s,cost2_s) 
  #define complete case variables on standardised scale in each arm 
  eff1_cc_s<-eff2_cc_s<-cost1_cc_s<-cost2_cc_s<-c() 
  eff1_cc_s<-na.omit(eff1) 
  eff2_cc_s<-na.omit(eff2) 
  eff_cc_s<-list(eff1_cc_s,eff2_cc_s) 
  cost1_cc_s<-na.omit(cost1) 
  cost2_cc_s<-na.omit(cost2) 
  cost_cc_s<-list(cost1_cc_s,cost2_cc_s) 
  #define length of complete cases in each arm 
  N1_cc<-N2_cc<-N1_mis<-N2_mis<-c() 
  N1_cc[1]<-length(eff1_cc_s) 
  N1_cc[2]<-length(cost1_cc_s) 
  N2_cc[1]<-length(eff2_cc_s) 
  N2_cc[2]<-length(cost2_cc_s) 
  N_cc<-cbind(N1_cc,N2_cc) 
  #define length of missing data in each arm 
  N1_mis<-N1-N1_cc 
  N2_mis<-N2-N2_cc 
  N_mis<-cbind(N1_mis,N2_mis) 
  #define variabes on natrual scale for each arm 
  effects<-list(eff1,eff2) 
  costs<-list(cost1,cost2) 
  #define complete case variables on natural scale in each arm 
  eff1_cc<-eff2_cc<-cost1_cc<-cost2_cc<-c() 
  eff1_cc<-na.omit(eff1) 
  eff2_cc<-na.omit(eff2) 
  eff_cc<-list(eff1_cc,eff2_cc) 
  cost1_cc<-na.omit(cost1) 
  cost2_cc<-na.omit(cost2) 
  cost_cc<-list(cost1_cc,cost2_cc)
  #check indicator for binary variables in X (apply standardisation only to continuous) 
  #other categorical variables must be provided as factors to avoid treat them as continuous
  #convert factors into series of dummy variables 
  cov1_e<-as.data.frame(x_e[t1_index,])
  names(cov1_e)<-colnames(x_e)
  cov2_e<-as.data.frame(x_e[t2_index,])
  names(cov2_e)<-colnames(x_e)
  cov_e<-list(cov1_e,cov2_e)
  x_c_hold<-x_c
  if("e"%in%colnames(x_c_hold)){
    x_c<-subset(x_c_hold,select = -c(e))
  }
  cov1_c<-as.data.frame(x_c[t1_index,])
  names(cov1_c)<-colnames(x_c)
  cov2_c<-as.data.frame(x_c[t2_index,])
  names(cov2_c)<-colnames(x_c)
  cov_c<-list(cov1_c,cov2_c)
  binary_check<-function(x){ 
    length(unique(x))<=2 
  } 
  #obtain matrix of standardised predictors for effects
  if(any(unlist(lapply(mf_e,is.factor)))==TRUE){
  #check column index for those binary 
  check1_e<-which(apply(cov1_e,2,binary_check)) 
  check1_e<-check1_e[-1]
  check2_e<-which(apply(cov2_e,2,binary_check)) 
  check2_e<-check2_e[-1]
  #standardise those not binary 
  cov1_stand_cont_e<-as.matrix(cov1_e[,-c(1,check1_e)])
  cov2_stand_cont_e<-as.matrix(cov2_e[,-c(1,check2_e)])
  if(length(cov1_stand_cont_e)!=0){
  cov1_stand_cont_e<-apply(as.matrix(cov1_stand_cont_e),2,stand_05) 
  cov2_stand_cont_e<-apply(as.matrix(cov2_stand_cont_e),2,stand_05)
  }
  #keep those binary and constant
  cov1_stand_bin_e<-as.matrix(cov1_e[,c(1,check1_e)])
  cov2_stand_bin_e<-as.matrix(cov2_e[,c(1,check2_e)])
  #obtain full covariate matrix standardised
  if(length(cov1_stand_cont_e)==0){
    cov1_s_e<-cov1_e
    }else if(length(cov1_stand_cont_e)!=0){
    cov1_s_e<-data.frame(cov1_stand_bin_e[,1],cov1_stand_cont_e,cov1_stand_bin_e[,-1])
    names(cov1_s_e)<-names(cov1_e)
    }
  if(length(cov2_stand_cont_e)==0){
    cov2_s_e<-cov2_e
  }else if(length(cov2_stand_cont_e)!=0){
    cov2_s_e<-data.frame(cov2_stand_bin_e[,1],cov2_stand_cont_e,cov2_stand_bin_e[,-1])
    names(cov2_s_e)<-names(cov2_e)
  }
 }else if(any(unlist(lapply(mf_e,is.factor)))==FALSE){
    cov1_stand_cont_e<-cov1_e[,-1]
    cov2_stand_cont_e<-cov2_e[,-1]
    if(length(cov1_stand_cont_e)!=0){
    cov1_stand_cont_e<-apply(as.matrix(cov1_stand_cont_e),2,stand_05)
    cov2_stand_cont_e<-apply(as.matrix(cov2_stand_cont_e),2,stand_05)
    }
    #obtain full covariate matrix standardised
    if(ncol(cov1_e)==1){
      cov1_s_e<-cov1_e
      cov2_s_e<-cov2_e
    }else if(ncol(cov1_e)>1){
    cov1_s_e<-data.frame(cov1_e$`(Intercept)`,cov1_stand_cont_e)
    names(cov1_s_e)<-names(cov1_e)
    cov2_s_e<-data.frame(cov2_e$`(Intercept)`,cov2_stand_cont_e)
    names(cov2_s_e)<-names(cov2_e)
    }
  }
  #obtain matrix of standardised predictors for costs
  if(any(unlist(lapply(mf_c,is.factor)))==TRUE){
    #check column index for those binary 
    check1_c<-which(apply(cov1_c,2,binary_check)) 
    check1_c<-check1_c[-1]
    check2_c<-which(apply(cov2_c,2,binary_check)) 
    check2_c<-check2_c[-1]
    #standardise those not binary 
    cov1_stand_cont_c<-as.matrix(cov1_c[,-c(1,check1_c)])
    cov2_stand_cont_c<-as.matrix(cov2_c[,-c(1,check2_c)])
    if(length(cov1_stand_cont_c)!=0){
    cov1_stand_cont_c<-apply(as.matrix(cov1_stand_cont_c),2,stand_05) 
    cov2_stand_cont_c<-apply(as.matrix(cov2_stand_cont_c),2,stand_05)
    }
    #keep those binary and constant
    cov1_stand_bin_c<-as.matrix(cov1_c[,c(1,check1_c)])
    cov2_stand_bin_c<-as.matrix(cov2_c[,c(1,check2_c)])
    #obtain full covariate matrix standardised
    if(length(cov1_stand_cont_c)==0){
      cov1_s_c<-cov1_c
    }else if(length(cov1_stand_cont_c)!=0){
      cov1_s_c<-data.frame(cov1_stand_bin_c[,1],cov1_stand_cont_c,cov1_stand_bin_c[,-1])
      names(cov1_s_c)<-names(cov1_c)
    }
    if(length(cov2_stand_cont_c)==0){
      cov2_s_c<-cov2_c
    }else if(length(cov2_stand_cont_c)!=0){
      cov2_s_c<-data.frame(cov2_stand_bin_c[,1],cov2_stand_cont_c,cov2_stand_bin_c[,-1])
      names(cov2_s_c)<-names(cov2_c)
    }
  }else if(any(unlist(lapply(mf_c,is.factor)))==FALSE){
    cov1_stand_cont_c<-cov1_c[,-1]
    cov2_stand_cont_c<-cov2_c[,-1]
    if(length(cov1_stand_cont_c)!=0){
    cov1_stand_cont_c<-apply(as.matrix(cov1_stand_cont_c),2,stand_05) 
    cov2_stand_cont_c<-apply(as.matrix(cov2_stand_cont_c),2,stand_05)
    }
    #obtain full covariate matrix standardised
    if(ncol(cov1_c)==1){
      cov1_s_c<-cov1_c
      cov2_s_c<-cov2_c
    }else if(ncol(cov1_c)>1){
    cov1_s_c<-data.frame(cov1_c$`(Intercept)`,cov1_stand_cont_c)
    names(cov1_s_c)<-names(cov1_c)
    cov2_s_c<-data.frame(cov2_c$`(Intercept)`,cov2_stand_cont_c)
    names(cov2_s_c)<-names(cov2_c)
    }
  }
  cov_es<-list(cov1_s_e,cov2_s_e) 
  cov_cs<-list(cov1_s_c,cov2_s_c) 
  cove<-list(cov1_e,cov2_e) 
  mean_cov_e<-list(apply(as.matrix(cov1_e),2,mean),apply(as.matrix(cov2_e),2,mean))
  mean_cov_e_t<-list(apply(as.matrix(cov1_s_e),2,mean),apply(as.matrix(cov2_s_e),2,mean))
  sd_cov_e_t<-list(apply(as.matrix(cov1_s_e),2,sd),apply(as.matrix(cov2_s_e),2,sd))
  sd_cov_e_t[[1]]<-sd_cov_e_t[[2]]<-1
  covc<-list(cov1_c,cov2_c) 
  mean_cov_c<-list(apply(as.matrix(cov1_c),2,mean),apply(as.matrix(cov2_c),2,mean))
  mean_cov_c_t<-list(apply(as.matrix(cov1_s_c),2,mean),apply(as.matrix(cov2_s_c),2,mean))
  sd_cov_c_t<-list(apply(as.matrix(cov1_s_c),2,sd),apply(as.matrix(cov2_s_c),2,sd))
  sd_cov_c_t[[1]]<-sd_cov_c_t[[2]]<-1
  #replace mean with 1 for categorical variables
  if(any(unlist(lapply(mf_e,is.factor)))==TRUE){
    mean_cov_e[[1]][check1_e]<-1
    mean_cov_e[[2]][check2_e]<-1
    mean_cov_e_t[[1]][check1_e]<-1
    sd_cov_e_t[[1]][check1_e]<-1
    mean_cov_e_t[[2]][check2_e]<-1
    sd_cov_e_t[[2]][check2_e]<-1
  }
  if(any(unlist(lapply(mf_c,is.factor)))==TRUE){
    mean_cov_c[[1]][check1_c]<-1
    mean_cov_c[[2]][check2_c]<-1
    mean_cov_c_t[[1]][check1_c]<-1
    sd_cov_c_t[[1]][check1_c]<-1
    mean_cov_c_t[[2]][check2_c]<-1
    sd_cov_c_t[[2]][check2_c]<-1
  }
  #give labels to output
  names(cov_es)<-names(cov_cs)<-names(cov_e)<-names(cov_c)<-names(mean_cov_e)<-names(mean_cov_c)<-names(mean_cov_e_t)<-names(mean_cov_c_t)<-names(sd_cov_e_t)<-names(sd_cov_c_t)<-c("Control","Intervention")
  names(eff_s)<-names(cost_s)<-names(eff_cc_s)<-names(cost_cc_s)<-names(m_eff)<-names(m_cost)<-c("Control","Intervention")
  names(effects)<-names(costs)<-names(eff_cc)<-names(cost_cc)<-c("Control","Intervention")
  #create list containing all variables on standardised scale to be called by the model 
  data_stand<-list("stand_effects"=eff_s,"stand_costs"=cost_s,"stand_effects_cc"=eff_cc_s, 
        "stand_costs_cc"=cost_cc_s,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis, 
        "mean_effects"=mean_eff,"mean_costs"=mean_cost,"sd_effects"=sd_eff,"sd_costs"=sd_cost, 
        "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates_effects"=cov_es,"covariates_costs"=cov_cs,
        "mean_cov_effects"=mean_cov_e_t, "mean_cov_costs"=mean_cov_c_t,"sd_cov_effects"=sd_cov_e_t,"sd_cov_costs"=sd_cov_c_t)
  #create list containing all variables on natural scale to be called by the model 
  data_raw<-list("raw_effects"=effects,"raw_costs"=costs,"raw_effects_cc"=eff_cc,"raw_costs_cc"=cost_cc,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis, 
            "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates_effects"=cov_e,"covariates_costs"=cov_c,
            "mean_cov_effects"=mean_cov_e, "mean_cov_costs"=mean_cov_c) 
  #return list outputs 
  return(list(data_raw=data_raw,data_stand=data_stand)) 
}


