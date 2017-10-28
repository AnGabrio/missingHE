#' A function to read and re-arrange the data in different ways for the hurdle model
#'
#' This internal function imports the data and outputs only those variables that are needed to run the hurdle model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.se}, \code{model.sc} (model formulas for the structural effect and cost models) . Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model.
#' @param model.se A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural effects are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the strcutural effects through a logistic-linear model.
#' @param model.sc A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural costs are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the strcutural costs through a logistic-linear model.
#' @param se Structural value to be found in the effect data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param type Type of structural value mechanism assumed, either 'SCAR' (Structural Completely At Random) or 'SAR' (Strcutural At Random).
#' @keywords read data hurdle models
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response
#' @export
#' @examples
#' \dontrun{
#' #create a data set which respects the requirements specified in "data" (see Arguments)
#' N1<-150
#' N2<-100
#' eff1<-rnorm(N1,0.5,0.5)
#' eff2<-rnorm(N2,0.5,0.5)
#' cost1<-rnorm(N1,90,20)
#' cost2<-rnorm(N2,90,20)
#' #introduce structural values
#' #ones for the effects
#' se=1
#' eff1[1:10]<-1
#' eff2[1:10]<-1
#' #zeros for the costs
#' sc=0
#' cost1[1:10]<-0
#' cost2[1:10]<-0
#' #arrange data frame
#' e<-c(eff1,eff2)
#' c<-c(cost1,cost2)
#' t<-c(t1,t2)
#' data<-data.frame(e,c,t)
#' #run the function
#' date_rearranged<-data_read_hurdle(data=data,model.eff=e~1,model.cost=c~1,
#' model.se=e~1,model.sc=c~1,se=1,sc=0)
#' }
#' #
#' #


data_read_hurdle<-function(data,model.eff,model.cost,model.se,model.sc,se,sc,type=type){
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
  index_mis_e<-index_mis2_e<-which(is.na(data$e))
  index_mis_c<-index_mis2_c<-which(is.na(data$c))
  data2<-data
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
  #define length of complete cases in each arm 
  N1_cc<-N2_cc<-N1_mis<-N2_mis<-c() 
  N1_cc[1]<-length(na.omit(eff1)) 
  N1_cc[2]<-length(na.omit(cost1)) 
  N2_cc[1]<-length(na.omit(eff2)) 
  N2_cc[2]<-length(na.omit(cost2)) 
  N_cc<-cbind(N1_cc,N2_cc) 
  #define length of missing data in each arm 
  N1_mis<-N1-N1_cc 
  N2_mis<-N2-N2_cc 
  N_mis<-cbind(N1_mis,N2_mis) 
  #define variabes for each arm 
  effects<-list(eff1,eff2) 
  costs<-list(cost1,cost2) 
  #define complete case variables in each arm 
  eff1_cc<-eff2_cc<-cost1_cc<-cost2_cc<-c() 
  eff1_cc<-na.omit(eff1) 
  eff2_cc<-na.omit(eff2) 
  eff_cc<-list(eff1_cc,eff2_cc) 
  cost1_cc<-na.omit(cost1) 
  cost2_cc<-na.omit(cost2) 
  cost_cc<-list(cost1_cc,cost2_cc)
  #check indicator for binary variables in X 
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
  cove<-list(cov1_e,cov2_e) 
  mean_cov_e<-list(apply(as.matrix(cov1_e),2,mean),apply(as.matrix(cov2_e),2,mean))
  covc<-list(cov1_c,cov2_c) 
  mean_cov_c<-list(apply(as.matrix(cov1_c),2,mean),apply(as.matrix(cov2_c),2,mean))
  #check column index for those binary 
  check1_e<-which(apply(cov1_e,2,binary_check)) 
  check1_e<-check1_e[-1]
  check2_e<-which(apply(cov2_e,2,binary_check)) 
  check2_e<-check2_e[-1]
  check1_c<-which(apply(cov1_c,2,binary_check)) 
  check1_c<-check1_c[-1]
  check2_c<-which(apply(cov2_c,2,binary_check)) 
  check2_c<-check2_c[-1]
  #replace mean with 1 for categorical variables
  if(any(unlist(lapply(mf_e,is.factor)))==TRUE){
    mean_cov_e[[1]][check1_e]<-1
    mean_cov_e[[2]][check2_e]<-1
  }
  if(any(unlist(lapply(mf_c,is.factor)))==TRUE){
    mean_cov_c[[1]][check1_c]<-1
    mean_cov_c[[2]][check2_c]<-1
  }
  #do the same for the structural values mechanism
  data2$e[is.na(data2$e)==TRUE]<--999999
  data2$c[is.na(data2$c)==TRUE]<--999999
  #create fake vectors to allow model.se and model.sc to depend on e and c when SAR selected
  if(type=="SCAR"|type=="SAR"){
    data2$se<-c(m_eff1,m_eff2)
    data2$sc<-c(m_cost1,m_cost2)
  }
  zf_e <- model.frame(formula=model.se, data=data2)
  if("e"%in%names(zf_e)|"c"%in%names(zf_e)){stop("only dependence on covariates is allowed. Please remove 'e' or 'c' from 'model.se'")}
  if(!any(names(zf_e)%in% names(data2))==TRUE){stop("you must provide names in the formula that correspond to those in the data")}
  zf_c <- model.frame(formula=model.sc, data=data2)
  if("c"%in%names(zf_c)|"e"%in%names(zf_c)){stop("only dependence on covariates is allowed. Please remove 'c' or 'e' from 'model.sc'")}
  if(!any(names(zf_c)%in% names(data2))==TRUE){stop("you must provide names in the formula that correspond to those in the data")}
  if("t"%in%names(zf_e)|"t"%in%names(zf_c)){
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  z_e <- model.matrix(attr(zf_e, "terms"), data=zf_e)
  z_c <- model.matrix(attr(zf_c, "terms"), data=zf_c)
  y_e <- model.response(mf_e)
  y_c <- model.response(mf_c)
  y_e[index_mis2_e]<-NA
  y_c[index_mis2_c]<-NA
  data2$e[index_mis2_e]<-NA
  data2$c[index_mis2_c]<-NA
  #create indicators for structural values
  if(is.null(se)==TRUE & is.null(sc)==TRUE){stop("Structural values in at least one outcome variable are required, 
  please provide the structural value")}
  if(is.null(se)==FALSE){
  index_str_e<-ifelse(y_e==se,1,0)
  if(any(na.omit(index_str_e)==1)==FALSE){stop("Provide structural values that are present in the data")}
  d_eff1<-d_eff2<-c()
  d_eff1<-index_str_e[t1_index]
  d_eff2<-index_str_e[t2_index]
  d_eff<-list(d_eff1,d_eff2) 
  #check indicator for binary variables in X 
  #other categorical variables must be provided as factors to avoid treat them as continuous
  #convert factors into series of dummy variables 
  covz1_e<-as.data.frame(z_e[t1_index,])
  names(covz1_e)<-colnames(z_e)
  covz2_e<-as.data.frame(z_e[t2_index,])
  names(covz2_e)<-colnames(z_e)
  covz_e<-list(covz1_e,covz2_e)
  covze<-list(covz1_e,covz2_e) 
  mean_covz_e<-list(apply(as.matrix(covz1_e),2,mean),apply(as.matrix(covz2_e),2,mean))
  #check column index for those binary 
  checkz1_e<-which(apply(covz1_e,2,binary_check)) 
  checkz1_e<-checkz1_e[-1]
  checkz2_e<-which(apply(covz2_e,2,binary_check)) 
  checkz2_e<-checkz2_e[-1]
  #replace mean with 1 for categorical variables
  if(any(unlist(lapply(zf_e,is.factor)))==TRUE){
    mean_covz_e[[1]][checkz1_e]<-1
    mean_covz_e[[2]][checkz2_e]<-1
  }
  #give labels
  names(covz_e)<-names(mean_covz_e)<-c("Control","Intervention")
  names(d_eff)<-c("Control","Intervention")
  }
  if(is.null(sc)==FALSE){
  index_str_c<-ifelse(y_c==sc,1,0)
  if(any(na.omit(index_str_c)==1)==FALSE){stop("Provide structural values that are present in the data")}
  d_cost1<-d_cost2<-c() 
  d_cost1<-index_str_c[t1_index]
  d_cost2<-index_str_c[t2_index]
  d_cost<-list(d_cost1,d_cost2)  
  #check indicator for binary variables in X 
  #other categorical variables must be provided as factors to avoid treat them as continuous
  #convert factors into series of dummy variables 
  covz1_c<-as.data.frame(z_c[t1_index,])
  names(covz1_c)<-colnames(z_c)
  covz2_c<-as.data.frame(z_c[t2_index,])
  names(covz2_c)<-colnames(z_c)
  covz_c<-list(covz1_c,covz2_c)
  covzc<-list(covz1_c,covz2_c) 
  mean_covz_c<-list(apply(as.matrix(covz1_c),2,mean),apply(as.matrix(covz2_c),2,mean))
  #check column index for those binary 
  checkz1_c<-which(apply(covz1_c,2,binary_check)) 
  checkz1_c<-checkz1_c[-1]
  checkz2_c<-which(apply(covz2_c,2,binary_check)) 
  checkz2_c<-checkz2_c[-1]
  #replace mean with 1 for categorical variables
  if(any(unlist(lapply(zf_c,is.factor)))==TRUE){
    mean_covz_c[[1]][checkz1_c]<-1
    mean_covz_c[[2]][checkz2_c]<-1
  }
  #give labels
  names(covz_c)<-names(mean_covz_c)<-c("Control","Intervention")
  names(d_cost)<-c("Control","Intervention")
  }
  #give labels to output
  names(cov_e)<-names(cov_c)<-names(mean_cov_e)<-names(mean_cov_c)<-c("Control","Intervention")
  names(m_eff)<-names(m_cost)<-c("Control","Intervention")
  names(effects)<-names(costs)<-names(eff_cc)<-names(cost_cc)<-c("Control","Intervention")
  #create list containing all variables to be called by the model according to different mechanisms for se and sc
  if(is.null(se)==FALSE & is.null(sc)==TRUE){
    data_raw<-list("raw_effects"=effects,"raw_costs"=costs,"raw_effects_cc"=eff_cc,"raw_costs_cc"=cost_cc,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis, 
                      "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates_effects"=cov_e,"covariates_costs"=cov_c,
                      "mean_cov_effects"=mean_cov_e, "mean_cov_costs"=mean_cov_c,"covariates_structural_effects"=covz_e,
                      "mean_cov_structural_effects"=mean_covz_e,"structural_effects"=d_eff,"data_ind"=data2) 
  }else if(is.null(se)==TRUE & is.null(sc)==FALSE){
    data_raw<-list("raw_effects"=effects,"raw_costs"=costs,"raw_effects_cc"=eff_cc,"raw_costs_cc"=cost_cc,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis, 
                   "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates_effects"=cov_e,"covariates_costs"=cov_c,
                   "mean_cov_effects"=mean_cov_e, "mean_cov_costs"=mean_cov_c,"covariates_structural_costs"=covz_c,
                   "mean_cov_structural_costs"=mean_covz_c,"structural_costs"=d_cost,"data_ind"=data2) 
  }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
    data_raw<-list("raw_effects"=effects,"raw_costs"=costs,"raw_effects_cc"=eff_cc,"raw_costs_cc"=cost_cc,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis, 
                   "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates_effects"=cov_e,"covariates_costs"=cov_c,
                   "mean_cov_effects"=mean_cov_e, "mean_cov_costs"=mean_cov_c,"covariates_structural_effects"=covz_e,
                   "mean_cov_structural_effects"=mean_covz_e,"covariates_structural_costs"=covz_c,"mean_cov_structural_costs"=mean_covz_c,
                   "structural_effects"=d_eff,"structural_costs"=d_cost,"data_ind"=data2) 
  }
  #return list outputs 
  return(data_raw) 
}


