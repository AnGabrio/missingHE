#' An internal function to read and re-arrange the data in different ways
#'
#' This function imports the data and outputs all types of information used to run the model.
#' @keywords read data
#' @param data a data frame with rows representing the individuals and with columns that must included 
#' the following elements: a column named 'e' containing the effectiveness data for the individuals, a column named 'c'
#' containing the cost data for the individuals, a column 't' representing the tretment arm allocation for each individual 
#' (only two arms are supported with 1 being the control and 2 being the new intervention). Additional elements that could
#' be provided are: a number of columns named 'X1,X2,...' each containing different covariate data for the individuals. 
#' Covariate data must be fully observed and if categorical they must be provided as factor variables. 
#' @importFrom stats na.omit sd as.formula model.matrix
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

read_data<-function(data){
  #define length of treatment arms
  N1<-N2<-c()
  N1<-sum(data$t==1)
  N2<-length(data$t)-N1
  N<-c(N1,N2)
  #create and define missing data indicators for each variable and arm
  indic<-function(x){
    w=NULL
    if (is.na(x)==TRUE)
      w=1
    else w=0
    return(w)
  }
  m_eff1<-m_eff2<-m_cost1<-m_cost2<-c()
  eff1<-data$e[data$t==1]
  eff2<-data$e[data$t==2]
  eff<-list(eff1,eff2)
  cost1<-data$c[data$t==1]
  cost2<-data$c[data$t==2]
  cost<-list(cost1,cost2)
  m_eff1<-sapply(eff1,indic)
  m_eff2<-sapply(eff2,indic)
  m_eff<-list(m_eff1,m_eff2)
  m_cost1<-sapply(cost1,indic)
  m_cost2<-sapply(cost2,indic)
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
  cov<-"no covariate data provided"
  cov_s<-"no covaraite data provided"
  mean_cov<-NULL
  sd_cov<-NULL
  #check indicator for binary variables in X (apply standardisation only to continuous)
  binary_check<-function(x){
    length(unique(x))<=2
  }
  #do the same for covariates if provided
  x_numb<-seq(1:1000)
  if(any(colnames(data)%in% paste("X",x_numb,sep=""))==TRUE){
    cov_data_total<-data[, -which(colnames(data) %in% c("e","c","t"))]
    #if only one predictor need to manually set colname of vector
    if(is.vector(cov_data_total)==TRUE){
      cov_data_total<-as.data.frame(cov_data_total)
      colnames(cov_data_total)<-"X1"
      }
    #convert factors into series of dummy variables
    formula<-as.formula(paste("", paste(c(0,colnames(cov_data_total)),collapse="+"), sep=" ~ "))
    cov_data<-model.matrix(formula,data=cov_data_total)
    cov1_data<-cov_data[data$t==1,]
    cov2_data<-cov_data[data$t==2,]
    #check if factors variables included
    if(any(unlist(lapply(data,is.factor)))==TRUE){
      #check column index for those binary
      check1<-which(apply(cov1_data,2,binary_check))
      check2<-which(apply(cov2_data,2,binary_check))
      #standardise those not binary
      cov1_stand_cont<-as.matrix(cov1_data[,-check1])
      cov1_stand_cont<-apply(cov1_stand_cont,2,stand_05)
      cov2_stand_cont<-as.matrix(cov2_data[,-check2])
      cov2_stand_cont<-apply(cov2_stand_cont,2,stand_05)
      #re-create matrix with all predictors and dummies
      cov1_stand<-cbind(cov1_stand_cont,cov1_data[,check1])
      cov2_stand<-cbind(cov2_stand_cont,cov2_data[,check2])
      #obtain mean and sd for each predictor (for dummy mean=mode and sd=1)
      mean_cov<-sd_cov<-matrix(NA,2,ncol(cov1_data))
      mean_cov1_cont<-apply(as.matrix(cov1_data[,-check1]),2,mean)
      mean_cov2_cont<-apply(as.matrix(cov2_data[,-check2]),2,mean)
      sd_cov1_cont<-apply(as.matrix(cov1_data[,-check1]),2,sd)
      sd_cov2_cont<-apply(as.matrix(cov2_data[,-check2]),2,sd)      
      #obtain mode and set sd=1 for binary
      mode <- function(v) {
        uniqv <- unique(v)
        uniqv[which.max(tabulate(match(v, uniqv)))]
      }
      mean_cov1_bin<-apply(as.matrix(cov1_data[,check1]),2,mode)
      mean_cov2_bin<-apply(as.matrix(cov2_data[,check2]),2,mode)
      mean_cov_bin<-rbind(mean_cov1_bin,mean_cov2_bin)
      sd_cov1_bin<-rep(1,ncol(cov1_data)-ncol(cov1_stand_cont))
      sd_cov2_bin<-rep(1,ncol(cov2_data)-ncol(cov2_stand_cont))
      #combine mean and sd for the two types of variables into one object
      mean_cov<-rbind(c(mean_cov1_cont,mean_cov1_bin),c(mean_cov2_cont,mean_cov2_bin))
      sd_cov<-rbind(c(sd_cov1_cont,sd_cov1_bin),c(sd_cov2_cont,sd_cov2_bin))  
    }else if(any(unlist(lapply(data,is.factor)))==FALSE){
      #standardise
      cov1_stand_cont<-cov1_data
      cov2_stand_cont<-cov2_data
      #if only one predictor need to manually compute standardised values and not use apply
      #single vector
      if(is.vector(cov1_stand_cont)==TRUE & is.vector(cov2_stand_cont)==TRUE){
        cov1_stand<-(cov1_stand_cont-mean(cov1_stand_cont))/sd(cov1_stand_cont)/2
        cov2_stand<-(cov2_stand_cont-mean(cov2_stand_cont))/sd(cov2_stand_cont)/2
        #obtain mean and sd for the predictor
        mean_cov<-sd_cov<-matrix(NA,2,1)
        mean_cov<-rbind(mean(cov1_data),mean(cov2_data))
        sd_cov<-rbind(sd(cov1_data),sd(cov2_data))
      } else{
        #matrix
        cov1_stand<-apply(cov1_stand_cont,2,stand_05)
        cov2_stand<-apply(cov2_stand_cont,2,stand_05)
        #obtain mean and sd for each predictor
        mean_cov<-sd_cov<-matrix(NA,2,ncol(cov1_data))
        mean_cov<-rbind(apply(cov1_data,2,mean),apply(cov2_data,2,mean))
        sd_cov<-rbind(apply(cov1_data,2,sd),apply(cov2_data,2,sd))
      }
    }    
    cov1_s<-matrix(NA,N[1],ncol(as.matrix(cov1_data)))
    cov2_s<-matrix(NA,N[2],ncol(as.matrix(cov2_data)))
    cov1_s<-cov1_stand
    cov2_s<-cov2_stand
    cov_s<-list(cov1_s,cov2_s)
    cov1<-matrix(NA,N[1],ncol(as.matrix(cov1_data)))
    cov2<-matrix(NA,N[2],ncol(as.matrix(cov2_data)))
    cov<-list(cov1_data,cov2_data)
  }
  #create list containing all variables on standardised scale to be called by the model
  data_stand<-list("stand_effects"=eff_s,"stand_costs"=cost_s,"stand_effects_cc"=eff_cc_s,
                   "stand_costs_cc"=cost_cc_s,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis,
                   "mean_effects"=mean_eff,"mean_costs"=mean_cost,"sd_effects"=sd_eff,"sd_costs"=sd_cost,
                   "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates"=cov_s,"mean_cov"=mean_cov,"sd_cov"=sd_cov)
  #create list containing all variables on natural scale to be called by the model
  data_raw<-list("raw_effects"=effects,"raw_costs"=costs,"raw_effects_cc"=eff_cc,"raw_costs_cc"=cost_cc,"arm_lengths"=N,"arm_lengths_cc"=N_cc,"arm_missing_data"=N_mis,
                 "missing_effects"=m_eff,"missing_costs"=m_cost,"covariates"=cov)
  #return list outputs
  return(list(data_raw=data_raw,data_stand=data_stand))
  
}




