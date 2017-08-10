#' An internal function to execute a BUGS model and get posterior results
#'
#' This function fits a BUGS using the \code{\link[R2OpenBUGS]{bugs}} funciton and obtain posterior inferences.
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
#' @keywords BUGS Bayesian model
#' @importFrom R2OpenBUGS bugs
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


run_bugs<-function(type,dist_e,dist_c,forward,inits)eval.parent(substitute({
  #load namespace of R2OpenBUGS else ask to install the package
  if(!isTRUE(requireNamespace("R2OpenBUGS",quietly=TRUE))) {
    stop("You need to install the R package 'R2OpenBUGS'. Please run in your R terminal:\n install.packages('R2OpenBUGS')")
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
  #adjust initialised values for BUGS so it will not stuck too often
  #for all combinations of models
  if(is.null(inits)==TRUE){
  if(dist_e=="beta" & dist_c=="gamma"){
    if(type=="MCAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=runif(2,0,1),mu_c=runif(2,0,100),s_c=runif(2,0,10),s_e=runif(2,0,0.001))
        inits <- function(){list_init}
        if(any(transf=="logit")==TRUE){#logit transformation
          list_init[[3]]<-rnorm(2,0,1)
          names(list_init)[3]<-"nu_e"
          inits <- function(){list_init}}
        if(any(transf=="log")==TRUE){#log transformation
          list_init[[4]]<-rnorm(2,0,1)
          names(list_init)[4]<-"nu_c"
          inits <- function(){list_init}}
    }else if(type=="MAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),beta0_e=runif(2,0,1),beta0_c=runif(2,0,100),s_c=runif(2,0,10),s_e=runif(2,0,0.01))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init[[7]]<-rnorm(2,0,1)
        inits <- function(){list_init}}
      if(any(transf=="log")==TRUE){#log transformation
        list_init[[8]]<-rnorm(2,0,1)
        inits <- function(){list_init}}
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=runif(2,0,1)
                      ,mu_c=runif(2,0,100),s_c=runif(2,0,10),s_e=runif(2,0,0.01),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init[[3]]<-rnorm(2,0,1)
        names(list_init)[3]<-"nu_e"}
      if(any(transf=="log")==TRUE){#log transformation
        list_init[[4]]<-rnorm(2,0,1)
        names(list_init)[4]<-"nu_c"}
      inits <- function(){list_init}
      if(type=="MNAR_eff"){inits<-function(){list_init[-8]}
      } else if(type=="MNAR_cost"){inits<-function(){list_init[-7]}
     }
    }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=runif(2,0,1),beta0_c=rnorm(2,0,1),s_c=runif(2,0,10),s_e=runif(2,0,0.01),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init$beta0_e<-rnorm(2,0,1)}
      if(any(transf=="log")==TRUE){#log transformation
        list_init$beta0_c<-rnorm(2,0,1)}
      inits <- function(){list_init}
      if(type=="MNAR_eff_cov"){inits<-function(){list_init[-12]}
      } else if(type=="MNAR_cost_cov"){inits<-function(){list_init[-11]}
      }
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }
  }else if(dist_e=="beta" & dist_c=="norm"){
    if(type=="MCAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=runif(2,0,1),mu_c=rnorm(2,0,1),ls_c=runif(2,0,1),s_e=runif(2,0,0.01))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init[[3]]<-rnorm(2,0,1)
        names(list_init)[3]<-"nu_e"
        inits <- function(){list_init}}
    }else if(type=="MAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),beta0_e=runif(2,0,1),beta0_c=rnorm(2,0,1),ls_c=runif(2,0,1),s_e=runif(2,0,0.01))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init$beta0_e<-rnorm(2,0,1)
        inits <- function(){list_init}}
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      list_init<-list(gamma0_e=rlogis(1,0,1),gamma0_c=rlogis(1,0,1),mu_e=runif(2,0,1)
                      ,mu_c=rnorm(2,0,1),ls_c=runif(2,0,1),s_e=runif(2,0,0.01),
                      delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init[[3]]<-rnorm(2,0,1)
        names(list_init)[3]<-"nu_e"}
        inits <- function(){list_init}
      if(type=="MNAR_eff"){inits<-function(){list_init[-8]}
       } else if(type=="MNAR_cost"){inits<-function(){list_init[-7]}
      }
    }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=runif(2,0,1),beta0_c=rnorm(2,0,1),
                      ls_c=runif(2,0,1),s_e=runif(2,0,0.01),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(any(transf=="logit")==TRUE){#logit transformation
        list_init$beta0_e<-rnorm(2,0,1)
        inits <- function(){list_init}}
      if(type=="MNAR_eff_cov"){inits<-function(){list_init[-12]}
       } else if(type=="MNAR_cost_cov"){inits<-function(){list_init[-11]}
       }
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
     }
  }else if(dist_e=="norm" & dist_c=="gamma"){
    if(type=="MCAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=rnorm(2,0,1),mu_c=runif(2,0,10),s_c=runif(2,0,10),ls_e=runif(2,0,1))
      inits <- function(){list_init}
    if(any(transf=="log")==TRUE){#log transformation
      list_init[[4]]<-rnorm(2,0,1)
      names(list_init)[4]<-"nu_c"
      inits <- function(){list_init}}
    }else if(type=="MAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=rnorm(2,0,1),beta0_c=runif(2,0,100),s_c=runif(2,0,10),ls_e=runif(2,0,1))
      inits <- function(){list_init}
    if(any(transf=="log")==TRUE){#log transformation
      list_init$beta0_c<-rnorm(2,0,1)
      inits <- function(){list_init}}
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=rnorm(2,0,1)
                      ,mu_c=runif(2,0,100),s_c=runif(2,0,10),ls_e=runif(2,0,1),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,0.001))
      inits <- function(){list_init}
      if(any(transf=="log")==TRUE){#log transformation
        list_init[[4]]<-rnorm(2,0,1)
        names(list_init)[4]<-"nu_c"
        inits <- function(){list_init}}
      if(type=="MNAR_eff"){inits<-function(){list_init[-8]}
      } else if(type=="MNAR_cost"){inits<-function(){list_init[-7]}
      }
    }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=rnorm(2,0,1),beta0_c=runif(2,0,100),s_c=runif(2,0,10),ls_e=runif(2,0,1),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(any(transf=="log")==TRUE){#log transformation
        list_init$beta0_c<-rnorm(2,0,1)
        inits <- function(){list_init}}
      if(type=="MNAR_eff_cov"){inits<-function(){list_init[-12]}
       } else if(type=="MNAR_cost_cov"){inits<-function(){list_init[-11]}
       }
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
     }
  }else if(dist_e=="norm" & dist_c=="norm"){
    if(type=="MCAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=rnorm(2,0,1)
                      ,mu_c=rnorm(2,0,1),ls_c=runif(2,-1,1),ls_e=runif(2,-1,1))
      inits <- function(){list_init}
      if(stand==TRUE){#scaled model
        names(list_init)[3]<-"mu_e_t"
        names(list_init)[4]<-"mu_c_t"
        names(list_init)[5]<-"ls_c_t"
        names(list_init)[6]<-"ls_e_t"
        inits <- function(){list_init}
      }
      if(ind==FALSE){#joint model
        list_init$theta<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MAR"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=rnorm(2,0,1),beta0_c=rnorm(2,0,1),
                      ls_c=runif(2,-1,1),ls_e=runif(2,-1,1))
      inits <- function(){list_init}
      if(stand==TRUE){#scaled model
        names(list_init)[9]<-"ls_c_t"
        names(list_init)[10]<-"ls_e_t"
        inits <- function(){list_init}
      }
      if(ind==FALSE){#joint model
        list_init$theta<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),mu_e=rnorm(2,0,1)
                      ,mu_c=rnorm(2,0,1),ls_c=runif(2,-1,1),ls_e=runif(2,-1,1),
                      delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(type=="MNAR_eff"){inits<-function(){list_init[-8]}
      } else if(type=="MNAR_cost"){inits<-function(){list_init[-7]}
      }
      if(stand==TRUE){#scaled model
        names(list_init)[3]<-"mu_e_t"
        names(list_init)[4]<-"mu_c_t"
        names(list_init)[5]<-"ls_c_t"
        names(list_init)[6]<-"ls_e_t"
        inits <- function(){list_init}
      }
      if(ind==FALSE){#joint model
        list_init$theta<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
    }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      list_init<-list(gamma0_e=rnorm(1,0,1),gamma0_c=rnorm(1,0,1),
                      beta_e=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),beta_c=cbind(rnorm(p,0,0.001),rnorm(p,0,0.001)),
                      gamma_e=rnorm(p,0,0.001),gamma_c=rnorm(p,0,0.001),
                      beta0_e=rnorm(2,0,1),beta0_c=rnorm(2,0,1),
                      ls_c=runif(2,-1,1),ls_e=runif(2,-1,1),delta_e=rnorm(1,0,1),delta_c=rnorm(1,0,1))
      inits <- function(){list_init}
      if(type=="MNAR_eff_cov"){inits<-function(){list_init[-12]}
       } else if(type=="MNAR_cost_cov"){inits<-function(){list_init[-11]}
       }
      if(stand==TRUE){#scaled model
        names(list_init)[9]<-"ls_c_t"
        names(list_init)[10]<-"ls_e_t"
        inits <- function(){list_init}
      }
      if(ind==FALSE){#joint model
        list_init$theta<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
      if(p==1){#one covariate
        list_init$beta_e<-rnorm(2,0,0.001)
        list_init$beta_c<-rnorm(2,0,0.001)
        inits <- function(){list_init}
      }
     }
    }
  }
  #write model and create the txt file in current WD
  model<-write_model(type = type ,dist_e = dist_e,dist_c = dist_c,program = "BUGS")
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
                       "p","X1_s","X2_s","mean_cov1_t","mean_cov2_t")
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1_s","eff2_s","cost1_s","cost2_s","X1_s","X2_s","p",
                       "mean_eff","sd_eff","mean_cost","sd_cost","m_eff1","m_eff2","m_cost1","m_cost2","mean_cov1_t","mean_cov2_t")
      }
    }else if(stand==FALSE){
      if(forward==TRUE){
        datalist<-list("N1","N2","p","X1","X2","mean_cov1","mean_cov2")
      }else if(forward==FALSE){
        datalist<-list("N1","N2","eff1","eff2","cost1","cost2","m_eff1","m_eff2","m_cost1","m_cost2","X1","X2",
                       "p","mean_cov1","mean_cov2")
      }
    }
  }
  #remove p if p=1
  if(p==1){
    pos_p<-which(datalist[]=="p")
    datalist[pos_p]<-NULL
  }
  #DIC is set to FALSE as no data provided
  DIC<-TRUE
  if(forward==FALSE){DIC=TRUE}else{DIC=FALSE}
  #define all parameters to monitor depending on MoM type selected
  if(type=="MCAR"){params<-c("mu_e","mu_c","s_e","s_c")}
  if(type=="MNAR_eff"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_e","delta_e")}
  if(type=="MNAR_cost"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_c","delta_c")}
  if(type=="MNAR"){params<-c("mu_e","mu_c","s_e","s_c","gamma0_e","gamma0_c","delta_e","delta_c")}
  if(type=="MAR"){params<-c("mu_e","mu_c","s_e","s_c","beta0_e","beta0_c","beta_e","beta_c")}
  if(type=="MNAR_eff_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta0_e","beta_e","gamma0_e","gamma_e","delta_e")}
  if(type=="MNAR_cost_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta0_c","beta_c","gamma0_c","gamma_c","delta_c")}
  if(type=="MNAR_cov"){params<-c("mu_e","mu_c","s_e","s_c","beta0_e","beta0_c","beta_e","beta_c","gamma0_e","gamma0_c","gamma_e","gamma_c","delta_e","delta_c")}
  if(ind==FALSE & dist_e=="norm" & dist_c=="norm"){params<-c(params,"theta")}
  if(forward==FALSE){params<-c(params,"eff1","cost1","eff2","cost2")}
  #index for data vectors
  #remove variable from parameters if no missing values
  index_eff1<-match("eff1",params)
  if(any(is.na(eff1))==FALSE){params<-params[-index_eff1]}
  index_eff2<-match("eff2",params)
  if(any(is.na(eff2))==FALSE){params<-params[-index_eff2]}
  index_cost1<-match("cost1",params)
  if(any(is.na(cost1))==FALSE){params<-params[-index_cost1]}
  index_cost2<-match("cost2",params)
  if(any(is.na(cost2))==FALSE){params<-params[-index_cost2]}
  #run model
  modelN1<-R2OpenBUGS::bugs(data=datalist,inits=inits,parameters.to.save=params,model.file=filein,n.chains=n.chains,
                        n.iter=n.iter,n.burnin = n.burnin,DIC=DIC,n.thin=n.thin)
  #call bugs function to perform sampling using all inputs previously defined
  #save parameters simulations into list object to be returned
  #parameters common to all models 
  if(forward==TRUE){
    mu_e<-modelN1$sims.list$mu_e
    mu_c<-modelN1$sims.list$mu_c
    s_e<-modelN1$sims.list$s_e
    s_c<-modelN1$sims.list$s_c
  }else if(forward==FALSE){
    mu_e<-modelN1$sims.list$mu_e
    mu_c<-modelN1$sims.list$mu_c
    s_e<-modelN1$sims.list$s_e
    s_c<-modelN1$sims.list$s_c
    eff1_pos<-matrix(eff1,N1,3)
    cost1_pos<-matrix(cost1,N1,3)
    eff2_pos<-matrix(eff2,N2,3)
    cost2_pos<-matrix(cost2,N2,3)
    if(any(is.na(eff1))==TRUE){
      eff1_pos[,1][is.na(eff1_pos[,1])]<-apply(modelN1$sims.list$eff1,2,mean)
      eff1_pos[,2][is.na(eff1_pos[,2])]<-apply(modelN1$sims.list$eff1,2,quantile,probs=prob[1])
      eff1_pos[,3][is.na(eff1_pos[,3])]<-apply(modelN1$sims.list$eff1,2,quantile,probs=prob[2])
    }
    if(any(is.na(eff2))==TRUE){
      eff2_pos[,1][is.na(eff2_pos[,1])]<-apply(modelN1$sims.list$eff2,2,mean) 
      eff2_pos[,2][is.na(eff2_pos[,2])]<-apply(modelN1$sims.list$eff2,2,quantile,probs=prob[1]) 
      eff2_pos[,3][is.na(eff2_pos[,3])]<-apply(modelN1$sims.list$eff2,2,quantile,probs=prob[2]) 
    }
    if(any(is.na(cost1))==TRUE){
      cost1_pos[,1][is.na(cost1_pos[,1])]<-apply(modelN1$sims.list$cost1,2,mean) 
      cost1_pos[,2][is.na(cost1_pos[,2])]<-apply(modelN1$sims.list$cost1,2,quantile,probs=prob[1]) 
      cost1_pos[,3][is.na(cost1_pos[,3])]<-apply(modelN1$sims.list$cost1,2,quantile,probs=prob[2]) 
    }
    if(any(is.na(cost2))==TRUE){
      cost2_pos[,1][is.na(cost2_pos[,1])]<-apply(modelN1$sims.list$cost2,2,mean) 
      cost2_pos[,2][is.na(cost2_pos[,2])]<-apply(modelN1$sims.list$cost2,2,quantile,probs=prob[1]) 
      cost2_pos[,3][is.na(cost2_pos[,3])]<-apply(modelN1$sims.list$cost2,2,quantile,probs=prob[2]) 
    }
  }
  if(ind==FALSE & dist_e=="norm" & dist_c=="norm"){theta<-modelN1$sims.list$theta}
  if(type=="MAR"){
    beta0_e<-modelN1$sims.list$beta0_e
    beta0_c<-modelN1$sims.list$beta0_c
    beta_e<-modelN1$sims.list$beta_e
    beta_c<-modelN1$sims.list$beta_c
  }
  if(type=="MNAR"){
    gamma0_e<-modelN1$sims.list$gamma0_e
    delta_e<-modelN1$sims.list$delta_e
    gamma0_c<-modelN1$sims.list$gamma0_c
    delta_c<-modelN1$sims.list$delta_c
  }
  if(type=="MNAR_eff"){
    gamma0_e<-modelN1$sims.list$gamma0_e
    delta_e<-modelN1$sims.list$delta_e
  }
  if(type=="MNAR_cost"){
    gamma0_c<-modelN1$sims.list$gamma0_c
    delta_c<-modelN1$sims.list$delta_c
  }
  if(type=="MNAR_cov"){
    beta0_e<-modelN1$sims.list$beta0_e
    beta0_c<-modelN1$sims.list$beta0_c
    gamma0_e<-modelN1$sims.list$gamma0_e
    beta_e<-modelN1$sims.list$beta_e
    gamma_e<-modelN1$sims.list$gamma_e
    delta_e<-modelN1$sims.list$delta_e
    gamma0_c<-modelN1$sims.list$gamma0_c
    beta_c<-modelN1$sims.list$beta_c
    gamma_c<-modelN1$sims.list$gamma_c
    delta_c<-modelN1$sims.list$delta_c
  }
  if(type=="MNAR_eff_cov"){
    beta0_e<-modelN1$sims.list$beta0_e
    gamma0_e<-modelN1$sims.list$gamma0_e
    beta_e<-modelN1$sims.list$beta_e
    gamma_e<-modelN1$sims.list$gamma_e
    delta_e<-modelN1$sims.list$delta_e
  }
  if(type=="MNAR_cost_cov"){
    beta0_c<-modelN1$sims.list$beta0_c
    gamma0_c<-modelN1$sims.list$gamma0_c
    beta_c<-modelN1$sims.list$beta_c
    gamma_c<-modelN1$sims.list$gamma_c
    delta_c<-modelN1$sims.list$delta_c
  }
  #hide constant variables data if standard sampling
  #adapt summary of results based on missingness in variables
  if(any(params=="eff1")==TRUE){
    index_eff1<-match("eff1",params)
    params<-params[-index_eff1]}
  if(any(params=="eff2")==TRUE){
    index_eff2<-match("eff2",params)
    params<-params[-index_eff2]}
  if(any(params=="cost1")==TRUE){
    index_cost1<-match("cost1",params)
    params<-params[-index_cost1]}
  if(any(params=="cost2")==TRUE){
    index_cost2<-match("cost2",params)
    params<-params[-index_cost2]}
  if(forward==FALSE){
    if(n.chains>1){model_sum<-round(bugsresults(x=modelN1, params=params, invert=FALSE),digits = 3)
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
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,"imputed"=imputed,"type"="BUGS")
  }
  if(type=="MAR"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_effects"=beta0_e,"baseline_parameter_costs"=beta0_c,
                            "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,"imputed"=imputed,"type"="BUGS")
  }
  if(type=="MNAR"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_effects"=gamma0_e,"baseline_parameter_costs"=gamma0_c,
                            "MNAR_parameter_effects"=delta_e,"MNAR_parameter_costs"=delta_c,"type"="BUGS")
  }
  if(type=="MNAR_eff"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_effects"=gamma0_e,"MNAR_parameter_effects"=delta_e,"type"="BUGS")
  }
  if(type=="MNAR_cost"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_costs"=gamma0_c,"MNAR_parameter_costs"=delta_c,"type"="BUGS")
  }
  if(type=="MNAR_cov"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_effects"=beta0_e,"baseline_parameter_costs"=beta0_c,
                            "baseline_parameter_effects"=gamma0_e,"baseline_parameter_costs"=gamma0_c,
                            "covariate_parameter_miss_effects"=gamma_e,"covariate_parameter_miss_costs"=gamma_c,
                            "covariate_parameter_effects"=beta_e,"covariate_parameter_costs"=beta_c,
                            "MNAR_parameter_effects"=delta_e,"MNAR_parameter_costs"=delta_c,"type"="BUGS")
  }
  if(type=="MNAR_eff_cov"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_effects"=beta0_e,
                            "baseline_parameter_effects"=gamma0_e,"covariate_parameter_miss_effects"=gamma_e,
                            "covariate_parameter_effects"=beta_e,"MNAR_parameter_effects"=delta_e,"type"="BUGS")
  }
  if(type=="MNAR_cost_cov"){
    model_output_bugs<-list("summary"=model_sum,"model summary"=modelN1,"mean_effects"=mu_e,"mean_costs"=mu_c,"sd_effects"=s_e,"sd_costs"=s_c,
                            "baseline_parameter_costs"=beta0_c,
                            "baseline_parameter_costs"=gamma0_c,"covariate_parameter_miss_costs"=gamma_c,
                            "covariate_parameter_costs"=beta_c,"MNAR_parameter_costs"=delta_c,"type"="BUGS")
  }
  if(forward==TRUE){model_output_bugs$summary<-NULL}
  if(n.chains==1){model_output_bugs<-model_output_bugs[-1]}
  return(model_output_bugs=model_output_bugs)
  #model output list
}))