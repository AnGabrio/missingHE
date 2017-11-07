#' An internal function to change the hyperprior parameters in the hurdle model provided by the user depending on the type of
#' structural value mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of hurdle model selected according 
#' to the type of structural value mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions hurdle models
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#' and Structural At Random (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_hurdle<-function(type,dist_e,dist_c)eval.parent(substitute({
    #For each of the parameters for the models check whether or not the values is provided
    #if provided, check whether or not it can be found in the model script created by the function write_model
    #match the distribution provided to the code and make the change in the model script
    #print out the updated model script to be called again by write_model
  if(type=="SCAR"){
    if(is.null(se)==FALSE){
      if(is.null(gamma0.prior.e)==FALSE & grepl("gamma_e[t]~",model_string_jags,fixed = T)==TRUE){
        if(length(gamma0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_pe<-gamma0.prior.e
        prior_pe_str<-paste("gamma_e[t]~dlogis(",prior_pe[1],",",prior_pe[2])
        model_string_jags<-gsub("gamma_e[t]~dlogis(0,1", prior_pe_str, model_string_jags,fixed=TRUE)}
    }
    if(is.null(sc)==FALSE){
    if(is.null(gamma0.prior.c)==FALSE & grepl("gamma_c[t]~",model_string_jags,fixed = T)==TRUE){
      if(length(gamma0.prior.c)!=2){stop("provide correct hyper prior values")}
      prior_pc<-gamma0.prior.c
      prior_pc_str<-paste("gamma_c[t]~dlogis(",prior_pc[1],",",prior_pc[2])
      model_string_jags<-gsub("gamma_c[t]~dlogis(0,1", prior_pc_str, model_string_jags,fixed=TRUE)} 
    }
  }
  if(type=="SAR"){
    if(is.null(se)==FALSE){
      if(is.null(ze)==FALSE){
        if(ze==1){
      if(is.null(gamma0.prior.e)==FALSE & grepl("gamma_e[1]~",model_string_jags,fixed = T)==TRUE & grepl("gamma_e[2]~",model_string_jags,fixed = T)==TRUE){
        if(length(gamma0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_pe<-gamma0.prior.e
        prior_pe_str<-paste("gamma_e[1]~dlogis(",prior_pe[1],",",prior_pe[2])
        model_string_jags<-gsub("gamma_e[1]~dlogis(0,1", prior_pe_str, model_string_jags,fixed=TRUE)
        prior_pe_str<-paste("gamma_e[2]~dlogis(",prior_pe[1],",",prior_pe[2])
        model_string_jags<-gsub("gamma_e[2]~dlogis(0,1", prior_pe_str, model_string_jags,fixed=TRUE)}
        }else if(ze>1){
          if(is.null(gamma0.prior.e)==FALSE & grepl("gamma_e[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("gamma_e[1,2]~",model_string_jags,fixed = T)==TRUE){
            if(length(gamma0.prior.e)!=2){stop("provide correct hyper prior values")}
            prior_pe<-gamma0.prior.e
            prior_pe_str<-paste("gamma_e[1,1]~dlogis(",prior_pe[1],",",prior_pe[2])
            model_string_jags<-gsub("gamma_e[1,1]~dlogis(0,1", prior_pe_str, model_string_jags,fixed=TRUE)
            prior_pe_str<-paste("gamma_e[1,2]~dlogis(",prior_pe[1],",",prior_pe[2])
            model_string_jags<-gsub("gamma_e[1,2]~dlogis(0,1", prior_pe_str, model_string_jags,fixed=TRUE)}
          if(is.null(gamma.prior.e)==FALSE & grepl("gamma_e[j,t]~",model_string_jags,fixed = T)==TRUE){
            if(length(gamma.prior.e)!=2){stop("provide correct hyper prior values")}
            prior_gammae<-gamma.prior.e
            prior_gammae_str<-paste("gamma_e[j,t]~dnorm(",prior_gammae[1],",",prior_gammae[2])
            model_string_jags<-gsub("gamma_e[j,t]~dnorm(0,0.01", prior_gammae_str, model_string_jags,fixed=TRUE)}
        }
      }
    }
    if(is.null(sc)==FALSE){
      if(is.null(zc)==FALSE){
        if(zc==1){
          if(is.null(gamma0.prior.c)==FALSE & grepl("gamma_c[1]~",model_string_jags,fixed = T)==TRUE & grepl("gamma_c[2]~",model_string_jags,fixed = T)==TRUE){
            if(length(gamma0.prior.c)!=2){stop("provide correct hyper prior values")}
            prior_pc<-gamma0.prior.c
            prior_pc_str<-paste("gamma_c[1]~dlogis(",prior_pc[1],",",prior_pc[2])
            model_string_jags<-gsub("gamma_c[1]~dlogis(0,1", prior_pc_str, model_string_jags,fixed=TRUE)
            prior_pc_str<-paste("gamma_c[2]~dlogis(",prior_pc[1],",",prior_pc[2])
            model_string_jags<-gsub("gamma_c[2]~dlogis(0,1", prior_pc_str, model_string_jags,fixed=TRUE)}
        }else if(zc>1){
          if(is.null(gamma0.prior.c)==FALSE & grepl("gamma_c[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("gamma_c[1,2]~",model_string_jags,fixed = T)==TRUE){
            if(length(gamma0.prior.c)!=2){stop("provide correct hyper prior values")}
            prior_pc<-gamma0.prior.c
            prior_pc_str<-paste("gamma_c[1,1]~dlogis(",prior_pc[1],",",prior_pc[2])
            model_string_jags<-gsub("gamma_c[1,1]~dlogis(0,1", prior_pc_str, model_string_jags,fixed=TRUE)
            prior_pc_str<-paste("gamma_c[1,2]~dlogis(",prior_pc[1],",",prior_pc[2])
            model_string_jags<-gsub("gamma_c[1,2]~dlogis(0,1", prior_pc_str, model_string_jags,fixed=TRUE)}
          if(is.null(gamma.prior.c)==FALSE & grepl("gamma_c[j,t]~",model_string_jags,fixed = T)==TRUE){
            if(length(gamma.prior.c)!=2){stop("provide correct hyper prior values")}
            prior_gammac<-gamma.prior.c
            prior_gammac_str<-paste("gamma_c[j,t]~dnorm(",prior_gammac[1],",",prior_gammac[2])
            model_string_jags<-gsub("gamma_c[j,t]~dnorm(0,0.01", prior_gammac_str, model_string_jags,fixed=TRUE)}
        }
      }
    }
  }
  if(type=="SCAR"|type=="SAR"){
    if(pe==1){
      if(is.null(beta0.prior.e)==FALSE & grepl("beta_e1[1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e2[1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_mue<-beta0.prior.e
        prior_mue_str<-paste("beta_e1[1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e1[1]~dnorm(0,0.0001", prior_mue_str, model_string_jags,fixed=TRUE)
        prior_mue_str<-paste("beta_e2[1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e2[1]~dnorm(0,0.0001", prior_mue_str, model_string_jags,fixed=TRUE)}
    }else if(pe>1){
      if(is.null(beta0.prior.e)==FALSE & grepl("beta_e1[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e2[1,1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_mue<-beta0.prior.e
        prior_mue_str<-paste("beta_e1[1,1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e1[1,1]~dnorm(0,0.001", prior_mue_str, model_string_jags,fixed=TRUE)
        prior_mue_str<-paste("beta_e2[1,1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e2[1,1]~dnorm(0,0.001", prior_mue_str, model_string_jags,fixed=TRUE)}
        if(is.null(beta.prior.e)==FALSE & grepl("beta_e1[j,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e2[j,1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_betae<-beta.prior.e
        prior_betae_str<-paste("beta_e1[j,1]~dnorm(",prior_betae[1],",",prior_betae[2])
        model_string_jags<-gsub("beta_e1[j,1]~dnorm(0,0.01", prior_betae_str, model_string_jags,fixed=TRUE)
        prior_betae_str<-paste("beta_e2[j,1]~dnorm(",prior_betae[1],",",prior_betae[2])
        model_string_jags<-gsub("beta_e2[j,1]~dnorm(0,0.01", prior_betae_str, model_string_jags,fixed=TRUE)}
    }
    if(pc==1){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_c1[1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_c2[1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_muc<-beta0.prior.c
        prior_muc_str<-paste("beta_c1[1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c1[1]~dnorm(0,0.0001", prior_muc_str, model_string_jags,fixed=TRUE)
        prior_muc_str<-paste("beta_c2[1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c2[1]~dnorm(0,0.0001", prior_muc_str, model_string_jags,fixed=TRUE)}
    }else if(pc>1){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_e1[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e2[1,1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_muc<-beta0.prior.c
        prior_muc_str<-paste("beta_c1[1,1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c1[1,1]~dnorm(0,0.001", prior_muc_str, model_string_jags,fixed=TRUE)
        prior_muc_str<-paste("beta_c2[1,1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c2[1,1]~dnorm(0,0.001", prior_muc_str, model_string_jags,fixed=TRUE)}
        if(is.null(beta.prior.c)==FALSE & grepl("beta_c1[j,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_c2[j,1]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_betac<-beta.prior.c
        prior_betac_str<-paste("beta_c1[j,1]~dnorm(",prior_betac[1],",",prior_betac[2])
        model_string_jags<-gsub("beta_c1[j,1]~dnorm(0,0.01", prior_betac_str, model_string_jags,fixed=TRUE)
        prior_betac_str<-paste("beta_c2[j,1]~dnorm(",prior_betac[1],",",prior_betac[2])
        model_string_jags<-gsub("beta_c2[j,1]~dnorm(0,0.01", prior_betac_str, model_string_jags,fixed=TRUE)}
    }
  }
  if(dist_e=="norm"){
  if(is.null(se)==FALSE){
  if(is.null(sigma.prior.e)==FALSE & grepl("ls_e1[1]~",model_string_jags,fixed = T)==TRUE & grepl("ls_e2[1]~",model_string_jags,fixed = T)==TRUE){
    if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
    prior_alphae<-sigma.prior.e
    prior_alphae_str<-paste("ls_e1[1]~dunif(",prior_alphae[1],",",prior_alphae[2])
    model_string_jags<-gsub("ls_e1[1]~dunif(-5,10", prior_alphae_str, model_string_jags,fixed=TRUE)
    prior_alphae_str<-paste("ls_e2[1]~dunif(",prior_alphae[1],",",prior_alphae[2])
    model_string_jags<-gsub("ls_e2[1]~dunif(-5,10", prior_alphae_str, model_string_jags,fixed=TRUE)}
  }else if(is.null(se)==TRUE){
    if(is.null(sigma.prior.e)==FALSE & grepl("ls_e1~",model_string_jags,fixed = T)==TRUE & grepl("ls_e2~",model_string_jags,fixed = T)==TRUE){
      if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
      prior_alphae<-sigma.prior.e
      prior_alphae_str<-paste("ls_e1~dunif(",prior_alphae[1],",",prior_alphae[2])
      model_string_jags<-gsub("ls_e1~dunif(-5,10", prior_alphae_str, model_string_jags,fixed=TRUE)
      prior_alphae_str<-paste("ls_e2~dunif(",prior_alphae[1],",",prior_alphae[2])
      model_string_jags<-gsub("ls_e2~dunif(-5,10", prior_alphae_str, model_string_jags,fixed=TRUE)}
     }
  }else if(dist_e=="beta"){
    if(is.null(se)==FALSE){
      if(is.null(sigma.prior.e)==FALSE & grepl("s_e1[1]~",model_string_jags,fixed = T)==TRUE & grepl("s_e2[1]~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_alphae<-sigma.prior.e
        prior_alphae_str<-paste("s_e1[1]~dunif(",prior_alphae[1],",",prior_alphae[2])
        model_string_jags<-gsub("s_e1[1]~dunif(0,sqrt(nu_e[1]*(1-nu_e[1]))", prior_alphae_str, model_string_jags,fixed=TRUE)
        prior_alphae_str<-paste("s_e2[1]~dunif(",prior_alphae[1],",",prior_alphae[2])
        model_string_jags<-gsub("s_e2[1]~dunif(0,sqrt(nu_e[2]*(1-nu_e[2]))", prior_alphae_str, model_string_jags,fixed=TRUE)}
    }else if(is.null(se)==TRUE){
      if(is.null(sigma.prior.e)==FALSE & grepl("s_e1~",model_string_jags,fixed = T)==TRUE & grepl("s_e2~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_alphae<-sigma.prior.e
        prior_alphae_str<-paste("s_e1~dunif(",prior_alphae[1],",",prior_alphae[2])
        model_string_jags<-gsub("s_e1~dunif(0,sqrt(nu_e[1]*(1-nu_e[1]))", prior_alphae_str, model_string_jags,fixed=TRUE)
        prior_alphae_str<-paste("s_e2~dunif(",prior_alphae[1],",",prior_alphae[2])
        model_string_jags<-gsub("s_e2~dunif(0,sqrt(nu_e[2]*(1-nu_e[2]))", prior_alphae_str, model_string_jags,fixed=TRUE)}
    }
  }
  if(dist_c=="norm"){
  if(is.null(sc)==FALSE){
    if(is.null(sigma.prior.c)==FALSE & grepl("ls_c1[1]~",model_string_jags,fixed = T)==TRUE & grepl("ls_c2[1]~",model_string_jags,fixed = T)==TRUE){
      if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
      prior_alphac<-sigma.prior.c
      prior_alphac_str<-paste("ls_c1[1]~dunif(",prior_alphac[1],",",prior_alphac[2])
      model_string_jags<-gsub("ls_c1[1]~dunif(-5,10", prior_alphac_str, model_string_jags,fixed=TRUE)
      prior_alphac_str<-paste("ls_c2[1]~dunif(",prior_alphac[1],",",prior_alphac[2])
      model_string_jags<-gsub("ls_c2[1]~dunif(-5,10", prior_alphac_str, model_string_jags,fixed=TRUE)}
  }else if(is.null(sc)==TRUE){
    if(is.null(sigma.prior.c)==FALSE & grepl("ls_c1~",model_string_jags,fixed = T)==TRUE & grepl("ls_c2~",model_string_jags,fixed = T)==TRUE){
      if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
      prior_alphac<-sigma.prior.c
      prior_alphac_str<-paste("ls_c1~dunif(",prior_alphac[1],",",prior_alphac[2])
      model_string_jags<-gsub("ls_c1~dunif(-5,10", prior_alphac_str, model_string_jags,fixed=TRUE)
      prior_alphac_str<-paste("ls_c2~dunif(",prior_alphac[1],",",prior_alphac[2])
      model_string_jags<-gsub("ls_c2~dunif(-5,10", prior_alphac_str, model_string_jags,fixed=TRUE)}
    }
  }else if(dist_c=="gamma"){
    if(is.null(sc)==FALSE){
      if(is.null(sigma.prior.c)==FALSE & grepl("s_c1[1]~",model_string_jags,fixed = T)==TRUE & grepl("s_c2[1]~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_alphac<-sigma.prior.c
        prior_alphac_str<-paste("s_c1[1]~dunif(",prior_alphac[1],",",prior_alphac[2])
        model_string_jags<-gsub("s_c1[1]~dunif(0,1000", prior_alphac_str, model_string_jags,fixed=TRUE)
        prior_alphac_str<-paste("s_c2[1]~dunif(",prior_alphac[1],",",prior_alphac[2])
        model_string_jags<-gsub("s_c2[1]~dunif(0,1000", prior_alphac_str, model_string_jags,fixed=TRUE)}
    }else if(is.null(sc)==TRUE){
      if(is.null(sigma.prior.c)==FALSE & grepl("s_c1~",model_string_jags,fixed = T)==TRUE & grepl("s_c2~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_alphac<-sigma.prior.c
        prior_alphac_str<-paste("s_c1~dunif(",prior_alphac[1],",",prior_alphac[2])
        model_string_jags<-gsub("s_c1~dunif(0,1000", prior_alphac_str, model_string_jags,fixed=TRUE)
        prior_alphac_str<-paste("s_c2~dunif(",prior_alphac[1],",",prior_alphac[2])
        model_string_jags<-gsub("s_c2~dunif(0,1000", prior_alphac_str, model_string_jags,fixed=TRUE)}
    }
  }
    #correlation only for joint models
      if(exists("rho.prior")==TRUE){
      if(is.null(rho.prior)==FALSE & grepl("rho",model_string_jags,fixed = T)==TRUE){
        if(length(rho.prior)!=2){stop("provide correct hyper prior values")}
        if(is.null(se)==FALSE & is.null(sc)==TRUE){
          prior_rho<-rho.prior
          prior_rho_str<-paste("rho[1]~dnorm(",prior_rho[1],",",prior_rho[2])
          model_string_jags<-gsub("rho[1]~dnorm(0,0.001", prior_rho_str, model_string_jags,fixed=TRUE)
          prior_rho_str<-paste("rho[2]~dnorm(",prior_rho[1],",",prior_rho[2])
          model_string_jags<-gsub("rho[2]~dnorm(0,0.001", prior_rho_str, model_string_jags,fixed=TRUE)}
        if(is.null(sc)==FALSE){
          prior_rho<-rho.prior
          prior_rho_str<-paste("rho1[1]~dnorm(",prior_rho[1],",",prior_rho[2])
          model_string_jags<-gsub("rho1[1]~dnorm(0,0.001", prior_rho_str, model_string_jags,fixed=TRUE)
          prior_rho_str<-paste("rho2[1]~dnorm(",prior_rho[1],",",prior_rho[2])
          model_string_jags<-gsub("rho2[1]~dnorm(0,0.001", prior_rho_str, model_string_jags,fixed=TRUE)}
        }
      }
  #save updated model file after prior change
  model_string_prior<-model_string_jags
    return(model_string_prior)
}))