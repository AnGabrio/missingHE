#' An internal function to change the hyperprior parameters in the selection model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of selection model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_selection<-function(type,dist_e,dist_c)eval.parent(substitute({
    #For each of the parameters for the models check whether or not the values is provided
    #if provided, check whether or not it can be found in the model script created by the function write_model
    #match the distribution provided to the code and make the change in the model script
    #print out the updated model script to be called again by write_model
  if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      if(is.null(delta.prior.e)==FALSE & grepl("delta_e[t]~",model_string_jags,fixed = T)==TRUE){
        if(length(delta.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_deltae<-delta.prior.e
        prior_deltae_str<-paste("delta_e[t]~dnorm(",prior_deltae[1],",",prior_deltae[2])
        model_string_jags<-gsub("delta_e[t]~dnorm(0,1", prior_deltae_str, model_string_jags,fixed=TRUE)}
    if(is.null(delta.prior.c)==FALSE & grepl("delta_c[t]~",model_string_jags,fixed = T)==TRUE){
      if(length(delta.prior.c)!=2){stop("provide correct hyper prior values")}
      prior_deltac<-delta.prior.c
      prior_deltac_str<-paste("delta_c[t]~dnorm(",prior_deltac[1],",",prior_deltac[2])
      model_string_jags<-gsub("delta_c[t]~dnorm(0,1", prior_deltac_str, model_string_jags,fixed=TRUE)} 
  }
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
    if(pe==1){
      if(is.null(beta0.prior.e)==FALSE & grepl("beta_e[1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e[2]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_mue<-beta0.prior.e
        prior_mue_str<-paste("beta_e[1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e[1]~dnorm(0,0.0001", prior_mue_str, model_string_jags,fixed=TRUE)
        prior_mue_str<-paste("beta_e[2]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e[2]~dnorm(0,0.0001", prior_mue_str, model_string_jags,fixed=TRUE)}
    }else if(pe>1){
      if(is.null(beta0.prior.e)==FALSE & grepl("beta_e[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e[1,2]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_mue<-beta0.prior.e
        prior_mue_str<-paste("beta_e[1,1]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e[1,1]~dnorm(0,0.001", prior_mue_str, model_string_jags,fixed=TRUE)
        prior_mue_str<-paste("beta_e[1,2]~dnorm(",prior_mue[1],",",prior_mue[2])
        model_string_jags<-gsub("beta_e[1,2]~dnorm(0,0.001", prior_mue_str, model_string_jags,fixed=TRUE)}
        if(is.null(beta.prior.e)==FALSE & grepl("beta_e[j,t]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_betae<-beta.prior.e
        prior_betae_str<-paste("beta_e[j,t]~dnorm(",prior_betae[1],",",prior_betae[2])
        model_string_jags<-gsub("beta_e[j,t]~dnorm(0,0.01", prior_betae_str, model_string_jags,fixed=TRUE)}
    }
    if(pc==1){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_c[1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_c[2]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_muc<-beta0.prior.c
        prior_muc_str<-paste("beta_c[1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c[1]~dnorm(0,0.0001", prior_muc_str, model_string_jags,fixed=TRUE)
        prior_muc_str<-paste("beta_c[2]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c[2]~dnorm(0,0.0001", prior_muc_str, model_string_jags,fixed=TRUE)}
    }else if(pc>1){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_e[1,1]~",model_string_jags,fixed = T)==TRUE & grepl("beta_e[1,2]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta0.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_muc<-beta0.prior.c
        prior_muc_str<-paste("beta_c[1,1]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c[1,1]~dnorm(0,0.001", prior_muc_str, model_string_jags,fixed=TRUE)
        prior_muc_str<-paste("beta_c[1,2]~dnorm(",prior_muc[1],",",prior_muc[2])
        model_string_jags<-gsub("beta_c[1,2]~dnorm(0,0.001", prior_muc_str, model_string_jags,fixed=TRUE)}
        if(is.null(beta.prior.c)==FALSE & grepl("beta_c[j,t]~",model_string_jags,fixed = T)==TRUE){
        if(length(beta.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_betac<-beta.prior.c
        prior_betac_str<-paste("beta_c[j,t]~dnorm(",prior_betac[1],",",prior_betac[2])
        model_string_jags<-gsub("beta_c[j,t]~dnorm(0,0.01", prior_betac_str, model_string_jags,fixed=TRUE)}
    }
  if(dist_e=="norm"){
  if(is.null(sigma.prior.e)==FALSE & grepl("ls_e[t]~",model_string_jags,fixed = T)==TRUE){
    if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
    prior_alphae<-sigma.prior.e
    prior_alphae_str<-paste("ls_e[t]~dunif(",prior_alphae[1],",",prior_alphae[2])
    model_string_jags<-gsub("ls_e[t]~dunif(-5,10", prior_alphae_str, model_string_jags,fixed=TRUE)}
  }else if(dist_e=="beta"){
      if(is.null(sigma.prior.e)==FALSE & grepl("s_e[t]~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.e)!=2){stop("provide correct hyper prior values")}
        prior_alphae<-sigma.prior.e
        prior_alphae_str<-paste("s_e[t]~dunif(",prior_alphae[1],",",prior_alphae[2])
        model_string_jags<-gsub("s_e[t]~dunif(0,sqrt(mu_e[t]*(1-mu_e[t]))", prior_alphae_str, model_string_jags,fixed=TRUE)}
  }
  if(dist_c=="norm"){
    if(is.null(sigma.prior.c)==FALSE & grepl("ls_c[t]~",model_string_jags,fixed = T)==TRUE){
      if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
      prior_alphac<-sigma.prior.c
      prior_alphac_str<-paste("ls_c[t]~dunif(",prior_alphac[1],",",prior_alphac[2])
      model_string_jags<-gsub("ls_c[t]~dunif(-5,10", prior_alphac_str, model_string_jags,fixed=TRUE)}
  }else if(dist_c=="gamma"){
      if(is.null(sigma.prior.c)==FALSE & grepl("s_c[t]~",model_string_jags,fixed = T)==TRUE){
        if(length(sigma.prior.c)!=2){stop("provide correct hyper prior values")}
        prior_alphac<-sigma.prior.c
        prior_alphac_str<-paste("s_c[t]~dunif(",prior_alphac[1],",",prior_alphac[2])
        model_string_jags<-gsub("s_c[t]~dunif(0,1000", prior_alphac_str, model_string_jags,fixed=TRUE)}
  }
    #correlation only for joint models
      if(exists("rho.prior")==TRUE){
      if(is.null(rho.prior)==FALSE & grepl("rho",model_string_jags,fixed = T)==TRUE){
        if(length(rho.prior)!=2){stop("provide correct hyper prior values")}
          prior_rho<-rho.prior
          prior_rho_str<-paste("rho[t]~dnorm(",prior_rho[1],",",prior_rho[2])
          model_string_jags<-gsub("rho[t]~dnorm(0,0.001", prior_rho_str, model_string_jags,fixed=TRUE)
        }
      }
  #save updated model file after prior change
  model_string_prior<-model_string_jags
    return(model_string_prior)
}))