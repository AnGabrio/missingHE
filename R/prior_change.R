#' An internal function to change the hyperprior parameters in the model provided by the user depending on the type of
#' missing data mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of model selected according 
#' to the type of missingness mechanism and distributions for the outcome variables assumed.
#' @keywords prior distributions
#' @param type type of missingness mechanism assumed. Choices are Missing Completely At Random (MCAR),
#'  Missing At Random (MAR), Missing Not At Random (MNAR). For 'MNAR' alternative versions are available 
#'  depending on whether the mechanism for only one variable is condiered, that is for the effects (MNAR_eff)
#'  or the costs (MNAR_cost), or also covariates are included either for the effects (MNAR_eff_cov),
#'  the costs (MNAR_cost), or both (MNAR_cov). For a complete list of all available hyper parameters 
#'  and types of models see the manual.
#' @param dist_e assumed distribution for the effects in the MoA
#' @param dist_c assumed distribution for the costs in the MoA
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_change<-function(type,dist_e,dist_c)eval.parent(substitute({
    #For each of the parameters for the models check whether or not the values is provided
    #if provided, check whether or not it can be found in the model script created by the function write_model
    #match the distribution provided to the BUGS/STAN code and make the change in the model script
    #print out the updated model script to be called again by write_model
  #MoM parameters for all types of  MoA
  #baseline parameters
  if(type=="MCAR"|type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
    if(is.null(gamma0.prior.e)==FALSE & grepl("gamma0_e",model_string_jags)==TRUE){
      prior_gamma0e<-gamma0.prior.e
      prior_gamma0e_str<-paste("gamma0_e~dlogis(",prior_gamma0e[1],",",prior_gamma0e[2])
      model_string_jags<-gsub("gamma0_e~dlogis(0,1", prior_gamma0e_str, model_string_jags,fixed=TRUE)}
    if(is.null(gamma0.prior.c)==FALSE & grepl("gamma0_c",model_string_jags)==TRUE){
      prior_gamma0c<-gamma0.prior.c
      prior_gamma0c_str<-paste("gamma0_c~dlogis(",prior_gamma0c[1],",",prior_gamma0c[2])
      model_string_jags<-gsub("gamma0_c~dlogis(0,1", prior_gamma0c_str, model_string_jags,fixed=TRUE)}  
  }
  #parameters specific to MNAR without covariates
  #MNAR parameters
  if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cov"|type=="MNAR_eff_cov"){
    if(is.null(delta.prior.e)==FALSE & grepl("delta_e",model_string_jags)==TRUE){
      prior_delta_e<-delta.prior.e
      prior_delta_e_str<-paste("delta_e~dnorm(",prior_delta_e[1],",",prior_delta_e[2])
      model_string_jags<-gsub("delta_e~dnorm(0,1", prior_delta_e_str, model_string_jags,fixed=TRUE)}
  }
  if(type=="MNAR"|type=="MNAR_cost"|type=="MNAR_cov"|type=="MNAR_cost_cov"){
    if(is.null(delta.prior.c)==FALSE & grepl("delta_c",model_string_jags)==TRUE){
      prior_delta_c<-delta.prior.c
      prior_delta_c_str<-paste("delta_c~dnorm(",prior_delta_c[1],",",prior_delta_c[2])
      model_string_jags<-gsub("delta_c~dnorm(0,1", prior_delta_c_str, model_string_jags,fixed=TRUE)}
  }
  #covariate parameters
  if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
    if(is.null(beta.prior.e)==FALSE & grepl("beta_e",model_string_jags)==TRUE){
      prior_mu_e_betas<-beta.prior.e
      prior_mu_e_betas_str<-paste("beta_e[j,t]~dnorm(",prior_mu_e_betas[1],",",prior_mu_e_betas[2])
      model_string_jags<-gsub("beta_e[j,t]~dnorm(0,0.01", prior_mu_e_betas_str, model_string_jags,fixed=TRUE)}
    if(is.null(beta.prior.c)==FALSE & grepl("beta_c",model_string_jags)==TRUE){
      prior_mu_c_betas<-beta.prior.c
      prior_mu_c_betas_str<-paste("beta_c[j,t]~dnorm(",prior_mu_c_betas[1],",",prior_mu_c_betas[2])
      model_string_jags<-gsub("beta_c[j,t]~dnorm(0,0.01", prior_mu_c_betas_str, model_string_jags,fixed=TRUE)}
    if(is.null(gamma.prior.e)==FALSE & grepl("gamma_e",model_string_jags)==TRUE){
      prior_gammae<-gamma.prior.e
      prior_gammae_str<-paste("gamma_e[j]~dnorm(",prior_gammae[1],",",prior_gammae[2])
      model_string_jags<-gsub("gamma_e[j]~dnorm(0,1", prior_gammae_str, model_string_jags,fixed=TRUE)}
    if(is.null(gamma.prior.c)==FALSE & grepl("gamma_c",model_string_jags)==TRUE){
      prior_gammac<-gamma.prior.c
      prior_gammac_str<-paste("gamma_c[j]~dnorm(",prior_gammac[1],",",prior_gammac[2])
      model_string_jags<-gsub("gamma_c[j]~dnorm(0,1", prior_gammac_str, model_string_jags,fixed=TRUE)}
  }
  #distribution-specific parameters of the MoA
  if(dist_e=="norm"){
    #standard deviations
    if(is.null(alpha.prior.e)==FALSE & grepl("ls_e",model_string_jags)==TRUE & grepl("ls_e_t",model_string_jags)==FALSE ){
      prior_ls_e<-alpha.prior.e
      prior_ls_e_str<-paste("ls_e[t]~dnorm(",prior_ls_e[1],",",prior_ls_e[2])
      model_string_jags<-gsub("ls_e[t]~dunif(-5,10", prior_ls_e_str, model_string_jags,fixed=TRUE)}
    if(is.null(alpha.prior.e)==FALSE & grepl("ls_e_t",model_string_jags)==TRUE){
      prior_ls_e<-alpha.prior.e
      prior_ls_e_str<-paste("ls_e_t[t]~dnorm(",prior_ls_e[1],",",prior_ls_e[2])
      model_string_jags<-gsub("ls_e_t[t]~dunif(-5,10", prior_ls_e_str, model_string_jags,fixed=TRUE)}
    if(type=="MCAR"|type=="MNAR"){
      #parameters in common between MCAR and MNAR
      #means
      if(is.null(mu.prior.e)==FALSE & grepl("mu_e",model_string_jags)==TRUE & grepl("mu_e_t",model_string_jags)==FALSE ){
        prior_mu_e<-mu.prior.e
        prior_mu_e_str<-paste("mu_e[t]~dnorm(",prior_mu_e[1],",",prior_mu_e[2])
        model_string_jags<-gsub("mu_e[t]~dnorm(0,0.0001", prior_mu_e_str, model_string_jags,fixed=TRUE)}
      if(is.null(mu.prior.e)==FALSE & grepl("mu_e_t",model_string_jags)==TRUE){
        prior_mu_e<-mu.prior.e
        prior_mu_e_str<-paste("mu_e_t[t]~dnorm(",prior_mu_e[1],",",prior_mu_e[2])
        model_string_jags<-gsub("mu_e_t[t]~dnorm(0,0.0001", prior_mu_e_str, model_string_jags,fixed=TRUE)}
    }
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      if(is.null(beta0.prior.e)==FALSE & grepl("beta_e",model_string_jags)==TRUE){
        prior_beta0_e<-beta0.prior.e
        prior_beta0_e_str1<-paste("beta_e[1,1]~dnorm(",prior_beta0_e[1],",",prior_beta0_e[2])
        prior_beta0_e_str2<-paste("beta_e[1,2]~dnorm(",prior_beta0_e[1],",",prior_beta0_e[2])
        model_string_jags<-gsub("beta_e[1,1]~dnorm(0,0.001", prior_beta0_e_str1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("beta_e[1,2]~dnorm(0,0.001", prior_beta0_e_str2, model_string_jags,fixed=TRUE)}
    }
     }
  if(dist_c=="norm"){
    if(is.null(alpha.prior.c)==FALSE & grepl("ls_c",model_string_jags)==TRUE & grepl("ls_c_t",model_string_jags)==FALSE){
      prior_ls_c<-alpha.prior.c
      prior_ls_c_str<-paste("ls_c[t]~dnorm(",prior_ls_c[1],",",prior_ls_c[2])
      model_string_jags<-gsub("ls_c[t]~dunif(-5,10", prior_ls_c_str, model_string_jags,fixed=TRUE)}
    if(is.null(alpha.prior.c)==FALSE & grepl("ls_c_t",model_string_jags)==TRUE){
      prior_ls_c<-alpha.prior.c
      prior_ls_c_str<-paste("ls_c_t[t]~dnorm(",prior_ls_c[1],",",prior_ls_c[2])
      model_string_jags<-gsub("ls_c_t[t]~dunif(-5,10", prior_ls_c_str, model_string_jags,fixed=TRUE)}
    if(type=="MCAR"|type=="MNAR"){
      #parameters in common between MCAR and MNAR
      #means
    if(is.null(mu.prior.c)==FALSE & grepl("mu_c",model_string_jags)==TRUE & grepl("mu_c_t",model_string_jags)==FALSE){
      prior_mu_c<-mu.prior.c
      prior_mu_c_str<-paste("mu_c[t]~dnorm(",prior_mu_c[1],",",prior_mu_c[2])
      model_string_jags<-gsub("mu_c[t]~dnorm(0,0.0001", prior_mu_c_str, model_string_jags,fixed=TRUE)}
    if(is.null(mu.prior.c)==FALSE & grepl("mu_c_t",model_string_jags)==TRUE){
      prior_mu_c<-mu.prior.c
      prior_mu_c_str<-paste("mu_c_t[t]~dnorm(",prior_mu_c[1],",",prior_mu_c[2])
      model_string_jags<-gsub("mu_c_t[t]~dnorm(0,0.0001", prior_mu_c_str, model_string_jags,fixed=TRUE)}
    }
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_c",model_string_jags)==TRUE){
        prior_beta0_c<-beta0.prior.c
        prior_beta0_c_str1<-paste("beta_c[1,1]~dnorm(",prior_beta0_c[1],",",prior_beta0_c[2])
        prior_beta0_c_str2<-paste("beta_c[1,2]~dnorm(",prior_beta0_c[1],",",prior_beta0_c[2])
        model_string_jags<-gsub("beta_c[1,1]~dnorm(0,0.001", prior_beta0_c_str1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("beta_c[1,2]~dnorm(0,0.001", prior_beta0_c_str2, model_string_jags,fixed=TRUE)}
    }
  }
    #correlation only for joitn models
      if(exists("theta.prior")==TRUE){
      if(is.null(theta.prior)==FALSE & grepl("theta",model_string_jags)==TRUE){
        prior_theta<-theta.prior
        prior_theta_str<-paste("theta[t]~dnorm(",prior_theta[1],",",prior_theta[2])
        model_string_jags<-gsub("theta[t]~dnorm(0,0.0001", prior_theta_str, model_string_jags,fixed=TRUE)}
      model_string_prior<-model_string_jags
      }
   if(dist_e=="beta"){
     #sd
     if(is.null(alpha.prior.e)==FALSE & grepl("s_e",model_string_jags)==TRUE){
       prior_scale_e<-alpha.prior.e
       prior_scale_e_str<-paste("s_e[t]~dunif(",prior_scale_e[1],",",prior_scale_e[2])
       model_string_jags<-gsub("s_e[t]~dunif(0,se.limit[t]", prior_scale_e_str, model_string_jags,fixed=TRUE)}
     if(type=="MCAR"|type=="MNAR"){
       #parameters in common between MCAR and MNAR
       #means
       if(is.null(mu.prior.e)==FALSE & grepl("mu_e",model_string_jags)==TRUE & any(transf=="default")==TRUE){
         prior_mu_e<-mu.prior.e
         prior_mu_e_str<-paste("mu_e[t]~dunif(",prior_mu_e[1],",",prior_mu_e[2])
         model_string_jags<-gsub("mu_e[t]~dunif(0,1", prior_mu_e_str, model_string_jags,fixed=TRUE)}
       if(is.null(mu.prior.e)==FALSE & grepl("nu_e",model_string_jags)==TRUE & any(transf=="logit")==TRUE){
         prior_mu_e<-mu.prior.e
         prior_mu_e_str<-paste("nu_e[t]~dnorm(",prior_mu_e[1],",",prior_mu_e[2])
         model_string_jags<-gsub("nu_e[t]~dnorm(0,0.001", prior_mu_e_str, model_string_jags,fixed=TRUE)}
     }
     if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
       if(is.null(beta0.prior.e)==FALSE & grepl("beta_e",model_string_jags)==TRUE & any(transf=="default")==TRUE){
         prior_beta0_e<-beta0.prior.e
         prior_beta0_e_str1<-paste("beta_e[1,1]~dunif(",prior_beta0_e[1],",",prior_beta0_e[2])
         prior_beta0_e_str2<-paste("beta_e[1,2]~dunif(",prior_beta0_e[1],",",prior_beta0_e[2])
         model_string_jags<-gsub("beta_e[1,1]~dunif(0,1", prior_beta0_e_str1, model_string_jags,fixed=TRUE)
         model_string_jags<-gsub("beta_e[1,2]~dunif(0,1", prior_beta0_e_str2, model_string_jags,fixed=TRUE)}
       if(is.null(beta0.prior.e)==FALSE & grepl("beta_e",model_string_jags)==TRUE & any(transf=="logit")==TRUE){
         prior_beta0_e<-beta0.prior.e
         prior_beta0_e_str1<-paste("beta_e[1,1]~dnorm(",prior_beta0_e[1],",",prior_beta0_e[2])
         prior_beta0_e_str2<-paste("beta_e[1,2]~dnorm(",prior_beta0_e[1],",",prior_beta0_e[2])
         model_string_jags<-gsub("beta_e[1,1]~dnorm(0,0.001", prior_beta0_e_str1, model_string_jags,fixed=TRUE)
         model_string_jags<-gsub("beta_e[1,2]~dnorm(0,0.001", prior_beta0_e_str2, model_string_jags,fixed=TRUE)}
     }
   }
  if(dist_c=="gamma"){
    #standard deviations
    if(is.null(alpha.prior.c)==FALSE & grepl("s_c",model_string_jags)==TRUE){
      prior_s_c<-alpha.prior.c
      prior_s_c_str<-paste("s_c[t]~dt(",prior_s_c[1],",",prior_s_c[2])
      model_string_jags<-gsub("s_c[t]~dt(0,0.16", prior_s_c_str, model_string_jags,fixed=TRUE)}
    if(type=="MCAR"|type=="MNAR"){
      #parameters in common between MCAR and MNAR
      #means
      if(is.null(mu.prior.c)==FALSE & grepl("mu_c",model_string_jags)==TRUE & any(transf=="default")==TRUE){
        prior_mu_c<-mu.prior.c
        prior_mu_c_str<-paste("mu_c[t]~dunif(",prior_mu_c[1],",",prior_mu_c[2])
        model_string_jags<-gsub("mu_c[t]~dunif(0,10000", prior_mu_c_str, model_string_jags,fixed=TRUE)}
      if(is.null(mu.prior.c)==FALSE & grepl("nu_c",model_string_jags)==TRUE & any(transf=="log")==TRUE){
        prior_mu_c<-mu.prior.c
        prior_mu_c_str<-paste("nu_c[t]~dnorm(",prior_mu_c[1],",",prior_mu_c[2])
        model_string_jags<-gsub("nu_c[t]~dnorm(0,0.001", prior_mu_c_str, model_string_jags,fixed=TRUE)}
    }
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_c",model_string_jags)==TRUE & any(transf=="default")==TRUE){
        prior_beta0_c<-beta0.prior.e
        prior_beta0_c_str1<-paste("beta_c[1,1]~dunif(",prior_beta0_c[1],",",prior_beta0_c[2])
        prior_beta0_c_str2<-paste("beta_c[1,2]~dunif(",prior_beta0_c[1],",",prior_beta0_c[2])
        model_string_jags<-gsub("beta_c[1,1]~dunif(0,1", prior_beta0_c_str1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("beta_c[1,2]~dunif(0,1", prior_beta0_c_str2, model_string_jags,fixed=TRUE)}
      if(is.null(beta0.prior.c)==FALSE & grepl("beta_c",model_string_jags)==TRUE & any(transf=="log")==TRUE){
        prior_beta0_c<-beta0.prior.c
        prior_beta0_c_str1<-paste("beta_c[1,1]~dnorm(",prior_beta0_c[1],",",prior_beta0_c[2])
        prior_beta0_c_str2<-paste("beta_c[1,2]~dnorm(",prior_beta0_c[1],",",prior_beta0_c[2])
        model_string_jags<-gsub("beta_c[1,1]~dnorm(0,0.001", prior_beta0_c_str1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("beta_c[1,2]~dnorm(0,0.001", prior_beta0_c_str2, model_string_jags,fixed=TRUE)}
    }
  }
  #save updated model file after prior change
  model_string_prior<-model_string_jags
    return(model_string_prior)
}))