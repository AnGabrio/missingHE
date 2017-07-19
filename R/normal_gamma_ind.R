#' An internal function to write BUGS code for normal-gamma model with natural-scaled variables (independence)
#'
#' This function writes in the current WD a txt file 
#' for the bivariate independent normal-gamma model with natural-scaled outcome variables.
#' @keywords JAGS BUGS models
#' @param type Type of missingness mechanism assumed. Choices are Missing Completely At Random (MCAR),
#'  Missing At Random (MAR), Missing Not At Random (MNAR). For 'MNAR' alternative versions are available 
#'  depending on whether the mechanism for only one variable is condiered, that is for the effects (MNAR_eff)
#'  or the costs (MNAR_cost), or also covariates are included either for the effects (MNAR_eff_cov),
#'  the costs (MNAR_cost), or both (MNAR_cov) 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


normal_gamma_ind<-function(type)eval.parent(substitute({
  #write model file in JAGS/BUGS code for MCAR
  if(type=="MCAR"){
    model_string_jags<-    "
    model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(alpha_c[1],beta_c[1])
    eff1[i]~dnorm(mu_e[1],tau_e[1])
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-gamma0_e
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-gamma0_c
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(alpha_c[2],beta_c[2])
    eff2[i]~dnorm(mu_e[2],tau_e[2])
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-gamma0_e
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-gamma0_c
    }
    
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])
    #from shape and rate to mean and standard deviation
    alpha_c[t]<-pow(mu_c[t],2)/pow(s_c[t],2)
    beta_c[t]<-mu_c[t]/pow(s_c[t],2)

    #missingness probability
    p_e[t]<-ilogit(gamma0_e)
    p_c[t]<-ilogit(gamma0_c)
    
    #priors
    #model for costs and effects
    s_c[t]~dt(0,0.16,1)T(0,)         
    ls_e[t]~dunif(-5,10)
    mu_c[t]~dunif(0,10000)
    mu_e[t]~dnorm(0,0.0001)
    }
    
    #missing data mechanism
    gamma0_e~dlogis(0,1)
    gamma0_c~dlogis(0,1)
    }
    "
    #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
    model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="gamma")
  } else if(type=="MAR"){
    #write model file in JAGS/BUGS code for MAR
    model_string_jags<-"
    model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(alpha_c1[i],beta_c1[i])
    eff1[i]~dnorm(mu_e1[i],tau_e[1])
    #mean regression
    alpha_c1[i]<-pow(mu_c1[i],2)/pow(s_c[1],2)
    beta_c1[i]<-mu_c1[i]/pow(s_c[1],2)
    mu_c1[i]<-beta0_c[1]+inprod(X1[i,],beta_c[,1])
    mu_e1[i]<-beta0_e[1]+inprod(X1[i,],beta_e[,1])
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-gamma0_e+inprod(X1[i,],gamma_e[])
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-gamma0_c+inprod(X1[i,],gamma_c[])
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(alpha_c2[i],beta_c2[i])
    eff2[i]~dnorm(mu_e2[i],tau_e[2])
    #mean regression
    alpha_c2[i]<-pow(mu_c2[i],2)/pow(s_c[2],2)
    beta_c2[i]<-mu_c2[i]/pow(s_c[2],2)
    mu_c2[i]<-beta0_c[2]+inprod(X2[i,],beta_c[,2])
    mu_e2[i]<-beta0_e[2]+inprod(X2[i,],beta_e[,2])
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-gamma0_e+inprod(X2[i,],gamma_e[])
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-gamma0_c+inprod(X2[i,],gamma_c[])
    }
    
    #obtain mean values for eff and cost
    xp_c1<-inprod(mean_cov1[],beta_c[,1])
    xp_c2<-inprod(mean_cov2[],beta_c[,2])
    xp_e1<-inprod(mean_cov1[],beta_e[,1])
    xp_e2<-inprod(mean_cov2[],beta_e[,2])
    mu_c[1]<-beta0_c[1]+xp_c1
    mu_c[2]<-beta0_c[2]+xp_c2
    mu_e[1]<-beta0_e[1]+xp_e1
    mu_e[2]<-beta0_e[2]+xp_e2
    
    #missingness probability
    xp_mis_c1<-inprod(mean_cov1[],gamma_c[])
    xp_mis_c2<-inprod(mean_cov2[],gamma_c[])
    xp_mis_e1<-inprod(mean_cov1[],gamma_e[])
    xp_mis_e2<-inprod(mean_cov2[],gamma_e[])
    p_e[1]<-ilogit(gamma0_e+xp_mis_e1)
    p_e[2]<-ilogit(gamma0_e+xp_mis_e2)
    p_c[1]<-ilogit(gamma0_c+xp_mis_c1)
    p_c[2]<-ilogit(gamma0_c+xp_mis_c2)
    
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])

    #priors
    #model for costs and effects
    s_c[t]~dt(0,0.16,1)T(0,)         
    ls_e[t]~dunif(-5,10)
    }
    
    #priors for mean regression coefficients
    for(t in 1:2){
    beta0_c[t]~dunif(0,10000)
    beta0_e[t]~dnorm(0,0.0001)
    }
    for (j in 1:p) {
    for(t in 1:2){
    beta_c[j,t]~dnorm(0,0.01)
    beta_e[j,t]~dnorm(0,0.01)
    }
    }
    
    #missing data mechanism
    #intercepts
    gamma0_c~dlogis(0,1)
    gamma0_e~dlogis(0,1)
    #logistic regression coefficients
    for(j in 1:p){
    gamma_c[j]~dnorm(0,1)
    gamma_e[j]~dnorm(0,1)
    }
    }
    "
    #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
    model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="gamma")
  } else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
    #write model file in JAGS/BUGS code for MNAR
    model_string_jags<-  "
    model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(alpha_c[1],beta_c[1])
    eff1[i]~dnorm(mu_e[1],tau_e[1])
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-gamma0_e+delta_e*(eff1[i])
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-gamma0_c+delta_c*(cost1[i])
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(alpha_c[2],beta_c[2])
    eff2[i]~dnorm(mu_e[2],tau_e[2])
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-gamma0_e+delta_e*(eff2[i])
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-gamma0_c+delta_c*(cost2[i])
    }
    
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])
    #from shape and rate to mean and standard deviation
    alpha_c[t]<-pow(mu_c[t],2)/pow(s_c[t],2)
    beta_c[t]<-mu_c[t]/pow(s_c[t],2)
    
    #priors
    #model for costs and effects
    s_c[t]~dt(0,0.16,1)T(0,)         
    ls_e[t]~dunif(-5,10)
    mu_c[t]~dunif(0,10000)
    mu_e[t]~dnorm(0,0.0001)
    }

    #missingness probability
    p_e[1]<-ilogit(gamma0_e+delta_e*mean(eff1[]))
    p_e[2]<-ilogit(gamma0_e+delta_e*mean(eff2[]))
    p_c[1]<-ilogit(gamma0_c+delta_c*mean(cost1[]))
    p_c[2]<-ilogit(gamma0_c+delta_c*mean(cost2[]))
    
    #missing data mechanism
    gamma0_e~dlogis(0,1)
    gamma0_c~dlogis(0,1)
    #MNAR parameters
    delta_e~dnorm(0,1)
    delta_c~dnorm(0,1)
    }
    "
    if(type=="MNAR_eff"){
      model_string_jags<-gsub("delta_c~dnorm(0,1)","delta_c<-0",model_string_jags,fixed=TRUE)
    } else if(type=="MNAR_cost"){
      model_string_jags<-gsub("delta_e~dnorm(0,1)","delta_e<-0",model_string_jags,fixed=TRUE)
    }
    #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
    model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="gamma")
  }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
    #write model file in JAGS/BUGS code for MCAR
    model_string_jags<-  "
    model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(alpha_c1[i],beta_c1[i])
    eff1[i]~dnorm(mu_e1[i],tau_e[1])
    #mean regression
    alpha_c1[i]<-pow(mu_c1[i],2)/pow(s_c[1],2)
    beta_c1[i]<-mu_c1[i]/pow(s_c[1],2)
    mu_c1[i]<-beta0_c[1]+inprod(X1[i,],beta_c[,1])
    mu_e1[i]<-beta0_e[1]+inprod(X1[i,],beta_e[,1])
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-gamma0_e+inprod(X1[i,],gamma_e[])+delta_e*(eff1[i])
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-gamma0_c+inprod(X1[i,],gamma_c[])+delta_c*(cost1[i])
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(alpha_c2[i],beta_c2[i])
    eff2[i]~dnorm(mu_e2[i],tau_e[2])
    #mean regression
    alpha_c2[i]<-pow(mu_c2[i],2)/pow(s_c[2],2)
    beta_c2[i]<-mu_c2[i]/pow(s_c[2],2)
    mu_c2[i]<-beta0_c[2]+inprod(X2[i,],beta_c[,2])
    mu_e2[i]<-beta0_e[2]+inprod(X2[i,],beta_e[,2])
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-gamma0_e+inprod(X2[i,],gamma_e[])+delta_e*(eff2[i])
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-gamma0_c+inprod(X2[i,],gamma_c[])+delta_c*(cost2[i])
    }
    
    #obtain mean values for eff and cost
    xp_c1<-inprod(mean_cov1[],beta_c[,1])
    xp_c2<-inprod(mean_cov2[],beta_c[,2])
    xp_e1<-inprod(mean_cov1[],beta_e[,1])
    xp_e2<-inprod(mean_cov2[],beta_e[,2])
    mu_c[1]<-beta0_c[1]+xp_c1
    mu_c[2]<-beta0_c[2]+xp_c2
    mu_e[1]<-beta0_e[1]+xp_e1
    mu_e[2]<-beta0_e[2]+xp_e2
    
    
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])
    
    #priors
    #model for costs and effects
    s_c[t]~dt(0,0.16,1)T(0,)         
    ls_e[t]~dunif(-5,10)
    }

    #missingness probability
    xp_mis_c1<-inprod(mean_cov1[],gamma_c[])
    xp_mis_c2<-inprod(mean_cov2[],gamma_c[])
    xp_mis_e1<-inprod(mean_cov1[],gamma_e[])
    xp_mis_e2<-inprod(mean_cov2[],gamma_e[])
    p_e[1]<-ilogit(gamma0_e+xp_mis_e1+delta_e*mean(eff1[]))
    p_e[2]<-ilogit(gamma0_e+xp_mis_e2+delta_e*mean(eff2[]))
    p_c[1]<-ilogit(gamma0_c+xp_mis_c1+delta_c*mean(cost1[]))
    p_c[2]<-ilogit(gamma0_c+xp_mis_c2+delta_c*mean(cost2[]))
    
    #priors for mean regression coefficients
    for(t in 1:2){
    beta0_c[t]~dunif(0,10000)
    beta0_e[t]~dnorm(0,0.0001)
    }
    for (j in 1:p) {
    for(t in 1:2){
    beta_c[j,t]~dnorm(0,0.01)
    beta_e[j,t]~dnorm(0,0.01)
    }
    }
    
    #missing data mechanism
    #intercepts
    gamma0_c~dlogis(0,1)
    gamma0_e~dlogis(0,1)
    #MNAR parameters
    delta_e~dnorm(0,1)
    delta_c~dnorm(0,1)
    
    #logistic regression coefficients
    for(j in 1:p){
    gamma_c[j]~dnorm(0,1)
    gamma_e[j]~dnorm(0,1)
    }
    }
    "
    if(type=="MNAR_eff_cov"){
      model_string_jags<-gsub("delta_c~dnorm(0,1)","delta_c<-0",model_string_jags,fixed=TRUE)
    } else if(type=="MNAR_cost_cov"){
      model_string_jags<-gsub("delta_e~dnorm(0,1)","delta_e<-0",model_string_jags,fixed=TRUE)
    }
  }
  #include transformations if indiciated
  if(any(transf=="log")==TRUE){
    if(type=="MCAR"|type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      nu_c<-"nu_c[t]~dnorm(0,0.001)"
      log_mu_c<-"log(mu_c[t])<-nu_c[t]"
      transf_paste<-paste(nu_c,log_mu_c,sep = "\n")
      model_string_jags<-gsub("mu_c[t]~dunif(0,10000)", transf_paste, model_string_jags,fixed=TRUE)
    }
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      beta0_c_transf<-"beta0_c[t]~dnorm(0,0.001)"
      mu_c1_transf<-"mu_c[1]<-exp(beta0_c[1]+xp_c1)"
      mu_c2_transf<-"mu_c[2]<-exp(beta0_c[2]+xp_c2)"
      log_mu_c1<-"log(mu_c1[i])<-beta0_c[1]+inprod(X1[i,],beta_c[,1])"
      log_mu_c2<-"log(mu_c2[i])<-beta0_c[2]+inprod(X2[i,],beta_c[,2])"
      model_string_jags<-gsub("mu_c1[i]<-beta0_c[1]+inprod(X1[i,],beta_c[,1])", log_mu_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("mu_c2[i]<-beta0_c[2]+inprod(X2[i,],beta_c[,2])", log_mu_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta0_c[t]~dunif(0,10000)", beta0_c_transf, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("mu_c[1]<-beta0_c[1]+xp_c1", mu_c1_transf, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("mu_c[2]<-beta0_c[2]+xp_c2", mu_c2_transf, model_string_jags,fixed=TRUE)
    }
  }
  #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="gamma")
      #assign name to model file to be called by run_jags or run_bugs
      writeLines(model_string_jags,"normal_gamma_ind.txt")
      model_string<-"normal_gamma_ind.txt"
  return(model_string)
}))