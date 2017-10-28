#' An internal function to write BUGS code for bivariate normal with an hurdle model for the effects (independence)
#'
#' This function writes in the current WD a txt file 
#' for the bivariate independent normal with an hurdle model for the effects
#' @keywords JAGS Hurdle models
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#'  and Structural At Random (MNAR) 
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #

normal_hurdle_e_ind<-function(type)eval.parent(substitute({
    #write model file in JAGS/BUGS code for SCAR
  if(type=="SCAR"){
    model_string_jags<-    "
                      model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dnorm(mu_c1[i],tau_c1)
    eff1[i]~dnorm(mu_e1[i],tau_e1[d_eff1[i]+1])
    #mean regression
    mu_c1[i]<-inprod(X1_c[i,],beta_c1[])
    mu_e1[i]<-inprod(X1_e[i,],beta_e1[,d_eff1[i]+1])    
    #structural values mechanism
    d_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-gamma_e[1]
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dnorm(mu_c2[i],tau_c2)
    eff2[i]~dnorm(mu_e2[i],tau_e2[d_eff2[i]+1])
    #mean regression
    mu_c2[i]<-inprod(X2_c[i,],beta_c2[])
    mu_e2[i]<-inprod(X2_e[i,],beta_e2[,d_eff2[i]+1])    
    #structural values mechanism
    d_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-gamma_e[2]
    }
    
    #transformation of parameters
    tau_c1<-1/ss_c1
    ss_c1<-s_c1*s_c1
    s_c1<-exp(ls_c1)
    tau_c2<-1/ss_c2
    ss_c2<-s_c2*s_c2
    s_c2<-exp(ls_c2)

    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e1[t]<-1/ss_e1[t]
    ss_e1[t]<-s_e1[t]*s_e1[t]
    s_e1[t]<-exp(ls_e1[t])
    
    tau_e2[t]<-1/ss_e2[t]
    ss_e2[t]<-s_e2[t]*s_e2[t]
    s_e2[t]<-exp(ls_e2[t])
    
    #structural values probability
    p_e[t]<-ilogit(gamma_e[t])
    }
    
    #calculate means at mean of covariates for non-structural values
    nu_c[1]<-inprod(mean_cov_c1[],beta_c1[])
    nu_c[2]<-inprod(mean_cov_c2[],beta_c2[])
    nu_e[1]<-inprod(mean_cov_e1[],beta_e1[,1])
    nu_e[2]<-inprod(mean_cov_e2[],beta_e2[,1])

    #obtain coefficients only for non-structural values
    for(j in 1:pe){#beta non-structural effects
    beta_e[j,1]<-beta_e1[j,1]
    beta_e[j,2]<-beta_e2[j,1]}
    for(j in 1:pc){#beta non-structural costs
    beta_c[j,1]<-beta_c1[j]
    beta_c[j,2]<-beta_c2[j]}

    #weighted mean parameters
    mu_e[1]<-nu_e[1]*(1-p_e[1])+se*p_e[1]
    mu_e[2]<-nu_e[2]*(1-p_e[2])+se*p_e[2]
    
    mu_c[1]<-nu_c[1]
    mu_c[2]<-nu_c[2]
    
    #sd parameters for non-structural values
    s_e[1]<-s_e1[1]
    s_e[2]<-s_e2[1]
    
    s_c[1]<-s_c1
    s_c[2]<-s_c2
    
    #priors

    #priors for mean regression coefficients
    for (j in 2:pe) {#begin beta priors effects
    beta_e1[j,1]~dnorm(0,0.01)
    beta_e2[j,1]~dnorm(0,0.01)
    beta_e1[j,2]<-0
    beta_e2[j,2]<-0
    }#end beta priors effects
    beta_e1[1,1]~dnorm(0,0.001)
    beta_e2[1,1]~dnorm(0,0.001)
    beta_e1[1,2]<-se
    beta_e2[1,2]<-se
                
    for (j in 2:pc) {#begin beta priors costs
    beta_c1[j]~dnorm(0,0.01)
    beta_c2[j]~dnorm(0,0.01)
    }#end beta priors costs
    beta_c1[1]~dnorm(0,0.001)
    beta_c2[1]~dnorm(0,0.001)

    #standard deviation priors
    ls_c1~dunif(-5,10)         
    ls_e1[1]~dunif(-5,10)
    ls_c2~dunif(-5,10)         
    ls_e2[1]~dunif(-5,10)
    ls_e1[2]<-sde
    ls_e2[2]<-sde
    
    #priors on structural values mechanism
    for(t in 1:2){
    gamma_e[t]~dlogis(0,1)
    }
  }
      "
      #call prior_change function to make changes to prior values and distributions if prior inputs in hurdle are provided
    #model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
  } else if(type=="SAR"){
    #write model file in JAGS/BUGS code for SAR
    model_string_jags<-"
                      model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dnorm(mu_c1[i],tau_c1)
    eff1[i]~dnorm(mu_e1[i],tau_e1[d_eff1[i]+1])
    #mean regression
    mu_c1[i]<-inprod(X1_c[i,],beta_c1[])
    mu_e1[i]<-inprod(X1_e[i,],beta_e1[,d_eff1[i]+1])    
    #structural values mechanism
    d_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-inprod(Z1_e[i,],gamma_e[,1])
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dnorm(mu_c2[i],tau_c2)
    eff2[i]~dnorm(mu_e2[i],tau_e2[d_eff2[i]+1])
    #mean regression
    mu_c2[i]<-inprod(X2_c[i,],beta_c2[])
    mu_e2[i]<-inprod(X2_e[i,],beta_e2[,d_eff2[i]+1])    
    #structural values mechanism
    d_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-inprod(Z2_e[i,],gamma_e[,2])
    }
    
    #transformation of parameters
    tau_c1<-1/ss_c1
    ss_c1<-s_c1*s_c1
    s_c1<-exp(ls_c1)
    tau_c2<-1/ss_c2
    ss_c2<-s_c2*s_c2
    s_c2<-exp(ls_c2)

    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e1[t]<-1/ss_e1[t]
    ss_e1[t]<-s_e1[t]*s_e1[t]
    s_e1[t]<-exp(ls_e1[t])
    
    tau_e2[t]<-1/ss_e2[t]
    ss_e2[t]<-s_e2[t]*s_e2[t]
    s_e2[t]<-exp(ls_e2[t])
    
    }
    
    #calculate means at mean of covariates for non-structural values
    nu_c[1]<-inprod(mean_cov_c1[],beta_c1[])
    nu_c[2]<-inprod(mean_cov_c2[],beta_c2[])
    nu_e[1]<-inprod(mean_cov_e1[],beta_e1[,1])
    nu_e[2]<-inprod(mean_cov_e2[],beta_e2[,1])

    #obtain coefficients only for non-structural values
    for(j in 1:pe){#beta non-structural effects
    beta_e[j,1]<-beta_e1[j,1]
    beta_e[j,2]<-beta_e2[j,1]}
    for(j in 1:pc){#beta non-structural costs
    beta_c[j,1]<-beta_c1[j]
    beta_c[j,2]<-beta_c2[j]}
    
    #weighted mean parameters
    mu_e[1]<-nu_e[1]*(1-p_e[1])+se*p_e[1]
    mu_e[2]<-nu_e[2]*(1-p_e[2])+se*p_e[2]
    
    mu_c[1]<-nu_c[1]
    mu_c[2]<-nu_c[2]
    
    #sd parameters for non-structural values
    s_e[1]<-s_e1[1]
    s_e[2]<-s_e2[1]
    
    s_c[1]<-s_c1
    s_c[2]<-s_c2

    #structural values probability
    p_e[1]<-ilogit(inprod(mean_z_e1[],gamma_e[,1]))
    p_e[2]<-ilogit(inprod(mean_z_e2[],gamma_e[,2]))
    
    #priors
    
    #priors for mean regression coefficients
    for (j in 2:pe) {#begin beta priors effects
    beta_e1[j,1]~dnorm(0,0.01)
    beta_e2[j,1]~dnorm(0,0.01)
    beta_e1[j,2]<-0
    beta_e2[j,2]<-0
    }#end beta priors effects
    beta_e1[1,1]~dnorm(0,0.001)
    beta_e2[1,1]~dnorm(0,0.001)
    beta_e1[1,2]<-se
    beta_e2[1,2]<-se
    
    for (j in 2:pc) {#begin beta priors costs
    beta_c1[j]~dnorm(0,0.01)
    beta_c2[j]~dnorm(0,0.01)
    }#end beta priors costs
    beta_c1[1]~dnorm(0,0.001)
    beta_c2[1]~dnorm(0,0.001)
    
    #standard deviation priors
    ls_c1~dunif(-5,10)         
    ls_e1[1]~dunif(-5,10)
    ls_c2~dunif(-5,10)         
    ls_e2[1]~dunif(-5,10)
    ls_e1[2]<-sde
    ls_e2[2]<-sde
    
    #priors on structural values mechanism
      for (j in 2:ze) {#begin gamma priors effects
      for(t in 1:2){gamma_e[j,t]~dnorm(0,0.01)}
    }#end gamma priors effects
    gamma_e[1,1]~dlogis(0,1)
    gamma_e[1,2]~dlogis(0,1)

    }
      "
  }
  #if predictor is one BUGS needs a different specification for inpprod function
   if(pe==1){
     if(type=="SCAR"|type=="SAR"){
       inprod_e1<-"X1_e[i]*beta_e1[d_eff1[i]+1]"
       inprod_e2<-"X2_e[i]*beta_e2[d_eff2[i]+1]"
       inprod_mean_e1<-"mean_cov_e1*beta_e1[1]"
       inprod_mean_e2<-"mean_cov_e2*beta_e2[1]"
       begin_beta_nons_e<-prior_beta_e1j<-prior_beta_e2j<-prior_beta_e10<-prior_beta_e20<-"#"
       beta_nons_e1<-"beta_e[1]<-beta_e1[1]"
       beta_nons_e2<-"beta_e[2]<-beta_e2[1]"
       begin_prior_beta<-"#begin beta priors effects"
       prior_beta_e1<-"beta_e1[1]~dnorm(0,0.0001)"
       prior_beta_e2<-"beta_e2[1]~dnorm(0,0.0001)"
       prior_beta_e1s<-"beta_e1[2]<-se"
       prior_beta_e2s<-"beta_e2[2]<-se"
       end_prior_beta<-"#end beta priors effects"
       model_string_jags<-gsub("inprod(X1_e[i,],beta_e1[,d_eff1[i]+1])", inprod_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(X2_e[i,],beta_e2[,d_eff2[i]+1])", inprod_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(mean_cov_e1[],beta_e1[,1])", inprod_mean_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(mean_cov_e2[],beta_e2[,1])", inprod_mean_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("for(j in 1:pe){#beta non-structural effects", begin_beta_nons_e, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e[j,1]<-beta_e1[j,1]", beta_nons_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e[j,2]<-beta_e2[j,1]}", beta_nons_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("for (j in 2:pe) {#begin beta priors effects", begin_prior_beta, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e1[j,1]~dnorm(0,0.01)", prior_beta_e1j, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e2[j,1]~dnorm(0,0.01)", prior_beta_e2j, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e1[j,2]<-0", prior_beta_e10, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e2[j,2]<-0", prior_beta_e20, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e1[1,1]~dnorm(0,0.001)", prior_beta_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e2[1,1]~dnorm(0,0.001)", prior_beta_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e1[1,2]<-se", prior_beta_e1s, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e2[1,2]<-se", prior_beta_e2s, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("}#end beta priors effects", end_prior_beta, model_string_jags,fixed=TRUE)
     }
   }
  if(pc==1){
    if(type=="SCAR"|type=="SAR"){
      inprod_c1<-"X1_c[i]*beta_c1[]"
      inprod_c2<-"X2_c[i]*beta_c2[]"
      inprod_mean_c1<-"mean_cov_c1*beta_c1[1]"
      inprod_mean_c2<-"mean_cov_c2*beta_c2[1]"
      begin_beta_nons_c<-prior_beta_c1j<-prior_beta_c2j<-prior_beta_c10<-prior_beta_c20<-"#"
      beta_nons_c1<-"beta_c[1]<-beta_c1[1]"
      beta_nons_c2<-"beta_c[2]<-beta_c2[1]"
      begin_prior_beta<-"#begin beta priors costs"
      prior_beta_c1<-"beta_c1[1]~dnorm(0,0.0001)"
      prior_beta_c2<-"beta_c2[1]~dnorm(0,0.0001)"
      end_prior_beta<-"#end beta priors costs"
      model_string_jags<-gsub("inprod(X1_c[i,],beta_c1[])", inprod_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2_c[i,],beta_c2[])", inprod_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c1[],beta_c1[])", inprod_mean_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c2[],beta_c2[])", inprod_mean_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(j in 1:pc){#beta non-structural costs", begin_beta_nons_c, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c[j,1]<-beta_c1[j]", beta_nons_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c[j,2]<-beta_c2[j]}", beta_nons_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c1[j]~dnorm(0,0.01)", prior_beta_c1j, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c2[j]~dnorm(0,0.01)", prior_beta_c2j, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c1[j]<-0", prior_beta_c10, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c2[j]<-0", prior_beta_c20, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c1[1]~dnorm(0,0.001)", prior_beta_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c2[1]~dnorm(0,0.001)", prior_beta_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end beta priors costs", end_prior_beta, model_string_jags,fixed=TRUE)
    }
  }
  if(ze==1){
    if(type=="SAR"){
      inprod_e1<-"Z1_e[i]*gamma_e[1]"
      inprod_e2<-"Z2_e[i]*gamma_e[2]"
      inprod_mean_e1<-"ilogit(mean_z_e1*gamma_e[1])"
      inprod_mean_e2<-"ilogit(mean_z_e2*gamma_e[2])"
      begin_prior_gamma<-"#begin gamma priors effects"
      begin_prior_gamma2<-"#"
      prior_gamma_e1<-"gamma_e[1]~dlogis(0,1)"
      prior_gamma_e2<-"gamma_e[2]~dlogis(0,1)"
      end_prior_gamma<-"#end gamma priors effects"
      model_string_jags<-gsub("inprod(Z1_e[i,],gamma_e[,1])", inprod_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(Z2_e[i,],gamma_e[,2])", inprod_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("ilogit(inprod(mean_z_e1[],gamma_e[,1]))", inprod_mean_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("ilogit(inprod(mean_z_e2[],gamma_e[,2]))", inprod_mean_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(t in 1:2){gamma_e[j,t]~dnorm(0,0.01)}", begin_prior_gamma2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[1,1]~dlogis(0,1)", prior_gamma_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[1,2]~dlogis(0,1)", prior_gamma_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags,fixed=TRUE)
    }
  }
  #call prior_change function to make changes to prior values and distributions if prior inputs in hurdle are provided
  model_string_jags<-prior_hurdle(type=type,dist_e ="norm",dist_c="norm")
      #assign name to model file to be called by run_jags or run_bugs
      writeLines(model_string_jags,"normal_hurdle_e_ind.txt")
      model_string<-"normal_hurdle_e_ind.txt"
  return(model_string)
}))
