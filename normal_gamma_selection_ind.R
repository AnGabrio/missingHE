#' An internal function to write BUGS code for bivariate normal-gamma selection model (independence)
#'
#' This function writes in the current WD a txt file 
#' for the bivariate independent normal-gamma selection model for both effects and costs.
#' @keywords JAGS Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #

normal_gamma_selection_ind<-function(type)eval.parent(substitute({
    #write model file in JAGS/BUGS code for MCAR
  if(type=="MAR"){
    model_string_jags<-    "
                      model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(mu_c1[i]*tau_c1[i],tau_c1[i])
    eff1[i]~dnorm(mu_e1[i],tau_e[1])
    #obtain mean and sd
    tau_c1[i]<-mu_c1[i]/pow(s_c[1],2)
    #mean regression
    log(mu_c1[i])<-inprod(X1_c[i,],beta_c[,1])
    mu_e1[i]<-inprod(X1_e[i,],beta_e[,1])    
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-inprod(Z1_e[i,],gamma_e[,1])
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-inprod(Z1_c[i,],gamma_c[,1])
    }
                      
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(mu_c2[i]*tau_c2[i],tau_c2[i])
    eff2[i]~dnorm(mu_e2[i],tau_e[2])
    #obtain mean and sd
    tau_c2[i]<-mu_c2[i]/pow(s_c[2],2)
    #mean regression
    log(mu_c2[i])<-inprod(X2_c[i,],beta_c[,2])
    mu_e2[i]<-inprod(X2_e[i,],beta_e[,2])    
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-inprod(Z2_e[i,],gamma_e[,2])
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-inprod(Z2_c[i,],gamma_c[,2])
    }
                      
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])
    }

    #missingness probability
    p_c[1]<-ilogit(inprod(mean_z_c1[],gamma_c[,1]))
    p_c[2]<-ilogit(inprod(mean_z_c2[],gamma_c[,2]))
    p_e[1]<-ilogit(inprod(mean_z_e1[],gamma_e[,1]))
    p_e[2]<-ilogit(inprod(mean_z_e2[],gamma_e[,2]))
    
    #calculate means at mean of covariates for non-structural values
    mu_c[1]<-exp(inprod(mean_cov_c1[],beta_c[,1]))
    mu_c[2]<-exp(inprod(mean_cov_c2[],beta_c[,2]))
    mu_e[1]<-inprod(mean_cov_e1[],beta_e[,1])
    mu_e[2]<-inprod(mean_cov_e2[],beta_e[,2])
    
    #priors
    
    #priors for mean regression coefficients
    for (j in 2:pe) {#begin beta priors effects
    for(t in 1:2){beta_e[j,t]~dnorm(0,0.01)}
    }#end beta priors effects
    beta_e[1,1]~dnorm(0,0.001)
    beta_e[1,2]~dnorm(0,0.001)
    
    for (j in 2:pc) {#begin beta priors costs
    for(t in 1:2){beta_c[j,t]~dnorm(0,0.01)}
    }#end beta priors costs
    beta_c[1,1]~dnorm(0,0.001)
    beta_c[1,2]~dnorm(0,0.001)

    #standard deviation priors
    for(t in 1:2){
    s_c[t]~dunif(0,1000)
    ls_e[t]~dunif(-5,10)
    }
    
    #priors on missing data mechanism
    for (j in 2:ze) {#begin gamma priors effects
    for(t in 1:2){gamma_e[j,t]~dnorm(0,0.01)}
    }#end gamma priors effects
    gamma_e[1,1]~dlogis(0,1)
    gamma_e[1,2]~dlogis(0,1)
    
    for (j in 2:zc) {#begin gamma priors costs
    for(t in 1:2){gamma_c[j,t]~dnorm(0,0.01)}
    }#end gamma priors costs
    gamma_c[1,1]~dlogis(0,1)
    gamma_c[1,2]~dlogis(0,1)

    }
      "
      #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  } else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
    #write model file in JAGS/BUGS code for MNAR
    model_string_jags<-  "
                      model{
    #control
    for(i in 1:N1){
    #costs and effects model
    cost1[i]~dgamma(mu_c1[i]*tau_c1[i],tau_c1[i])
    eff1[i]~dnorm(mu_e1[i],tau_e[1])
    #obtain mean and sd
    tau_c1[i]<-mu_c1[i]/pow(s_c[1],2)
    #mean regression
    log(mu_c1[i])<-inprod(X1_c[i,],beta_c[,1])
    mu_e1[i]<-inprod(X1_e[i,],beta_e[,1])    
    #missing data mechanism
    m_eff1[i]~dbern(pq_1[i])
    logit(pq_1[i])<-inprod(Z1_e[i,],gamma_e[,1])+delta_e[1]*eff1[i]
    m_cost1[i]~dbern(pc_1[i])
    logit(pc_1[i])<-inprod(Z1_c[i,],gamma_c[,1])+delta_c[1]*cost1[i]
    }
    
    #intervention
    for(i in 1:N2){
    #costs and effects model
    cost2[i]~dgamma(mu_c2[i]*tau_c2[i],tau_c2[i])
    eff2[i]~dnorm(mu_e2[i],tau_e[2])
    #obtain mean and sd
    tau_c2[i]<-mu_c2[i]/pow(s_c[2],2)
    #mean regression
    log(mu_c2[i])<-inprod(X2_c[i,],beta_c[,2])
    mu_e2[i]<-inprod(X2_e[i,],beta_e[,2])    
    #missing data mechanism
    m_eff2[i]~dbern(pq_2[i])
    logit(pq_2[i])<-inprod(Z2_e[i,],gamma_e[,2])+delta_e[2]*eff2[i]
    m_cost2[i]~dbern(pc_2[i])
    logit(pc_2[i])<-inprod(Z2_c[i,],gamma_c[,2])+delta_c[2]*cost2[i]
    }
    
    #transformation of parameters
    for (t in 1:2){
    #from precision to variance and standard deviation
    tau_e[t]<-1/ss_e[t]
    ss_e[t]<-s_e[t]*s_e[t]
    s_e[t]<-exp(ls_e[t])
    }

    #missingness probability
    p_c[1]<-ilogit(inprod(mean_z_c1[],gamma_c[,1])+delta_c[1]*mean(cost1[]))
    p_c[2]<-ilogit(inprod(mean_z_c2[],gamma_c[,2])+delta_c[2]*mean(cost2[]))
    p_e[1]<-ilogit(inprod(mean_z_e1[],gamma_e[,1])+delta_e[1]*mean(eff1[]))
    p_e[2]<-ilogit(inprod(mean_z_e2[],gamma_e[,2])+delta_e[2]*mean(eff2[]))
    
    #calculate means at mean of covariates for non-structural values
    mu_c[1]<-exp(inprod(mean_cov_c1[],beta_c[,1]))
    mu_c[2]<-exp(inprod(mean_cov_c2[],beta_c[,2]))
    mu_e[1]<-inprod(mean_cov_e1[],beta_e[,1])
    mu_e[2]<-inprod(mean_cov_e2[],beta_e[,2])
    
    #priors
    
    #priors for mean regression coefficients
    for (j in 2:pe) {#begin beta priors effects
    for(t in 1:2){beta_e[j,t]~dnorm(0,0.01)}
    }#end beta priors effects
    beta_e[1,1]~dnorm(0,0.001)
    beta_e[1,2]~dnorm(0,0.001)
    
    for (j in 2:pc) {#begin beta priors costs
    for(t in 1:2){beta_c[j,t]~dnorm(0,0.01)}
    }#end beta priors costs
    beta_c[1,1]~dnorm(0,0.001)
    beta_c[1,2]~dnorm(0,0.001)
    
    #standard deviation priors
    for(t in 1:2){
    s_c[t]~dunif(0,1000)
    ls_e[t]~dunif(-5,10)
    }
    
    #priors on missing data mechanism
    for (j in 2:ze) {#begin gamma priors effects
    for(t in 1:2){gamma_e[j,t]~dnorm(0,0.01)}
    }#end gamma priors effects
    gamma_e[1,1]~dlogis(0,1)
    gamma_e[1,2]~dlogis(0,1)
    
    for (j in 2:zc) {#begin gamma priors costs
    for(t in 1:2){gamma_c[j,t]~dnorm(0,0.01)}
    }#end gamma priors costs
    gamma_c[1,1]~dlogis(0,1)
    gamma_c[1,2]~dlogis(0,1)

    #mnar parameters
    for(t in 1:2){
    delta_e[t]~dnorm(0,1)
    delta_c[t]~dnorm(0,1)
    }
    
    }
       "
    if(type=="MNAR_eff"){
      model_string_jags<-gsub("+delta_c[1]*mean(cost1[])","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_c[2]*mean(cost2[])","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_c[1]*cost1[i]","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_c[2]*cost2[i]","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("delta_c[t]~dnorm(0,1)","",model_string_jags,fixed=TRUE)
    } else if(type=="MNAR_cost"){
      model_string_jags<-gsub("+delta_e[1]*mean(eff1[])","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_e[2]*mean(eff2[])","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_e[1]*eff1[i]","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("+delta_e[2]*eff2[i]","",model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("delta_e[t]~dnorm(0,1)","",model_string_jags,fixed=TRUE)
    }
    #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  }
  #if predictor is one BUGS needs a different specification for inpprod function
   if(pe==1){
       inprod_e1<-"X1_e[i]*beta_e[1]"
       inprod_e2<-"X2_e[i]*beta_e[2]"
       inprod_mean_e1<-"mean_cov_e1*beta_e[1]"
       inprod_mean_e2<-"mean_cov_e2*beta_e[2]"
       begin_prior_beta<-"#begin beta priors effects"
       prior_beta<-"#"
       end_prior_beta<-"#end beta priors effects"
       prior_beta_e1<-"beta_e[1]~dnorm(0,0.0001)"
       prior_beta_e2<-"beta_e[2]~dnorm(0,0.0001)"
       model_string_jags<-gsub("inprod(X1_e[i,],beta_e[,1])", inprod_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(X2_e[i,],beta_e[,2])", inprod_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(mean_cov_e1[],beta_e[,1])", inprod_mean_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("inprod(mean_cov_e2[],beta_e[,2])", inprod_mean_e2, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("for (j in 2:pe) {#begin beta priors effects", begin_prior_beta, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("for(t in 1:2){beta_e[j,t]~dnorm(0,0.01)}",prior_beta,model_string_jags,fixed = TRUE)
       model_string_jags<-gsub("}#end beta priors effects",end_prior_beta,model_string_jags,fixed = TRUE)
       model_string_jags<-gsub("beta_e[1,1]~dnorm(0,0.001)", prior_beta_e1, model_string_jags,fixed=TRUE)
       model_string_jags<-gsub("beta_e[1,2]~dnorm(0,0.001)", prior_beta_e2, model_string_jags,fixed=TRUE)
       
   }
  if(pc==1){
    inprod_c1<-"X1_c[i]*beta_c[1]"
    inprod_c2<-"X2_c[i]*beta_c[2]"
    inprod_mean_c1<-"mean_cov_c1*beta_c[1]"
    inprod_mean_c2<-"mean_cov_c2*beta_c[2]"
    begin_prior_beta<-"#begin beta priors costs"
    prior_beta<-"#"
    end_prior_beta<-"#end beta priors costs"
    prior_beta_c1<-"beta_c[1]~dnorm(0,0.0001)"
    prior_beta_c2<-"beta_c[2]~dnorm(0,0.0001)"
    model_string_jags<-gsub("inprod(X1_c[i,],beta_c[,1])", inprod_c1, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("inprod(X2_c[i,],beta_c[,2])", inprod_c2, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("inprod(mean_cov_c1[],beta_c[,1])", inprod_mean_c1, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("inprod(mean_cov_c2[],beta_c[,2])", inprod_mean_c2, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("for(t in 1:2){beta_c[j,t]~dnorm(0,0.01)}",prior_beta,model_string_jags,fixed = TRUE)
    model_string_jags<-gsub("}#end beta priors costs",end_prior_beta,model_string_jags,fixed = TRUE)
    model_string_jags<-gsub("beta_c[1,1]~dnorm(0,0.001)", prior_beta_c1, model_string_jags,fixed=TRUE)
    model_string_jags<-gsub("beta_c[1,2]~dnorm(0,0.001)", prior_beta_c2, model_string_jags,fixed=TRUE)
    
  }
  if(ze==1){
      inprod_e1<-"Z1_e[i]*gamma_e[1]"
      inprod_e2<-"Z2_e[i]*gamma_e[2]"
      if(type=="MAR"|type=="MNAR_cost"){
        inprod_mean_e1<-"ilogit(mean_z_e1*gamma_e[1])"
        inprod_mean_e2<-"ilogit(mean_z_e2*gamma_e[2])"
        model_string_jags<-gsub("ilogit(inprod(mean_z_e1[],gamma_e[,1]))", inprod_mean_e1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("ilogit(inprod(mean_z_e2[],gamma_e[,2]))", inprod_mean_e2, model_string_jags,fixed=TRUE)
      }
      if(type=="MNAR_eff"|type=="MNAR"){
        inprod_mean_e1<-"ilogit(mean_z_e1*gamma_e[1]+delta_e[1]*mean(eff1[]))"
        inprod_mean_e2<-"ilogit(mean_z_e2*gamma_e[2]+delta_e[2]*mean(eff2[]))"
        model_string_jags<-gsub("ilogit(inprod(mean_z_e1[],gamma_e[,1])+delta_e[1]*mean(eff1[]))", inprod_mean_e1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("ilogit(inprod(mean_z_e2[],gamma_e[,2])+delta_e[2]*mean(eff2[]))", inprod_mean_e2, model_string_jags,fixed=TRUE)
      }
      begin_prior_gamma<-"#begin gamma priors effects"
      begin_prior_gamma2<-"#"
      prior_gamma_e1<-"gamma_e[1]~dlogis(0,1)"
      prior_gamma_e2<-"gamma_e[2]~dlogis(0,1)"
      end_prior_gamma<-"#end gamma priors effects"
      model_string_jags<-gsub("inprod(Z1_e[i,],gamma_e[,1])", inprod_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(Z2_e[i,],gamma_e[,2])", inprod_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:ze) {#begin gamma priors effects", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(t in 1:2){gamma_e[j,t]~dnorm(0,0.01)}", begin_prior_gamma2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[1,1]~dlogis(0,1)", prior_gamma_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[1,2]~dlogis(0,1)", prior_gamma_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags,fixed=TRUE)
  }
  if(zc==1){
      inprod_c1<-"Z1_c[i]*gamma_c[1]"
      inprod_c2<-"Z2_c[i]*gamma_c[2]"
      if(type=="MAR"|type=="MNAR_eff"){
        inprod_mean_c1<-"ilogit(mean_z_c1*gamma_c[1])"
        inprod_mean_c2<-"ilogit(mean_z_c2*gamma_c[2])"
        model_string_jags<-gsub("ilogit(inprod(mean_z_c1[],gamma_c[,1]))", inprod_mean_c1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("ilogit(inprod(mean_z_c2[],gamma_c[,2]))", inprod_mean_c2, model_string_jags,fixed=TRUE)
      }
      if(type=="MNAR_cost"|type=="MNAR"){
        inprod_mean_c1<-"ilogit(mean_z_c1*gamma_c[1]+delta_c[1]*mean(cost1[]))"
        inprod_mean_c2<-"ilogit(mean_z_c2*gamma_c[2]+delta_c[2]*mean(cost2[]))"
        model_string_jags<-gsub("ilogit(inprod(mean_z_c1[],gamma_c[,1])+delta_c[1]*mean(cost1[]))", inprod_mean_c1, model_string_jags,fixed=TRUE)
        model_string_jags<-gsub("ilogit(inprod(mean_z_c2[],gamma_c[,2])+delta_c[2]*mean(cost2[]))", inprod_mean_c2, model_string_jags,fixed=TRUE)
      }
      begin_prior_gamma<-"#begin gamma priors costs"
      begin_prior_gamma2<-"#"
      prior_gamma_c1<-"gamma_c[1]~dlogis(0,1)"
      prior_gamma_c2<-"gamma_c[2]~dlogis(0,1)"
      end_prior_gamma<-"#end gamma priors costs"
      model_string_jags<-gsub("inprod(Z1_c[i,],gamma_c[,1])", inprod_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(Z2_c[i,],gamma_c[,2])", inprod_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:zc) {#begin gamma priors costs", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(t in 1:2){gamma_c[j,t]~dnorm(0,0.01)}", begin_prior_gamma2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_c[1,1]~dlogis(0,1)", prior_gamma_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_c[1,2]~dlogis(0,1)", prior_gamma_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags,fixed=TRUE)
  }
  
  #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  model_string_jags<-prior_selection(type=type,dist_e ="norm",dist_c="gamma")
      #assign name to model file to be called by run_jags or run_bugs
      writeLines(model_string_jags,"normal_gamma_selection_ind.txt")
      model_string<-"normal_gamma_selection_ind.txt"
  return(model_string)
}))
