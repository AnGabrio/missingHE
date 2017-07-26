#' An internal function to write BUGS code for bivariate normal model with natural-scaled variables (joint)
#'
#' This function writes in the current WD a txt file 
#' for the bivariate joint normal model with natural-scaled outcome variables.
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

normal_model_joint<-function(type=type)eval.parent(substitute({
    if(type=="MCAR"){
    #write model file in JAGS/BUGS code for MCAR
    model_string_jags<-    "
                      model{
                      #control
                      for(i in 1:N1){
                      #costs and effects model
                      cost1[i]~dnorm(phi_c1[i],lambda[1])
                      phi_c1[i]<-mu_c[1]+theta[1]*(eff1[i])
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
                      cost2[i]~dnorm(mu_c[2],lambda[2])
                      phi_c2[i]<-mu_c[2]+theta[2]*(eff2[i])
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
                      ss_c[t]<-s_c[t]*s_c[t]
                      s_c[t]<-exp(ls_c[t])
                      tau_e[t]<-1/ss_e[t]
                      ss_e[t]<-s_e[t]*s_e[t]
                      s_e[t]<-exp(ls_e[t])
                      
                      #conditional variance relation
                      lambda[t]<-1/psi[t]                 #conditional precision
                      psi[t]<-ss_c[t]-ss_e[t]*pow(theta[t],2)  #conditional variance
                      
                      #missingness probability
                      p_e[t]<-ilogit(gamma0_e)
                      p_c[t]<-ilogit(gamma0_c)
                      
                      #priors
                      #model for costs and effects
                      theta[t]~dnorm(0,0.0001)
                      ls_c[t]~dunif(-5,10)         
                      ls_e[t]~dunif(-5,10)
                      mu_c[t]~dnorm(0,0.0001)
                      mu_e[t]~dnorm(0,0.0001)
                      }
                      
                      #missing data mechanism
                      gamma0_e~dlogis(0,1)
                      gamma0_c~dlogis(0,1)
                      }
                      "
    
                #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
    model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
    } else if(type=="MAR"){
      #write model file in JAGS/BUGS code for MAR
      model_string_jags<-    "
      model{
      #control
      for(i in 1:N1){
      #costs and effects model
      cost1[i]~dnorm(phi_c1[i],lambda[1])
      phi_c1[i]<-mu_c1[i]+theta[1]*(eff1[i])
      eff1[i]~dnorm(mu_e1[i],tau_e[1])
      #mean regression
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
      cost2[i]~dnorm(phi_c2[i],lambda[2])
      phi_c2[i]<-mu_c2[i]+theta[2]*(eff2[i])
      eff2[i]~dnorm(mu_e2[i],tau_e[2])
      #mean regression
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
      tau_c[t]<-1/ss_c[t]
      ss_c[t]<-s_c[t]*s_c[t]
      s_c[t]<-exp(ls_c[t])
      tau_e[t]<-1/ss_e[t]
      ss_e[t]<-s_e[t]*s_e[t]
      s_e[t]<-exp(ls_e[t])
      
      #conditional precision and variance
      lambda[t]<-1/psi[t]
      psi[t]<-ss_c[t]-ss_e[t]*pow(theta[t],2)
      
      #priors
      #model for costs and effects
      ls_c[t]~dunif(-5,10)         
      ls_e[t]~dunif(-5,10)
      theta[t]~dnorm(0,0.0001)
      }
      
      
      #priors for mean regression coefficients
      for(t in 1:2){
      beta0_c[t]~dnorm(0,0.0001)
      beta0_e[t]~dnorm(0,0.0001)
      }
      for (j in 1:p) {#begin beta priors
      for(t in 1:2){
      beta_c[j,t]~dnorm(0,0.01)
      beta_e[j,t]~dnorm(0,0.01)
        }
      }#end beta priors
      
      #missing data mechanism
      #intercepts
      gamma0_c~dlogis(0,1)
      gamma0_e~dlogis(0,1)
      #logistic regression coefficients
      for (j in 1:p){#begin gamma priors
      gamma_c[j]~dnorm(0,1)
      gamma_e[j]~dnorm(0,1)
      }#end gamma priors
      }
      "
      
      #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
      model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
    } else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      #write model file in JAGS/BUGS code for MNAR
      model_string_jags<-  "
      model{
      #control
      for(i in 1:N1){
      #costs and effects model
      cost1[i]~dnorm(phi_c1[i],lambda[1])
      phi_c1[i]<-mu_c[1]+theta[1]*(eff1[i])
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
      cost2[i]~dnorm(phi_c2[i],lambda[2])
      phi_c2[i]<-mu_c[2]+theta[2]*(eff2[i])
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
      tau_c[t]<-1/ss_c[t]
      ss_c[t]<-s_c[t]*s_c[t]
      s_c[t]<-exp(ls_c[t])
      tau_e[t]<-1/ss_e[t]
      ss_e[t]<-s_e[t]*s_e[t]
      s_e[t]<-exp(ls_e[t])
      
      #conditional variance relation
      lambda[t]<-1/psi[t]                 #conditional precision
      psi[t]<-ss_c[t]-ss_e[t]*pow(theta[t],2)  #conditional variance
      
      #priors
      #model for costs and effects
      theta[t]~dnorm(0,0.0001)
      ls_c[t]~dunif(-5,10)         
      ls_e[t]~dunif(-5,10)
      mu_c[t]~dnorm(0,0.0001)
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
      model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
    }else if(type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      #write model file in JAGS/BUGS code for MCAR
      model_string_jags<-  "
      model{
      #control
      for(i in 1:N1){
      #costs and effects model
      cost1[i]~dnorm(phi_c1[i],lambda[1])
      phi_c1[i]<-mu_c1[i]+theta[1]*(eff1[i])
      eff1[i]~dnorm(mu_e1[i],tau_e[1])
      #mean regression
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
      cost2[i]~dnorm(phi_c2[i],lambda[2])
      phi_c2[i]<-mu_c2[i]+theta[2]*(eff2[i])
      eff2[i]~dnorm(mu_e2[i],tau_e[2])
      #mean regression
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
      
      #missingness probability
      xp_mis_c1<-inprod(mean_cov1[],gamma_c[])
      xp_mis_c2<-inprod(mean_cov2[],gamma_c[])
      xp_mis_e1<-inprod(mean_cov1[],gamma_e[])
      xp_mis_e2<-inprod(mean_cov2[],gamma_e[])
      p_e[1]<-ilogit(gamma0_e+xp_mis_e1+delta_e*mean(eff1[]))
      p_e[2]<-ilogit(gamma0_e+xp_mis_e2+delta_e*mean(eff2[]))
      p_c[1]<-ilogit(gamma0_c+xp_mis_c1+delta_c*mean(cost1[]))
      p_c[2]<-ilogit(gamma0_c+xp_mis_c2+delta_c*mean(cost2[]))
      
      #transformation of parameters
      for (t in 1:2){
      #from precision to variance and standard deviation
      tau_c[t]<-1/ss_c[t]
      ss_c[t]<-s_c[t]*s_c[t]
      s_c[t]<-exp(ls_c[t])
      tau_e[t]<-1/ss_e[t]
      ss_e[t]<-s_e[t]*s_e[t]
      s_e[t]<-exp(ls_e[t])
      
      #conditional variance relation
      lambda[t]<-1/psi[t]                 #conditional precision
      psi[t]<-ss_c[t]-ss_e[t]*pow(theta[t],2)  #conditional variance

      #priors
      #model for costs and effects
      theta[t]~dnorm(0,0.0001)
      ls_c[t]~dunif(-5,10)         
      ls_e[t]~dunif(-5,10)
      }
      
      #priors for mean regression coefficients
      for(t in 1:2){
      beta0_c[t]~dnorm(0,0.0001)
      beta0_e[t]~dnorm(0,0.0001)
      }
      for (j in 1:p) {#begin beta priors
      for(t in 1:2){
      beta_c[j,t]~dnorm(0,0.01)
      beta_e[j,t]~dnorm(0,0.01)
        }
      }#end beta priors
      
      #missing data mechanism
      #intercepts
      gamma0_c~dlogis(0,1)
      gamma0_e~dlogis(0,1)
      #MNAR parameters
      delta_e~dnorm(0,1)
      delta_c~dnorm(0,1)
      
      #logistic regression coefficients
      for (j in 1:p){#begin gamma priors
      gamma_c[j]~dnorm(0,1)
      gamma_e[j]~dnorm(0,1)
      }#end gamma priors
      }
      "
      if(type=="MNAR_eff_cov"){
        model_string_jags<-gsub("delta_c~dnorm(0,1)","delta_c<-0",model_string_jags,fixed=TRUE)
      } else if(type=="MNAR_cost_cov"){
        model_string_jags<-gsub("delta_e~dnorm(0,1)","delta_e<-0",model_string_jags,fixed=TRUE)
      }
    }
  #if predictor is one BUGS need a different specification for inpprod function
  if(p==1){
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      X1<-as.vector(X1)
      X2<-as.vector(X2)
      inprod_e1<-"X1[i]*beta_e[1]"
      inprod_e2<-"X2[i]*beta_e[2]"
      inprod_c1<-"X1[i]*beta_c[1]"
      inprod_c2<-"X2[i]*beta_c[2]"
      inprod_pq1<-"X1[i]*gamma_e"
      inprod_pq2<-"X2[i]*gamma_e"
      inprod_pc1<-"X1[i]*gamma_c"
      inprod_pc2<-"X2[i]*gamma_c"
      inprod_mean_e1<-"mean_cov1*beta_e[1]"
      inprod_mean_e2<-"mean_cov2*beta_e[2]"
      inprod_mean_c1<-"mean_cov1*beta_c[1]"
      inprod_mean_c2<-"mean_cov2*beta_c[2]"
      inprod_mean_pq1<-"mean_cov1*gamma_e"
      inprod_mean_pq2<-"mean_cov2*gamma_e"
      inprod_mean_pc1<-"mean_cov1*gamma_c"
      inprod_mean_pc2<-"mean_cov2*gamma_c"
      begin_prior_beta<-"#begin beta priors"
      prior_beta_e<-"beta_e[t]~dnorm(0,0.01)"
      prior_beta_c<-"beta_c[t]~dnorm(0,0.01)"
      end_prior_beta<-"#end beta priors"
      begin_prior_gamma<-"#begin gamma priors"
      prior_gamma_e<-"gamma_e~dnorm(0,1)"
      prior_gamma_c<-"gamma_c~dnorm(0,1)"
      end_prior_gamma<-"#end gamma priors"
      model_string_jags<-gsub("inprod(X1[i,],beta_e[,1])", inprod_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2[i,],beta_e[,2])", inprod_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X1[i,],beta_c[,1])", inprod_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2[i,],beta_c[,2])", inprod_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X1[i,],gamma_e[])", inprod_pq1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2[i,],gamma_e[])", inprod_pq2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X1[i,],gamma_c[])", inprod_pc1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2[i,],gamma_c[])", inprod_pc2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov1[],beta_e[,1])", inprod_mean_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov2[],beta_e[,2])", inprod_mean_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov1[],beta_c[,1])", inprod_mean_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov2[],beta_c[,2])", inprod_mean_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov1[],gamma_e[])", inprod_mean_pq1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov2[],gamma_e[])", inprod_mean_pq2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov1[],gamma_c[])", inprod_mean_pc1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov2[],gamma_c[])", inprod_mean_pc2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 1:p) {#begin beta priors", begin_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_e[j,t]~dnorm(0,0.01)", prior_beta_e, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c[j,t]~dnorm(0,0.01)", prior_beta_c, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end beta priors", end_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 1:p){#begin gamma priors", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[j]~dnorm(0,1)", prior_gamma_e, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_c[j]~dnorm(0,1)", prior_gamma_c, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors", end_prior_gamma, model_string_jags,fixed=TRUE)
    }
  }
  #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
                #assign name to model file to be called by run_jags or run_bugs
                 writeLines(model_string_jags,"normal_model_joint.txt")
                 model_string<-"normal_model_joint.txt"
  return(model_string)
  }))
