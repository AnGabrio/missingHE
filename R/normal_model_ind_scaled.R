#' An internal function to write BUGS code for bivariate normal model with scaled variables (independence)
#'
#' This function writes in the current WD a txt file 
#' for the bivariate independent normal model with scaled outcome variables.
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

normal_model_ind_scaled<-function(type)eval.parent(substitute({
    if(type=="MCAR"){
    #write model file in JAGS/BUGS code for MCAR
    model_string_jags<-  "
                       model{
                       #variables denoted by _s are related to their standardised version 
                       #variables denoted by _t are are related to their standardised version that need to be back tranformed on their natural scale for conducting analysis
                       
                       #sd_ and mean_ are vectors containing the values, per treatment arm, of sd and mean necessary to back transform variables
                       
                       #control
                       for(i in 1:N1){
                       #costs and effects model
                       cost1_s[i]~dnorm(mu_c_t[1],tau_c_t[1])
                       eff1_s[i]~dnorm(mu_e_t[1],tau_e_t[1])
                       #missing data mechanism
                       m_eff1[i]~dbern(pq_1[i])
                       logit(pq_1[i])<-gamma0_e
                       m_cost1[i]~dbern(pc_1[i])
                       logit(pc_1[i])<-gamma0_c
                       }
                       
                       #intervention
                       for(i in 1:N2){
                       #costs and effects model
                       cost2_s[i]~dnorm(mu_c_t[2],tau_c_t[2])
                       eff2_s[i]~dnorm(mu_e_t[2],tau_e_t[2])
                       #missing data mechanism
                       m_eff2[i]~dbern(pq_2[i])
                       logit(pq_2[i])<-gamma0_e
                       m_cost2[i]~dbern(pc_2[i])
                       logit(pc_2[i])<-gamma0_c
                       }
                       
                       #transformation of parameters
                       for (t in 1:2){
                       #retrieve original variable scales
                       
                       mu_e[t]<-mu_e_t[t]*sd_eff[t]*2+mean_eff[t]
                       mu_c[t]<-mu_c_t[t]*sd_cost[t]*2+mean_cost[t]
                       
                       s_e[t]<-s_e_t[t]*sd_eff[t]*2+mean_eff[t]
                       s_c[t]<-s_c_t[t]*sd_cost[t]*2+mean_cost[t]
                       
                       #from precision to variance and standard deviation
                       tau_c_t[t]<-1/ss_c_t[t]
                       ss_c_t[t]<-s_c_t[t]*s_c_t[t]
                       s_c_t[t]<-exp(ls_c_t[t])
                       tau_e_t[t]<-1/ss_e_t[t]
                       ss_e_t[t]<-s_e_t[t]*s_e_t[t]
                       s_e_t[t]<-exp(ls_e_t[t])
                       
                       #missingness probability
                       p_e[t]<-ilogit(gamma0_e)
                       p_c[t]<-ilogit(gamma0_c)
                       
                       #priors
                       #model for costs and effects
                       ls_c_t[t]~dunif(-5,10)         
                       ls_e_t[t]~dunif(-5,10)
                       mu_c_t[t]~dnorm(0,0.0001)
                       mu_e_t[t]~dnorm(0,0.0001)
                       }
                       
                       #missing data mechanism
                       gamma0_e~dlogis(0,1)
                       gamma0_c~dlogis(0,1)

                       #retrieve original cost and effect data
                       for(i in 1:N1){
                       eff1[i]<-eff1_s[i]*sd_eff[1]*2+mean_eff[1]
                       cost1[i]<-cost1_s[i]*sd_cost[1]*2+mean_cost[1]
                       }
                       for(i in 1:N2){
                       eff2[i]<-eff2_s[i]*sd_eff[2]*2+mean_eff[2]
                       cost2[i]<-cost2_s[i]*sd_cost[2]*2+mean_cost[2]
                       }

                       }
                       "
  
                       #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
    model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
    }else if(type=="MAR"){
      #write model file in JAGS/BUGS code for MAR
      model_string_jags<-  "
      model{
      #control
      for(i in 1:N1){
      #costs and effects model
      cost1_s[i]~dnorm(mu_c1_t[i],tau_c_t[1])
      eff1_s[i]~dnorm(mu_e1_t[i],tau_e_t[1])
      #mean regression
      mu_c1_t[i]<-inprod(X1_cs[i,],beta_c[,1])
      mu_e1_t[i]<-inprod(X1_es[i,],beta_e[,1])
      #missing data mechanism
      m_eff1[i]~dbern(pq_1[i])
      logit(pq_1[i])<-inprod(X1_es[i,],gamma_e[])
      m_cost1[i]~dbern(pc_1[i])
      logit(pc_1[i])<-inprod(X1_cs[i,],gamma_c[])
      }
      
      #intervention
      for(i in 1:N2){
      #costs and effects model
      cost2_s[i]~dnorm(mu_c2_t[i],tau_c_t[2])
      eff2_s[i]~dnorm(mu_e2_t[i],tau_e_t[2])
      #mean regression
      mu_c2_t[i]<-inprod(X2_cs[i,],beta_c[,2])
      mu_e2_t[i]<-inprod(X2_es[i,],beta_e[,2])
      #missing data mechanism
      m_eff2[i]~dbern(pq_2[i])
      logit(pq_2[i])<-inprod(X2_es[i,],gamma_e[])
      m_cost2[i]~dbern(pc_2[i])
      logit(pc_2[i])<-inprod(X2_cs[i,],gamma_c[])
      }
      
      #obtain mean values for eff and cost
      mu_c_t[1]<-inprod(mean_cov_c1_t[],beta_c[,1])
      mu_c_t[2]<-inprod(mean_cov_c2_t[],beta_c[,2])
      mu_e_t[1]<-inprod(mean_cov_e1_t[],beta_e[,1])
      mu_e_t[2]<-inprod(mean_cov_e2_t[],beta_e[,2])
      #missingness probability
      #p_c[1]<-ilogit(inprod(mean_cov_c1_t[],gamma_c[]))
      #p_c[2]<-ilogit(inprod(mean_cov_c2_t[],gamma_c[]))
      #p_e[1]<-ilogit(inprod(mean_cov_e1_t[],gamma_e[]))
      #p_e[2]<-ilogit(inprod(mean_cov_e2_t[],gamma_e[]))
      
      #transformation of parameters
      for (t in 1:2){
      #from precision to variance and standard deviation
      tau_c_t[t]<-1/ss_c_t[t]
      ss_c_t[t]<-s_c_t[t]*s_c_t[t]
      s_c_t[t]<-exp(ls_c_t[t])
      tau_e_t[t]<-1/ss_e_t[t]
      ss_e_t[t]<-s_e_t[t]*s_e_t[t]
      s_e_t[t]<-exp(ls_e_t[t])

      #obtain parameters on natural scale
      mu_e[t]<-mu_e_t[t]*sd_eff[t]*2+mean_eff[t]
      mu_c[t]<-mu_c_t[t]*sd_cost[t]*2+mean_cost[t]

      s_e[t]<-s_e_t[t]*sd_eff[t]*2+mean_eff[t]
      s_c[t]<-s_c_t[t]*sd_cost[t]*2+mean_cost[t]
      
      #priors
      #model for costs and effects
      ls_c_t[t]~dunif(-5,10)         
      ls_e_t[t]~dunif(-5,10)
      }
      
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
      
      #missing data mechanism
      #logistic regression coefficients
      for (j in 1:pe){#begin gamma priors effects
      gamma_e[j]~dnorm(0,1)
      }#end gamma priors effects
      
      for (j in 1:pc){#begin gamma priors costs
      gamma_c[j]~dnorm(0,1)
      }#end gamma priors costs

      #retrieve original cost and effect data
      for(i in 1:N1){
      eff1[i]<-eff1_s[i]*sd_eff[1]*2+mean_eff[1]
      cost1[i]<-cost1_s[i]*sd_cost[1]*2+mean_cost[1]
      }
      for(i in 1:N2){
      eff2[i]<-eff2_s[i]*sd_eff[2]*2+mean_eff[2]
      cost2[i]<-cost2_s[i]*sd_cost[2]*2+mean_cost[2]
      }
    }
      "
      
      #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
      model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
    } else if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      #write model file in JAGS/BUGS code for MNAR
      model_string_jags<-  "
      model{
      #variables denoted by _s are related to their standardised version 
      #variables denoted by _t are are related to their standardised version that need to be back tranformed on their natural scale for conducting analysis
      
      #sd_ and mean_ are vectors containing the values, per treatment arm, of sd and mean necessary to back transform variables
      
      #control
      for(i in 1:N1){
      #costs and effects model
      cost1_s[i]~dnorm(mu_c_t[1],tau_c_t[1])
      eff1_s[i]~dnorm(mu_e_t[1],tau_e_t[1])
      #missing data mechanism
      m_eff1[i]~dbern(pq_1[i])
      logit(pq_1[i])<-gamma0_e+delta_e*(eff1_s[i])
      m_cost1[i]~dbern(pc_1[i])
      logit(pc_1[i])<-gamma0_c+delta_c*(cost1_s[i])
      }
      
      #intervention
      for(i in 1:N2){
      #costs and effects model
      cost2_s[i]~dnorm(mu_c_t[2],tau_c_t[2])
      eff2_s[i]~dnorm(mu_e_t[2],tau_e_t[2])
      #missing data mechanism
      m_eff2[i]~dbern(pq_2[i])
      logit(pq_2[i])<-gamma0_e+delta_e*(eff2_s[i])
      m_cost2[i]~dbern(pc_2[i])
      logit(pc_2[i])<-gamma0_c+delta_c*(cost2_s[i])
      }
      
      #transformation of parameters
      for (t in 1:2){
      #retrieve original variable scales
      
      mu_e[t]<-mu_e_t[t]*sd_eff[t]*2+mean_eff[t]
      mu_c[t]<-mu_c_t[t]*sd_cost[t]*2+mean_cost[t]
      
      s_e[t]<-s_e_t[t]*sd_eff[t]*2+mean_eff[t]
      s_c[t]<-s_c_t[t]*sd_cost[t]*2+mean_cost[t]
      
      #from precision to variance and standard deviation
      tau_c_t[t]<-1/ss_c_t[t]
      ss_c_t[t]<-s_c_t[t]*s_c_t[t]
      s_c_t[t]<-exp(ls_c_t[t])
      tau_e_t[t]<-1/ss_e_t[t]
      ss_e_t[t]<-s_e_t[t]*s_e_t[t]
      s_e_t[t]<-exp(ls_e_t[t])
      
      #missingness probability
      #p_e[1]<-ilogit(gamma0_e+delta_e*mean(eff1_s[]))
      #p_e[2]<-ilogit(gamma0_e+delta_e*mean(eff2_s[]))
      #p_c[1]<-ilogit(gamma0_c+delta_c*mean(cost1_s[]))
      #p_c[2]<-ilogit(gamma0_c+delta_c*mean(cost2_s[]))
      
      #priors
      #model for costs and effects
      ls_c_t[t]~dunif(-5,10)         
      ls_e_t[t]~dunif(-5,10)
      mu_c_t[t]~dnorm(0,0.0001)
      mu_e_t[t]~dnorm(0,0.0001)
      }
      
      #missing data mechanism
      gamma0_e~dlogis(0,1)
      gamma0_c~dlogis(0,1)
      #MNAR parameters
      delta_e~dnorm(0,1)
      delta_c~dnorm(0,1)
      
      #retrieve original cost and effect data
      for(i in 1:N1){
      eff1[i]<-eff1_s[i]*sd_eff[1]*2+mean_eff[1]
      cost1[i]<-cost1_s[i]*sd_cost[1]*2+mean_cost[1]
      }
      for(i in 1:N2){
      eff2[i]<-eff2_s[i]*sd_eff[2]*2+mean_eff[2]
      cost2[i]<-cost2_s[i]*sd_cost[2]*2+mean_cost[2]
      }
      
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
      cost1_s[i]~dnorm(mu_c1_t[i],tau_c_t[1])
      eff1_s[i]~dnorm(mu_e1_t[i],tau_e_t[1])
      #mean regression
      mu_c1_t[i]<-inprod(X1_cs[i,],beta_c[,1])
      mu_e1_t[i]<-inprod(X1_es[i,],beta_e[,1])
      #missing data mechanism
      m_eff1[i]~dbern(pq_1[i])
      logit(pq_1[i])<-inprod(X1_es[i,],gamma_e[])+delta_e*(eff1_s[i])
      m_cost1[i]~dbern(pc_1[i])
      logit(pc_1[i])<-inprod(X1_cs[i,],gamma_c[])+delta_c*(cost1_s[i])
      }
      
      #intervention
      for(i in 1:N2){
      #costs and effects model
      cost2_s[i]~dnorm(mu_c2_t[i],tau_c_t[2])
      eff2_s[i]~dnorm(mu_e2_t[i],tau_e_t[2])
      #mean regression
      mu_c2_t[i]<-inprod(X2_cs[i,],beta_c[,2])
      mu_e2_t[i]<-inprod(X2_es[i,],beta_e[,2])
      #missing data mechanism
      m_eff2[i]~dbern(pq_2[i])
      logit(pq_2[i])<-inprod(X2_es[i,],gamma_e[])+delta_e*(eff2_s[i])
      m_cost2[i]~dbern(pc_2[i])
      logit(pc_2[i])<-inprod(X2_cs[i,],gamma_c[])+delta_c*(cost2_s[i])
      }
      
      #obtain mean values for eff and cost
      mu_c_t[1]<-inprod(mean_cov_c1_t[],beta_c[,1])
      mu_c_t[2]<-inprod(mean_cov_c2_t[],beta_c[,2])
      mu_e_t[1]<-inprod(mean_cov_e1_t[],beta_e[,1])
      mu_e_t[2]<-inprod(mean_cov_e2_t[],beta_e[,2])
      
      #missingness probability
      #p_e[1]<-ilogit(inprod(mean_cov_e1[],gamma_e[])+delta_e*mean(eff1_s[]))
      #p_e[2]<-ilogit(inprod(mean_cov_e2[],gamma_e[])+delta_e*mean(eff2_s[]))
      #p_c[1]<-ilogit(inprod(mean_cov_c1[],gamma_c[])+delta_c*mean(cost1_s[]))
      #p_c[2]<-ilogit(inprod(mean_cov_c2[],gamma_c[])+delta_c*mean(cost2_s[]))
      
      #transformation of parameters
      for (t in 1:2){
      #from precision to variance and standard deviation
      tau_c_t[t]<-1/ss_c_t[t]
      ss_c_t[t]<-s_c_t[t]*s_c_t[t]
      s_c_t[t]<-exp(ls_c_t[t])
      tau_e_t[t]<-1/ss_e_t[t]
      ss_e_t[t]<-s_e_t[t]*s_e_t[t]
      s_e_t[t]<-exp(ls_e_t[t])
      
      #obtain parameters on natural scale
      mu_e[t]<-mu_e_t[t]*sd_eff[t]*2+mean_eff[t]
      mu_c[t]<-mu_c_t[t]*sd_cost[t]*2+mean_cost[t]
      
      s_e[t]<-s_e_t[t]*sd_eff[t]*2+mean_eff[t]
      s_c[t]<-s_c_t[t]*sd_cost[t]*2+mean_cost[t]
      
      #priors
      #model for costs and effects
      ls_c_t[t]~dunif(-5,10)         
      ls_e_t[t]~dunif(-5,10)
      }
      
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
      
      #missing data mechanism
      #MNAR parameters
      delta_e~dnorm(0,1)
      delta_c~dnorm(0,1)
      
      #logistic regression coefficients
      for (j in 1:pe){#begin gamma priors effects
      gamma_e[j]~dnorm(0,1)
      }#end gamma priors effects
      
      for (j in 1:pc){#begin gamma priors costs
      gamma_c[j]~dnorm(0,1)
      }#end gamma priors costs

      #retrieve original cost and effect data
      for(i in 1:N1){
      eff1[i]<-eff1_s[i]*sd_eff[1]*2+mean_eff[1]
      cost1[i]<-cost1_s[i]*sd_cost[1]*2+mean_cost[1]
      }
      for(i in 1:N2){
      eff2[i]<-eff2_s[i]*sd_eff[2]*2+mean_eff[2]
      cost2[i]<-cost2_s[i]*sd_cost[2]*2+mean_cost[2]
      }
      }
      "
      if(type=="MNAR_eff_cov"){
        model_string_jags<-gsub("delta_c~dnorm(0,1)","delta_c<-0",model_string_jags,fixed=TRUE)
      } else if(type=="MNAR_cost_cov"){
        model_string_jags<-gsub("delta_e~dnorm(0,1)","delta_e<-0",model_string_jags,fixed=TRUE)
      }
    }
  #if predictor is one BUGS needs a different specification for inpprod function
  if(pe==1){
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      inprod_e1<-"X1_es[i]*beta_e[1]"
      inprod_e2<-"X2_es[i]*beta_e[2]"
      inprod_pq1<-"X1_es[i]*gamma_e"
      inprod_pq2<-"X2_es[i]*gamma_e"
      inprod_mean_e1<-"mean_cov_e1_t*beta_e[1]"
      inprod_mean_e2<-"mean_cov_e2_t*beta_e[2]"
      inprod_mean_pq1<-"mean_cov_e1_t*gamma_e"
      inprod_mean_pq2<-"mean_cov_e2_t*gamma_e"
      begin_prior_beta<-"#begin beta priors effects"
      prior_beta_e<-"#"
      prior_beta_e1_inter<-"beta_e[1]~dnorm(0,0.001)"
      prior_beta_e2_inter<-"beta_e[2]~dnorm(0,0.001)"
      end_prior_beta<-"#end beta priors effects"
      begin_prior_gamma<-"#begin gamma priors effetcs"
      prior_gamma_e<-"gamma_e~dnorm(0,1)"
      end_prior_gamma<-"#end gamma priors effects"
      model_string_jags<-gsub("inprod(X1_es[i,],beta_e[,1])", inprod_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2_es[i,],beta_e[,2])", inprod_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X1_es[i,],gamma_e[])", inprod_pq1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2_es[i,],gamma_e[])", inprod_pq2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_e1_t[],beta_e[,1])", inprod_mean_e1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_e2_t[],beta_e[,2])", inprod_mean_e2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_e1_t[],gamma_e[])", inprod_mean_pq1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_e2_t[],gamma_e[])", inprod_mean_pq2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:pe) {#begin beta priors effects", begin_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(t in 1:2){beta_e[j,t]~dnorm(0,0.01)}", prior_beta_e, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_e[1,1]~dnorm(0,0.001)", prior_beta_e1_inter, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_e[1,2]~dnorm(0,0.001)", prior_beta_e2_inter, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end beta priors effects", end_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 1:pe){#begin gamma priors effects", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_e[j]~dnorm(0,1)", prior_gamma_e, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors effects", end_prior_gamma, model_string_jags,fixed=TRUE)
    }
  }
  if(pc==1){
    if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
      inprod_c1<-"X1_cs[i]*beta_c[1]"
      inprod_c2<-"X2_cs[i]*beta_c[2]"
      inprod_pc1<-"X1_cs[i]*gamma_c"
      inprod_pc2<-"X2_cs[i]*gamma_c"
      inprod_mean_c1<-"mean_cov_c1_t*beta_c[1]"
      inprod_mean_c2<-"mean_cov_c2_t*beta_c[2]"
      inprod_mean_pc1<-"mean_cov_c1_t*gamma_c"
      inprod_mean_pc2<-"mean_cov_c2_t*gamma_c"
      begin_prior_beta<-"#begin beta priors costs"
      prior_beta_c<-"#"
      prior_beta_c1_inter<-"beta_c[1]~dnorm(0,0.001)"
      prior_beta_c2_inter<-"beta_c[2]~dnorm(0,0.001)"
      end_prior_beta<-"#end beta priors costs"
      begin_prior_gamma<-"#begin gamma priors costs"
      prior_gamma_c<-"gamma_c~dnorm(0,1)"
      end_prior_gamma<-"#end gamma priors costs"
      model_string_jags<-gsub("inprod(X1_cs[i,],beta_c[,1])", inprod_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2_cs[i,],beta_c[,2])", inprod_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X1_cs[i,],gamma_c[])", inprod_pc1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(X2_cs[i,],gamma_c[])", inprod_pc2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c1_t[],beta_c[,1])", inprod_mean_c1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c2_t[],beta_c[,2])", inprod_mean_c2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c1_t[],gamma_c[])", inprod_mean_pc1, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("inprod(mean_cov_c2_t[],gamma_c[])", inprod_mean_pc2, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 2:pc) {#begin beta priors costs", begin_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for(t in 1:2){beta_c[j,t]~dnorm(0,0.01)}", prior_beta_c, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c[1,1]~dnorm(0,0.001)", prior_beta_c1_inter, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("beta_c[1,2]~dnorm(0,0.001)", prior_beta_c2_inter, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end beta priors costs", end_prior_beta, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("for (j in 1:pc){#begin gamma priors costs", begin_prior_gamma, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("gamma_c[j]~dnorm(0,1)", prior_gamma_c, model_string_jags,fixed=TRUE)
      model_string_jags<-gsub("}#end gamma priors costs", end_prior_gamma, model_string_jags,fixed=TRUE)
    }
  }
  #call prior_change function to make changes to prior values and distributions if prior inputs in run_model are provided
  model_string_jags<-prior_change(type=type,dist_e ="norm",dist_c="norm")
     #assign name to model file to be called by run_jags or run_bugs
     writeLines(model_string_jags,"normal_model_ind_scaled.txt")
     model_string<-"normal_model_ind_scaled.txt"                         
  return(model_string)
}))
