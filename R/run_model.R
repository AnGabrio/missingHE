#' Full Bayesian Models to handle missingness in Health Economic Evaluations
#' 
#' Full Bayesian cost-effectiveness models to handle missing data in the outcomes under different missing data 
#' mechanism assumptions and using a variatey of alternative parametric distributions for the effect and cost variables and 
#' using a selection model to specify the missingness mechanism. The analysis is performed using the \code{BUGS} language, 
#' which is implemented either in \code{JAGS} or \code{OpenBUGS} using the functions \code{\link[R2jags]{jags}} and 
#' \code{\link[R2OpenBUGS]{bugs}}, respectively. The output is stored in an object of class 'missingHE'.
#' 
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively.
#' @param model.eff A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#' @param model.cost A formula expression in conventional R linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model.
#' @param type Type of missingness mechanism assumed. Choices are: Missing Completely At Random (MCAR),
#'  Missing At Random (MAR), Missing Not At Random (MNAR). Different 'MNAR' alternative versions are available 
#'  depending on whether the mechanism for only one or both outcomes is considered. Specifically it is possible
#'  to select MNAR for the effect variables ('MNAR_eff'), for the cost variables ('MNAR_cost') or both ('MNAR').
#'  If covariate data are provided in \code{data}, \code{run_model} automatically includes them in the model.
#' @param dist_e Distribution assumed for the effects. Current available chocies are: Normal ("norm") and Beta ("beta").
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ("norm") and Gamma ("gamma").
#' @param program type of software used to run the model. Current alternatives are 'OpenBUGS'('BUGS') and 'JAGS' ('JAGS').
#' @param forward Logical. If \code{forward} is \code{TRUE}, the model is run in forward sampling mode
#'  without using any data, else if \code{forward} is \code{FALSE} standard sampling mode is selected.
#' @param save_model Logical. If \code{save_model} is \code{TRUE} a \code{txt} file containing the model code is printed
#'  in the current working directory.
#' @param prob A numeric vector of probabilities within the range (0,1), representing the upper and lower
#'  CI sample quantiles to be calculated and returned for the imputed values.
#' @param n.chains Number of chains.
#' @param n.burnin Number of warmup iterations.
#' @param n.iter Number of iterations.
#' @param n.thin Thinning interval.
#' @param inits A list with elements equal to the number of chains selected; each element of the list is itself a list of starting values for the
#' \code{BUGS/JAGS} model, or a function creating (possibly random) initial values. If \code{inits} is \code{NULL}, \code{BUGS/JAGS}
#'  will generate initial values for all the model parameters.
#' @param prior A list containing the hyperprior values provided by the user. Each element of this list must be a vector of length two
#' containing the user-provided hyperprior values and must be named with the name of the corresponding parameter. For example, the hyperprior
#' values for the mean effect parameter, when a normal distribution is assumed, can be provided using \code{prior=list('mu.prior.e'=c(0,1))}.
#' For more information about how to provide prior hypervalues for each model parameter see details. If \code{prior} is 'default' the default values will be used.  
#' @param ... additional input parameters. Examples are \code{stand} (logical) and \code{transf}. When \code{stand} is \code{TRUE} outcomes are scaled so to have mean \code{0} and
#' standard deviation \code{0.5} (only for bivariate normal). If \code{transf} is set to 'logit', mean effect variables are modelled on the log-odds scale 
#' (only for beta distribution), while if 'log' is chosen, mean cost variables are modelled on the log scale (only for gamma distribution). If a beta-gamma
#' model is selected both transformations can be carried out by setting \code{transf} equal to a vector that contains both transformation string names.
#' Other additional arguments contained in the function \code{\link[BCEA]{bcea}} can be provided. 
#' @return An object of the class 'missingHE' containing the following elements
#' \describe{
#'   \item{data_set}{A list containing the original data set provided in \code{data} (see Arguments), the number of observed and missing individuals 
#'   and the total number of individuals by treatment arm}
#'   \item{model_output}{A list containing the output of a \code{BUGS/JAGS} model generated from the functions \code{\link[R2jags]{jags}} 
#'   or \code{\link[R2OpenBUGS]{bugs}}}
#'   \item{mean_effects}{A matrix with \code{nsim} rows and \code{2} columns containing the posterior samples for the mean effect parameters in the
#'   treatment arms}
#'   \item{mean_costs}{A matrix with \code{nsim} rows and \code{2} columns containing the posterior samples for the mean cost parameters in the
#'   treatment arms}
#'   \item{sd_effects}{A matrix with \code{nsim} rows and \code{2} columns containing the posterior samples for the standard 
#'   deviation effect parameters in the treatment arms}
#'   \item{sd_costs}{A matrix with \code{nsim} rows and \code{2} columns containing the posterior samples for the standard 
#'   deviation cost parameters in the treatment arms}
#'   \item{imputed}{A list containing the the posterior samples for the imputed individuals in each arm and for each outcome.
#'   The stored imputed values represent the mean and the upper and lower quantiles of the posterior distribution for the
#'   missing individuals according to the values defined in \code{prob} (see Arguments)}
#'   \item{type}{A character variable that indicate which type of program has been used to run the model, either \code{JAGS} or \code{BUGS}}
#' }
#' @seealso \code{\link[R2OpenBUGS]{bugs}}, \code{\link[R2jags]{jags}}, \code{\link[BCEA]{bcea}}
#' @keywords CEA OpenBUGS JAGS Missingness
#' @importFrom stats model.frame 
#' @details Depending on the distributional assumptions specified for the outcome variables in the arguments \code{dist_e} and
#' \code{dist_c} (Model of Analysis, MoA) and the type of missingness mechanism assumed in the argument \code{type} (Model of Missingness, MoM)
#' , different types of models are built and run in the background by \code{run_model}.
#' Specifically, the MoM module is handled using a Selection Model (SM) appraoch, which can be summarised with a simple example.
#' Consider a data set comprising a partially observed response variable \eqn{y}, the corresponding missing data indicator vector \eqn{m}
#' , and a fully-observed covariate \eqn{x}. Under the SM approach, the joint distribution \eqn{p(y,m)} is factored as the product of 
#' the marginal distribution \eqn{p(y)} and the conditional distribution \eqn{p(m|y)}.
#' \deqn{p(y,m|x,\theta^MoA,\theta^MoM)=p(y|x,\theta^MoA)p(m|y,x,\theta^MoM)}
#' where, \eqn{\theta^MoA},\eqn{\theta^MoM} are respectively the main parameters of interest (in the MoA) 
#' and the parameters associated with the missing data mechanism (in the MoM).
#' While the MoA specification depends on the assumed distributional for for \eqn{y}, the MoM can have a more generalised structure.
#' \code{run_model} assigns a Bernoulli probability distribution to the missing data indicator \eqn{m ~ Bernoulli(\pi)}, 
#' where \eqn{\pi} is the parameter specifying the probability of \eqn{y} being missing. Such probability is modelled using
#' a logistic regression transformation. More specifically:
#' \deqn{logit(\pi)=\gamma_0+\gamma_1 x+\delta y}
#' where
#' \itemize{
#' \item \eqn{\gamma_0} represents the baseline probability of missingness in \eqn{y} that does not depend on any variable. When \eqn{\gamma_1=\delta=0}, the model assumes MCAR.
#' \item \eqn{\gamma_1} represents the impact on the probability of missingness in \eqn{y} of the fully observed covariate \eqn{x}. When \eqn{\delta=0}, the model assumes MAR.
#' \item \eqn{\delta} represents the impact on the probability of missingness in \eqn{y} of the possibly unobserved values in \eqn{y} (MNAR).
#' }
#' 
#' It is important to keep in mind both the model structure of the MoA and MoM when providing user-defined hyperprior values for the parameters indexing such components
#' using the argument \code{prior} in the function \code{run_model}. Specifically, for the MoA module the default prior distributions that can be overwritten by the user
#' are the following:
#' \itemize{
#' \item Normal. Priors can be supplied for the location (mean) and auxiliary (log-sd) parameter with default values: \eqn{\mu ~ Normal(0,0.00001), \alpha~Uniform(-5,10)}
#' \item Beta.  Priors can be supplied for the location (mean) and auxiliary (sd) parameter with default values: \eqn{\mu ~ Uniform(0,1), \alpha~Uniform(0,\sqrt\mu(1-\mu))}
#' \item Gamma.  Priors can be supplied for the location (mean) and auxiliary (sd) parameter with default values: \eqn{\mu ~ Uniform(0,10000), \alpha~Cauchy(0,0.16)T(0,)}
#' }
#' When covariate data are included in the model with a linear regression to the mean of the outcome variable, then instead of specifying priors on \eqn{\mu} 
#' it is necessary to define the priors on the coefficient parameters of the linear regression such as
#' \deqn{\mu=\beta_0+\sum\beta_j X_j}
#' where \eqn{\beta_0} is the marginal mean to which the default prior is the same as that specified on \eqn{\mu} if no covariate is included, while
#' \eqn{\beta_j} are the covariate coefficients whose defualt prior is the same foe all j: \eqn{\beta_j~Normal(0,0.000001)}. 
#' For the MoM parameters the default prior distributions assumed are the following:
#' \itemize{
#' \item \eqn{gamma_0~Logisitc(0,1)}
#' \item \eqn{gamma_j~Normal(0,1)}
#' \item \eqn{delta~Normal(0,1)}
#' }
#' When user-defined hyperprior values are supplied via the argument \code{prior} in the function \code{run_model}, the list elements inside such an object (see Arguments)
#' must be vectors of length 2 containing the user-provided hyperprior values and must take specific names according to the parameters whose priors the user wants to modify. 
#' Specifically, the names accepted by \strong{missingHE} are the following:
#' \itemize{
#' \item location parameters in the MoA: "mean.prior.e"(effects) and/or "mean.prior.c"(costs)
#' \item auxiliary parameters in the MoA: "alpha.prior.e"(effects) and/or "alpha.prior.c"(costs)
#' \item marginal mean parameter (if covariate data provided) in the MoA: "beta0.prior.e"(effects) and/or "beta0.prior.c"(costs)
#' \item covariate parameters (if covariate data provided) in the MoA: "beta.prior.e"(effects) and/or "beta.prior.c"(costs)
#' \item baseline parameter in the MoM: "gamma0.prior.e"(effects) and/or "gamma0.prior.c"(costs)
#' \item covariate parameters in the MoM (if covariate data provided) in the MoM: "gamma.prior.e"(effects) and/or "gamma.prior.c"(costs)
#' \item mnar parameter in the MoM: "delta.prior.e"(effects) and/or "delta.prior.c"(costs)
#' } 
#' To notice that the priors for the marginal mean and covariate parameters can be provided only if covariate data are inlcuded in the model, with the
#' marginal mean parameter being a substitute of the location parameter under such model framework.
#' 
#' 
#' 
#' @author Andrea Gabrio
#' @references  
#' Daniels, MJ. Hogan, JW. \emph{Missing Data in Longitudinal Studies: strategies for Bayesian modelling and sensitivity analysis}, CRC/Chapman Hall.
#' 
#' Baio, G.(2012). \emph{Bayesian Methods in Health Economics}. CRC/Chapman Hall, London.
#' 
#' Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). \emph{Bayesian Data Analysis, 2nd edition}, CRC Press.
#' 
#' Plummer, M. \emph{JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling.} (2003).
#' @export
#'
#' @examples
#'#Simple example to simulate and analyse a data set
#'#Define the number of individuals fer treatment arm
#'N1<-150
#'N2<-100
#'#Create the misisngness indicators totally random (MCAR mechanism)
#'m_eff1<-m_cost1<-rbinom(N1,1,0.25)
#'m_eff2<-m_cost2<-rbinom(N2,1,0.25)
#'#Simulate data from normal distributions for both arms
#'eff1<-rnorm(N1,0.5,0.5)
#'eff2<-rnorm(N2,0.5,0.5)
#'cost1<-rnorm(N1,90,20)
#'cost2<-rnorm(N2,90,20)
#'#Set value missing if indicator is 1
#'eff1[m_eff1==1]<-NA
#'eff2[m_eff2==1]<-NA
#'cost1[m_cost1==1]<-NA
#'cost2[m_cost2==1]<-NA
#'#Create treatment arm indicators
#'t1<-rep(1,length(eff1))
#'t2<-rep(2,length(eff2))
#'#Combine variables and define a data set
#'e<-c(eff1,eff2)
#'c<-c(cost1,cost2)
#'m_eff<-c(m_eff1,m_eff2)
#'m_cost<-c(m_cost1,m_cost2)
#'t<-c(t1,t2)
#'data<-data.frame(e,c,t)
#'
#'#Run the model using run_model with JAGS assuming a MCAR missingness mechanism
#'x<-run_model(data=data,model.eff=e~1,model.cost=c~1,
#'dist_e="norm",dist_c="norm",type="MCAR",program="JAGS")
#'#print the results of the JAGS/BUGS model
#'print(x)
#'#
#'#Assess model convergence using graphical tools
#'#Produce histograms of the posterior samples for the mean effect
#'#parameters in the two treatment arms. 
#'diagnostic_checks(x,type = "histogram",param = "mu.e")
#'#
#'#Compare observed outcome data with imputations from the model
#'# (posteiror means and credible intervals)
#'plot(x,class="scatter",outcome="all")
#'#
#'#Summarise the CEA information from model results
#'summary(x)
#'#
#'#


run_model<-function(data,model.eff,model.cost,dist_e,dist_c,type,program="JAGS",forward=FALSE,prob=c(0.05,0.95),
                    n.chains=2,n.iter=20000,n.burnin=floor(n.iter/2),inits=NULL,n.thin=1,
                    save_model=FALSE,prior="default",...){
  #prevent global error message
  filein<-NULL
  if(is.data.frame(data)==FALSE){
    stop("data must be in data frame format")
  }
  if(!any(c("e","c","t")%in% names(data))==TRUE){
    stop("Please rename or provide variables in the data as 'e', 'c' and 't' for the effectiveness, cost and treatment indicator")
  }
  if(any(names(data)=="e")==TRUE & any(names(data)=="c")==TRUE){
    e<-as.name("e")
    c<-as.name("c")
  }
  cov_matrix<-subset(data, select = -c(e,c))
  if(any(is.na(cov_matrix))==TRUE){
    stop("no missing covariate or treatment indicator is allowed")
  }
  if(any(levels(as.factor(cov_matrix$t))!=c("1","2"))==TRUE){
    stop("A two arm indicator variable must be provided")
  }
  #need to specify type of missingness mechanism among those available
  if(!type %in% c("MCAR","MAR","MNAR_eff","MNAR_cost","MNAR")) {
    stop("Types available for use are 'MCAR','MAR','MNAR_eff','MNAR_cost','MNAR'")
  }
  #need to specify distributions for effects and costs among those available
  if(!dist_e %in% c("norm","beta")|!dist_c %in% c("norm","gamma")) {
    stop("Distributions available for use are 'norm','beta' for the effects and 'norm','gamma' for the costs")
  }
  #return value for standard or forward sampling 
  model_class<-NULL
  if(forward==TRUE){
    model_class<-"forward"
  } else{model_class<-"standard"}
  #need to specify software among those available
  if(!program %in% c("JAGS","BUGS")) {
    stop("Programs available for use are 'JAGS','BUGS'")
  }
  #quantile bounds must be valid
  if(length(prob)!=2|is.numeric(prob)==FALSE|any(prob<0)!=FALSE|any(prob>1)!=FALSE){
    stop("You must provide valid lower/upper quantiles for the imputed data distribution")}
  #call read_data to read and output the variables to be modelled
  #and define each variable to be used in the model
  data_read<-read_data(data=data,model.eff=model.eff,model.cost=model.cost)
  N1<-data_read$data_raw$arm_lengths[1]
  N2<-data_read$data_raw$arm_lengths[2]
  #number of predictors
  pe<-ncol(data_read$data_raw$covariates_effects$Intervention)
  pc<-ncol(data_read$data_raw$covariates_costs$Intervention)
  if(pe==1 & pc==1){
   if(type=="MCAR"|type=="MAR"){
     type="MCAR"
   }
  }else if(pe>1 | pc>1){
    if(type=="MNAR"|type=="MNAR_eff"|type=="MNAR_cost"){
      type=paste(type,"cov",sep = "_")
    }
  }
  #missing data indicators
  m_eff1<-data_read$data_raw$missing_effects$Control
  m_eff2<-data_read$data_raw$missing_effects$Intervention
  m_cost1<-data_read$data_raw$missing_costs$Control
  m_cost2<-data_read$data_raw$missing_costs$Intervention
  #outcomes
  eff1<-data_read$data_raw$raw_effects$Control
  eff2<-data_read$data_raw$raw_effects$Intervention
  cost1<-data_read$data_raw$raw_costs$Control
  cost2<-data_read$data_raw$raw_costs$Intervention
  #standardised outcome
  eff1_s<-data_read$data_stand$stand_effects$Control
  eff2_s<-data_read$data_stand$stand_effects$Intervention
  cost1_s<-data_read$data_stand$stand_costs$Control
  cost2_s<-data_read$data_stand$stand_costs$Intervention
  mean_eff = data_read$data_stand$mean_effects
  mean_cost = data_read$data_stand$mean_costs
  sd_eff<-data_read$data_stand$sd_effects
  sd_cost<-data_read$data_stand$sd_costs
  eff1_cc_s<-data_read$data_stand$stand_effects_cc$Control
  eff2_cc_s<-data_read$data_stand$stand_effects_cc$Intervention
  cost1_cc_s<-data_read$data_stand$stand_costs_cc$Control
  cost2_cc_s<-data_read$data_stand$stand_costs_cc$Intervention
  #number of observations and missing data
  N1_cc<-data_read$data_raw$arm_lengths_cc[,1]
  N2_cc<-data_read$data_raw$arm_lengths_cc[,2]
  N1_mis<-data_read$data_raw$arm_missing_data[,1]
  N2_mis<-data_read$data_raw$arm_missing_data[,2]
  #covariates (vectors of 1s if no variables)
  #standardised
  X1_es<-as.matrix(data_read$data_stand$covariates_effects$Control)
  X2_es<-as.matrix(data_read$data_stand$covariates_effects$Intervention)
  X1_cs<-as.matrix(data_read$data_stand$covariates_costs$Control)
  X2_cs<-as.matrix(data_read$data_stand$covariates_costs$Intervention)
  if(pe==1){
    X1_es<-as.vector(X1_es)
    X2_es<-as.vector(X2_es)
  }
  if(pc==1){
    X1_cs<-as.vector(X1_cs)
    X2_cs<-as.vector(X2_cs)
  }
  mean_cov_e1_t<-as.vector(data_read$data_stand$mean_cov_effects$Control)
  mean_cov_e2_t<-as.vector(data_read$data_stand$mean_cov_effects$Intervention)
  mean_cov_c1_t<-as.vector(data_read$data_stand$mean_cov_costs$Control)
  mean_cov_c2_t<-as.vector(data_read$data_stand$mean_cov_costs$Intervention)
  #original
  X1_e<-as.matrix(data_read$data_raw$covariates_effects$Control)
  X2_e<-as.matrix(data_read$data_raw$covariates_effects$Intervention)
  X1_c<-as.matrix(data_read$data_raw$covariates_costs$Control)
  X2_c<-as.matrix(data_read$data_raw$covariates_costs$Intervention)
  if(pe==1){
    X1_e<-as.vector(X1_e)
    X2_e<-as.vector(X2_e)
  }
  if(pc==1){
    X1_c<-as.vector(X1_c)
    X2_c<-as.vector(X2_c)
  }
  mean_cov_e1<-as.vector(data_read$data_raw$mean_cov_effects$Control)
  mean_cov_e2<-as.vector(data_read$data_raw$mean_cov_effects$Intervention)
  mean_cov_c1<-as.vector(data_read$data_raw$mean_cov_costs$Control)
  mean_cov_c2<-as.vector(data_read$data_raw$mean_cov_costs$Intervention)
  #define additional inputs (prior changes) as a list
  exArgs <- list(...)
  #define additional inputs (prior changes) as a list
  #if additional inputs are provided
  if(any(prior=="default")==TRUE){
    prior<-list(default="default")
    }else if(any(prior=="default")==FALSE){
    #check that they are vectors 
    list_check_vector<-lapply(prior,is.vector)
    #stop if any of the additional input provided is not a vector
    if(all(as.logical(list_check_vector))==FALSE){
      stop("all user-supplied priors should be in vector format")
    }
    #set of all possible parameters to which priors can be assigned
    par_prior<-c("mu.prior.e","mu.prior.c","alpha.prior.e","alpha.prior.c","gamma0.prior.e","gamma0.prior.c","gamma.prior.e","gamma.prior.c",
                 "beta0.prior.e","beta0.prior.c","beta.prior.e","beta.prior.c","theta.prior","delta.prior.e","delta.prior.c")
    #stop message
    stop_mes<-"priors can be assigned only using specific string parameter names depending on the type of model assumed. Type ''help(run_model)'' for more details"
    if(type=="MCAR"){if(!any(names(list_check_vector) %in% par_prior[c(1:6)]==TRUE))stop(stop_mes)}
    if(type=="MNAR"){if(!any(names(list_check_vector) %in% par_prior[c(1:6,14:15)]==TRUE))stop(stop_mes)}
    if(type=="MNAR_eff"){if(!any(names(list_check_vector) %in% par_prior[c(1:6,14)]==TRUE))stop(stop_mes)}
    if(type=="MNAR_cost"){if(!any(names(list_check_vector) %in% par_prior[c(1:6,15)]==TRUE))stop(stop_mes)}
    if(type=="MAR"){if(!any(names(list_check_vector) %in% par_prior[c(3:12)]==TRUE))stop(stop_mes)}
    if(type=="MNAR_cov"){if(!any(names(list_check_vector) %in% par_prior[c(3:12,14:15)]==TRUE))stop(stop_mes)}
    if(type=="MNAR_eff_cov"){if(!any(names(list_check_vector) %in% par_prior[c(3:12,14)]==TRUE))stop(stop_mes)}
    if(type=="MNAR_cost_cov"){if(!any(names(list_check_vector) %in% par_prior[c(3:12,15)]==TRUE))stop(stop_mes)}
      if(length(exArgs$ind)!=0){
        if(exists("ind",where=exArgs) & exArgs$ind!=TRUE){
          if("theta.prior" %in% names(list_check_vector)){stop(stop_mes)}}}
    }
  #include transformation if selected
  if(exists("transf",where=exArgs)) {transf=exArgs$transf} else {transf="default"}
  if(!any(transf %in% c("logit","log","default")==TRUE)){
    stop("Only 'log' or 'logit' transformations available")
  }
  if(any(transf=="logit")==TRUE & dist_e!="beta"){
    stop("Logit transformation available only for beta distribution for the effects")
  }
  if(any(transf=="log")==TRUE & dist_c!="gamma"){
    stop("Log transformation available only for gamma distribution for the costs")
  }
  #set prior input names and if not provided set the values to null (MCAR,MAR,MNAR)
  #MoA parameters for different types of distributions
  #all
  if(exists("alpha.prior.e",where=prior)) {alpha.prior.e=prior$alpha.prior.e} else {alpha.prior.e=NULL}
  if(exists("alpha.prior.c",where=prior)) {alpha.prior.c=prior$alpha.prior.c} else {alpha.prior.c=NULL}
  if(exists("gamma0.prior.e",where=prior)) {gamma0.prior.e=prior$gamma0.prior.e} else {gamma0.prior.e=NULL}
  if(exists("gamma0.prior.c",where=prior)) {gamma0.prior.c=prior$gamma0.prior.c} else {gamma0.prior.c=NULL}
  #MNAR specific parameters (both mechanisms) 
  if(type=="MCAR"|type=="MNAR"){
    if(exists("mu.prior.e",where=prior)) {mu.prior.e=prior$mu.prior.e} else {mu.prior.e=NULL}
    if(exists("mu.prior.c",where=prior)) {mu.prior.c=prior$mu.prior.c} else {mu.prior.c=NULL}
  }
  #MNAR specific parameters (effect mechanism)
  if(type=="MNAR"|type=="MNAR_cov"|type=="MNAR_eff"|type=="MNAR_eff_cov"){
    if(exists("delta.prior.e",where=prior)) {delta.prior.e=prior$delta.prior.e} else {delta.prior.e=NULL}
  }
  #MNAR specific parameters (cost mechanism)
  if(type=="MNAR"|type=="MNAR_cov"|type=="MNAR_cost"|type=="MNAR_cost_cov"){
    if(exists("delta.prior.c",where=prior)) {delta.prior.c=prior$delta.prior.c} else {delta.prior.c=NULL}
  }
  #covariate specfic parameters
  if(type=="MAR"|type=="MNAR_cov"|type=="MNAR_eff_cov"|type=="MNAR_cost_cov"){
    if(exists("beta0.prior.e",where=prior)) {beta0.prior.e=prior$beta0.prior.e} else {beta0.prior.e=NULL}
    if(exists("beta0.prior.c",where=prior)) {beta0.prior.c=prior$beta0.prior.c} else {beta0.prior.c=NULL}
    if(exists("beta.prior.e",where=prior)) {beta.prior.e=prior$beta.prior.e} else {beta.prior.e=NULL}
    if(exists("beta.prior.c",where=prior)) {beta.prior.c=prior$beta.prior.c} else {beta.prior.c=NULL}
    if(exists("gamma.prior.e",where=prior)) {gamma.prior.e=prior$gamma.prior.e} else {gamma.prior.e=NULL}
    if(exists("gamma.prior.c",where=prior)) {gamma.prior.c=prior$gamma.prior.c} else {gamma.prior.c=NULL}
  }
  #standardisation and joint assumption available for normal only
  #default values for stand and ind in the case of normal assumptions
  stand=FALSE
  corr_assumption <- model.frame(formula=model.cost, data=data)
  if("e"%in%names(corr_assumption)){
  ind=FALSE  
  } else{ind=TRUE}
  if(dist_e!="norm" | dist_c!="norm"){
    if(length(exArgs$stand)!=0){
    if(exists("stand",where=exArgs) & exArgs$stand!=FALSE) {
    stop("Standardised variables only available for bivariate normal distribution")
    }}
    if(ind==FALSE){
      stop("Joint assumption only available for bivariate normal distribution")
    }
   }
  #bivariate normal is assumed
  if(dist_e=="norm" & dist_c=="norm"){
    #either standardise or not all variables
    if(exists("stand",where=exArgs)) {stand=as.logical(exArgs$stand)} else {stand=FALSE}
    #correlation parameter only if joint model
    if(exists("theta.prior",where=prior)) {theta.prior=prior$theta.prior} else {theta.prior=NULL}
    #variables are standardised before being modelled
    if(stand==TRUE){
      #create list of natural and standardised data output with total number of observed, missing and complete data
      data_set<-list("effects"=data_read$data_raw$raw_effects,"effects_stand"=data_read$data_stand$stand_effects ,
                     "costs"=data_read$data_raw$raw_costs,"costs_stand"=data_read$data_stand$stand_costs,"N in reference arm"=N1,"N in comparator arm"=N2,
                     "N observed in reference arm"=N1_cc,"N observed in comparator arm"=N2_cc,"N missing in reference arm"=N1_mis,"N missing in comparator arm"=N2_mis,
                     "covariates_stand_effects"=data_read$data_stand$covariates_effects,"covariates_stand_costs"=data_read$data_stand$covariates_costs,
                     "covariates_effects"=data_read$data_raw$covariates_effects,"covariates_costs"=data_read$data_raw$covariates_costs)
        #call run_jags to run the model in JAGS for each possible model setup 
        #call run_bugs to run the model in BUGS for each possible model setup 
        if(program=="JAGS"){
            model_output<-run_jags(type=type,dist_e=dist_e,dist_c=dist_c,forward=forward,inits=inits)
        } else if (program=="BUGS"){
            model_output<-run_bugs(type=type,dist_e=dist_e,dist_c=dist_c,forward=forward,inits=inits)
        }
    } else if(stand==FALSE){
      #variables are used on their natural scale
      #create list of natural scaled data output with total number of observed, missing and complete data
      data_set<-list("effects"=data_read$data_raw$raw_effects,"costs"=data_read$data_raw$raw_costs,"N in reference arm"=N1,"N in comparator arm"=N2,
                     "N observed in reference arm"=N1_cc,"N observed in comparator arm"=N2_cc,"N missing in reference arm"=N1_mis,"N missing in comparator arm"=N2_mis,
                     "covariates_effects"=data_read$data_raw$covariates_effects,"covariates_costs"=data_read$data_raw$covariates_costs)
        #call run_jags to run the model in JAGS for each possible model setup 
        #call run_bugs to run the model in BUGS for each possible model setup 
        if(program=="JAGS"){
            model_output<-run_jags(type=type,dist_e="norm",dist_c="norm",forward=forward,inits=inits)
        } else if(program=="BUGS"){
            model_output<-run_bugs(type=type,dist_e="norm",dist_c="norm",forward=forward,inits=inits)
        } 
      }
  } 
  if(dist_e!="norm" | dist_c!="norm"){
    #variables are used on their natural scale
    #create list of natural scaled data output with total number of observed, missing and complete data
    data_set<-list("effects"=data_read$data_raw$raw_effects,"costs"=data_read$data_raw$raw_costs,"N in reference arm"=N1,"N in comparator arm"=N2,
                   "N observed in reference arm"=N1_cc,"N observed in comparator arm"=N2_cc,"N missing in reference arm"=N1_mis,"N missing in comparator arm"=N2_mis,
                   "covariates_effects"=data_read$data_raw$covariates_effects,"covariates_costs"=data_read$data_raw$covariates_costs)
    #independence is assumed between effects and costs
    if(program=="JAGS"){
      #call run_jags to run the model in JAGS for each possible model setup 
      #call run_bugs to run the model in BUGS for each possible model setup 
        model_output<-run_jags(type=type,dist_e=dist_e,dist_c=dist_c,forward=forward,inits=inits)
    } else if(program=="BUGS"){
        model_output<-run_bugs(type=type,dist_e=dist_e,dist_c=dist_c,forward=forward,inits=inits)
    }
  }
  #delete from wd model file if save_model=FALSE (default)
  #unlink/remove file in directory for JAGS/BUGS
  if(save_model==FALSE & program=="JAGS"){
    unlink(filein)
  } else if(save_model==FALSE & program=="BUGS"){ 
    unlink(filein)
  }
  #perform CEA calling BCEA (only if forward sampling not selected)
  if(forward==FALSE){
    #set default argument values if not provided by user
    if(exists("ref",where=exArgs)) {ref=exArgs$ref} else {ref=2}
    if(exists("interventions",where=exArgs)) {interventions=exArgs$interventions} else {interventions=NULL}
    if(exists("Kmax",where=exArgs)) {Kmax=exArgs$Kmax} else {Kmax=50000}
    if(exists("wtp",where=exArgs)) {wtp=exArgs$wtp} else {wtp=NULL}
    if(exists("plot",where=exArgs)) {plot=exArgs$plot} else {plot=FALSE}
    cea<-BCEA::bcea(e=model_output$mean_effects,c=,model_output$mean_costs,ref=ref,interventions=interventions,Kmax=Kmax,wtp=wtp,plot=plot)
    #create list containing data list output and model list output
    res<-list(data_set=data_set,model_output=model_output,cea=cea,type=type,model_class=model_class)
  } else if(forward==TRUE){
    res<-list(data_set=data_set,model_output=model_output,type=type,model_class=model_class)
  }
  #assign new class name to the object
  class(res)<-"missingHE"
  #output.
  #data list: for each arm, natural or scaled variables, number of observed, missing and complete observations
  #model output list: model summary, posterior samples for main parameter of interests (means, sd, missingness probabilities)
  return(res)
}
