#' Diagnostic checks for assessing MCMC convergence 
#'
#' Convergence issues need to be assessed every time a  probabilistic model is implemented (via MCMC). The focus 
#' is restricted to full Bayesian models in cost-effectiveness analyses run using the function \code{\link{run_model}}, 
#' with convergence of the chains that is assessed through Visual checks of the posterior distribution of the parameters of interest,
#' such as density plots, trace plots, autocrrelation plots, etc. Other types of posterior checks are related to some summary MCMC statistics 
#' able to detect possible issues in the convergence of the algorithm, such as the potential scale reduction factor or the effective sample size.
#' \code{diagnostic_checks} uses many different types of diagnostic tools and statistics to assess model convergence from MCMC algorithms
#' using functions contained in the package \code{ggmcmc} and \code{mcmcplots}. Graphics display is managed and tailored using functions contained 
#' in the package \code{ggplot2} and \code{ggthemes}.
#' @keywords diagnostics MCMC convergence checks
#' @param x An object of class "missingHE" containing the posterior results of a full Bayesian model implemented usgin the function \code{\link{run_model}}
#' @param class type of diagnostic check to be plotted for model parameters, mostly taken from the packages \strong{ggmcmc} and \strong{mcmcplots}. Available choices include: 'histogram' for histrogram plots,
#' 'denplot' for density plots, 'traceplot' for trace plots and 'autocorrplot' for autocorrelation plots, 'running' for running mean plots,
#' 'compare' for comparing the distribution of the whole chain with only its last part, 'crosscorplot' for crosscorrelation plots, 'Rhat' for the potential scale reduction factor, 'geweke' for the geweke diagnostic,
#' 'pairs' for posterior correlation among the parameters,'caterpillar' for caterpillar plots. In addition the class 'summary' provides an overview of some of the most popular
#' diagnostic checks for each parameter selected.
#' @param parm Name of the family of parameters to process, as given by a regular expression. For example the mean parameters 
#' for the effect and cost variables can be specified using 'mu.e' and 'mu.c', respectively. Different types of
#' models may have different parameters depending on the assumed distributions and missing data mechanisms. 
#' To see a complete list of all possible paraemters by types of models assumed see details.
#' @param theme Type of ggplot theme among some pre-defined themes, mostly taken from the package \strong{ggthemes}. For a full list of available themes see details.
#' @param ... additional parameters that can be provided to manage the output of \code{diagnostic_checks}. 
#' These are mostly taken from the package \strong{ggmcmc} and are for graphical purposes.
#' @return A \code{ggplot} object containing the plots specified in the argument \code{class}
#' @seealso \code{\link[ggmcmc]{ggs}} \code{\link{run_model}}
#' @details Depending on the types of plots specified in the argument \code{class} the output of \code{diagnostic_checks} can produce
#' different combinations of MCMC visual posterior checks for the parameters of interest indicated in the argument \code{parm}.
#' For a full list of the available plots see the description of the argument \code{class} or see the corresponding plots in the package \strong{ggmcmc}.
#' 
#' The parameters that can be assessed through \code{diagnostic_checks} are only those inlcuded in the object \code{x} (see Arguments)
#' and are related to different model aspects in the setting of cost-effectiveness analyses. Specifically, in order to assess each parameter,
#' a different character name must be specified in the argument \code{parm}. The available names and the parameters associated with them are:
#' \itemize{
#' \item "mu.e" the mean parameters of the effect variables in the two treatment arms.
#' \item "mu.c" the mean parameters of the cost variables in the two treatment arms.
#' \item "sd.e" the standard deviation parameters of the effect variables in the two treatment arms.
#' \item "sd.c" the standard deviation parameters of the cost variables in the two treatment arms.
#' \item "beta0.e" the marginal mean parameters for the effect variables in the two treatment arms (only when covariates were included in the model).
#' \item "beta.e" the covariate coefficient parameters for the effect variables in the two treatment arms (only when covariates were included in the model).
#' \item "beta0.c" the marginal mean parameters for the cost variables in the two treatment arms (only when covariates were included in the model).
#' \item "beta.c" the covariate coefficient parameters for the cost variables in the two treatment arms (only when covariates were included in the model).
#' \item "gamma0.e" the baseline parameters of the missingness mechanism for the effect variables in the two treatment arms.
#' \item "gamma0.c" the baseline parameters of the missingness mechanism for the cost variables in the two treatment arms.
#' \item "gamma.e" the covariate coefficient parameters of the missingness mechanism for the effect variables in the two treatment arms.
#' \item "gamma.c" the covariate coefficient parameters of the missingness mechanism for the cost variables in the two treatment arms.
#' \item "delta.e" the mnar parameters of the missingness mechanism for the effect variables in the two treatment arms.
#' \item "delta.c" the mnar parameters of the missingness mechanism for the cost variables in the two treatment arms.
#' \item "all" all available parameters stored in the object \code{x}.
#' }
#' To notice that the marginal mean and the covariate coefficient parameters can be accessed only if covariate data were included in the model
#' and stored in the object \code{x} of class \strong{MissingHE}. Specifically, the mean parameters "beta0.e" and/or "beta0.c" are mere substitute of the 
#' effect and cost mean parameters in the case covariates are included in the model (must be used as substitutes of "mu.e" and/or "mu.c").
#' Finally, the parameters associated with the missingness model can be accessed only if an \emph{informative} or MNAR missing data mechanism was assumed
#' when runinng the model using \code{\link{run_model}}.
#' 
#' The argument \code{theme} allows to customise the graphical aspect of the plots generated by \code{diagnostic_checks} and
#' allows to choose among a set of possible pre-defined themes taken form the package \strong{ggtheme}. Those available can
#' be indicated using the following character names: "base","calc","economist","excel","few","538","gdocs","hc","par","pander","solarized","stata","tufte","wsj".
#' 
#' @author Andrea Gabrio
#' @references 
#' Gelman, A. Carlin, JB., Stern, HS. Rubin, DB.(2003). \emph{Bayesian Data Analysis, 2nd edition}, CRC Press.
#'
#'Brooks, S. Gelman, A. Jones, JL. Meng, XL. (2011). \emph{Handbook of Markov Chain Monte Carlo}, CRC/Chapman and Hall.
#' @import ggplot2 coda
#' @importFrom stats quantile
#' @export 
#' @examples 
#' #For examples see the function run_model
#' #
#' #

diagnostic_checks<-function(x,class="histogram",parm="all",theme=NULL,...){
  #diagnostic posterior checks to assess convergence of chains
  #define additional inputs as a list
  exArgs <- list(...)
  #x can only be object of class missing
  if(class(x)!="missingHE"){
    stop("Only objects of class 'missing' can be used")
  }
  #load namespace of ggmcmc and mcmcplots else ask to install the packages
  if(!isTRUE(requireNamespace("ggmcmc"))|!isTRUE(requireNamespace("mcmcplots"))) {
    stop("You need to install the R packages 'ggmcmc' and 'mcmcplots'. Please run in your R terminal:\n install.packages('ggmcmc','mcmcplots')")
  }
  #if theme specified load namespace of ggmcmc and mcmcplots else ask to install the packages
  if(length(theme)!=0){
    theme_names=c("base","calc","economist","excel","few","538","gdocs","hc","par","pander","solarized","stata","tufte","wsj")
    if(theme %in% theme_names){
      if(!isTRUE(requireNamespace("ggmcmc"))){
        stop("You need to install the R packages 'ggmcmc' and 'mcmcplots'. Please run in your R terminal:\n install.packages('ggmcmc','mcmcplots')")
      }
    } else if(!theme %in% theme_names){
      stop("You must provide one of the available theme styles")
    }}
  #need to specify parameters available in the model
  if(all(parm %in% c("all","mu.e","mu.c","sd.e","sd.c","corr",
                     "beta0.e","beta0.c","beta.e","beta.c",
                     "gammao.e","gamma0.c","gamma.e","gamma.c","delta.e","delta.c"))==FALSE ){
    stop("You must provide valid parameter names contained in the output of run_model")
  }
  if(length(parm)!=1){
    stop("You can only visualise diagnostic checks for one family of parameters at a time or for all parameters together \n setting the default value parm='all'")
  }
  #need to specify type of checks among those available
  if(!class %in% c("summary","histogram","running","denplot","compare","traceplot","autocorrplot","crosscorplot","Rhat","geweke","caterpillar","pairs")) {
    stop("Classes available for use are 'summary','histogram','running','denplot','compare','traceplot','autocorrplot','crosscorplot','Rhat','geweke','caterpillar','pairs'")
  }
  labs <- parm
  labs[pmatch("corr",labs)] <- "theta"
  labs[pmatch("mu.e",labs)] <- "mu_e"
  labs[pmatch("mu.c",labs)] <- "mu_c"
  labs[pmatch("sd.e",labs)] <- "s_e"
  labs[pmatch("sd.c",labs)] <- "s_c"
  labs[pmatch("gamma0.e",labs)] <- "gamma0_e"
  labs[pmatch("gamma0.c",labs)] <- "gamma0_c"
  labs[pmatch("beta0.e",labs)] <- "beta0_e"
  labs[pmatch("beta0.c",labs)] <- "beta0_c"
  labs[pmatch("beta.e",labs)] <- "beta_e"
  labs[pmatch("beta.c",labs)] <- "beta_c"
  labs[pmatch("gamma.e",labs)] <- "gamma_e"
  labs[pmatch("gamma.c",labs)] <- "gamma_c"
  labs[pmatch("delta.e",labs)] <- "delta_e"
  labs[pmatch("delta.c",labs)] <- "delta_c"
  #warning for forward sampling 
  if(x$model_class=="forward"){
    stop("you cannot run diagnostic checks for forward sampling output")
  }
  #create mcmc object from model results
  mcmc_object<-coda::as.mcmc(x$model_output$`model summary`)
  #set type of parameters to call in the plots (either all or a family of parameters)
  if(parm=="all"& class!="summary"){
    #exclude missing values from parameter plotting
    #get names of variables
    v_name<-coda::varnames(mcmc_object[,,drop=FALSE])
    check_name_eff<-grepl("eff",v_name)
    check_index_eff<-which(check_name_eff,TRUE)
    check_name_cost<-grepl("cost",v_name)
    check_index_cost<-which(check_name_cost,TRUE)
    check_index<-c(check_index_cost,check_index_eff)
    #exclude all missing cost and effect data
    param<-v_name[-check_index]
    #assign a common pre-name to the parameters for plotting using family argument
    param_all<-paste("model.",param,sep="")
    v_name[-check_index]<-param_all
    coda::varnames(mcmc_object)<-v_name
    #define data frame for change name displayed in the plots
    #P<-data.frame(
    #Parameter<-v_name,
    #Label<-c(v_name[check_index],param))
    #colnames(P)<-c("Parameter","Label")
    ggmcmc_object<-ggmcmc::ggs(mcmc_object)
    family=c("model")
  } else {
  family=labs
  ggmcmc_object<-ggmcmc::ggs(mcmc_object)
  }
  #check whether additional inputs are specified or use default values
  if(class=="summary"){
    #summary results html page 
    if(parm=="all"){
      v_name<-coda::varnames(mcmc_object[,,drop=FALSE])
      check_name_eff<-grepl("eff",v_name)
      check_index_eff<-which(check_name_eff,TRUE)
      check_name_cost<-grepl("cost",v_name)
      check_index_cost<-which(check_name_cost,TRUE)
      check_index<-c(check_index_cost,check_index_eff)
      param<-v_name[-check_index]
      check_index2<-unique(substr(param,1,4))
      check_index3<-gsub("\\[|\\]", "", check_index2)
      param<-gsub("devi", "deviance", check_index3,fixed=TRUE)
      if(x$type=="MNAR_eff"|x$type=="MNAR_eff_cov"){
      param<-gsub("delt", "delta_e", param,fixed=TRUE)
      param<-gsub("gamm", "gamma0_e", param,fixed=TRUE)
      } else if(x$type=="MNAR_cost"|x$type=="MNAR_cost_cov"){
        param<-gsub("delt", "delta_c", param,fixed=TRUE)
        param<-gsub("gamm", "gamma0_c", param,fixed=TRUE)
      } else if(x$type=="MNAR"|x$type=="MNAR_cov"){
        param<-c(gsub("delt", "delta_e", param,fixed=TRUE),"delta_c")
        param<-c(gsub("gamm", "gamma0_e", param,fixed=TRUE),"gamma0_c")
      }
       } else {param=labs}
    if(exists("regex",where=exArgs)) {regex=exArgs$regex} else {regex=NULL}
    if(exists("leaf.marker",where=exArgs)) {leaf.marker=exArgs$leaf.marker} else {leaf.marker= "[\\[_]"}
    if(exists("random",where=exArgs)) {random=exArgs$random} else {random=NULL}
    if(exists("dir",where=exArgs)) {dir=exArgs$dir} else {dir=tempdir()}
    if(exists("filename",where=exArgs)) {filename=exArgs$filename} else {filename="MCMCoutput"}
    if(exists("extension",where=exArgs)) {extension=exArgs$extension} else {extension="html"}
    if(exists("title",where=exArgs)) {title=exArgs$title} else {title=NULL}
    if(exists("col",where=exArgs)) {col=exArgs$col} else {col=NULL}
    if(exists("lty",where=exArgs)) {lty=exArgs$lty} else {lty=1}
    if(exists("xlim",where=exArgs)) {xlim=exArgs$xlim} else {xlim=NULL}
    if(exists("ylim",where=exArgs)) {ylim=exArgs$ylim} else {ylim=NULL}
    if(exists("style",where=exArgs)) {style=exArgs$style} else {style=c("gray", "plain")}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=TRUE}
    ggmcmc_out<-mcmcplots::mcmcplot(mcmc_object,parms = param,regex = regex,leaf.marker = leaf.marker,random = random,dir = dir,filename = filename,
                                    extension = extension,title = title,col = col,lty = lty,xlim = xlim,ylim = ylim,style = style,greek = greek)
  } else if(class=="histogram"){
    #histogram plots
    if(exists("bins",where=exArgs)) {bins=exArgs$bins} else {bins=30}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_histogram(ggmcmc_object,family=family,bins = bins,greek = greek)
  } else if(class=="denplot"){
    #density plots
    if(exists("rug",where=exArgs)) {rug=exArgs$rug} else {rug=FALSE}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_density(ggmcmc_object,family=family,rug = rug,greek = greek)
  } else if(class=="running"){
    ##running mean plots
    if(exists("original_burnin",where=exArgs)) {original_burnin=exArgs$original_burnin} else {original_burnin=TRUE}
    if(exists("original_thin",where=exArgs)) {original_thin=exArgs$original_thin} else {original_thin=TRUE}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_running(ggmcmc_object,family=family,original_burnin = original_burnin,original_thin = original_thin,greek = greek)
  } else if(class=="compare"){
    #compare all chain with a fraction
    if(exists("partial",where=exArgs)) {partial=exArgs$partial} else {partial=0.1}
    if(exists("rug",where=exArgs)) {rug=exArgs$rug} else {rug=FALSE}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_compare_partial(ggmcmc_object,family=family,partial = partial,rug = rug,greek = greek)
  } else if(class=="traceplot"){
    #trace plots
    if(exists("original_burnin",where=exArgs)) {original_burnin=exArgs$original_burnin} else {original_burnin=TRUE}
    if(exists("original_thin",where=exArgs)) {original_thin=exArgs$original_thin} else {original_thin=TRUE}
    if(exists("simplify",where=exArgs)) {simplify=exArgs$simplify} else {simplify=NULL}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_traceplot(ggmcmc_object,family=family,original_burnin = original_burnin,original_thin = original_thin,simplify = simplify,greek = greek)
  } else if(class=="autocorrplot"){
    #autocorrelation plots
    if(exists("nLags",where=exArgs)) {nLags=exArgs$nLags} else {nLags=50}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_autocorrelation(ggmcmc_object,family=family,nLags = nLags,greek = greek)
  } else if(class=="crosscorplot"){
    #crosscorrelation plots
    if(exists("absolute_scale",where=exArgs)) {absolute_scale=exArgs$absolute_scale} else {absolute_scale=TRUE}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_crosscorrelation(ggmcmc_object,family=family,absolute_scale = absolute_scale,greek = greek)
  } else if(class=="Rhat"){
    #scale factor statistic plot
    if(exists("scaling",where=exArgs)) {scaling=exArgs$scaling} else {scaling=1.5}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_Rhat(ggmcmc_object,family=family, scaling = scaling,greek = greek)+ ggplot2::xlab("R_hat")
  } else if(class=="geweke"){
    #geweke stat plot
    if(exists("frac1",where=exArgs)) {frac1=exArgs$frac1} else {frac1=0.1}
    if(exists("frac2",where=exArgs)) {frac2=exArgs$frac2} else {frac2=0.5}
    if(exists("shadow_limit",where=exArgs)) {shadow_limit=exArgs$shadow_limit} else {shadow_limit=TRUE}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_geweke(ggmcmc_object,family=family, frac1 = frac1,frac2 = frac2,shadow_limit = shadow_limit,greek = greek)
  } else if(class=="caterpillar"){
    #caterpillar plots
    if(exists("X",where=exArgs)) {X=exArgs$X} else {X=NA}
    if(exists("thick_ci",where=exArgs)) {thick_ci=exArgs$thick_ci} else {thick_ci=c(0.05, 0.95)}
    if(exists("thin_ci",where=exArgs)) {thin_ci=exArgs$thin_ci} else {thin_ci=c(0.025, 0.975)}
    if(exists("line",where=exArgs)) {line=exArgs$line} else {line=NA}
    if(exists("horizontal",where=exArgs)) {horizontal=exArgs$horizontal} else {horizontal=TRUE}
    if(exists("model_labels",where=exArgs)) {model_labels=exArgs$model_labels} else {model_labels=NULL}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_caterpillar(ggmcmc_object,family=family, X=X,thick_ci = thick_ci,thin_ci = thin_ci,line = line,horizontal = horizontal,model_labels = model_labels,greek = greek)
  } else if(class=="pairs"){
    #plot parameter summary checks by pairs
    if(exists("title",where=exArgs)) {title=exArgs$title} else {title=NULL}
    if(exists("upper",where=exArgs)) {upper=exArgs$upper} else {upper=list(continuous = "cor", combo = "box_no_facet", discrete = "facetbar", na = "na")}
    if(exists("lower",where=exArgs)) {lower=exArgs$lower} else {lower=list(continuous = "points", combo = "facethist", discrete = "facetbar", na = "na")}
    if(exists("diag",where=exArgs)) {diag=exArgs$diag} else {diag= list(continuous = "densityDiag", discrete = "barDiag", na = "naDiag")}
    if(exists("xlab",where=exArgs)) {xlab=exArgs$xlab} else {xlab=NULL}
    if(exists("ylab",where=exArgs)) {ylab=exArgs$ylab} else {ylab=NULL}
    if(exists("axisLabels",where=exArgs)) {axisLabels=exArgs$axisLabels} else {axisLabels=c("show", "internal", "none")}
    if(exists("labeller",where=exArgs)) {labeller=exArgs$labeller} else {labeller="label_value"}
    if(exists("showStrips",where=exArgs)) {showStrips=exArgs$showStrips} else {showStrips=NULL}
    if(exists("legend",where=exArgs)) {legend=exArgs$legend} else {legend=NULL}
    if(exists("greek",where=exArgs)) {greek=exArgs$greek} else {greek=FALSE}
    ggmcmc_out<-ggmcmc::ggs_pairs(ggmcmc_object,family=family,greek = greek,title=title,upper=upper,lower=lower,diag=diag,xlab=xlab,ylab=ylab,
                                  axisLabels=axisLabels,labeller=labeller,showStrips=showStrips,legend=legend)
  }
  if(length(theme)!=0 & class!="summary"){
    #call theme for ggplots
    if(theme=="base") ggmcmc_out<-ggmcmc_out+ggthemes::theme_base()
    if(theme=="calc") ggmcmc_out<-ggmcmc_out+ggthemes::theme_calc()
    if(theme=="economist") ggmcmc_out<-ggmcmc_out+ggthemes::theme_economist()
    if(theme=="excel") ggmcmc_out<-ggmcmc_out+ggthemes::theme_excel()
    if(theme=="few") ggmcmc_out<-ggmcmc_out+ggthemes::theme_few()
    if(theme=="538") ggmcmc_out<-ggmcmc_out+ggthemes::theme_fivethirtyeight()
    if(theme=="gdocs") ggmcmc_out<-ggmcmc_out+ggthemes::theme_gdocs()
    if(theme=="hc") ggmcmc_out<-ggmcmc_out+ggthemes::theme_hc()
    if(theme=="par") ggmcmc_out<-ggmcmc_out+ggthemes::theme_par()
    if(theme=="solarized") ggmcmc_out<-ggmcmc_out+ggthemes::theme_solarized()
    if(theme=="pander") ggmcmc_out<-ggmcmc_out+ggthemes::theme_pander()
    if(theme=="stata") ggmcmc_out<-ggmcmc_out+ggthemes::theme_stata()
    if(theme=="tufte") ggmcmc_out<-ggmcmc_out+ggthemes::theme_tufte()
    if(theme=="wsj") ggmcmc_out<-ggmcmc_out+ggthemes::theme_wsj()
  }
  return(print(ggmcmc_out))
}




