#An internal function to select which type of model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables, choice of natural-scaled or scaled outcome
#variables, class of missingness mechanism assumed, independence or joint modelling  and choice
#between forward or standard sampling mode.

#' This function selects which type of model to execute.
#' @keywords BUGS JAGS models
#' @param type Type of missingness mechanism assumed. Choices are Missing Completely At Random (MCAR),
#'  Missing At Random (MAR), Missing Not At Random (MNAR). For 'MNAR' alternative versions are available 
#'  depending on whether the mechanism for only one variable is condiered, that is for the effects (MNAR_eff)
#'  or the costs (MNAR_cost), or also covariates are included either for the effects (MNAR_eff_cov),
#'  the costs (MNAR_cost), or both (MNAR_cov). 
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal.
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal.
#' @param program type of software used to run the model. Current alternatives are 'OpenBUGS' and 'JAGS'.
#' @param forward Logical. If \code{forward} is \code{TRUE}, the model is run in forward sampling mode
#'  without providing any data, else if \code{forward} is \code{FALSE} standard sampling mode is selected.
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_model<-function(type,dist_e,dist_c,program,forward)eval.parent(substitute({
  #if bivariate normal distribution is selected
  if(dist_e=="norm" & dist_c=="norm"){
      #if variables are standardised
      if(stand==TRUE){
        #if independence is assumed
        if(ind==TRUE){
          #if JAGS or BUGS program is selected
          if(program=="JAGS"|program=="BUGS"){
            #call function to write independent bivariate normal (standardised variables) model file/script in JAGS code for both forward/standard sampling 
            model_string<-normal_model_ind_scaled(type=type)
          } else {
            #stop if program name given is different from those available
            stop(paste("only BUGS or JAGS supported"))
          }
        } else if(ind==FALSE){
          #same as before but writing joint bivariate normal (standardised variables) model file 
          if(program=="JAGS"|program=="BUGS"){
            model_string<-normal_model_joint_scaled(type=type)
          } else {
            stop(paste("only BUGS or JAGS supported"))
          }
        }
      } else if(stand==FALSE){
        #same as before but writing independent bivariate normal (natural scale variables) model file 
        if(ind==TRUE){
          if(program=="JAGS"|program=="BUGS"){
            model_string<-normal_model_ind(type=type)
          } else {
            stop(paste("only BUGS or JAGS supported"))
          }
        } else if(ind==FALSE){
          #same as before but writing joint bivariate normal (natural scale variables) model file 
          if(program=="JAGS"|program=="BUGS"){
            model_string<-normal_model_joint(type=type)
          } else {
            stop(paste("only BUGS or JAGS supported"))
          }
        }
      } 
    #if normal-gamma model selected
  } else if(dist_e=="norm"& dist_c=="gamma"){
    #if JAGS or BUGS program is selected
    if(program=="JAGS"|program=="BUGS"){
      #call function to write independent normal-gamma model file/script in JAGS code for both forward/standard sampling 
      model_string<-normal_gamma_ind(type=type)
    } else {
      #stop if program name given is different from those available
      stop(paste("only BUGS or JAGS supported"))
    }
    #if beta-normal model selected
  }else if(dist_e=="beta"& dist_c=="norm"){
    #if JAGS or BUGS program is selected
    if(program=="JAGS"|program=="BUGS"){
      #call function to write independent normal-gamma model file/script in JAGS code for both forward/standard sampling 
      model_string<-beta_normal_ind(type=type)
    } else {
      #stop if program name given is different from those available
      stop(paste("only BUGS or JAGS supported"))
    }
    #if beta-gamma model selected
  } else if(dist_e=="beta" & dist_c=="gamma"){
    #if JAGS or BUGS program is selected
    if(program=="JAGS"|program=="BUGS"){
      #call function to write independent normal-gamma model file/script in JAGS code for both forward/standard sampling 
      model_string<-beta_gamma_ind(type=type)
    } else {
      #stop if program name given is different from those available
      stop(paste("only BUGS or JAGS supported"))
    }
  }
  #output model file written
  return(list(model_string=model_string))
}))