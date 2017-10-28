#'An internal function to select which type of selection model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of missingness mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR)
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_selection<-function(type,dist_e,dist_c)eval.parent(substitute({
  #if bivariate normal distribution is selected
  if(dist_e=="norm" & dist_c=="norm"){
        if(ind==TRUE){
            model_string<-normal_selection_ind(type=type)
        }else if(ind==FALSE){
           model_string<-normal_selection_joint(type=type)
        }
  }else if(dist_e=="norm" & dist_c=="gamma"){
    if(ind==TRUE){
      model_string<-normal_gamma_selection_ind(type=type)
    }else if(ind==FALSE){
      model_string<-normal_gamma_selection_joint(type=type)
    }
  }else if(dist_e=="beta" & dist_c=="norm"){
    if(ind==TRUE){
      model_string<-beta_normal_selection_ind(type=type)
    }else if(ind==FALSE){
      model_string<-beta_normal_selection_joint(type=type)
    }
  }else if(dist_e=="beta" & dist_c=="gamma"){
    if(ind==TRUE){
      model_string<-beta_gamma_selection_ind(type=type)
    }else if(ind==FALSE){
      model_string<-beta_gamma_selection_joint(type=type)
    }
  }
  #output model file written
  return(list(model_string=model_string))
}))