#'An internal function to select which type of hurdle model to execute. 

#'Alternatives vary depending on the type of distribution assumed for the effect and cost variables,
#'type of strcutural value mechanism assumed and independence or joint modelling

#' This function selects which type of model to execute.
#' @keywords JAGS hurdle models
#' @param type Type of structural value mechanism assumed. Choices are Structural Completely At Random (SCAR),
#'  and Structural At Random (MNAR) 
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm') or Beta ('beta').
#' @param dist_c distribution assumed for the costs. Current available chocies are: Normal ('norm') or Gamma ('gamma').
#' @param se Structural value to be found in the effect data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @examples
#' #Internal function only
#' #No examples
#' #
#' #


write_hurdle<-function(type,dist_e,dist_c,se=se,sc=sc)eval.parent(substitute({
  #if bivariate normal distribution is selected
  if(dist_e=="norm" & dist_c=="norm"){
    #distinguish model if structurals in e, c or both and based on independence assumption
    if(is.null(se)==TRUE & is.null(sc)==FALSE){
        if(ind==TRUE){
            model_string<-normal_hurdle_c_ind(type=type)
        } else if(ind==FALSE){
            model_string<-normal_hurdle_c_joint(type=type)
        }
    }else if(is.null(se)==FALSE & is.null(sc)==TRUE){
      if(ind==TRUE){
        model_string<-normal_hurdle_e_ind(type=type)
      } else if(ind==FALSE){
        model_string<-normal_hurdle_e_joint(type=type)
      }
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-normal_hurdle_ec_ind(type=type)
      } else if(ind==FALSE){
        model_string<-normal_hurdle_ec_joint(type=type)
      }
    }
    #if normal-gamma model selected
  } else if(dist_e=="norm"& dist_c=="gamma"){
    if(is.null(se)==TRUE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-normal_gamma_hurdle_c_ind(type=type)
      } else if(ind==FALSE){
        model_string<-normal_gamma_hurdle_c_joint(type=type)
      }
    }else if(is.null(se)==FALSE & is.null(sc)==TRUE){
      if(ind==TRUE){
        model_string<-normal_gamma_hurdle_e_ind(type=type)
      } else if(ind==FALSE){
        model_string<-normal_gamma_hurdle_e_joint(type=type)
      }
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-normal_gamma_hurdle_ec_ind(type=type)
      } else if(ind==FALSE){
        model_string<-normal_gamma_hurdle_ec_joint(type=type)
      }
    }
    #if beta-normal model selected
  } else if(dist_e=="beta" & dist_c=="norm"){
    if(is.null(se)==TRUE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-beta_normal_hurdle_c_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_normal_hurdle_c_joint(type=type)
      }
    }else if(is.na(se)==FALSE & is.null(sc)==TRUE){
      if(ind==TRUE){
        model_string<-beta_normal_hurdle_e_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_normal_hurdle_e_joint(type=type)
      }
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-beta_normal_hurdle_ec_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_normal_hurdle_ec_joint(type=type)
      }
    }      
    #if beta-gamma model selected
  } else if(dist_e=="beta" & dist_c=="gamma"){
    if(is.null(se)==TRUE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-beta_gamma_hurdle_c_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_gamma_hurdle_c_joint(type=type)
      }
    }else if( is.null(se)==FALSE & is.null(sc)==TRUE){
      if(ind==TRUE){
        model_string<-beta_gamma_hurdle_e_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_gamma_hurdle_e_joint(type=type)
      }
    }else if(is.null(se)==FALSE & is.null(sc)==FALSE){
      if(ind==TRUE){
        model_string<-beta_gamma_hurdle_ec_ind(type=type)
      } else if(ind==FALSE){
        model_string<-beta_gamma_hurdle_ec_joint(type=type)
      }
    }      
  }
  #output model file written
  return(list(model_string=model_string))
}))