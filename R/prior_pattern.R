#' An internal function to change the hyperprior parameters in the selection model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of selection model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm') or Beta ('beta').
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param pe Number of covariates for the effectiveness model
#' @param pc Number of covariates for the cost model
#' @param d_list a list of the number and types of patterns in the data
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_pattern <- function(type, dist_e, dist_c, pe, pc, d_list) eval.parent( substitute( {
    if(pe == 1) {
      if(is.null(alpha0.prior) == FALSE) {
        if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_mue <- alpha0.prior 
         if(d_list$n_patterns[1] == 4) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p1[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[2] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         }
        if(d_list$n_patterns[2] == 4) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[2] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        }
      }
    } else if(pe > 1){
      if(is.null(alpha0.prior) == FALSE | is.null(alpha.prior) == FALSE) {
        if(is.null(alpha0.prior) == FALSE) {
         if(length(alpha0.prior) != 2) {stop("provide correct hyper prior values") }
        }
        if(is.null(alpha.prior) == FALSE) {
          if(length(alpha.prior) != 2) {stop("provide correct hyper prior values") }
        }
           if(d_list$n_patterns[1] == 4) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
               prior_alphae_str <- paste("alpha_p1[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
               prior_alphae_str <- paste("alpha_p1[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
               prior_alphae_str <- paste("alpha_p1[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 2] <- ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 2] ~ dnorm(", prior_mue[1], ",", prior_mue[2], ")")
               model_string_jags <- gsub("alpha_p1[1, 2] <- alpha_p1[1, 1]", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
               prior_alphae_str <- paste("alpha_p1[j, 2] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 2] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           }
        if(d_list$n_patterns[2] == 4) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
            prior_alphae_str <- paste("alpha_p2[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
            prior_alphae_str <- paste("alpha_p2[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
            prior_alphae_str <- paste("alpha_p2[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 2] <- ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 2] ~ dnorm(", prior_mue[1], ",", prior_mue[2], ")")
            model_string_jags <- gsub("alpha_p2[1, 2] <- alpha_p2[1, 1]", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) 
            prior_alphae_str <- paste("alpha_p2[j, 2] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 2] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        }
      }
    }
  if(pc == 1) {
    if(is.null(beta0.prior) == FALSE) {
      if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
      prior_muc <- beta0.prior 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(pc > 1){
    if(is.null(beta0.prior) == FALSE | is.null(beta.prior) == FALSE) {
      if(is.null(beta0.prior) == FALSE) {
        if(length(beta0.prior) != 2) {stop("provide correct hyper prior values") }
      }
      if(is.null(beta.prior) == FALSE) {
        if(length(beta.prior) != 2) {stop("provide correct hyper prior values") }
      }
      if(d_list$n_patterns[1] == 4) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p1[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p1[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p1[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p1[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p2[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p2[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p2[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) 
          prior_betac_str <- paste("beta_p2[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      }
    }
  }
  if(dist_e == "norm") {
    if(is.null(sigma.prior.e) == FALSE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[2] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[2] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_e == "beta") {
    if(is.null(sigma.prior.e) == FALSE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[2] ~ dunif(0, sqrt(meane_p1[2] * (1 - meane_p1[2]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[2] ~ dunif(0, sqrt(meane_p2[2] * (1 - meane_p2[2]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
    }
  }
  if(dist_c == "norm") {
    if(is.null(sigma.prior.c) == FALSE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_c == "lnorm") {
    if(is.null(sigma.prior.c) == FALSE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_c == "gamma") {
    if(is.null(sigma.prior.c) == FALSE) {
      if(length(sigma.prior.c) != 2) {stop("provide correct hyper prior values") }
      prior_alphac <- sigma.prior.c 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 1000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      }
    }
  }
    if(exists("beta_f.prior") == TRUE) {
      if(is.null(beta_f.prior) == FALSE) {
        if(length(beta_f.prior) != 2) {stop("provide correct hyper prior values") }
        prior_beta_f <- beta_f.prior 
          if(grepl("beta_f_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_beta_f_str <- paste("beta_f_p1[1] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
            model_string_jags <- gsub("beta_f_p1[1] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags,fixed = TRUE) }
          if(grepl("beta_f_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_beta_f_str <- paste("beta_f_p2[1] ~ dnorm(", prior_beta_f[1], ",", prior_beta_f[2])
            model_string_jags <- gsub("beta_f_p2[1] ~ dnorm(0, 0.0000001", prior_beta_f_str, model_string_jags,fixed = TRUE) }
        }
      }
  model_string_prior <- model_string_jags
    return(model_string_prior)
}))