#' An internal function to change the hyperprior parameters in the selection model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of selection model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Selection models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available chocies are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weibull'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('nbinom') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available chocies are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param pe_fixed Number of fixed effects for the effectiveness model
#' @param pc_fixed Number of fixed effects for the cost model
#' @param pe_random Number of random effects for the effectiveness model
#' @param pc_random Number of random effects for the cost model
#' @param model_e_random Random effects formula for the effectiveness model
#' @param model_c_random Random effects formula for the costs model
#' @param d_list a list of the number and types of patterns in the data
#' @param restriction type of identifying restriction to be imposed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_pattern <- function(type, dist_e, dist_c, pe_fixed, pc_fixed, model_e_random, model_c_random, pe_random, pc_random, 
                          d_list, restriction) eval.parent( substitute( {
    if(pe_fixed == 1) {
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
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
           if(grepl("alpha_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
         } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
           if(grepl("alpha_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[2] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
           if(grepl("alpha_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
             prior_mue_str <- paste("alpha_p1[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
             model_string_jags <- gsub("alpha_p1[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
         } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
          if(grepl("alpha_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(grepl("alpha_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[2] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[2] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
          if(grepl("alpha_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue_str <- paste("alpha_p2[1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
        }
      }
    } else if(pe_fixed > 1){
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
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
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
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
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
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
             if(is.null(alpha0.prior) == FALSE & grepl("alpha_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_mue <- alpha0.prior
               prior_mue_str <- paste("alpha_p1[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
               model_string_jags <- gsub("alpha_p1[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
             if(is.null(alpha.prior) == FALSE & grepl("alpha_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
               prior_alphae <- alpha.prior
               prior_alphae_str <- paste("alpha_p1[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
               model_string_jags <- gsub("alpha_p1[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
           } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
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
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 3] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 3] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) } 
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 3] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 3] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
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
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
          if(is.null(alpha0.prior) == FALSE & grepl("alpha_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_mue <- alpha0.prior
            prior_mue_str <- paste("alpha_p2[1, 1] ~ dnorm(", prior_mue[1], ",", prior_mue[2])
            model_string_jags <- gsub("alpha_p2[1, 1] ~ dnorm(0, 0.0000001", prior_mue_str, model_string_jags,fixed = TRUE) }
          if(is.null(alpha.prior) == FALSE & grepl("alpha_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
            prior_alphae <- alpha.prior
            prior_alphae_str <- paste("alpha_p2[j, 1] ~ dnorm(", prior_alphae[1], ",", prior_alphae[2])
            model_string_jags <- gsub("alpha_p2[j, 1] ~ dnorm(0, 0.0000001", prior_alphae_str, model_string_jags, fixed = TRUE) }
        } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
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
  if(pc_fixed == 1) {
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("beta_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p1[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("beta_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(grepl("beta_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc_str <- paste("beta_p2[2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(pc_fixed > 1){
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p1[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p1[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p1[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p1[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p1[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p1[j, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 2] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 2] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 1] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 1] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
          prior_betac_str <- paste("beta_p2[j, 1] ~ dnorm(", prior_betac[1], ",", prior_betac[2])
          model_string_jags <- gsub("beta_p2[j, 1] ~ dnorm(0, 0.0000001", prior_betac_str, model_string_jags, fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(is.null(beta0.prior) == FALSE & grepl("beta_p2[1, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_muc <- beta0.prior
          prior_muc_str <- paste("beta_p2[1, 2] ~ dnorm(", prior_muc[1], ",", prior_muc[2])
          model_string_jags <- gsub("beta_p2[1, 2] ~ dnorm(0, 0.0000001", prior_muc_str, model_string_jags,fixed = TRUE) }
        if(is.null(beta.prior) == FALSE & grepl("beta_p2[j, 2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_betac <- beta.prior
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[2] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p1[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[3] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[2] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("ls_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("ls_e_p2[1] ~ dunif(-5, 10", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, sqrt(meane_p1[3] * (1 - meane_p1[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[2] ~ dunif(0, sqrt(meane_p1[2] * (1 - meane_p1[2]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, sqrt(meane_p1[1] * (1 - meane_p1[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, sqrt(meane_p2[3] * (1 - meane_p2[3]))", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[2] ~ dunif(0, sqrt(meane_p2[2] * (1 - meane_p2[2]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, sqrt(meane_p2[1] * (1 - meane_p2[1]))", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_e == "gamma" | dist_e == "logis") {
    if(is.null(sigma.prior.e) == FALSE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[2] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[2] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 10000", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_e == "exp" | dist_e == "bern" | dist_e == "pois") {
    if(is.null(sigma.prior.e) == FALSE) {
      stop("no prior for sigma required for the effects under the 'exp', 'bern', 'pois' distributions") }
  } else if(dist_e == "weibull") {
    if(is.null(sigma.prior.e) == FALSE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[2] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
        if(grepl("s_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[2] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(grepl("s_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("s_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("s_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
    }
  } else if(dist_e == "nbinom") {
    if(is.null(sigma.prior.e) == FALSE) {
      if(length(sigma.prior.e) != 2) {stop("provide correct hyper prior values") }
      prior_alphae <- sigma.prior.e 
      if(d_list$n_patterns[1] == 4) {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("tau_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_mis == FALSE) {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p1[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[2] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
        if(grepl("tau_e_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p1[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p1[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      }
      if(d_list$n_patterns[2] == 4) {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("tau_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_mis == FALSE) {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p2[3] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[3] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[3] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) } 
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
        if(grepl("tau_e_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[2] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[2] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
        if(grepl("tau_e_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphae_str <- paste("tau_e_p2[1] ~ dunif(", prior_alphae[1], ",", prior_alphae[2])
          model_string_jags <- gsub("tau_e_p2[1] ~ dunif(0, 100", prior_alphae_str, model_string_jags,fixed = TRUE) }
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(-5, 10", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p1[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("ls_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[1] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("ls_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("ls_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("ls_c_p2[2] ~ dunif(0, 100", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 3 & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_c_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_e_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_c_p1[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p1[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p1[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p1[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[1] == 2 & d_list$d1$d1_ec_obs == FALSE & d_list$d1$d1_ec_mis == FALSE & restriction == "AC") {
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
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_c_obs == FALSE & restriction == "CC") {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_ec_obs == FALSE & restriction == "AC") {
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 3 & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
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
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_e_obs == FALSE & restriction == "CC") {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_c_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_e_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "CC") {
        if(grepl("s_c_p2[1] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[1] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[1] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
        if(grepl("s_c_p2[2] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          prior_alphac_str <- paste("s_c_p2[2] ~ dunif(", prior_alphac[1], ",", prior_alphac[2])
          model_string_jags <- gsub("s_c_p2[2] ~ dunif(0, 10000", prior_alphac_str, model_string_jags,fixed = TRUE) }
      } else if(d_list$n_patterns[2] == 2 & d_list$d2$d2_ec_obs == FALSE & d_list$d2$d2_ec_mis == FALSE & restriction == "AC") {
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
    if(length(model_e_random) != 0 & pe_random == 1) {
      if(is.null(mu.a0.prior) == FALSE & grepl("mu_a_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.a0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_a0 <- mu.a0.prior
        prior_a0_str <- paste("mu_a_hat[t] ~ dnorm(", prior_a0[1], ",", prior_a0[2])
        model_string_jags <- gsub("mu_a_hat[t] ~ dnorm(0, 0.001", prior_a0_str, model_string_jags, fixed = TRUE) }
        if(is.null(s.a0.prior) == FALSE & grepl("s_a_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          if(length(s.a0.prior) != 2) {stop("provide correct hyper prior values") }
          prior_a0 <- s.a0.prior
          prior_a0_str <- paste("s_a_hat[t] ~ dunif(", prior_a0[1], ",", prior_a0[2])
          model_string_jags <- gsub("s_a_hat[t] ~ dunif(0, 100", prior_a0_str, model_string_jags, fixed = TRUE) }
    } else if(length(model_e_random) != 0 & pe_random > 1) {
      if(is.null(mu.a.prior) == FALSE & grepl("mu_a_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
       if(length(mu.a.prior) != 2) {stop("provide correct hyper prior values") }
       prior_a <- mu.a.prior
       prior_a_str <- paste("mu_a_hat[j, t] ~ dnorm(", prior_a[1], ",", prior_a[2])
       model_string_jags <- gsub("mu_a_hat[j, t] ~ dnorm(0, 0.001", prior_a_str, model_string_jags, fixed = TRUE) }
       if(is.null(s.a.prior) == FALSE & grepl("s_a_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
         if(length(s.a.prior) != 2) {stop("provide correct hyper prior values") }
         prior_a <- s.a.prior
         prior_a_str <- paste("s_a_hat[j, t] ~ dunif(", prior_a[1], ",", prior_a[2])
         model_string_jags <- gsub("s_a_hat[j, t] ~ dunif(0, 100", prior_a_str, model_string_jags, fixed = TRUE) }
    }
    if(length(model_c_random) != 0 & pc_random == 1) {
      if(is.null(mu.b0.prior) == FALSE & grepl("mu_b_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b0.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b0 <- mu.b0.prior
        prior_b0_str <- paste("mu_b_hat[t] ~ dnorm(", prior_b0[1], ",", prior_b0[2])
        model_string_jags <- gsub("mu_b_hat[t] ~ dnorm(0, 0.001", prior_b0_str, model_string_jags, fixed = TRUE) }
        if(is.null(s.b0.prior) == FALSE & grepl("s_b_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          if(length(s.b0.prior) != 2) {stop("provide correct hyper prior values") }
          prior_b0 <- s.b0.prior
          prior_b0_str <- paste("s_b_hat[t] ~ dunif(", prior_b0[1], ",", prior_b0[2])
          model_string_jags <- gsub("s_b_hat[t] ~ dunif(0, 100", prior_b0_str, model_string_jags, fixed = TRUE) }
    } else if(length(model_c_random) != 0 & pc_random > 1) {
      if(is.null(mu.b.prior) == FALSE & grepl("mu_b_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b <- mu.b.prior
        prior_b_str <- paste("mu_b_hat[j, t] ~ dnorm(", prior_b[1], ",", prior_b[2])
        model_string_jags <- gsub("mu_b_hat[j, t] ~ dnorm(0, 0.001", prior_b_str, model_string_jags, fixed = TRUE) }
        if(is.null(s.b.prior) == FALSE & grepl("s_b_hat[j, t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
          if(length(s.b.prior) != 2) {stop("provide correct hyper prior values") }
          prior_b <- s.b.prior
          prior_b_str <- paste("s_b_hat[j, t] ~ dunif(", prior_b[1], ",", prior_b[2])
          model_string_jags <- gsub("s_b_hat[j, t] ~ dunif(0, 100", prior_b_str, model_string_jags, fixed = TRUE) }
    }   
    if(exists("mu.b_f.prior") == TRUE) {
      if(is.null(mu.b_f.prior) == FALSE & grepl("mu_b_f_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(mu.b_f.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_f <- mu.b_f.prior
        prior_b_f_str <- paste("mu_b_f_hat[t] ~ dnorm(", prior_b_f[1], ",", prior_b_f[2])
        model_string_jags <- gsub("mu_b_f_hat[t] ~ dnorm(0, 0.001", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    }      
    if(exists("s.b_f.prior") == TRUE) {
      if(is.null(s.b_f.prior) == FALSE & grepl("s_b_f_hat[t] ~ ", model_string_jags, fixed = TRUE) == TRUE) {
        if(length(s.b_f.prior) != 2) {stop("provide correct hyper prior values") }
        prior_b_f <- s.b_f.prior
        prior_b_f_str <- paste("s_b_f_hat[t] ~ dunif(", prior_b_f[1], ",", prior_b_f[2])
        model_string_jags <- gsub("s_b_f_hat[t] ~ dunif(0, 100", prior_b_f_str, model_string_jags, fixed = TRUE)
      }
    } 
    model_string_prior <- model_string_jags
    return(model_string_prior)
}))