#' An internal function to change the hyperprior parameters in the pattern mixture model provided by the user depending on the type of
#' missingness mechanism and outcome distributions assumed
#' 
#' This function modifies default hyper prior parameter values in the type of pattern mixture model selected according 
#' to the type of missingness mechanism and distributions for the outcomes assumed.
#' @keywords priors distributions Pattern mixture models
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR), Missing Not At Random for the effects (MNAR_eff),
#' Missing Not At Random for the costs (MNAR_cost), and Missing Not At Random for both (MNAR). For a complete list of all available hyper parameters 
#' and types of models see the manual.
#' @param dist_e distribution assumed for the effects. Current available choices are: Normal ('norm'), Beta ('beta'), Gamma ('gamma'), Exponential ('exp'),
#' Weibull ('weib'), Logistic ('logis'), Poisson ('pois'), Negative Binomial ('negbin') or Bernoulli ('bern')
#' @param dist_c Distribution assumed for the costs. Current available choices are: Normal ('norm'), Gamma ('gamma') or LogNormal ('lnorm')
#' @param model_txt_info list containing model specification information used to write the txt file of the JAGS model.
#' @param model_string_jags text file of the model.
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

prior_pattern <- function(type, dist_e, dist_c, model_txt_info, model_string_jags) {
  prior.list <- model_txt_info$prior
  prior.dist <- sapply(prior.list, "[[", 1)
  if(!any(is.character(prior.dist))) { stop("Please provide prior distribution names as characters")}
  prior.dist <- tolower(prior.dist)
  prior.dist.valid <- c("norm", "bern", "beta", "bin", "chisqr", "exp", "f", "gamma", "gen.gamma", "dexp", "unif", "nt", "dirch",
                        "hyper", "logis", "lnorm", "negbin", "nchisqr", "par", "pois", "t", "weib", "cauchy", "cat", "half-norm", "half-cauchy", "default")  
  if(any(prior.dist %in% c("normal","gaussian"))) { 
    prior.dist[which(prior.dist %in% c("normal","gaussian"))] <- "norm"}
  if(any(prior.dist %in% c("exponential"))) { 
    prior.dist[which(prior.dist %in% c("exponential"))] <- "exp"}
  if(any(prior.dist %in% c("weib"))) { 
    prior.dist[which(prior.dist %in% c("weibull"))] <- "weib"}
  if(any(prior.dist %in% c("logistic"))) { 
    prior.dist[which(prior.dist %in% c("logistic"))] <- "logis"}
  if(any(prior.dist %in% c("bernoulli"))) { 
    prior.dist[which(prior.dist %in% c("bernoulli"))] <- "bern"}
  if(any(prior.dist %in% c("poisson"))) { 
    prior.dist[which(prior.dist %in% c("poisson"))] <- "pois"}
  if(any(prior.dist %in% c("negative binomial", "negbinom", "negbin", "negbinomial"))) { 
    prior.dist[which(prior.dist %in% c("negative binomial", "negbinom", "negbin", "negbinomial"))] <- "negbin"}
  if(any(prior.dist %in% c("lognormal", "lognorm", "lnormal"))) { 
    prior.dist[which(prior.dist %in% c("lognormal", "lognorm", "lnormal"))] <- "lnorm"}
  if(any(prior.dist %in% c("binom", "binomial"))) { 
    prior.dist[which(prior.dist %in% c("binom", "binomial"))] <- "bin"}
  if(any(prior.dist %in% c("binom", "binomial"))) { 
    prior.dist[which(prior.dist %in% c("chi2", "chisq", "chi-squared", "chisquared"))] <- "chisqr"}
  if(any(prior.dist %in% c("double exponential", "dexponential"))) { 
    prior.dist[which(prior.dist %in% c("double exponential", "dexponential"))] <- "dexp"}
  if(any(prior.dist %in% c("gengamma", "generalised gamma", "ggamma"))) { 
    prior.dist[which(prior.dist %in% c("gengamma", "generalised gamma", "ggamma"))] <- "gen.gamma"}
  if(any(prior.dist %in% c("hypergeometric", "hypergeom", "hgeometric", "hgeom"))) { 
    prior.dist[which(prior.dist %in% c("hypergeometric", "hypergeom", "hgeometric", "hgeom"))] <- "hyper"}
  if(any(prior.dist %in% c("hypergeometric", "hypergeom", "hgeometric", "hgeom"))) { 
    prior.dist[which(prior.dist %in% c("nchi2", "nchisq", "nchi-squared", "nchisquared"))] <- "nchisqr"}
  if(any(prior.dist %in% c("pareto"))) { 
    prior.dist[which(prior.dist %in% c("pareto"))] <- "par"}
  if(any(prior.dist %in% c("t-student", "t-stud"))) { 
    prior.dist[which(prior.dist %in% c("t-student", "t-stud"))] <- "t"}
  if(any(prior.dist %in% c("t-student", "t-stud"))) { 
    prior.dist[which(prior.dist %in% c("nt-student", "nt-stud"))] <- "nt"}
  if(any(prior.dist %in% c("uniform"))) { 
    prior.dist[which(prior.dist %in% c("uniform"))] <- "unif"}
  if(any(prior.dist %in% c("categorical", "categ"))) { 
    prior.dist[which(prior.dist %in% c("categorical", "categ"))] <- "cat"}
  if(any(prior.dist %in% c("cauchy"))) { 
    prior.dist[which(prior.dist %in% c("cauchy"))] <- "cauchy"}
  if(any(prior.dist %in% c("half-normal", "half-gaussian"))) { 
    prior.dist[which(prior.dist %in% c("half-normal", "half-gaussian"))] <- "half-norm"}
  if(any(prior.dist %in% c("half-cauchy"))) { 
    prior.dist[which(prior.dist %in% c("half-cauchy"))] <- "half-cauchy"}
  if(any(prior.dist %in% c("dirch"))) { 
    prior.dist[which(prior.dist %in% c("dirichlet", "dir"))] <- "dirch"}
  if(!any(prior.dist %in% prior.dist.valid)) { stop("Please provide valid names for prior distributions")}
  prior.dist.d <- paste("d", prior.dist, sep = "")
  names(prior.dist.d) <- names(prior.dist)
  stop_dist_val <- "Please provide correct number and values for prior parameters."  
  if(type %in% c("MNAR", "MNAR_eff", "MNAR_cost")) {
    if(!is.null(prior.list$delta.prior.e) & grepl("delta_e ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["delta.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$delta.prior.e[-1]) != 1) { stop(stop_dist_val)}
        delta.prior.e <- as.numeric(prior.list$delta.prior.e[-1])
        prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e, ")", sep = "")
      }
      if(prior.dist["delta.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$delta.prior.e[-1]) != 2) { stop(stop_dist_val)}
        delta.prior.e <- as.numeric(prior.list$delta.prior.e[-1])
        if(prior.dist.d["delta.prior.e"] == "dhalf-norm") { 
          prior.dist.d["delta.prior.e"] <- "dnorm"
          prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["delta.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["delta.prior.e"] <- "dt"
          prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["delta.prior.e"] == "dcauchy") { 
          prior.dist.d["delta.prior.e"] <- "dt"
          prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ")", sep = "")
      }
      if(prior.dist["delta.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$delta.prior.e[-1]) != 3) { stop(stop_dist_val)}
        delta.prior.e <- as.numeric(prior.list$delta.prior.e[-1])
        prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ", ", delta.prior.e[3], ")", sep = "")
      }
      if(prior.dist["delta.prior.e"] %in% c("hyper")) {
        if(length(prior.list$delta.prior.e[-1]) != 4) { stop(stop_dist_val)}
        delta.prior.e <- as.numeric(prior.list$delta.prior.e[-1])
        prior_deltae_str <- paste("delta_e ~ ", paste(prior.dist.d["delta.prior.e"]), "(", delta.prior.e[1], ", ", delta.prior.e[2], ", ", delta.prior.e[3], ", ", delta.prior.e[4], ")", sep = "")
      }
      model_string_jags <- gsub("delta_e ~ dunif(0, 1)", prior_deltae_str, model_string_jags, fixed = TRUE)
    }    
    if(!is.null(prior.list$delta.prior.c) & grepl("delta_c ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["delta.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$delta.prior.c[-1]) != 1) { stop(stop_dist_val)}
        delta.prior.c <- as.numeric(prior.list$delta.prior.c[-1])
        prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c, ")", sep = "")
      }
      if(prior.dist["delta.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$delta.prior.c[-1]) != 2) { stop(stop_dist_val)}
        delta.prior.c <- as.numeric(prior.list$delta.prior.c[-1])
        if(prior.dist.d["delta.prior.c"] == "dhalf-norm") { 
          prior.dist.d["delta.prior.c"] <- "dnorm"
          prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["delta.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["delta.prior.c"] <- "dt"
          prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["delta.prior.c"] == "dcauchy") { 
          prior.dist.d["delta.prior.c"] <- "dt"
          prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ")", sep = "")
      }
      if(prior.dist["delta.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$delta.prior.c[-1]) != 3) { stop(stop_dist_val)}
        delta.prior.c <- as.numeric(prior.list$delta.prior.c[-1])
        prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ", ", delta.prior.c[3], ")", sep = "")
      }
      if(prior.dist["delta.prior.c"] %in% c("hyper")) {
        if(length(prior.list$delta.prior.c[-1]) != 4) { stop(stop_dist_val)}
        delta.prior.c <- as.numeric(prior.list$delta.prior.c[-1])
        prior_deltae_str <- paste("delta_c ~ ", paste(prior.dist.d["delta.prior.c"]), "(", delta.prior.c[1], ", ", delta.prior.c[2], ", ", delta.prior.c[3], ", ", delta.prior.c[4], ")", sep = "")
      }
      model_string_jags <- gsub("delta_c ~ dunif(0, 1)", prior_deltae_str, model_string_jags, fixed = TRUE)
    }    
  }
  if(model_txt_info$n_patterns == 4) {
   if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
    if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
      alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
      prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
    }
    if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                        "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
      alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
      if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
        prior.dist.d["alpha.prior"] <- "dnorm"
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
        prior.dist.d["alpha.prior"] <- "dt"
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["alpha.prior"] == "dcauchy") { 
        prior.dist.d["alpha.prior"] <- "dt"
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
      prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
    }
    if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
      alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
      prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
    }
    if(prior.dist["alpha.prior"] %in% c("hyper")) {
      if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
      alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
      prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
    }
    model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
    model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", prior_deltae_str3, model_string_jags, fixed = TRUE)
   }  
  }
  if(model_txt_info$n_patterns == 3 & !2 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }    
  }
  if(model_txt_info$n_patterns == 3 & !1 %in% model_txt_info$d_or & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 3] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }    
  }
  if(model_txt_info$n_patterns == 3 & !3 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 3 & !4 %in% model_txt_info$d_or) {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")
          prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
        prior_deltae_str3 <- paste("alpha_p[j, 3] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p[j, 3] ~ dnorm(0, 0.0000001)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 2 & !all(c(2, 3) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
    }    
  }   
  if(model_txt_info$n_patterns == 2 & !all(c(2, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
        prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
        prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
        prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
        prior_deltae_str2 <- paste("alpha_p[j, 2] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("alpha_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltae_str2, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 2 & !all(c(3, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
    }    
  }     
  if(model_txt_info$n_patterns == 2 & !all(c(1, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$alpha.prior) & grepl("alpha_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["alpha.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$alpha.prior[-1]) != 1) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior, ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$alpha.prior[-1]) != 2) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        if(prior.dist.d["alpha.prior"] == "dhalf-norm") { 
          prior.dist.d["alpha.prior"] <- "dnorm"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dhalf-cauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["alpha.prior"] == "dcauchy") { 
          prior.dist.d["alpha.prior"] <- "dt"
          prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$alpha.prior[-1]) != 3) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ")", sep = "")
      }
      if(prior.dist["alpha.prior"] %in% c("hyper")) {
        if(length(prior.list$alpha.prior[-1]) != 4) { stop(stop_dist_val)}
        alpha.prior <- as.numeric(prior.list$alpha.prior[-1])
        prior_deltae_str1 <- paste("alpha_p[j, 1] ~ ", paste(prior.dist.d["alpha.prior"]), "(", alpha.prior[1], ", ", alpha.prior[2], ", ", alpha.prior[3], ", ", alpha.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("alpha_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltae_str1, model_string_jags, fixed = TRUE)
    }    
  } 
  if(length(model_txt_info$model_e_random) != 0) { 
    if(model_txt_info$pe_random == 1) { 
      mu_a_hat_text <- paste("mu_a_hat ~ ")
      s_a_hat_text <- paste("s_a_hat ~ ")} 
    if(model_txt_info$pe_random > 1) { 
      mu_a_hat_text <- paste("mu_a_hat[j] ~ ")
      s_a_hat_text <- paste("s_a_hat[j] ~ ")}
    if(!is.null(prior.list$mu.a.prior) & grepl(mu_a_hat_text, model_string_jags, fixed = TRUE)) {
      if(prior.dist["mu.a.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$mu.a.prior[-1]) != 1) { stop(stop_dist_val)}
        mu.a.prior <- as.numeric(prior.list$mu.a.prior[-1])
        prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior, ")", sep = "")
      }
      if(prior.dist["mu.a.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$mu.a.prior[-1]) != 2) { stop(stop_dist_val)}
        mu.a.prior <- as.numeric(prior.list$mu.a.prior[-1])
        if(prior.dist.d["mu.a.prior"] == "dhalf-norm") { 
          prior.dist.d["mu.a.prior"] <- "dnorm"
          prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["mu.a.prior"] == "dhalf-cauchy") { 
          prior.dist.d["mu.a.prior"] <- "dt"
          prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["mu.a.prior"] == "dcauchy") { 
          prior.dist.d["mu.a.prior"] <- "dt"
          prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ")", sep = "")
      }
      if(prior.dist["mu.a.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$mu.a.prior[-1]) != 3) { stop(stop_dist_val)}
        mu.a.prior <- as.numeric(prior.list$mu.a.prior[-1])
        prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ", ", mu.a.prior[3], ")", sep = "")
      }
      if(prior.dist["mu.a.prior"] %in% c("hyper")) {
        if(length(prior.list$mu.a.prior[-1]) != 4) { stop(stop_dist_val)}
        mu.a.prior <- as.numeric(prior.list$mu.a.prior[-1])
        prior_deltae_str <- paste(mu_a_hat_text, paste(prior.dist.d["mu.a.prior"]), "(", mu.a.prior[1], ", ", mu.a.prior[2], ", ", mu.a.prior[3], ", ", mu.a.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub(paste(mu_a_hat_text, "dnorm(0, 0.001)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
    }
    if(!is.null(prior.list$s.a.prior) & grepl(s_a_hat_text, model_string_jags, fixed = TRUE)) {
      if(prior.dist["s.a.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$s.a.prior[-1]) != 1) { stop(stop_dist_val)}
        s.a.prior <- as.numeric(prior.list$s.a.prior[-1])
        prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior, ")", sep = "")
      }
      if(prior.dist["s.a.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                        "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$s.a.prior[-1]) != 2) { stop(stop_dist_val)}
        s.a.prior <- as.numeric(prior.list$s.a.prior[-1])
        if(prior.dist.d["s.a.prior"] == "dhalf-norm") { 
          prior.dist.d["s.a.prior"] <- "dnorm"
          prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["s.a.prior"] == "dhalf-cauchy") { 
          prior.dist.d["s.a.prior"] <- "dt"
          prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["s.a.prior"] == "dcauchy") { 
          prior.dist.d["s.a.prior"] <- "dt"
          prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ")", sep = "")
      }
      if(prior.dist["s.a.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$s.a.prior[-1]) != 3) { stop(stop_dist_val)}
        s.a.prior <- as.numeric(prior.list$s.a.prior[-1])
        prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ", ", s.a.prior[3], ")", sep = "")
      }
      if(prior.dist["s.a.prior"] %in% c("hyper")) {
        if(length(prior.list$s.a.prior[-1]) != 4) { stop(stop_dist_val)}
        s.a.prior <- as.numeric(prior.list$s.a.prior[-1])
        prior_deltae_str <- paste(s_a_hat_text, paste(prior.dist.d["s.a.prior"]), "(", s.a.prior[1], ", ", s.a.prior[2], ", ", s.a.prior[3], ", ", s.a.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub(paste(s_a_hat_text, "dunif(0, 100)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
    }
  }  
  if(model_txt_info$n_patterns == 4) {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str2, model_string_jags, fixed = TRUE)
    }  
  }
  if(model_txt_info$n_patterns == 3 & !2 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 3 & !1 %in% model_txt_info$d_or & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 2] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 3 & !3 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str2, model_string_jags, fixed = TRUE)
    } 
  } 
  if(model_txt_info$n_patterns == 3 & !4 %in% model_txt_info$d_or) {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str2, model_string_jags, fixed = TRUE)
    } 
  }        
  if(model_txt_info$n_patterns == 2 & !all(c(2, 3) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
    }    
  }   
  if(model_txt_info$n_patterns == 2 & !all(c(2, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
    }    
  }  
  if(model_txt_info$n_patterns == 2 & !all(c(3, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str1 <- paste("beta_p[j, 1] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 1] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str2, model_string_jags, fixed = TRUE)
    }    
  }     
  if(model_txt_info$n_patterns == 2 & !all(c(1, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$beta.prior) & grepl("beta_p[j, 2] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["beta.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$beta.prior[-1]) != 1) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior, ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$beta.prior[-1]) != 2) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        if(prior.dist.d["beta.prior"] == "dhalf-norm") { 
          prior.dist.d["beta.prior"] <- "dnorm"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dhalf-cauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["beta.prior"] == "dcauchy") { 
          prior.dist.d["beta.prior"] <- "dt"
          prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", 1", ")", sep = "")}
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$beta.prior[-1]) != 3) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ")", sep = "")
      }
      if(prior.dist["beta.prior"] %in% c("hyper")) {
        if(length(prior.list$beta.prior[-1]) != 4) { stop(stop_dist_val)}
        beta.prior <- as.numeric(prior.list$beta.prior[-1])
        prior_deltac_str2 <- paste("beta_p[j, 2] ~ ", paste(prior.dist.d["beta.prior"]), "(", beta.prior[1], ", ", beta.prior[2], ", ", beta.prior[3], ", ", beta.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub("beta_p[j, 2] ~ dnorm(0, 0.0000001)", prior_deltac_str1, model_string_jags, fixed = TRUE)
    }    
  }
  if(length(model_txt_info$model_c_random) != 0) { 
    if(model_txt_info$pc_random == 1) { 
      mu_b_hat_text <- paste("mu_b_hat ~ ")
      s_b_hat_text <- paste("s_b_hat ~ ")} 
    if(model_txt_info$pc_random > 1) { 
      mu_b_hat_text <- paste("mu_b_hat[j] ~ ")
      s_b_hat_text <- paste("s_b_hat[j] ~ ")}
    if(!is.null(prior.list$mu.b.prior) & grepl(mu_b_hat_text, model_string_jags, fixed = TRUE)) {
      if(prior.dist["mu.b.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$mu.b.prior[-1]) != 1) { stop(stop_dist_val)}
        mu.b.prior <- as.numeric(prior.list$mu.b.prior[-1])
        prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior, ")", sep = "")
      }
      if(prior.dist["mu.b.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$mu.b.prior[-1]) != 2) { stop(stop_dist_val)}
        mu.b.prior <- as.numeric(prior.list$mu.b.prior[-1])
        if(prior.dist.d["mu.b.prior"] == "dhalf-norm") { 
          prior.dist.d["mu.b.prior"] <- "dnorm"
          prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["mu.b.prior"] == "dhalf-cauchy") { 
          prior.dist.d["mu.b.prior"] <- "dt"
          prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["mu.b.prior"] == "dcauchy") { 
          prior.dist.d["mu.b.prior"] <- "dt"
          prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ")", sep = "")
      }
      if(prior.dist["mu.b.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$mu.b.prior[-1]) != 3) { stop(stop_dist_val)}
        mu.b.prior <- as.numeric(prior.list$mu.b.prior[-1])
        prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ", ", mu.b.prior[3], ")", sep = "")
      }
      if(prior.dist["mu.b.prior"] %in% c("hyper")) {
        if(length(prior.list$mu.b.prior[-1]) != 4) { stop(stop_dist_val)}
        mu.b.prior <- as.numeric(prior.list$mu.b.prior[-1])
        prior_deltae_str <- paste(mu_b_hat_text, paste(prior.dist.d["mu.b.prior"]), "(", mu.b.prior[1], ", ", mu.b.prior[2], ", ", mu.b.prior[3], ", ", mu.b.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub(paste(mu_b_hat_text, "dnorm(0, 0.001)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
    }
    if(!is.null(prior.list$s.b.prior) & grepl(s_b_hat_text, model_string_jags, fixed = TRUE)) {
      if(prior.dist["s.b.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$s.b.prior[-1]) != 1) { stop(stop_dist_val)}
        s.b.prior <- as.numeric(prior.list$s.b.prior[-1])
        prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior, ")", sep = "")
      }
      if(prior.dist["s.b.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                        "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$s.b.prior[-1]) != 2) { stop(stop_dist_val)}
        s.b.prior <- as.numeric(prior.list$s.b.prior[-1])
        if(prior.dist.d["s.b.prior"] == "dhalf-norm") { 
          prior.dist.d["s.b.prior"] <- "dnorm"
          prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["s.b.prior"] == "dhalf-cauchy") { 
          prior.dist.d["s.b.prior"] <- "dt"
          prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["s.b.prior"] == "dcauchy") { 
          prior.dist.d["s.b.prior"] <- "dt"
          prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ", 1", ")", sep = "")}
        prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ")", sep = "")
      }
      if(prior.dist["s.b.prior"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$s.b.prior[-1]) != 3) { stop(stop_dist_val)}
        s.b.prior <- as.numeric(prior.list$s.b.prior[-1])
        prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ", ", s.b.prior[3], ")", sep = "")
      }
      if(prior.dist["s.b.prior"] %in% c("hyper")) {
        if(length(prior.list$s.b.prior[-1]) != 4) { stop(stop_dist_val)}
        s.b.prior <- as.numeric(prior.list$s.b.prior[-1])
        prior_deltae_str <- paste(s_b_hat_text, paste(prior.dist.d["s.b.prior"]), "(", s.b.prior[1], ", ", s.b.prior[2], ", ", s.b.prior[3], ", ", s.b.prior[4], ")", sep = "")
      }
      model_string_jags <- gsub(paste(s_b_hat_text, "dunif(0, 100)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
    }
  }
  if(model_txt_info$n_patterns == 4) {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE) |
       !is.null(prior.list$sigma.prior.e) & grepl("tau_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
     if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
        prior.dist.d["sigma.prior.e"] <- "dnorm"
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
        prior.dist.d["sigma.prior.e"] <- "dt"
        prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("s_e[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
        prior.dist.d["sigma.prior.e"] <- "dt"
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
      if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
    }
    if(dist_e == "norm") {
      model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e == "beta") {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e %in% c("gamma", "logis")) {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e %in% c("exp", "bern", "pois", "weib")) {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }    
    if(dist_e == "negbin") {
      prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
      prior_deltae_str <- gsub("s_e_p[3]", "tau_e_p[3]", prior_deltae_str3, fixed = TRUE)
      model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
   }
  }
  if(model_txt_info$n_patterns == 3 & !2 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
     if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                          "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
        prior.dist.d["sigma.prior.e"] <- "dnorm"
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
        prior.dist.d["sigma.prior.e"] <- "dt"
        prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")
        prior_deltae_str3 <- paste("s_e[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
        prior.dist.d["sigma.prior.e"] <- "dt"
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
    }
    if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
      if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
      sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
      prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
    }
    if(dist_e == "norm") {
      model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e == "beta") {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e %in% c("gamma", "logis")) {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
    if(dist_e %in% c("exp", "bern", "pois", "weib")) {
      model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }    
    if(dist_e == "negbin") {
      prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
      prior_deltae_str <- gsub("s_e_p[3]", "tau_e_p[3]", prior_deltae_str3, fixed = TRUE)
      model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
    }
   }    
  }
  if(model_txt_info$n_patterns == 3 & !1 %in% model_txt_info$d_or & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[3] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str3 <- paste("s_e[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[3]", "tau_e_p[3]", prior_deltae_str3, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
    }       
  }
  if(model_txt_info$n_patterns == 3 & !3 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 3 & !4 %in% model_txt_info$d_or) {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str3 <- paste("s_e[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")
          prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
        prior_deltae_str3 <- paste("s_e_p[3] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[3] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, sqrt(mu_e_p[3] * (1 - mu_e_p[3])))", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 10000)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str <- gsub("s_e_p[3]", "tau_e_p[3]", prior_deltae_str3, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[3] ~ dunif(0, 100)", prior_deltae_str3, model_string_jags, fixed = TRUE)
      }
    }      
  }
  if(model_txt_info$n_patterns == 2 & !all(c(2, 3) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(2, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
        prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_e[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
        prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
        prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
        prior_deltae_str2 <- paste("s_e_p[2] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[2] ~ dunif(0, sqrt(mu_e_p[2] * (1 - mu_e_p[2])))", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_e_p[2] ~ dunif(0, 100)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str <- gsub("s_e_p[2]", "tau_e_p[2]", prior_deltae_str2, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[2] ~ dunif(0, 100)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(3, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }       
  }
  if(model_txt_info$n_patterns == 2 & !all(c(1, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$sigma.prior.e) & grepl("s_e_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.e"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.e[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e, ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.e[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        if(prior.dist.d["sigma.prior.e"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.e"] <- "dnorm"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.e"] == "dcauchy") { 
          prior.dist.d["sigma.prior.e"] <- "dt"
          prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.e[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.e"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.e[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.e <- as.numeric(prior.list$sigma.prior.e[-1])
        prior_deltae_str1 <- paste("s_e_p[1] ~ ", paste(prior.dist.d["sigma.prior.e"]), "(", sigma.prior.e[1], ", ", sigma.prior.e[2], ", ", sigma.prior.e[3], ", ", sigma.prior.e[4], ")", sep = "")
      }
      if(dist_e == "norm") {
        model_string_jags <- gsub("s_e_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e == "beta") {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, sqrt(mu_e_p[1] * (1 - mu_e_p[1])))", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("gamma", "logis")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_e %in% c("exp", "bern", "pois", "weib")) {
        model_string_jags <- gsub("s_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }    
      if(dist_e == "negbin") {
        prior_deltae_str <- gsub("s_e_p[1]", "tau_e_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("tau_e_p[1] ~ dunif(0, 100)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }       
  }
  if(model_txt_info$n_patterns == 4) {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }    
  }
  if(model_txt_info$n_patterns == 3 & !2 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }    
  }
  if(model_txt_info$n_patterns == 3 & !1 %in% model_txt_info$d_or & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[2] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str2, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 3 & !3 %in% model_txt_info$d_or & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 3 & !4 %in% model_txt_info$d_or) {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(2, 3) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(2, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(3, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "CC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[1] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str1 <- paste("s_c_p[1] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[1] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str1 <- gsub("s_c_p[1]", "ls_c_p[1]", prior_deltae_str1, fixed = TRUE)
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str1, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[1] ~ dunif(0, 10)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[1] ~ dunif(0, 10000)", prior_deltae_str1, model_string_jags, fixed = TRUE)
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }     
  }
  if(model_txt_info$n_patterns == 2 & !all(c(1, 4) %in% model_txt_info$d_or) & model_txt_info$restriction == "AC") {
    if(!is.null(prior.list$sigma.prior.c) & grepl("s_c_p[2] ~ ", model_string_jags, fixed = TRUE)) {
      if(prior.dist["sigma.prior.c"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
        if(length(prior.list$sigma.prior.c[-1]) != 1) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c, ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                            "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
        if(length(prior.list$sigma.prior.c[-1]) != 2) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        if(prior.dist.d["sigma.prior.c"] == "dhalf-norm") { 
          prior.dist.d["sigma.prior.c"] <- "dnorm"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dhalf-cauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", "T(0, )", sep = "")}
        if(prior.dist.d["sigma.prior.c"] == "dcauchy") { 
          prior.dist.d["sigma.prior.c"] <- "dt"
          prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", 1", ")", sep = "")}
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("gen.gamma", "nt", "t")) {
        if(length(prior.list$sigma.prior.c[-1]) != 3) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ")", sep = "")
      }
      if(prior.dist["sigma.prior.c"] %in% c("hyper")) {
        if(length(prior.list$sigma.prior.c[-1]) != 4) { stop(stop_dist_val)}
        sigma.prior.c <- as.numeric(prior.list$sigma.prior.c[-1])
        prior_deltae_str2 <- paste("s_c_p[2] ~ ", paste(prior.dist.d["sigma.prior.c"]), "(", sigma.prior.c[1], ", ", sigma.prior.c[2], ", ", sigma.prior.c[3], ", ", sigma.prior.c[4], ")", sep = "")
      }
      if(dist_c == "norm") {
        model_string_jags <- gsub("s_c_p[2] ~ dt(0, pow(2.5, -2), 1)T(0,)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "lnorm") {
        prior_deltae_str2 <- gsub("s_c_p[2]", "ls_c_p[2]", prior_deltae_str2, fixed = TRUE)
        model_string_jags <- gsub("ls_c_p[2] ~ dunif(0, 10)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
      if(dist_c == "gamma") {
        model_string_jags <- gsub("s_c_p[2] ~ dunif(0, 10000)", prior_deltae_str2, model_string_jags, fixed = TRUE)
      }
    }       
  }
  if(!is.null(prior.list$beta_f.prior) & grepl("beta_f_p[1] ~ ", model_string_jags, fixed = TRUE)) {
    if(prior.dist["beta_f.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$beta_f.prior[-1]) != 1) { stop(stop_dist_val)}
      beta_f.prior <- as.numeric(prior.list$beta_f.prior[-1])
      prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior, ")", sep = "")
    }
    if(prior.dist["beta_f.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$beta_f.prior[-1]) != 2) { stop(stop_dist_val)}
      beta_f.prior <- as.numeric(prior.list$beta_f.prior[-1])
      if(prior.dist.d["beta_f.prior"] == "dhalf-norm") { 
        prior.dist.d["beta_f.prior"] <- "dnorm"
        prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["beta_f.prior"] == "dhalf-cauchy") { 
        prior.dist.d["beta_f.prior"] <- "dt"
        prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["beta_f.prior"] == "dcauchy") { 
        prior.dist.d["beta_f.prior"] <- "dt"
        prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ", 1", ")", sep = "")}
      prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ")", sep = "")
    }
    if(prior.dist["beta_f.prior"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$beta_f.prior[-1]) != 3) { stop(stop_dist_val)}
      beta_f.prior <- as.numeric(prior.list$beta_f.prior[-1])
      prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ", ", beta_f.prior[3], ")", sep = "")
    }
    if(prior.dist["beta_f.prior"] %in% c("hyper")) {
      if(length(prior.list$beta_f.prior[-1]) != 4) { stop(stop_dist_val)}
      beta_f.prior <- as.numeric(prior.list$beta_f.prior[-1])
      prior_deltae_str <- paste("beta_f_p[1] ~ ", paste(prior.dist.d["beta_f.prior"]), "(", beta_f.prior[1], ", ", beta_f.prior[2], ", ", beta_f.prior[3], ", ", beta_f.prior[4], ")", sep = "")
    }
    model_string_jags <- gsub("beta_f_p[1] ~ dnorm(0, 0.0000001)", prior_deltae_str, model_string_jags, fixed = TRUE)
  }  
  if(!is.null(prior.list$mu.b_f.prior) & grepl("mu_b_f_hat ~", model_string_jags, fixed = TRUE)) {
    mu_b_f_hat_text <- paste("mu_b_f_hat ~ ")
    if(prior.dist["mu.b_f.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$mu.b_f.prior[-1]) != 1) { stop(stop_dist_val)}
      mu.b_f.prior <- as.numeric(prior.list$mu.b_f.prior[-1])
      prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior, ")", sep = "")
    }
    if(prior.dist["mu.b_f.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                         "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$mu.b_f.prior[-1]) != 2) { stop(stop_dist_val)}
      mu.b_f.prior <- as.numeric(prior.list$mu.b_f.prior[-1])
      if(prior.dist.d["mu.b_f.prior"] == "dhalf-norm") { 
        prior.dist.d["mu.b_f.prior"] <- "dnorm"
        prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["mu.b_f.prior"] == "dhalf-cauchy") { 
        prior.dist.d["mu.b_f.prior"] <- "dt"
        prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["mu.b_f.prior"] == "dcauchy") { 
        prior.dist.d["mu.b_f.prior"] <- "dt"
        prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ", 1", ")", sep = "")}
      prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ")", sep = "")
    }
    if(prior.dist["mu.b_f.prior"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$mu.b_f.prior[-1]) != 3) { stop(stop_dist_val)}
      mu.b_f.prior <- as.numeric(prior.list$mu.b_f.prior[-1])
      prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ", ", mu.b_f.prior[3], ")", sep = "")
    }
    if(prior.dist["mu.b_f.prior"] %in% c("hyper")) {
      if(length(prior.list$mu.b_f.prior[-1]) != 4) { stop(stop_dist_val)}
      mu.b_f.prior <- as.numeric(prior.list$mu.b_f.prior[-1])
      prior_deltae_str <- paste(mu_b_f_hat_text, paste(prior.dist.d["mu.b_f.prior"]), "(", mu.b_f.prior[1], ", ", mu.b_f.prior[2], ", ", mu.b_f.prior[3], ", ", mu.b_f.prior[4], ")", sep = "")
    }
    model_string_jags <- gsub(paste(mu_b_f_hat_text, "dnorm(0, 0.001)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
  }  
  if(!is.null(prior.list$s.b_f.prior) & grepl("s_b_f_hat ~", model_string_jags, fixed = TRUE)) {
    s_b_f_hat_text <- paste("s_b_f_hat ~ ")
    if(prior.dist["s.b_f.prior"] %in% c("chisqr", "exp", "bern", "cat", "pois")) {
      if(length(prior.list$s.b_f.prior[-1]) != 1) { stop(stop_dist_val)}
      s.b_f.prior <- as.numeric(prior.list$s.b_f.prior[-1])
      prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior, ")", sep = "")
    }
    if(prior.dist["s.b_f.prior"] %in% c("norm", "beta", "bin", "f", "gamma", "dexp", "unif", "logis", "lnorm", 
                                        "negbin", "nchisqr", "par", "weib", "cauchy", "half-norm", "half-cauchy")) {
      if(length(prior.list$s.b_f.prior[-1]) != 2) { stop(stop_dist_val)}
      s.b_f.prior <- as.numeric(prior.list$s.b_f.prior[-1])
      if(prior.dist.d["s.b_f.prior"] == "dhalf-norm") { 
        prior.dist.d["s.b_f.prior"] <- "dnorm"
        prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ")", "T(0, )", sep = "")}
      if(prior.dist.d["s.b_f.prior"] == "dhalf-cauchy") { 
        prior.dist.d["s.b_f.prior"] <- "dt"
        prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ", 1", ")", "T(0, )", sep = "")}
      if(prior.dist.d["s.b_f.prior"] == "dcauchy") { 
        prior.dist.d["s.b_f.prior"] <- "dt"
        prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ", 1", ")", sep = "")}
      prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ")", sep = "")
    }
    if(prior.dist["s.b_f.prior"] %in% c("gen.gamma", "nt", "t")) {
      if(length(prior.list$s.b_f.prior[-1]) != 3) { stop(stop_dist_val)}
      s.b_f.prior <- as.numeric(prior.list$s.b_f.prior[-1])
      prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ", ", s.b_f.prior[3], ")", sep = "")
    }
    if(prior.dist["s.b_f.prior"] %in% c("hyper")) {
      if(length(prior.list$s.b_f.prior[-1]) != 4) { stop(stop_dist_val)}
      s.b_f.prior <- as.numeric(prior.list$s.b_f.prior[-1])
      prior_deltae_str <- paste(s_b_f_hat_text, paste(prior.dist.d["s.b_f.prior"]), "(", s.b_f.prior[1], ", ", s.b_f.prior[2], ", ", s.b_f.prior[3], ", ", s.b_f.prior[4], ")", sep = "")
    }
    model_string_jags <- gsub(paste(s_b_f_hat_text, "dunif(0, 100)", sep = ""), prior_deltae_str, model_string_jags, fixed = TRUE)
  }  
  if(!is.null(prior.list$patterns.prior) & grepl("p_prob ~ ", model_string_jags, fixed = TRUE)) {
    if(prior.dist["patterns.prior"] %in% c("ddirch")) {
      if(length(prior.list$patterns.prior[-1]) != model_txt_info$n_patterns) { stop(stop_dist_val)}
      patterns.prior <- as.numeric(prior.list$patterns.prior[-1])
      if(model_txt_info$n_patterns == 2) {
        prior_deltae_str <- paste("p_prob ~ ", paste(prior.dist.d["patterns.prior"]), "(", patterns.prior[1], ", ", patterns.prior[2], ")", sep = "")
      }
      if(model_txt_info$n_patterns == 3) {
        prior_deltae_str <- paste("p_prob ~ ", paste(prior.dist.d["patterns.prior"]), "(", patterns.prior[1], ", ", patterns.prior[2], ", ", patterns.prior[3], ")", sep = "")
      }
      if(model_txt_info$n_patterns == 4) {
        prior_deltae_str <- paste("p_prob ~ ", paste(prior.dist.d["patterns.prior"]), "(", patterns.prior[1], ", ", patterns.prior[2], ", ", patterns.prior[3], ", ", patterns.prior[4], ")", sep = "")
      }      
    }
    model_string_jags <- gsub("p_prob ~ ddirch(pp[])", prior_deltae_str, model_string_jags, fixed = TRUE)
  }
  model_string_prior <- model_string_jags
  return(model_string_prior)
}