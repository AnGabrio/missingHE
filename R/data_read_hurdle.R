#' A function to read and re-arrange the data in different ways for the hurdle model
#'
#' This internal function imports the data and outputs only those variables that are needed to run the hurdle model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.se}, \code{model.sc} (model formulas for the structural effect and cost models) . Among these,
#' effectiveness, cost and treatment indicator variables must always be provided and named 'e', 'c' and 'trt' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side. 
#'  Random effects can also be specified for each model parameter.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side. 
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model. Random effects can also be specified for each model parameter.
#' @param model.se A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural effects are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the structural effects through a logistic-linear model. Random effects can also be specified for each model parameter.
#' @param model.sc A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural costs are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the structural costs through a logistic-linear model. Random effects can also be specified for each model parameter.
#' @param se Structural value to be found in the effect data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param cov_matrix Covariance matrix containing all model covariates. 
#' @param fixed_e Formula of the fixed effects specified in the effectiveness model.
#' @param fixed_c Formula of the fixed effects specified in the cost model.
#' @param random_e Formula of the random effects specified in the effectiveness model.
#' @param random_c Formula of the random effects specified in the cost model.
#' @param trt_lev Levels of the treatment indicator factor.
#' @param trt_pos Numeric position of the treatment indicator within the model formula.
#' @param type Type of structural value mechanism assumed, either 'SCAR' (Structural Completely At Random) or 'SAR' (Strcutural At Random).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @keywords read data hurdle models
#' @importFrom stats setNames na.omit sd as.formula model.matrix model.frame model.response terms
#' @export
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

data_read_hurdle <- function(data, model.eff, model.cost, 
                             model.se, model.sc, se, sc, cov_matrix, 
                             type, center, fixed_e, fixed_c, random_e, random_c,
                             trt_lev, trt_pos) {
  data2 <- data
  trt_index <- lapply(trt_lev, function(target) which(data$trt == target))
  names(trt_index) <- trt_lev
  data$e[is.na(data$e)] <- -999999
  data$c[is.na(data$c)] <- -999999
  mf_e_fixed <- model.frame(formula = fixed_e, data = data)
  mf_c_fixed <- model.frame(formula = fixed_c, data = data)
  trt_name_pos_e <- which(names(mf_e_fixed) == "trt")
  trt_name_pos_c <- which(names(mf_c_fixed) == "trt")
  terms <- e <- NULL
  x_e_fixed <- model.matrix(attr(mf_e_fixed, "terms"), data = mf_e_fixed)
  x_c_fixed <- model.matrix(attr(mf_c_fixed, "terms"), data = mf_c_fixed)
  t_m_names_pos_e <- which(colnames(x_e_fixed) %in% paste("trt", trt_lev, sep = ""))
  t_m_names_pos_c <- which(colnames(x_c_fixed) %in% paste("trt", trt_lev, sep = ""))
  fname_re_e_coeff <- as.formula(paste("e", "0", sep=" ~ "))
  fname_re_c_coeff <- as.formula(paste("c", "0", sep=" ~ "))
  fname_re_se_coeff <- as.formula(paste("se", "0", sep=" ~ "))
  fname_re_sc_coeff <- as.formula(paste("sc", "0", sep=" ~ "))  
  name_re_e_coeff <- name_re_c_coeff <- NULL
  if(!is.null(random_e)){
    name_re_e_coeff <- sub("\\|.*", "", random_e)
    if(grepl("0 + 1", name_re_e_coeff, fixed = TRUE)) { stop("Please either remove or add a random intercept")}
    name_clus_e <- sub('.*\\|', '', random_e)
    if(lengths(strsplit(name_clus_e, " ")) > 2) { stop("Please provide a single clustering variable for each model")}
    name_clus_e <- gsub(" ", "", name_clus_e, fixed = TRUE)
    if(!name_clus_e %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_e_coeff, "")[[1]][1] == 0) {
      no_random_int_e <- TRUE} else { no_random_int_e <- FALSE}
    if(no_random_int_e) { 
      name_re_e_coeff <- sub("[0]", "", name_re_e_coeff) 
      name_re_e_coeff <- sub("[+]", "", name_re_e_coeff)}
    if(name_re_e_coeff == "" | name_re_e_coeff == " ") { stop("Please state for which variables random effects are assumed")}
    fname_re_e_coeff <- as.formula(paste("e", name_re_e_coeff, sep = " ~ "))
    if(!all(names(model.frame(fname_re_e_coeff, data = data)) %in% c("0", "1", names(model.frame(fixed_e, data = data))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("e" %in% labels(terms(fname_re_e_coeff))) {
      stop("Please remove 'e' from the random effects of 'model.eff'")}
    if("c" %in% labels(terms(fname_re_e_coeff))) {
      stop("Dependence allowed only through the cost model; please remove 'c' from 'model.eff'")}
    mf_e_random <- model.frame(formula = fname_re_e_coeff, data = data)
    x_e_random <- model.matrix(attr(mf_e_random, "terms"), data = mf_e_random)
    if(no_random_int_e) {
      x_e_random <- as.matrix(x_e_random[, !colnames(x_e_random) == "(Intercept)"])
      if(is.null(colnames(x_e_random)) & dim(x_e_random)[2] == 1) {
        colnames(x_e_random) <- gsub(" ", "", name_re_e_coeff)}
    }
    clus_e <- data[, name_clus_e]
    if(!is.factor(clus_e)) { stop("Please define clustering variables as factors")}
    clus_e_lev <- levels(clus_e)
    clus_e_pos <- which(clus_e_lev %in% levels(clus_e))
    clusn_e <- as.numeric(clus_e)
    if(!all(diff(sort(unique(clusn_e))) == 1) | !min(clusn_e) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_e <- list()
    for(i in trt_lev) { clusnt_e[[i]] <- clusn_e[trt_index[[i]]]}
  }  
  if(!is.null(random_c)){
    name_re_c_coeff <- sub("\\|.*", "", random_c)
    if(grepl("0 + 1", name_re_c_coeff, fixed = TRUE)) { stop("Please either remove or add a random intercept")}
    name_clus_c <- sub('.*\\|', '', random_c)
    if(lengths(strsplit(name_clus_c, " ")) > 2) { stop("Please provide a single clustering variable for each model")}
    name_clus_c <- gsub(" ", "", name_clus_c, fixed = TRUE)
    if(!name_clus_c %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_c_coeff, "")[[1]][1] == 0) {
      no_random_int_c <- TRUE} else {no_random_int_c <- FALSE}
    if(no_random_int_c) { 
      name_re_c_coeff <- sub("[0]", "", name_re_c_coeff) 
      name_re_c_coeff <- sub("[+]", "", name_re_c_coeff)}
    if(name_re_c_coeff == "" | name_re_c_coeff == " ") { stop("Please state for which variables random effects are assumed")}
    if(gsub(" ", "", name_re_c_coeff) == "e" & !no_random_int_c) { name_re_c_coeff <- "1 + e"}
    fname_re_c_coeff <- as.formula(paste("c", name_re_c_coeff, sep = " ~ "))
    if(!all(names(model.frame(fname_re_c_coeff, data = data)) %in% c("0", "1", names(model.frame(fixed_c, data = data))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("c" %in% labels(terms(fname_re_c_coeff))) {
      stop("Please remove 'c' from the random effects of 'model.cost'")}
    if("e" %in% labels(terms(fname_re_c_coeff))) {
      if(length(grep(":e", labels(terms(fname_re_c_coeff)))) != 0 | length(grep("e:", labels(terms(fname_re_c_coeff)))) != 0) {
        stop("No interaction effects for 'e' are allowed")} 
    }
    mf_c_random <- model.frame(formula = fname_re_c_coeff, data = data)
    x_c_random <- model.matrix(attr(mf_c_random, "terms"), data = mf_c_random)
    if("e" %in% labels(terms(fname_re_c_coeff)) & length(labels(terms(fname_re_c_coeff))) == 1) {
      x_c_random <- subset(x_c_random, select = -c(e))}
    if(no_random_int_c) {
      x_c_random <- as.matrix(x_c_random[, !colnames(x_c_random) == "(Intercept)"])
      if(is.null(colnames(x_c_random)) & dim(x_c_random)[2] == 1) {
        colnames(x_c_random) <- gsub(" ", "", name_re_c_coeff)}
    }
    clus_c <- data[, name_clus_c]
    if(!is.factor(clus_c)) { stop("Please define clustering variables as factors")}
    clus_c_lev <- levels(clus_c)
    clus_c_pos <- which(clus_c_lev %in% levels(clus_c))
    clusn_c <- as.numeric(clus_c)
    if(!all(diff(sort(unique(clusn_c))) == 1) | !min(clusn_c) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_c <- list()
    for(i in trt_lev) { clusnt_c[[i]] <- clusn_c[trt_index[[i]]]}
  }  
  y_e <- model.response(mf_e_fixed)
  y_c <- model.response(mf_c_fixed)
  y_e[which(is.na(data2$e))] <- NA
  y_c[which(is.na(data2$c))] <- NA
  data$e[which(is.na(data2$e))] <- NA
  data$c[which(is.na(data2$c))] <- NA
  nt <- length(trt_lev)
  n <- setNames(table(data$trt), trt_lev)
  m_eff <- ifelse(is.na(data$e), 1, 0)
  m_cost <- ifelse(is.na(data$c), 1, 0)
  data$me <- m_eff
  data$mc <- m_cost
  if(!is.null(se)) { 
    s_eff <- ifelse(data$e == se, 1, 0)
    data$se <- s_eff
  } else {
    s_eff <- NULL
    }
  if(!is.null(sc)) { 
    s_cost <- ifelse(data$c == sc, 1, 0)
    data$sc <- s_cost
  } else {
    s_cost <- NULL
    }
  cov_e_fixed <- cov_c_fixed <- list()
  cov_e_center_fixed <- cov_c_center_fixed <- list()
  x_c_hold_fixed <- x_c_fixed
  if("e" %in% colnames(x_c_hold_fixed)) {
    x_c_fixed <- subset(x_c_hold_fixed, select = -c(e))
  }
  efft <- costt <- m_efft <- m_costt <- s_efft <- s_costt <- list()
  n_obs_eff <- n_obs_cost <- setNames(trt_pos, trt_lev) 
  if(!is.null(se)) { n_ns_eff <- setNames(trt_pos, trt_lev)} else { n_ns_eff <- NULL}
  if(!is.null(sc)) { n_ns_cost <- setNames(trt_pos, trt_lev)} else { n_ns_cost <- NULL}
  for(i in trt_lev) {
    efft[[i]] <- y_e[trt_index[[i]]]
    costt[[i]] <- y_c[trt_index[[i]]]
    cov_e_fixed[[i]] <- as.data.frame(x_e_fixed[trt_index[[i]], ])
    cov_c_fixed[[i]] <- as.data.frame(x_c_fixed[trt_index[[i]], ])
    cov_e_center_fixed[[i]] <- as.data.frame(scale(cov_e_fixed[[i]], scale = FALSE))
    cov_c_center_fixed[[i]] <- as.data.frame(scale(cov_c_fixed[[i]], scale = FALSE))
    cov_e_center_fixed[[i]][, 1] <- rep(1, n[i])
    cov_c_center_fixed[[i]][, 1] <- rep(1, n[i])
    m_efft[[i]] <- m_eff[trt_index[[i]]]
    m_costt[[i]] <- m_cost[trt_index[[i]]]
    n_obs_eff[i] <- length(na.omit(efft[[i]]))
    n_obs_cost[i] <- length(na.omit(costt[[i]]))
    if(!is.null(se)) {
      s_efft[[i]] <- s_eff[trt_index[[i]]]
      n_ns_eff[i] <- length(na.omit(efft[[i]])[na.omit(efft[[i]]) != se])
    }
    if(!is.null(sc)) {
      s_costt[[i]] <- s_cost[trt_index[[i]]]
      n_ns_cost[i] <- length(na.omit(costt[[i]])[na.omit(costt[[i]]) != sc])
    }
  }  
  n_mis_eff <- n - n_obs_eff
  n_mis_cost <- n - n_obs_cost
  if(!is.null(se)) { n_s_eff <- n_obs_eff - n_ns_eff} else { n_s_eff <- NULL}
  if(!is.null(sc)) { n_s_cost <- n_obs_cost - n_ns_cost} else { n_s_cost <- NULL}
  mean_cov_e_fixed <- lapply(cov_e_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_cov_c_fixed <- lapply(cov_c_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_cov_e_center_fixed <- lapply(cov_e_center_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_cov_c_center_fixed <- lapply(cov_c_center_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  if(center == TRUE) {
    cov_e_fixed <- cov_e_center_fixed
    cov_c_fixed <- cov_c_center_fixed
    mean_cov_e_fixed <- mean_cov_e_center_fixed
    mean_cov_c_fixed <- mean_cov_c_center_fixed
  }  
  if(!is.null(random_e)){
    n_clus_e <- setNames(table(clus_e), clus_e_lev)
    cov_e_random <- cov_e_center_random <- list()
    for(i in trt_lev) {
      cov_e_random[[i]] <- as.data.frame(x_e_random[trt_index[[i]], ])
      cov_e_center_random[[i]] <- as.data.frame(scale(cov_e_random[[i]], scale = FALSE))
      if(!no_random_int_e) {
        cov_e_center_random[[i]][, 1] <- rep(1, n[i])}
    }
    cov_e_random <- lapply(cov_e_random, setNames, colnames(x_e_random))
    cov_e_center_random <- lapply(cov_e_center_random, setNames, colnames(x_e_random))
    mean_cov_e_random <- lapply(cov_e_random, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    mean_cov_e_center_random <- lapply(cov_e_center_random, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    if(center) {
      cov_e_random <- cov_e_center_random
      mean_cov_e_random <- mean_cov_e_center_random
    }
  }
  if(!is.null(random_c)){
    x_c_hold_random <- x_c_random
    if("e" %in% colnames(x_c_hold_random)) {
      x_c_random <- subset(x_c_hold_random, select = -c(e))
    }
    n_clus_c <- setNames(table(clus_c), clus_c_lev)
    cov_c_random <- cov_c_center_random <- list()
    for(i in trt_lev) {
      cov_c_random[[i]] <- as.data.frame(x_c_random[trt_index[[i]], ])
      cov_c_center_random[[i]] <- as.data.frame(scale(cov_c_random[[i]], scale = FALSE))
      if(!no_random_int_c) {
        cov_c_center_random[[i]][, 1] <- rep(1, n[i])}
    }
    cov_c_random <- lapply(cov_c_random, setNames, colnames(x_c_random))
    cov_c_center_random <- lapply(cov_c_center_random, setNames, colnames(x_c_random))
    mean_cov_c_random <- lapply(cov_c_random, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    mean_cov_c_center_random <- lapply(cov_c_center_random, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    if(center) {
      cov_c_random <- cov_c_center_random
      mean_cov_c_random <- mean_cov_c_center_random
    }
  }  
  data2$e[is.na(data2$e)] <- -999999
  data2$c[is.na(data2$c)] <- -999999
  data2$me <- m_eff
  data2$mc <- m_cost
  if(!is.null(se)) { 
    data2$se <- s_eff
    data2$se[is.na(data2$se)] <- -999999
    }
  if(!is.null(sc)) { 
    data2$sc <- s_cost
    data2$sc[is.na(data2$sc)] <- -999999
    }
  if(!inherits(model.se, "formula") | !inherits(model.sc, "formula")) {
    stop("'model.se' and/or 'model.sc' must be formula objects")}
  fixed_se <- nobars_(model.se)
  fixed_sc <- nobars_(model.sc)
  if(!is.null(se)) {
    if(!all(names(model.frame(fixed_se, data = data2)) %in% c("se", "e", names(cov_matrix)))) {
      stop("Partially-observed covariates cannot be included in the model")}
    if(!all(names(model.frame(fixed_se, data = data2)) %in% names(data2))) {
      stop("Please provide names in the formula that correspond to those in the data")}
    if(names(model.frame(fixed_se, data = data2)[1]) != "se") {
      stop("Please set 'se' as the response in `model.se`")}
    if("c" %in% names(model.frame(fixed_se, data = data2))) {
      stop("Please remove 'c' from the right-hand side of 'model.se'")}
    if(!is.null(sc)){
      if("e" %in% names(model.frame(fixed_sc, data = data2))) {
        stop("Please remove 'e' from the right-hand side of 'model.sc'")} 
    }
    if("e" %in% labels(terms(fixed_se))) {
      if(length(grep(":e", labels(terms(fixed_se)))) != 0 | length(grep("e:", labels(terms(fixed_se)))) != 0) {
        stop("No interaction effects for 'e' are allowed")}
    }
  }
  if(!is.null(sc)) {
    if(!all(names(model.frame(fixed_sc, data = data2)) %in% c("sc", "c", names(cov_matrix)))) {
      stop("Partially-observed covariates cannot be included in the model")}
    if(!all(names(model.frame(fixed_sc, data = data2)) %in% names(data2))) {
      stop("Please provide names in the formula that correspond to those in the data")}
    if(names(model.frame(fixed_sc, data = data2)[1]) != "sc") {
      stop("Please set 'sc' as the response in `model.sc`")}
    if(!is.null(se)){
      if("c" %in% names(model.frame(fixed_se, data = data2))) {
        stop("Please remove 'c' from the right-hand side of 'model.se'")} 
    }
    if("c" %in% labels(terms(fixed_sc))) {
      if(length(grep(":c", labels(terms(fixed_sc)))) != 0 | length(grep("c:", labels(terms(fixed_sc)))) != 0) {
        stop("No interaction effects for 'c' are allowed")} 
    }
  }
  if(is.null(se)) { 
    fixed_se <- se ~ 1
    covz_e_fixed <- mean_covz_e_fixed <- NULL}
  if(is.null(sc)) { fixed_sc <- sc ~ 1
  covz_c_fixed <- mean_covz_c_fixed <- NULL}
  if(is.null(se) & is.null(sc)) { 
    stop("Structural values in at least one outcome variable are required, 
         please provide the structural value to be used in 'se' and/or in 'sc'")}
  random_se <- fb(model.se)
  random_sc <- fb(model.sc)
  if(is.null(se)) { random_se <- NULL}
  if(is.null(sc)) { random_sc <- NULL}
  if(!is.null(random_se) & length(random_se) > 1 | !is.null(random_sc) & length(random_sc) > 1) {
    stop("Random effects can be included in the formula only through a single expression within brackets")}
  name_re_se_coeff <- NULL
  name_re_sc_coeff <- NULL  
  if(!is.null(random_se)){
    name_re_se_coeff <- sub("\\|.*", "", random_se)
    if(grepl("0 + 1", name_re_se_coeff, fixed = TRUE)) { stop("Either remove or add a random intercept")}
    name_clus_se <- sub('.*\\|', '', random_se)
    if(lengths(strsplit(name_clus_se, " ")) > 2) {stop("Please provide a single clustering variable for each model")}
    name_clus_se <- gsub(" ", "", name_clus_se, fixed = TRUE)
    if(!name_clus_se %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_se_coeff, "")[[1]][1] == 0) {
      no_random_int_se <- TRUE} else {no_random_int_se <- FALSE}
    if(no_random_int_se) { 
      name_re_se_coeff <- sub("[0]", "", name_re_se_coeff) 
      name_re_se_coeff <- sub("[+]", "", name_re_se_coeff) 
    }
    if(name_re_se_coeff == "" | name_re_se_coeff == " ") { stop("Please state for which variables the random effects are assumed")}
    if(gsub(" ", "", name_re_se_coeff) == "e" & !no_random_int_se) {name_re_se_coeff <- "1 + e"}
    fname_re_se_coeff <- as.formula(paste("se", name_re_se_coeff, sep=" ~ "))
    if(!all(names(model.frame(fname_re_se_coeff, data = data2)) %in% c("0", "1", names(model.frame(fixed_se, data = data2))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("se" %in% labels(terms(fname_re_se_coeff))) {
      stop("Please remove 'se' from the random effects of 'model.se'")}
    if("e" %in% labels(terms(fname_re_se_coeff))) {
      if(length(grep(":e", labels(terms(fname_re_se_coeff)))) != 0 | length(grep("e:", labels(terms(fname_re_se_coeff)))) != 0) {
        stop("No interaction effects for 'e' are allowed")} 
    }
    clus_se <- data[, name_clus_se]
    clus_se_lev <- levels(clus_se)
    clus_se_pos <- which(clus_se_lev %in% levels(clus_se))
    if(!is.factor(clus_se)) { stop("Please define clustering variables as factors")}
    clusn_se <- as.numeric(clus_se)
    if(!all(diff(sort(unique(clusn_se))) == 1) | !min(clusn_se) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_se <- list()
    for(i in trt_lev) { clusnt_se[[i]] <- clusn_se[trt_index[[i]]]}
  }
  if(!is.null(random_sc)){
    name_re_sc_coeff <- sub("\\|.*", "", random_sc)
    if(grepl("0 + 1", name_re_sc_coeff, fixed = TRUE)) { stop("Either remove or add a random intercept")}
    name_clus_sc <- sub('.*\\|', '', random_sc)
    if(lengths(strsplit(name_clus_sc, " ")) > 2) {stop("Please provide a single clustering variable for each model")}
    name_clus_sc <- gsub(" ", "", name_clus_sc, fixed = TRUE)
    if(!name_clus_sc %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_sc_coeff, "")[[1]][1] == 0) {
      no_random_int_sc <- TRUE} else {no_random_int_sc <- FALSE}
    if(no_random_int_sc) { 
      name_re_sc_coeff <- sub("[0]", "", name_re_sc_coeff) 
      name_re_sc_coeff <- sub("[+]", "", name_re_sc_coeff) 
    }
    if(name_re_sc_coeff == "" | name_re_sc_coeff == " ") { stop("Please state for which variables the random effects are assumed")}
    if(gsub(" ", "", name_re_sc_coeff) == "c" & !no_random_int_sc) {name_re_sc_coeff <- "1 + c"}
    fname_re_sc_coeff <- as.formula(paste("sc", name_re_sc_coeff, sep=" ~ "))
    if(!all(names(model.frame(fname_re_sc_coeff, data = data2)) %in% c("0", "1", names(model.frame(fixed_sc, data = data2))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("sc" %in% labels(terms(fname_re_sc_coeff))) {
      stop("Please remove 'sc' from the random effects of 'model.sc'")}
    if("c" %in% labels(terms(fname_re_sc_coeff))) {
      if(length(grep(":c", labels(terms(fname_re_sc_coeff)))) != 0 | length(grep("c:", labels(terms(fname_re_sc_coeff)))) != 0) {
        stop("No interaction effects for 'c' are allowed")} 
    }
    clus_sc <- data[, name_clus_sc]
    clus_sc_lev <- levels(clus_sc)
    clus_sc_pos <- which(clus_sc_lev %in% levels(clus_sc))
    if(!is.factor(clus_sc)) { stop("Please define clustering variables as factors")}
    clusn_sc <- as.numeric(clus_sc)
    if(!all(diff(sort(unique(clusn_sc))) == 1) | !min(clusn_sc) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_sc <- list()
    for(i in trt_lev) { clusnt_sc[[i]] <- clusn_sc[trt_index[[i]]]}
  }
  if(!is.null(se)) { 
    mf_se_fixed <- model.frame(formula = fixed_se, data = data2)
    z_e_fixed <- model.matrix(attr(mf_se_fixed, "terms"), data = mf_se_fixed)
    } else { mf_se_fixed <- z_e_fixed <- NULL}
  if(!is.null(sc)) { 
    mf_sc_fixed <- model.frame(formula = fixed_sc, data = data2)
    z_c_fixed <- model.matrix(attr(mf_sc_fixed, "terms"), data = mf_sc_fixed)
  } else { mf_sc_fixed <- z_c_fixed <- NULL}
  if(!is.null(se)) {
    if(!any(na.omit(s_eff) == 1)) { 
      stop("Please select a structural value in 'se' that is present in the data")}
    if(dim(z_e_fixed)[2] > 1) {
      colnames(z_e_fixed) <- c(colnames(z_e_fixed)[1], names(mf_se_fixed)[-1])}
    z_e_hold_fixed <- z_e_fixed
    if("e" %in% colnames(z_e_hold_fixed)) {
      z_e_fixed <- subset(z_e_hold_fixed, select = -c(e))}
    covz_e_fixed <- covz_e_center_fixed <- list()
    for(i in trt_lev) {
      covz_e_fixed[[i]] <- as.data.frame(z_e_fixed[trt_index[[i]], ])
      names(covz_e_fixed[[i]]) <- colnames(z_e_fixed)
      covz_e_center_fixed[[i]] <- as.data.frame(scale(covz_e_fixed[[i]], scale = FALSE))
      covz_e_center_fixed[[i]][, 1] <- rep(1, n[i])
    }  
    mean_covz_e_fixed <- lapply(covz_e_fixed, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    mean_covz_e_center_fixed <- lapply(covz_e_center_fixed, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    if(center) {
      covz_e_fixed <- covz_e_center_fixed
      mean_covz_e_fixed <- mean_covz_e_center_fixed
    }
    if(!is.null(random_se)){
      mf_se_random <- model.frame(formula = fname_re_se_coeff, data = data2)
      z_e_random <- model.matrix(attr(mf_se_random, "terms"), data = mf_se_random)
      if(no_random_int_se) {
        z_e_random <- as.matrix(z_e_random[, !colnames(z_e_random) == "(Intercept)"])
        if(dim(z_e_random)[2] == 1) { colnames(z_e_random) <- colnames(model.matrix(attr(mf_se_random, "terms"), data = mf_se_random))[2]}
      }
      z_e_hold_random <- z_e_random
      if("e" %in% colnames(z_e_hold_random)) {
        z_e_random <- subset(z_e_hold_random, select = -c(e))
      }
      n_clus_se <- setNames(table(clus_se), clus_se_lev)
      covz_e_random <- covz_e_center_random <- list()
      for(i in trt_lev) {
        covz_e_random[[i]] <- as.data.frame(z_e_random[trt_index[[i]], ])
        covz_e_center_random[[i]] <- as.data.frame(scale(covz_e_random[[i]], scale = FALSE))
        if(!no_random_int_se) {
          covz_e_center_random[[i]][, 1] <- rep(1, n[i])}
      }
      covz_e_random <- lapply(covz_e_random, setNames, colnames(z_e_random))
      covz_e_center_random <- lapply(covz_e_center_random, setNames, colnames(z_e_random))
      mean_covz_e_random <- lapply(covz_e_random, function(m) {
        apply(m, 2, mean, na.rm = TRUE)})
      mean_covz_e_center_random <- lapply(covz_e_center_random, function(m) {
        apply(m, 2, mean, na.rm = TRUE)})
      if(center) {
        covz_e_random <- covz_e_center_random
        mean_covz_e_random <- mean_covz_e_center_random
      }
    }
  }
  if(!is.null(sc)) {
    if(!any(na.omit(s_cost) == 1)) { 
      stop("Please select a structural value in 'sc' that is present in the data")}
    if(dim(z_c_fixed)[2] > 1) {
      colnames(z_c_fixed) <- c(colnames(z_c_fixed)[1], names(mf_sc_fixed)[-1])}
    z_c_hold_fixed <- z_c_fixed
    if("c" %in% colnames(z_c_hold_fixed)) {
      z_c_fixed <- subset(z_c_hold_fixed, select = -c(c))}
    covz_c_fixed <- covz_c_center_fixed <- list()
    for(i in trt_lev) {
      covz_c_fixed[[i]] <- as.data.frame(z_c_fixed[trt_index[[i]], ])
      names(covz_c_fixed[[i]]) <- colnames(z_c_fixed)
      covz_c_center_fixed[[i]] <- as.data.frame(scale(covz_c_fixed[[i]], scale = FALSE))
      covz_c_center_fixed[[i]][, 1] <- rep(1, n[i])
    }  
    mean_covz_c_fixed <- lapply(covz_c_fixed, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    mean_covz_c_center_fixed <- lapply(covz_c_center_fixed, function(m) {
      apply(m, 2, mean, na.rm = TRUE)})
    if(center) {
      covz_c_fixed <- covz_c_center_fixed
      mean_covz_c_fixed <- mean_covz_c_center_fixed
    }
    if(!is.null(random_sc)){
      mf_sc_random <- model.frame(formula = fname_re_sc_coeff, data = data2)
      z_c_random <- model.matrix(attr(mf_sc_random, "terms"), data = mf_sc_random)
      if(no_random_int_sc) {
        z_c_random <- as.matrix(z_c_random[, !colnames(z_c_random) == "(Intercept)"])
        if(dim(z_c_random)[2] == 1) { colnames(z_c_random) <- colnames(model.matrix(attr(mf_sc_random, "terms"), data = mf_sc_random))[2]}
      }
      z_c_hold_random <- z_c_random
      if("c" %in% colnames(z_c_hold_random)) {
        z_c_random <- subset(z_c_hold_random, select = -c(c))
      }
      n_clus_sc <- setNames(table(clus_sc), clus_sc_lev)
      covz_c_random <- covz_c_center_random <- list()
      for(i in trt_lev) {
        covz_c_random[[i]] <- as.data.frame(z_c_random[trt_index[[i]], ])
        covz_c_center_random[[i]] <- as.data.frame(scale(covz_c_random[[i]], scale = FALSE))
        if(!no_random_int_sc) {
          covz_c_center_random[[i]][, 1] <- rep(1, n[i])}
      }
      covz_c_random <- lapply(covz_c_random, setNames, colnames(z_c_random))
      covz_c_center_random <- lapply(covz_c_center_random, setNames, colnames(z_c_random))
      mean_covz_c_random <- lapply(covz_c_random, function(m) {
        apply(m, 2, mean, na.rm = TRUE)})
      mean_covz_c_center_random <- lapply(covz_c_center_random, function(m) {
        apply(m, 2, mean, na.rm = TRUE)})
      if(center) {
        covz_c_random <- covz_c_center_random
        mean_covz_c_random <- mean_covz_c_center_random
      }
    }
  }
  if(is.null(random_c)) {
    cov_c_random <- mean_cov_c_random <- x_c_random <- NULL
    clusnt_c <- n_clus_c <- name_clus_c <- clusn_c <- clus_c <- clus_c_lev <- NULL}
  if(is.null(random_sc)) { covz_c_random <- mean_covz_c_random <- z_c_random <- NULL
  clusnt_sc <- n_clus_sc <- name_clus_sc <- clusn_sc <- clus_sc <- clus_sc_lev <- NULL}
  if(is.null(random_e)) { cov_e_random <- mean_cov_e_random <- x_e_random <- NULL
  clusnt_e <- n_clus_e <- name_clus_e <- clusn_e <- clus_e <- clus_e_lev <- NULL}
  if(is.null(random_se)) { covz_e_random <- mean_covz_e_random <- z_e_random <- NULL
  clusnt_se <- n_clus_se <- name_clus_se <- clusn_se <- clus_se <- clus_se_lev <- NULL}
  data_raw <- setNames(list(y_e, y_c, m_eff, m_cost, s_eff, s_cost,
                            n_obs_eff, n_obs_cost, n_mis_eff, n_mis_cost,
                            n_ns_eff, n_ns_cost, n_s_eff, n_s_cost,
                            n, cov_e_fixed, cov_c_fixed, 
                            mean_cov_e_fixed, mean_cov_c_fixed, covz_e_fixed,
                            covz_c_fixed, mean_covz_e_fixed, mean_covz_c_fixed,
                            cov_e_random, cov_c_random, mean_cov_e_random, 
                            mean_cov_c_random, covz_e_random, covz_c_random, 
                            mean_covz_e_random, mean_covz_c_random, 
                            clusnt_e, clusnt_c, clusnt_se, clusnt_sc,
                            n_clus_e, n_clus_c, n_clus_se, n_clus_sc, data,
                            trt_name_pos_e, trt_name_pos_c, t_m_names_pos_e, t_m_names_pos_c,
                            trt_index, efft, costt, m_efft, m_costt, s_efft, s_costt, 
                            name_clus_e, name_clus_c,
                            name_clus_se, name_clus_sc, x_e_fixed, x_c_fixed,
                            x_e_random, x_c_random, z_e_fixed, z_c_fixed,
                            z_e_random, z_c_random, clusn_e, clusn_c, clusn_se, clusn_sc,
                            clus_e_lev, clus_c_lev, clus_se_lev, clus_sc_lev), 
                       c("e", "c", "me", "mc", "se", "sc", "n_obs_e", "n_obs_c", 
                         "n_mis_e", "n_mis_c", "n_ns_eff", "n_ns_cost", 
                         "n_s_eff", "n_s_cost", "n", "cov_fixed_e", 
                         "cov_fixed_c", "avg_cov_fixed_e", "avg_cov_fixed_c", 
                         "cov_fixed_se", "cov_fixed_sc", "avg_cov_fixed_se",
                         "avg_cov_fixed_sc", "cov_random_e", "cov_random_c",
                         "avg_cov_random_e", "avg_cov_random_c", "cov_random_se",
                         "cov_random_sc", "avg_cov_random_se", "avg_cov_random_sc",
                         "clus_e", "clus_c", "clus_se", "clus_sc",
                         "n_clus_e", "n_clus_c", "n_clus_se", "n_clus_sc", "data",
                         "trt_pos_e", "trt_pos_c", "trt_pos_se", "trt_pos_sc", "trt_index",
                         "efft", "costt", "m_efft", "m_costt", "s_efft", "s_costt", 
                         "name_clus_e", "name_clus_c",
                         "name_clus_se", "name_clus_sc", "x_e_fixed", "x_c_fixed",
                         "x_e_random", "x_c_random", "z_e_fixed", "z_c_fixed",
                         "z_e_random", "z_c_random", "clusn_e", "clusn_c", "clusn_se", "clusn_sc",
                         "clus_e_lev", "clus_c_lev", "clus_se_lev", "clus_sc_lev"))
  model_formula <- list("mf_model.e_fixed" = fixed_e, "mf_model.c_fixed" = fixed_c, 
                        "mf_model.se_fixed" = fixed_se, "mf_model.sc_fixed" = fixed_sc,
                        "mf_model.e_random" = fname_re_e_coeff, "mf_model.c_random" = fname_re_c_coeff, 
                        "mf_model.se_random" = fname_re_se_coeff, "mf_model.sc_random" = fname_re_sc_coeff)
  data_list <- list("data_raw" = data_raw, "model_formula" = model_formula)  
  return(data_list) 
}