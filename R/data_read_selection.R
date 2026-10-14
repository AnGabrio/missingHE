#' A function to read and re-arrange the data in different ways
#'
#' This internal function imports the data and outputs only those variables that are needed to run the model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator variables must always be provided and named 'e', 'c' and 'trt' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('e') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#'  Random effects can also be specified for each model parameter.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model. Random effects can also be specified for each model parameter.
#' @param model.me A formula expression in conventional \code{R} linear modelling syntax.  The response must be indicated with the 
#' term 'me'(missing effects) and any covariates used to estimate the probability of missing effects are given on the right-hand side. 
#' If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "probability" parameter for the missing effects through a logistic-linear model.
#' Random effects can also be specified for each model parameter.
#' @param model.mc A formula expression in conventional R linear modelling syntax. The response must be indicated with the 
#' term 'mc'(missing costs) and any covariates used to estimate the probability of missing costs should be given on the right-hand side. 
#' If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "probability" parameter for the missing costs through a logistic-linear model.
#' Random effects can also be specified for each model parameter.
#' @param cov_matrix Data frame containing the covariate matrix of the model.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @param trt_lev Vector of names of each treatment factor level
#' @param trt_pos Vector of name indices of each treatment factor level
#' @param fixed_e Fixed effects variables to be included in the effects model
#' @param fixed_c Fixed effects variables to be included in the costs model
#' @param random_e Random effects variables to be included in the effects model
#' @param random_c Random effects variables to be included in the costs model
#' @keywords internal
#' @importFrom stats setNames na.omit sd as.formula model.matrix model.frame model.response terms
#' @export
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

data_read_selection <- function(data, model.eff, model.cost, 
                                model.me, model.mc, cov_matrix, 
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
  fname_re_me_coeff <- as.formula(paste("me", "0", sep=" ~ "))
  fname_re_mc_coeff <- as.formula(paste("mc", "0", sep=" ~ "))
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
  cov_e_fixed <- cov_c_fixed <- list()
  cov_e_center_fixed <- cov_c_center_fixed <- list()
  x_c_hold_fixed <- x_c_fixed
  if("e" %in% colnames(x_c_hold_fixed)) {
    x_c_fixed <- subset(x_c_hold_fixed, select = -c(e))
  }
  efft <- costt <- m_efft <- m_costt <- list()
  n_obs_eff <- n_obs_cost <- setNames(trt_pos, trt_lev)
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
  }
  n_mis_eff <- n - n_obs_eff
  n_mis_cost <- n - n_obs_cost
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
  if(!inherits(model.me, "formula") | !inherits(model.mc, "formula")) {
    stop("'model.me' and/or 'model.mc' must be formula objects")}
  fixed_me <- nobars_(model.me)
  fixed_mc <- nobars_(model.mc)
  random_me <- fb(model.me)
  random_mc <- fb(model.mc)
  if(!is.null(random_me) & length(random_me) > 1 | !is.null(random_mc) & length(random_mc) > 1) {
    stop("Random effects can be included in the formula only through a single expression within brackets")}
  if(!all(names(model.frame(fixed_me, data = data2)) %in% c("me", "e", names(cov_matrix))) | 
     !all(names(model.frame(fixed_mc, data = data2)) %in% c("mc", "c", names(cov_matrix)))) {
    stop("Partially-observed covariates cannot be included in the model")}
  if(!all(names(model.frame(fixed_me, data = data2)) %in% names(data2)) | 
     !all(names(model.frame(fixed_mc, data = data2)) %in% names(data2))) {
    stop("Please provide names in the formula that correspond to those in the data")}
  if(names(model.frame(fixed_me, data = data2)[1]) != "me") {
    stop("Please set 'me' as the response in `model.me`")}
  if(names(model.frame(fixed_mc, data = data2)[1]) != "mc") {
    stop("Please set 'mc' as the response in `model.mc`")}
  if("c" %in% names(model.frame(fixed_me, data = data2)) | "e" %in% names(model.frame(fixed_mc, data = data2))) {
    stop("Please remove 'e'/'c' from the right-hand side of 'model.me' and/or 'c'/'e' from the right-hand side of 'model.mc'")}
  if("e" %in% labels(terms(fixed_me))) {
    if(length(grep(":e", labels(terms(fixed_me)))) != 0 | length(grep("e:", labels(terms(fixed_me)))) != 0) {
      stop("No interaction effects for 'e' are allowed")}
  }
  if("c" %in% labels(terms(fixed_mc))) {
    if(length(grep(":c", labels(terms(fixed_mc)))) != 0 | length(grep("c:", labels(terms(fixed_mc)))) != 0) {
      stop("No interaction effects for 'c' are allowed")} 
  }
  name_re_me_coeff <- NULL
  name_re_mc_coeff <- NULL
  if(!is.null(random_me)){
    name_re_me_coeff <- sub("\\|.*", "", random_me)
    if(grepl("0 + 1", name_re_me_coeff, fixed = TRUE)) { stop("Either remove or add a random intercept")}
    name_clus_me <- sub('.*\\|', '', random_me)
    if(lengths(strsplit(name_clus_me, " ")) > 2) {stop("Please provide a single clustering variable for each model")}
    name_clus_me <- gsub(" ", "", name_clus_me, fixed = TRUE)
    if(!name_clus_me %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_me_coeff, "")[[1]][1] == 0) {
      no_random_int_me <- TRUE} else {no_random_int_me <- FALSE}
    if(no_random_int_me) { 
      name_re_me_coeff <- sub("[0]", "", name_re_me_coeff) 
      name_re_me_coeff <- sub("[+]", "", name_re_me_coeff) 
      }
    if(name_re_me_coeff == "" | name_re_me_coeff == " ") { stop("Please state for which variables the random effects are assumed")}
    if(gsub(" ", "", name_re_me_coeff) == "e" & !no_random_int_me) {name_re_me_coeff <- "1 + e"}
    fname_re_me_coeff <- as.formula(paste("me", name_re_me_coeff, sep=" ~ "))
    if(!all(names(model.frame(fname_re_me_coeff, data = data2)) %in% c("0", "1", names(model.frame(fixed_me, data = data2))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("me" %in% labels(terms(fname_re_me_coeff))) {
      stop("Please remove 'me' from the random effects of 'model.me'")}
    if("e" %in% labels(terms(fname_re_me_coeff))) {
      if(length(grep(":e", labels(terms(fname_re_me_coeff)))) != 0 | length(grep("e:", labels(terms(fname_re_me_coeff)))) != 0) {
        stop("No interaction effects for 'e' are allowed")} 
    }
    clus_me <- data[, name_clus_me]
    clus_me_lev <- levels(clus_me)
    clus_me_pos <- which(clus_me_lev %in% levels(clus_me))
    if(!is.factor(clus_me)) { stop("Please define clustering variables as factors")}
    clusn_me <- as.numeric(clus_me)
    if(!all(diff(sort(unique(clusn_me))) == 1) | !min(clusn_me) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_me <- list()
    for(i in trt_lev) { clusnt_me[[i]] <- clusn_me[trt_index[[i]]]}
  }
  if(!is.null(random_mc)){
    name_re_mc_coeff <- sub("\\|.*", "", random_mc)
    if(grepl("0 + 1", name_re_mc_coeff, fixed = TRUE)) { stop("Either remove or add a random intercept")}
    name_clus_mc <- sub('.*\\|', '', random_mc)
    if(lengths(strsplit(name_clus_mc, " ")) > 2) {stop("Please provide a single clustering variable for each model")}
    name_clus_mc <- gsub(" ", "", name_clus_mc, fixed = TRUE)
    if(!name_clus_mc %in% names(cov_matrix)) { stop("Please provide the clustering variable in the dataset")}
    if(strsplit(name_re_mc_coeff, "")[[1]][1] == 0) {
      no_random_int_mc <- TRUE} else {no_random_int_mc <- FALSE}
    if(no_random_int_mc) { 
      name_re_mc_coeff <- sub("[0]", "", name_re_mc_coeff) 
      name_re_mc_coeff <- sub("[+]", "", name_re_mc_coeff) 
      }
    if(name_re_mc_coeff == "" | name_re_mc_coeff == " ") { stop("Please state for which variables the random effects are assumed")}
    if(gsub(" ", "", name_re_mc_coeff) == "c" & !no_random_int_mc) {name_re_mc_coeff <- "1 + c"}
    fname_re_mc_coeff <- as.formula(paste("mc", name_re_mc_coeff, sep=" ~ "))
    if(!all(names(model.frame(fname_re_mc_coeff, data = data2)) %in% c("0", "1", names(model.frame(fixed_mc, data = data2))))) {
      stop("Only covariates defined as fixed effects can be included as random effects")}
    if("mc" %in% labels(terms(fname_re_mc_coeff))) {
      stop("Please remove 'mc' from the random effects of 'model.mc'")}
    if("c" %in% labels(terms(fname_re_mc_coeff))) {
      if(length(grep(":c", labels(terms(fname_re_mc_coeff)))) != 0 | length(grep("c:", labels(terms(fname_re_mc_coeff)))) != 0) {
        stop("No interaction effects for 'c' are allowed")} 
      }
    clus_mc <- data[, name_clus_mc]
    clus_mc_lev <- levels(clus_mc)
    clus_mc_pos <- which(clus_mc_lev %in% levels(clus_mc))
    if(!is.factor(clus_mc)) { stop("Please define clustering variables as factors")}
    clusn_mc <- as.numeric(clus_mc)
    if(!all(diff(sort(unique(clusn_mc))) == 1) | !min(clusn_mc) == 1) {
      stop("Please make sure ordered levels of clustering variables have no gaps and start from 1")}
    clusnt_mc <- list()
    for(i in trt_lev) { clusnt_mc[[i]] <- clusn_mc[trt_index[[i]]]}
  }
  mf_me_fixed <- model.frame(formula = fixed_me, data = data2)
  mf_mc_fixed <- model.frame(formula = fixed_mc, data = data2)
  z_e_fixed <- model.matrix(attr(mf_me_fixed, "terms"), data = mf_me_fixed)
  z_c_fixed <- model.matrix(attr(mf_mc_fixed, "terms"), data = mf_mc_fixed)
  if(dim(z_e_fixed)[2] > 1) {
  colnames(z_e_fixed) <- c(colnames(z_e_fixed)[1], names(mf_me_fixed)[-1])}
  if(dim(z_c_fixed)[2] > 1) {
    colnames(z_c_fixed) <- c(colnames(z_c_fixed)[1], names(mf_mc_fixed)[-1])}
  z_e_hold_fixed <- z_e_fixed
  if("e" %in% colnames(z_e_hold_fixed)) {
    z_e_fixed <- subset(z_e_hold_fixed, select = -c(e))}
  z_c_hold_fixed <- z_c_fixed
  if("c" %in% colnames(z_c_hold_fixed)) {
    z_c_fixed <- subset(z_c_hold_fixed, select = -c(c))}
  covz_e_fixed <- covz_c_fixed <- list()
  covz_e_center_fixed <- covz_c_center_fixed <- list()
  for(i in trt_lev) {
    covz_e_fixed[[i]] <- as.data.frame(z_e_fixed[trt_index[[i]], ])
    covz_c_fixed[[i]] <- as.data.frame(z_c_fixed[trt_index[[i]], ])
    names(covz_e_fixed[[i]]) <- colnames(z_e_fixed)
    names(covz_c_fixed[[i]]) <- colnames(z_c_fixed)
    covz_e_center_fixed[[i]] <- as.data.frame(scale(covz_e_fixed[[i]], scale = FALSE))
    covz_c_center_fixed[[i]] <- as.data.frame(scale(covz_c_fixed[[i]], scale = FALSE))
    covz_e_center_fixed[[i]][, 1] <- rep(1, n[i])
    covz_c_center_fixed[[i]][, 1] <- rep(1, n[i])
  }
  mean_covz_e_fixed <- lapply(covz_e_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_covz_c_fixed <- lapply(covz_c_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_covz_e_center_fixed <- lapply(covz_e_center_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
  mean_covz_c_center_fixed <- lapply(covz_c_center_fixed, function(m) {
    apply(m, 2, mean, na.rm = TRUE)})
   if(center) {
    covz_e_fixed <- covz_e_center_fixed
    covz_c_fixed <- covz_c_center_fixed
    mean_covz_e_fixed <- mean_covz_e_center_fixed
    mean_covz_c_fixed <- mean_covz_c_center_fixed
  }
  if(!is.null(random_me)){
    mf_me_random <- model.frame(formula = fname_re_me_coeff, data = data2)
    z_e_random <- model.matrix(attr(mf_me_random, "terms"), data = mf_me_random)
    if(no_random_int_me) {
      z_e_random <- as.matrix(z_e_random[, !colnames(z_e_random) == "(Intercept)"])
      if(dim(z_e_random)[2] == 1) { colnames(z_e_random) <- colnames(model.matrix(attr(mf_me_random, "terms"), data = mf_me_random))[2]}
    }
    z_e_hold_random <- z_e_random
    if("e" %in% colnames(z_e_hold_random)) {
      z_e_random <- subset(z_e_hold_random, select = -c(e))
    }
    n_clus_me <- setNames(table(clus_me), clus_me_lev)
    covz_e_random <- covz_e_center_random <- list()
    for(i in trt_lev) {
      covz_e_random[[i]] <- as.data.frame(z_e_random[trt_index[[i]], ])
      covz_e_center_random[[i]] <- as.data.frame(scale(covz_e_random[[i]], scale = FALSE))
      if(!no_random_int_me) {
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
  if(!is.null(random_mc)){
    mf_mc_random <- model.frame(formula = fname_re_mc_coeff, data = data2)
    z_c_random <- model.matrix(attr(mf_mc_random, "terms"), data = mf_mc_random)
    if(no_random_int_mc) {
      z_c_random <- as.matrix(z_c_random[, !colnames(z_c_random) == "(Intercept)"])
      if(dim(z_c_random)[2] == 1) { colnames(z_c_random) <- colnames(model.matrix(attr(mf_mc_random, "terms"), data = mf_mc_random))[2]}
    }
    z_c_hold_random <- z_c_random
    if("c" %in% colnames(z_c_hold_random)) {
      z_c_random <- subset(z_c_hold_random, select = -c(c))
    }
    n_clus_mc <- setNames(table(clus_mc), clus_mc_lev)
    covz_c_random <- covz_c_center_random <- list()
    for(i in trt_lev) {
      covz_c_random[[i]] <- as.data.frame(z_c_random[trt_index[[i]], ])
      covz_c_center_random[[i]] <- as.data.frame(scale(covz_c_random[[i]], scale = FALSE))
      if(!no_random_int_mc) {
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
  if(is.null(random_c)) {
    cov_c_random <- mean_cov_c_random <- x_c_random <- NULL
    clusnt_c <- n_clus_c <- name_clus_c <- clusn_c <- clus_c <- clus_c_lev <- NULL}
  if(is.null(random_mc)) { covz_c_random <- mean_covz_c_random <- z_c_random <- NULL
    clusnt_mc <- n_clus_mc <- name_clus_mc <- clusn_mc <- clus_mc <- clus_mc_lev <- NULL}
  if(is.null(random_e)) { cov_e_random <- mean_cov_e_random <- x_e_random <- NULL
    clusnt_e <- n_clus_e <- name_clus_e <- clusn_e <- clus_e <- clus_e_lev <- NULL}
  if(is.null(random_me)) { covz_e_random <- mean_covz_e_random <- z_e_random <- NULL
    clusnt_me <- n_clus_me <- name_clus_me <- clusn_me <- clus_me <- clus_me_lev <- NULL}
  data_raw <- setNames(list(y_e, y_c, m_eff, m_cost, 
                            n_obs_eff, n_obs_cost, n_mis_eff, n_mis_cost,
                            n, cov_e_fixed, cov_c_fixed, 
                            mean_cov_e_fixed, mean_cov_c_fixed, covz_e_fixed,
                            covz_c_fixed, mean_covz_e_fixed, mean_covz_c_fixed,
                            cov_e_random, cov_c_random, mean_cov_e_random, 
                            mean_cov_c_random, covz_e_random, covz_c_random, 
                            mean_covz_e_random, mean_covz_c_random, 
                            clusnt_e, clusnt_c, clusnt_me, clusnt_mc,
                            n_clus_e, n_clus_c, n_clus_me, n_clus_mc, data,
                            trt_name_pos_e, trt_name_pos_c, t_m_names_pos_e, t_m_names_pos_c,
                            trt_index, efft, costt, m_efft, m_costt, name_clus_e, name_clus_c,
                            name_clus_me, name_clus_mc, x_e_fixed, x_c_fixed,
                            x_e_random, x_c_random, z_e_fixed, z_c_fixed,
                            z_e_random, z_c_random, clusn_e, clusn_c, clusn_me, clusn_mc,
                            clus_e_lev, clus_c_lev, clus_me_lev, clus_mc_lev), 
                       c("e", "c", "me", "mc", "n_obs_e", "n_obs_c", 
                         "n_mis_e", "n_mis_c", "n", "cov_fixed_e", 
                         "cov_fixed_c", "avg_cov_fixed_e", "avg_cov_fixed_c", 
                         "cov_fixed_me", "cov_fixed_mc", "avg_cov_fixed_me",
                         "avg_cov_fixed_mc", "cov_random_e", "cov_random_c",
                         "avg_cov_random_e", "avg_cov_random_c", "cov_random_me",
                         "cov_random_mc", "avg_cov_random_me", "avg_cov_random_mc",
                         "clus_e", "clus_c", "clus_me", "clus_mc",
                         "n_clus_e", "n_clus_c", "n_clus_me", "n_clus_mc", "data",
                         "trt_pos_e", "trt_pos_c", "trt_pos_me", "trt_pos_mc", "trt_index",
                         "efft", "costt", "m_efft", "m_costt", "name_clus_e", "name_clus_c",
                         "name_clus_me", "name_clus_mc", "x_e_fixed", "x_c_fixed",
                         "x_e_random", "x_c_random", "z_e_fixed", "z_c_fixed",
                         "z_e_random", "z_c_random", "clusn_e", "clusn_c", "clusn_me", "clusn_mc",
                         "clus_e_lev", "clus_c_lev", "clus_me_lev", "clus_mc_lev"))
  model_formula <- list("mf_model.e_fixed" = fixed_e, "mf_model.c_fixed" = fixed_c, 
                        "mf_model.me_fixed" = fixed_me, "mf_model.mc_fixed" = fixed_mc,
                        "mf_model.e_random" = fname_re_e_coeff, "mf_model.c_random" = fname_re_c_coeff, 
                        "mf_model.me_random" = fname_re_me_coeff, "mf_model.mc_random" = fname_re_mc_coeff)
  data_list <- list("data_raw" = data_raw, "model_formula" = model_formula)
  return(data_list) 
}