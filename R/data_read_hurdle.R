#' A function to read and re-arrange the data in different ways for the hurdle model
#'
#' This internal function imports the data and outputs only those variables that are needed to run the hurdle model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff}, \code{model.cost} (model formulas for effects and costs) 
#' and \code{model.se}, \code{model.sc} (model formulas for the structural effect and cost models) . Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
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
#'  By default, covariates are placed on the "probability" parameter for the strcutural effects through a logistic-linear model. Random effects can also be specified for each model parameter.
#' @param model.sc A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and 
#'  any covariates used to estimate the probability of structural costs are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "probability" parameter for the strcutural costs through a logistic-linear model. Random effects can also be specified for each model parameter.
#' @param se Structural value to be found in the effect data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the effects is run.
#' @param sc Structural value to be found in the cost data defined in \code{data}. If set to \code{NULL}, 
#' no structural value is chosen and a standard model for the costs is run.
#' @param type Type of structural value mechanism assumed, either 'SCAR' (Structural Completely At Random) or 'SAR' (Strcutural At Random).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @keywords read data hurdle models
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response terms
#' @export
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

data_read_hurdle <- function(data, model.eff, model.cost, model.se, model.sc, se, sc, type, center) {
  if(is.data.frame(data) == FALSE) {
    stop("object data must be provided as data frame")
  }
  if(any(names(data) == "e") == TRUE & any(names(data) == "c") == TRUE) {
    e <- as.name("e")
    c <- as.name("c")
  }
  cov_matrix <- subset(data, select = -c(e, c))
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  is.formula<-function (x) { inherits(x, "formula") }
  if(is.formula(model.eff) == FALSE | is.formula(model.cost) == FALSE) {
    stop("model.eff and/or model.cost must be formula objects")
  }
  fixed_e <- nobars_(model.eff)
  fixed_c <- nobars_(model.cost)
  random_e <- fb(model.eff)
  random_c <- fb(model.cost)
  fname_re_e_coeff <- as.formula(paste("e", "0", sep=" ~ "))
  fname_re_c_coeff <- as.formula(paste("c", "0", sep=" ~ "))
  fname_re_se_coeff <- as.formula(paste("se", "0", sep=" ~ "))
  fname_re_sc_coeff <- as.formula(paste("sc", "0", sep=" ~ "))
  clusn_e <- clusn_c <- NULL
  clusn_se <- clusn_sc <- NULL
  if(!is.null(random_e) & length(random_e) > 1 | !is.null(random_c) & length(random_c) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(all(names(model.frame(fixed_e, data = data)) %in% c("e", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_c, data = data)) %in% c("c", "e", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the fixed effects model")
  }
  if("e" %in% labels(terms(fixed_e)) | "c" %in% labels(terms(fixed_c))) {
    stop("please remove 'e' from the right hand side of model.eff and/or 'c' from the right hand side of model.cost")
  }
  if(names(model.frame(fixed_e, data = data)[1]) != "e") {
    stop("you must set 'e' as the response in the formula model.eff")
  }
  if("c" %in% names(model.frame(fixed_e, data = data))) {
    stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
  }
  if(names(model.frame(fixed_c, data = data)[1]) != "c") {
    stop("you must set 'c' as the response in the formula model.cost")
  }
  if("e" %in% labels(terms(fixed_c))) {
    if(length(grep(":e", labels(terms(fixed_c)))) != 0 | length(grep("e:", labels(terms(fixed_c)))) != 0) {
      stop("no interaction effects for 'e' is allowed")
    } 
  }
  if("t" %in% names(model.frame(fixed_c, data = data)) | "t" %in% names(model.frame(fixed_e, data = data))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  index_mis_e <- index_mis2_e <- which(is.na(data$e))
  index_mis_c <- index_mis2_c <- which(is.na(data$c))
  data2 <- data
  data$e[is.na(data$e) == TRUE] <- -999999
  data$c[is.na(data$c) == TRUE] <- -999999
  mf_e_fixed <- model.frame(formula = fixed_e, data = data)
  mf_c_fixed <- model.frame(formula = fixed_c, data = data)
  terms <- NULL
  x_e_fixed <- model.matrix(attr(mf_e_fixed, "terms"), data = mf_e_fixed)
  x_c_fixed <- model.matrix(attr(mf_c_fixed, "terms"), data = mf_c_fixed)
  if("e" %in% names(mf_c_fixed)){
    mf_c_fixed$e[index_mis_e] <- NA
  }
  name_re_e_coeff <- NULL
  name_re_c_coeff <- NULL
  if(!is.null(random_e)){
    name_re_e_coeff <- sub("\\|.*", "", random_e)
    if(grepl("0 + 1", name_re_e_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_e <- sub('.*\\|', '', random_e)
    if(lengths(strsplit(name_clus_e, " ")) > 2) { stop("a single clustering variable must selected for each formula") }
    name_clus_e <- gsub(" ", "", name_clus_e, fixed = TRUE)
    if(!name_clus_e %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_e_coeff, "")[[1]][1] == 0) {
      no_random_int_e <- TRUE} else {no_random_int_e <- FALSE }
    if(no_random_int_e == TRUE) { 
      name_re_e_coeff <- sub("[0]", "", name_re_e_coeff) 
      name_re_e_coeff <- sub("[+]", "", name_re_e_coeff) 
    }
    if(name_re_e_coeff == "" | name_re_e_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    fname_re_e_coeff <- as.formula(paste("e", name_re_e_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_e_coeff, data = data)) %in% c("0", "1", names(model.frame(fixed_e, data = data)))) == FALSE) {
      stop("only covariates defined as fixed effects can be included in the random effects model")
    }
    if("e" %in% labels(terms(fname_re_e_coeff))) {
      stop("please remove 'e' from the random effects expression of model.eff")
    }
    if("c" %in% labels(terms(fname_re_e_coeff))) {
      stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
    }
    mf_e_random <- model.frame(formula = fname_re_e_coeff, data = data)
    x_e_random <- model.matrix(attr(mf_e_random, "terms"), data = mf_e_random)
    if(no_random_int_e == TRUE) {
      x_e_random <- as.matrix(x_e_random[, !colnames(x_e_random) == "(Intercept)"])
      if(is.null(colnames(x_e_random)) == TRUE & dim(x_e_random)[2] == 1) {
        colnames(x_e_random) <- gsub(" ", "", name_re_e_coeff)
      }
    }
    clus_e <- data[, name_clus_e]
    if(!is.factor(clus_e)) { stop("clustering variables must be defined as factors") }
    clusn_e <- as.numeric(clus_e)
    if(!all(diff(sort(unique(clusn_e))) == 1) | !min(clusn_e) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  if(!is.null(random_c)){
    name_re_c_coeff <- sub("\\|.*", "", random_c)
    if(grepl("0 + 1", name_re_c_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_c <- sub('.*\\|', '', random_c)
    if(lengths(strsplit(name_clus_c, " ")) > 2) { stop("a single clustering variable must selected for each formula") }
    name_clus_c <- gsub(" ", "", name_clus_c, fixed = TRUE)
    if(!name_clus_c %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_c_coeff, "")[[1]][1] == 0) {
      no_random_int_c <- TRUE} else {no_random_int_c <- FALSE }
    if(no_random_int_c == TRUE) { 
      name_re_c_coeff <- sub("[0]", "", name_re_c_coeff) 
      name_re_c_coeff <- sub("[+]", "", name_re_c_coeff) 
    }
    if(name_re_c_coeff == "" | name_re_c_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    if(gsub(" ", "", name_re_c_coeff) == "e" & no_random_int_c == FALSE) {name_re_c_coeff <- "1 + e" }
    fname_re_c_coeff <- as.formula(paste("c", name_re_c_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_c_coeff, data = data)) %in% c("0", "1", names(model.frame(fixed_c, data = data)))) == FALSE) {
      stop("only covariates defined as fixed effects can be included in the random effects model")
    }
    if("c" %in% labels(terms(fname_re_c_coeff))) {
      stop("please remove 'c' from the random effects expression of model.cost")
    }
    if("e" %in% labels(terms(fname_re_c_coeff))) {
      if(length(grep(":e", labels(terms(fname_re_c_coeff)))) != 0 | length(grep("e:", labels(terms(fname_re_c_coeff)))) != 0) {
        stop("no interaction effects for 'e' is allowed")
      } 
    }
    mf_c_random <- model.frame(formula = fname_re_c_coeff, data = data)
    x_c_random <- model.matrix(attr(mf_c_random, "terms"), data = mf_c_random)
    if("e" %in% labels(terms(fname_re_c_coeff)) & length(labels(terms(fname_re_c_coeff))) == 1) {
      x_c_random <- subset(x_c_random, select = -c(e))
    }
    if(no_random_int_c == TRUE) {
      x_c_random <- as.matrix(x_c_random[, !colnames(x_c_random) == "(Intercept)"])
      if(is.null(colnames(x_c_random)) == TRUE & dim(x_c_random)[2] == 1) {
        colnames(x_c_random) <- gsub(" ", "", name_re_c_coeff)
      }
    }
    clus_c <- data[, name_clus_c]
    if(!is.factor(clus_c)) { stop("clustering variables must be defined as factors") }
    clusn_c <- as.numeric(clus_c)
    if(!all(diff(sort(unique(clusn_c))) == 1) | !min(clusn_c) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  y_e <- model.response(mf_e_fixed)
  y_c <- model.response(mf_c_fixed)
  y_e[index_mis_e] <- NA
  y_c[index_mis_c] <- NA
  data$e[index_mis_e] <- NA
  data$c[index_mis_c] <- NA
  N1 <- N2 <- c() 
  N1 <- sum(data$t == 1) 
  N2 <- length(data$t) - N1 
  N <- c(N1, N2)
  m_eff <- rep(0, length(data$e))
  m_eff[index_mis_e] <- 1
  m_cost <- rep(0, length(data$c))
  m_cost[index_mis_c] <- 1
  m_eff1 <- m_eff2 <- m_cost1 <- m_cost2 <- c() 
  t1_index <- which(data$t == 1)
  t2_index <- which(data$t == 2)
  eff1 <- y_e[t1_index] 
  eff2 <- y_e[t2_index] 
  eff <- list(eff1, eff2)
  cost1 <- y_c[t1_index] 
  cost2 <- y_c[t2_index] 
  cost <- list(cost1, cost2)
  m_eff1 <- m_eff[t1_index]
  m_eff2 <- m_eff[t2_index]
  m_eff <- list(m_eff1, m_eff2) 
  m_cost1 <- m_cost[t1_index]
  m_cost2 <- m_cost[t2_index]
  m_cost <- list(m_cost1, m_cost2) 
  N1_cc <- N2_cc <- N1_mis <- N2_mis <- c() 
  N1_cc[1] <- length(na.omit(eff1)) 
  N1_cc[2] <- length(na.omit(cost1)) 
  N2_cc[1] <- length(na.omit(eff2)) 
  N2_cc[2] <- length(na.omit(cost2)) 
  N_cc <- cbind(N1_cc, N2_cc) 
  N1_mis <- N1 - N1_cc 
  N2_mis <- N2 - N2_cc 
  N_mis <- cbind(N1_mis, N2_mis) 
  effects <- list(eff1, eff2) 
  costs <- list(cost1, cost2) 
  eff1_cc <- eff2_cc <- cost1_cc <- cost2_cc <- c() 
  eff1_cc <- na.omit(eff1) 
  eff2_cc <- na.omit(eff2) 
  eff_cc <- list(eff1_cc, eff2_cc) 
  cost1_cc <- na.omit(cost1) 
  cost2_cc <- na.omit(cost2) 
  cost_cc <- list(cost1_cc, cost2_cc)
  cov1_e_fixed <- as.data.frame(x_e_fixed[t1_index, ])
  names(cov1_e_fixed) <- colnames(x_e_fixed)
  cov2_e_fixed <- as.data.frame(x_e_fixed[t2_index, ])
  names(cov2_e_fixed) <- colnames(x_e_fixed)
  cov_e_fixed <- list(cov1_e_fixed, cov2_e_fixed)
  x_c_hold_fixed <- x_c_fixed
  if("e" %in% colnames(x_c_hold_fixed)) {
    x_c_fixed <- subset(x_c_hold_fixed, select = -c(e))
  }
  cov1_c_fixed <- as.data.frame(x_c_fixed[t1_index, ])
  names(cov1_c_fixed) <- colnames(x_c_fixed)
  cov2_c_fixed <- as.data.frame(x_c_fixed[t2_index, ])
  names(cov2_c_fixed) <- colnames(x_c_fixed)
  cov_c_fixed <- list(cov1_c_fixed, cov2_c_fixed)
  cove_fixed <- list(cov1_e_fixed, cov2_e_fixed) 
  mean_cov_e_fixed <- list(apply(as.matrix(cov1_e_fixed), 2, mean), apply(as.matrix(cov2_e_fixed), 2, mean))
  covc_fixed <- list(cov1_c_fixed, cov2_c_fixed) 
  mean_cov_c_fixed <- list(apply(as.matrix(cov1_c_fixed), 2, mean), apply(as.matrix(cov2_c_fixed), 2, mean))
  cov1_e_center_fixed <- as.data.frame(scale(cov1_e_fixed, scale = FALSE))
  cov2_e_center_fixed <- as.data.frame(scale(cov2_e_fixed, scale = FALSE))
  cov1_e_center_fixed[, 1] <- rep(1, nrow(cov1_e_fixed))
  cov2_e_center_fixed[, 1] <- rep(1, nrow(cov2_e_fixed))
  cov_e_center_fixed <- list(cov1_e_center_fixed, cov2_e_center_fixed)
  mean_cov_e_center_fixed <- list(apply(as.matrix(cov1_e_center_fixed), 2, mean), apply(as.matrix(cov2_e_center_fixed), 2, mean))
  cov1_c_center_fixed <- as.data.frame(scale(cov1_c_fixed, scale = FALSE))
  cov2_c_center_fixed <- as.data.frame(scale(cov2_c_fixed, scale = FALSE))
  cov1_c_center_fixed[, 1] <- rep(1, nrow(cov1_c_fixed))
  cov2_c_center_fixed[, 1] <- rep(1, nrow(cov2_c_fixed))
  cov_c_center_fixed <- list(cov1_c_center_fixed, cov2_c_center_fixed)
  mean_cov_c_center_fixed <- list(apply(as.matrix(cov1_c_center_fixed), 2, mean), apply(as.matrix(cov2_c_center_fixed), 2, mean))
  if(center == TRUE) {
    cov_e_fixed <- cov_e_center_fixed
    cov_c_fixed <- cov_c_center_fixed
    mean_cov_e_fixed <- mean_cov_e_center_fixed
    mean_cov_c_fixed <- mean_cov_c_center_fixed
  }
  if(!is.null(random_e)){
    cov1_e_random <- as.data.frame(x_e_random[t1_index, ])
    names(cov1_e_random) <- colnames(x_e_random)
    cov2_e_random <- as.data.frame(x_e_random[t2_index, ])
    names(cov2_e_random) <- colnames(x_e_random)
    cov_e_random <- list(cov1_e_random, cov2_e_random)
    cove_random <- list(cov1_e_random, cov2_e_random) 
    mean_cov_e_random <- list(apply(as.matrix(cov1_e_random), 2, mean), apply(as.matrix(cov2_e_random), 2, mean))
    cov1_e_center_random <- as.data.frame(scale(cov1_e_random, scale = FALSE))
    cov2_e_center_random <- as.data.frame(scale(cov2_e_random, scale = FALSE))
    if(no_random_int_e == FALSE) {
      cov1_e_center_random[, 1] <- rep(1, nrow(cov1_e_random))
      cov2_e_center_random[, 1] <- rep(1, nrow(cov2_e_random))
    }
    cov_e_center_random <- list(cov1_e_center_random, cov2_e_center_random)
    mean_cov_e_center_random <- list(apply(as.matrix(cov1_e_center_random), 2, mean), apply(as.matrix(cov2_e_center_random), 2, mean))
    if(center == TRUE) {
      cov_e_random <- cov_e_center_random
      mean_cov_e_random <- mean_cov_e_center_random
    }
    clusn_e1 <- clusn_e[t1_index]
    clusn_e1 <- factor(clusn_e1, levels = unique(clusn_e1))
    clusn_e2 <- clusn_e[t2_index]
    clusn_e2 <- factor(clusn_e2, levels = unique(clusn_e2))
  }
  if(!is.null(random_c)){
    x_c_hold_random <- x_c_random
    if("e" %in% colnames(x_c_hold_random)) {
      x_c_random <- subset(x_c_hold_random, select = -c(e))
    }
    cov1_c_random <- as.data.frame(x_c_random[t1_index, ])
    names(cov1_c_random) <- colnames(x_c_random)
    cov2_c_random <- as.data.frame(x_c_random[t2_index, ])
    names(cov2_c_random) <- colnames(x_c_random)
    cov_c_random <- list(cov1_c_random, cov2_c_random)
    covc_random <- list(cov1_c_random, cov2_c_random) 
    mean_cov_c_random <- list(apply(as.matrix(cov1_c_random), 2, mean), apply(as.matrix(cov2_c_random), 2, mean))
    cov1_c_center_random <- as.data.frame(scale(cov1_c_random, scale = FALSE))
    cov2_c_center_random <- as.data.frame(scale(cov2_c_random, scale = FALSE))
    if(no_random_int_c == FALSE) {
      cov1_c_center_random[, 1] <- rep(1, nrow(cov1_c_random))
      cov2_c_center_random[, 1] <- rep(1, nrow(cov2_c_random))
    }
    cov_c_center_random <- list(cov1_c_center_random, cov2_c_center_random)
    mean_cov_c_center_random <- list(apply(as.matrix(cov1_c_center_random), 2, mean), apply(as.matrix(cov2_c_center_random), 2, mean))
    if(center == TRUE) {
      cov_c_random <- cov_c_center_random
      mean_cov_c_random <- mean_cov_c_center_random
    }
    clusn_c1 <- clusn_c[t1_index]
    clusn_c1 <- factor(clusn_c1, levels = unique(clusn_c1))
    clusn_c2 <- clusn_c[t2_index]
    clusn_c2 <- factor(clusn_c2, levels = unique(clusn_c2))
  }
  data2$e[is.na(data2$e) == TRUE] <- -999999
  data2$c[is.na(data2$c) == TRUE] <- -999999
  if(type == "SCAR" | type == "SAR") {
    data2$se <- c(m_eff1, m_eff2)
    data2$sc <- c(m_cost1, m_cost2)
  }
  if(!is.formula(model.se) | !is.formula(model.sc)) {
    stop("model.se and/or model.sc must be formula objects")
  }
  fixed_se <- nobars_(model.se)
  fixed_sc <- nobars_(model.sc)
  if(is.null(se) == TRUE) {fixed_se <- se ~ 1 }
  if(is.null(sc) == TRUE) {fixed_sc <- sc ~ 1 }
  random_se <- fb(model.se)
  random_sc <- fb(model.sc)
  if(is.null(se) == TRUE) {random_se <- NULL }
  if(is.null(sc) == TRUE) {random_sc <- NULL }
  if(!is.null(random_se) & length(random_se) > 1 | !is.null(random_sc) & length(random_sc) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(all(names(model.frame(fixed_se, data = data2)) %in% c("se", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_sc, data = data2)) %in% c("sc", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the model")
  }
  if(all(names(model.frame(fixed_se, data = data2)) %in% names(data2)) == FALSE | 
     all(names(model.frame(fixed_sc, data = data2)) %in% names(data2)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if(names(model.frame(fixed_se, data = data2)[1]) != "se") {
    stop("you must set 'se' as the response in the formula model.se")
  }
  if(names(model.frame(fixed_sc, data = data2)[1]) != "sc") {
    stop("you must set 'sc' as the response in the formula model.sc")
  }
  if("t" %in% names(model.frame(fixed_sc, data = data2)) | "t" %in% names(model.frame(fixed_se, data = data2))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.se' and/or 'model.sc'")
  }
  if("c" %in% labels(terms(fixed_se)) | "e" %in% labels(terms(fixed_se)) |
     "e" %in% labels(terms(fixed_sc)) | "c" %in% labels(terms(fixed_sc))) {
    stop("please remove 'e' and/or 'c' from model.se and/or model.sc")
  }
  name_re_se_coeff <- NULL
  name_re_sc_coeff <- NULL
  if(!is.null(random_se)){
    name_re_se_coeff <- sub("\\|.*", "", random_se)
    if(grepl("0 + 1", name_re_se_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_se <- sub('.*\\|', '', random_se)
    if(lengths(strsplit(name_clus_se, " ")) > 2) {stop("a single clustering variable must be selected for each formula") }
    name_clus_se <- gsub(" ", "", name_clus_se, fixed = TRUE)
    if(!name_clus_se %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_se_coeff, "")[[1]][1] == 0) {
      no_random_int_se <- TRUE} else {no_random_int_se <- FALSE }
    if(no_random_int_se == TRUE) { 
      name_re_se_coeff <- sub("[0]", "", name_re_se_coeff) 
      name_re_se_coeff <- sub("[+]", "", name_re_se_coeff) 
    }
    if(name_re_se_coeff == "" | name_re_se_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    fname_re_se_coeff <- as.formula(paste("se", name_re_se_coeff, sep=" ~ "))
    if(all(names(model.frame(fname_re_se_coeff, data = data2)) %in% c("0","1", names(model.frame(fixed_se, data = data2)))) == FALSE) {
      stop("only covariates inlcued as fixed effects can be included in the random effects model")
    }
    if("se" %in% labels(terms(fname_re_se_coeff))) {
      stop("please remove 'se' from the right hand side of model.se")
    }
    clus_se <- data[, name_clus_se]
    if(!is.factor(clus_se)) { stop("clustering variables must be defined as factors") }
    clusn_se <- as.numeric(clus_se)
    if(!all(diff(sort(unique(clusn_se))) == 1) | !min(clusn_se) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  if(!is.null(random_sc)){
    name_re_sc_coeff <- sub("\\|.*", "", random_sc)
    if(grepl("0 + 1", name_re_sc_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_sc <- sub('.*\\|', '', random_sc)
    if(lengths(strsplit(name_clus_sc, " ")) > 2) {stop("a single clustering variable must be selected for each formula") }
    name_clus_sc <- gsub(" ", "", name_clus_sc, fixed = TRUE)
    if(!name_clus_sc %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_sc_coeff, "")[[1]][1] == 0) {
      no_random_int_sc <- TRUE} else {no_random_int_sc <- FALSE }
    if(no_random_int_sc == TRUE) { 
      name_re_sc_coeff <- sub("[0]", "", name_re_sc_coeff) 
      name_re_sc_coeff <- sub("[+]", "", name_re_sc_coeff) 
    }
    if(name_re_sc_coeff == "" | name_re_sc_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    fname_re_sc_coeff <- as.formula(paste("sc", name_re_sc_coeff, sep=" ~ "))
    if(all(names(model.frame(fname_re_sc_coeff, data = data2)) %in% c("0", "1", names(model.frame(fixed_sc, data = data2)))) == FALSE) {
      stop("only covariates inlcued as fixed effects can be included in the random effects model")
    }
    if("sc" %in% labels(terms(fname_re_sc_coeff))) {
      stop("please remove 'sc' from the right hand side of model.sc")
    }
    clus_sc <- data[, name_clus_sc]
    if(!is.factor(clus_sc)) { stop("clustering variables must be defined as factors") }
    clusn_sc <- as.numeric(clus_sc)
    if(!all(diff(sort(unique(clusn_sc))) == 1) | !min(clusn_sc) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  mf_se_fixed <- model.frame(formula = fixed_se, data = data2)
  mf_sc_fixed <- model.frame(formula = fixed_sc, data = data2)
  z_e_fixed <- model.matrix(attr(mf_se_fixed, "terms"), data = mf_se_fixed)
  z_c_fixed <- model.matrix(attr(mf_sc_fixed, "terms"), data = mf_sc_fixed)
  y_e <- model.response(mf_e_fixed)
  y_c <- model.response(mf_c_fixed)
  y_e[index_mis2_e] <- NA
  y_c[index_mis2_c] <- NA
  data2$e[index_mis2_e] <- NA
  data2$c[index_mis2_c] <- NA
  if(is.null(se) == TRUE & is.null(sc) == TRUE) { stop("Structural values in at least one outcome variable are required, please provide the structural value") }
  if(is.null(se) == FALSE) {
    index_str_e <- ifelse(round(y_e, digits = 4) == se, 1, 0)
    if(any(na.omit(index_str_e) == 1) == FALSE) { stop("Provide structural values that are present in the data") }
    d_eff1 <- d_eff2 <- c()
    d_eff1 <- index_str_e[t1_index]
    d_eff2 <- index_str_e[t2_index]
    d_eff <- list(d_eff1, d_eff2) 
    covz1_e_fixed <- as.data.frame(z_e_fixed[t1_index, ])
    names(covz1_e_fixed) <- colnames(z_e_fixed)
    covz2_e_fixed <- as.data.frame(z_e_fixed[t2_index, ])
    names(covz2_e_fixed) <- colnames(z_e_fixed)
    covz_e_fixed <- list(covz1_e_fixed, covz2_e_fixed)
    covze_fixed <- list(covz1_e_fixed, covz2_e_fixed) 
    mean_covz_e_fixed <- list(apply(as.matrix(covz1_e_fixed), 2, mean), apply(as.matrix(covz2_e_fixed), 2, mean))
    names(covz_e_fixed) <- names(mean_covz_e_fixed) <- c("Control", "Intervention")
    covz1_e_center_fixed <- as.data.frame(scale(covz1_e_fixed, scale = FALSE))
    covz2_e_center_fixed <- as.data.frame(scale(covz2_e_fixed, scale = FALSE))
    covz1_e_center_fixed[, 1] <- rep(1, nrow(covz1_e_fixed))
    covz2_e_center_fixed[, 1] <- rep(1, nrow(covz2_e_fixed))
    covz_e_center_fixed <- list(covz1_e_center_fixed, covz2_e_center_fixed)
    mean_covz_e_center_fixed <- list(apply(as.matrix(covz1_e_center_fixed), 2, mean), apply(as.matrix(covz2_e_center_fixed), 2, mean))
    if(center == TRUE) {
      covz_e_fixed <- covz_e_center_fixed
      mean_covz_e_fixed <- mean_covz_e_center_fixed
    }
    names(covz_e_fixed) <- names(mean_covz_e_fixed) <- c("Control", "Intervention")
    names(d_eff) <- c("Control", "Intervention")
    if(!is.null(random_se)){
      mf_se_random <- model.frame(formula = fname_re_se_coeff, data = data2)
      z_e_random <- model.matrix(attr(mf_se_random, "terms"), data = mf_se_random)
      if(no_random_int_se == TRUE) {
        z_e_random <- as.matrix(z_e_random[, !colnames(z_e_random) == "(Intercept)"])
        if(dim(z_e_random)[2] == 1) { colnames(z_e_random) <- colnames(model.matrix(attr(mf_se_random, "terms"), data = mf_se_random))[2] }
      }
      covz1_e_random <- as.data.frame(z_e_random[t1_index, ])
      names(covz1_e_random) <- colnames(z_e_random)
      covz2_e_random <- as.data.frame(z_e_random[t2_index, ])
      names(covz2_e_random) <- colnames(z_e_random)
      covz_e_random <- list(covz1_e_random, covz2_e_random)
      covze_random <- list(covz1_e_random, covz2_e_random) 
      mean_covz_e_random <- list(apply(as.matrix(covz1_e_random), 2, mean), apply(as.matrix(covz2_e_random), 2, mean))
      names(covz_e_random) <- names(mean_covz_e_random) <- c("Control", "Intervention")
      covz1_e_center_random <- as.data.frame(scale(covz1_e_random, scale = FALSE))
      covz2_e_center_random <- as.data.frame(scale(covz2_e_random, scale = FALSE))
      if(no_random_int_se == FALSE) {
        covz1_e_center_random[, 1] <- rep(1, nrow(covz1_e_random))
        covz2_e_center_random[, 1] <- rep(1, nrow(covz2_e_random))
      }
      covz_e_center_random <- list(covz1_e_center_random, covz2_e_center_random)
      mean_covz_e_center_random <- list(apply(as.matrix(covz1_e_center_random), 2, mean), apply(as.matrix(covz2_e_center_random), 2, mean))
      if(center == TRUE) {
        covz_e_random <- covz_e_center_random
        mean_covz_e_random <- mean_covz_e_center_random
      }
      clusn_se1 <- clusn_se[t1_index]
      clusn_se1 <- factor(clusn_se1, levels = unique(clusn_se1))
      clusn_se2 <- clusn_se[t2_index]
      clusn_se2 <- factor(clusn_se2, levels = unique(clusn_se2))
    }
  }
  if(is.null(sc) == FALSE) {
    index_str_c <- ifelse(round(y_c, digits = 0) == sc, 1, 0)
    if(any(na.omit(index_str_c) == 1) == FALSE){ stop("Provide structural values that are present in the data") }
    d_cost1 <- d_cost2 <- c() 
    d_cost1 <- index_str_c[t1_index]
    d_cost2 <- index_str_c[t2_index]
    d_cost <- list(d_cost1, d_cost2)  
    covz1_c_fixed <- as.data.frame(z_c_fixed[t1_index, ])
    names(covz1_c_fixed) <- colnames(z_c_fixed)
    covz2_c_fixed <- as.data.frame(z_c_fixed[t2_index, ])
    names(covz2_c_fixed) <- colnames(z_c_fixed)
    covz_c_fixed <- list(covz1_c_fixed, covz2_c_fixed)
    covzc_fixed <- list(covz1_c_fixed, covz2_c_fixed) 
    mean_covz_c_fixed <- list(apply(as.matrix(covz1_c_fixed), 2, mean), apply(as.matrix(covz2_c_fixed), 2, mean))
    covz1_c_center_fixed <- as.data.frame(scale(covz1_c_fixed, scale = FALSE))
    covz2_c_center_fixed <- as.data.frame(scale(covz2_c_fixed, scale = FALSE))
    covz1_c_center_fixed[, 1] <- rep(1, nrow(covz1_c_fixed))
    covz2_c_center_fixed[, 1] <- rep(1, nrow(covz2_c_fixed))
    covz_c_center_fixed <- list(covz1_c_center_fixed, covz2_c_center_fixed)
    mean_covz_c_center_fixed <- list(apply(as.matrix(covz1_c_center_fixed), 2, mean), apply(as.matrix(covz2_c_center_fixed), 2, mean))
    if(center == TRUE) {
      covz_c_fixed <- covz_c_center_fixed
      mean_covz_c_fixed <- mean_covz_c_center_fixed
    }
    names(covz_c_fixed) <- names(mean_covz_c_fixed) <- c("Control", "Intervention")
    names(d_cost) <- c("Control", "Intervention")
    if(!is.null(random_sc)){
      mf_sc_random <- model.frame(formula = fname_re_sc_coeff, data = data2)
      z_c_random <- model.matrix(attr(mf_sc_random, "terms"), data = mf_sc_random)
      if(no_random_int_sc == TRUE) {
        z_c_random <- as.matrix(z_c_random[, !colnames(z_c_random) == "(Intercept)"])
        if(dim(z_c_random)[2] == 1) { colnames(z_c_random) <- colnames(model.matrix(attr(mf_sc_random, "terms"), data = mf_sc_random))[2] }
      }
      covz1_c_random <- as.data.frame(z_c_random[t1_index, ])
      names(covz1_c_random) <- colnames(z_c_random)
      covz2_c_random <- as.data.frame(z_c_random[t2_index, ])
      names(covz2_c_random) <- colnames(z_c_random)
      covz_c_random <- list(covz1_c_random, covz2_c_random)
      covzc_random <- list(covz1_c_random, covz2_c_random) 
      mean_covz_c_random <- list(apply(as.matrix(covz1_c_random), 2, mean), apply(as.matrix(covz2_c_random), 2, mean))
      names(covz_c_random) <- names(mean_covz_c_random) <- c("Control", "Intervention")
      covz1_c_center_random <- as.data.frame(scale(covz1_c_random, scale = FALSE))
      covz2_c_center_random <- as.data.frame(scale(covz2_c_random, scale = FALSE))
      if(no_random_int_sc == FALSE) {
        covz1_c_center_random[, 1] <- rep(1, nrow(covz1_c_random))
        covz2_c_center_random[, 1] <- rep(1, nrow(covz2_c_random))
      }
      covz_c_center_random <- list(covz1_c_center_random, covz2_c_center_random)
      mean_covz_c_center_random <- list(apply(as.matrix(covz1_c_center_random), 2, mean), apply(as.matrix(covz2_c_center_random), 2, mean))
      if(center == TRUE) {
        covz_c_random <- covz_c_center_random
        mean_covz_c_random <- mean_covz_c_center_random
      }
      clusn_sc1 <- clusn_sc[t1_index]
      clusn_sc1 <- factor(clusn_sc1, levels = unique(clusn_sc1))
      clusn_sc2 <- clusn_sc[t2_index]
      clusn_sc2 <- factor(clusn_sc2, levels = unique(clusn_sc2))
    }
  }
  if(is.null(se) == FALSE) { names(covz_e_fixed) <- names(mean_covz_e_fixed) <- c("Control", "Intervention") }
  if(is.null(sc) == FALSE) { names(covz_c_fixed) <- names(mean_covz_c_fixed) <- c("Control", "Intervention") }
  names(cov_e_fixed) <- names(cov_c_fixed) <- names(mean_cov_e_fixed) <- names(mean_cov_c_fixed) <- c("Control", "Intervention")
  if(!is.null(random_c)) {
    names(cov_c_random) <- names(mean_cov_c_random) <- c("Control", "Intervention")
    clusn_c <- list("Control" = clusn_c1, "Intervention" = clusn_c2)
  } else {cov_c_random <- mean_cov_c_random <- NULL}
  if(!is.null(random_sc)) {
    names(covz_c_random) <- names(mean_covz_c_random) <- c("Control", "Intervention")
    clusn_sc <- list("Control" = clusn_sc1, "Intervention" = clusn_sc2)
  } else {covz_c_random <- mean_covz_c_random <- NULL}
  if(!is.null(random_e)) {
    names(cov_e_random) <- names(mean_cov_e_random) <- c("Control", "Intervention")
    clusn_e <- list("Control" = clusn_e1, "Intervention" = clusn_e2)
  } else {cov_e_random <- mean_cov_e_random <- NULL}
  if(!is.null(random_se)) {
    names(covz_e_random) <- names(mean_covz_e_random) <- c("Control", "Intervention")
    clusn_se <- list("Control" = clusn_se1, "Intervention" = clusn_se2)
  } else {covz_e_random <- mean_covz_e_random <- NULL}
  names(m_eff) <- names(m_cost) <- c("Control", "Intervention")
  names(effects) <- names(costs) <- names(eff_cc) <- names(cost_cc) <- c("Control", "Intervention")
  if(is.null(se) == FALSE & is.null(sc) == TRUE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects_fixed" = cov_e_fixed, "covariates_costs_fixed" = cov_c_fixed, "mean_cov_effects_fixed" = mean_cov_e_fixed, "mean_cov_costs_fixed" = mean_cov_c_fixed, 
                     "covariates_structural_effects_fixed" = covz_e_fixed, "mean_cov_structural_effects_fixed" = mean_covz_e_fixed, 
                     "covariates_effects_random" = cov_e_random, "covariates_costs_random" = cov_c_random, "mean_cov_effects_random" = mean_cov_e_random, "mean_cov_costs_random" = mean_cov_c_random, 
                     "covariates_structural_effects_random" = covz_e_random, "mean_cov_structural_effects_random" = mean_covz_e_random, 
                     "clus_e" = clusn_e, "clus_c" = clusn_c, "clus_se" = clusn_se, "structural_effects" = d_eff, "data_ind" = data2)
  } else if(is.null(se) == TRUE & is.null(sc) == FALSE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects_fixed" = cov_e_fixed, "covariates_costs_fixed" = cov_c_fixed, "mean_cov_effects_fixed" = mean_cov_e_fixed, "mean_cov_costs_fixed" = mean_cov_c_fixed, 
                     "covariates_structural_costs_fixed" = covz_c_fixed, "mean_cov_structural_costs_fixed" = mean_covz_c_fixed, 
                     "covariates_costs_random" = cov_c_random, "covariates_effects_random" = cov_e_random, "mean_cov_effects_random" = mean_cov_e_random, "mean_cov_costs_random" = mean_cov_c_random, 
                     "covariates_structural_costs_random" = covz_c_random, "mean_cov_structural_costs_random" = mean_covz_c_random, 
                     "clus_e" = clusn_e, "clus_c" = clusn_c, "clus_sc" = clusn_sc, "structural_costs" = d_cost, "data_ind" = data2)
  }else if(is.null(se) == FALSE & is.null(sc) == FALSE) {
    data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                     "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, 
                     "covariates_effects_fixed" = cov_e_fixed, "covariates_costs_fixed" = cov_c_fixed, "mean_cov_effects_fixed" = mean_cov_e_fixed, "mean_cov_costs_fixed" = mean_cov_c_fixed, 
                     "covariates_structural_effects_fixed" = covz_e_fixed, "mean_cov_structural_effects_fixed" = mean_covz_e_fixed, "covariates_structural_costs_fixed" = covz_c_fixed, 
                     "mean_cov_structural_costs_fixed" = mean_covz_c_fixed, "covariates_effects_random" = cov_e_random, "covariates_costs_random" = cov_c_random, "mean_cov_effects_random" = mean_cov_e_random, "mean_cov_costs_random" = mean_cov_c_random, 
                     "covariates_structural_effects_random" = covz_e_random, "mean_cov_structural_effects_random" = mean_covz_e_random, "covariates_structural_costs_random" = covz_c_random, 
                     "mean_cov_structural_costs_random" = mean_covz_c_random, "clus_e" = clusn_e, "clus_c" = clusn_c, "clus_se" = clusn_se, "clus_sc" = clusn_sc, 
                     "structural_effects" = d_eff, "structural_costs" = d_cost, "data_ind" = data2)
  }
  model_formula <- list("mf_model.e_fixed" = fixed_e, "mf_model.c_fixed" = fixed_c, "mf_model.se_fixed" = fixed_se, "mf_model.sc_fixed" = fixed_sc,
                        "mf_model.e_random" = fname_re_e_coeff, "mf_model.c_random" = fname_re_c_coeff, "mf_model.se_random" = fname_re_se_coeff, "mf_model.sc_random" = fname_re_sc_coeff)
  data_list <- list("data_raw" = data_raw, "model_formula" = model_formula)
  return(data_list) 
}