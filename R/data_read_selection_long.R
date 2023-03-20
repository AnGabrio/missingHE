#' A function to read and re-arrange the data in different ways
#'
#' This internal function imports the data and outputs only those variables that are needed to run the model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'u', 'c' and 't' respectively. 
#' @param model.eff A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  effectiveness outcome ('u') whose name must correspond to that used in \code{data}, and 
#'  any covariates are given on the right-hand side. If there are no covariates, specify \code{1} on the right hand side.
#'  By default, covariates are placed on the "location" parameter of the distribution through a linear model.
#'  Random effects can also be specified for each model parameter.
#' @param model.cost A formula expression in conventional \code{R} linear modelling syntax. The response must be a health economics
#'  cost outcome ('c') whose name must correspond to that used in \code{data}, and any covariates are given on the right-hand side.
#'  If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "location" 
#'  parameter of the distribution through a linear model. Random effects can also be specified for each model parameter.
#' @param model.mu A formula expression in conventional \code{R} linear modelling syntax.  The response must be indicated with the 
#' term 'mu'(missing effects) and any covariates used to estimate the probability of missing effects are given on the right-hand side. 
#' If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "probability" parameter for the missing effects through a logistic-linear model.
#' Random effects can also be specified for each model parameter.
#' @param model.mc A formula expression in conventional \code{R} linear modelling syntax. The response must be indicated with the 
#' term 'mc'(missing costs) and any covariates used to estimate the probability of missing costs should be given on the right-hand side. 
#' If there are no covariates, specify \code{1} on the right hand side. By default, covariates are placed on the "probability" parameter for the missing costs through a logistic-linear model.
#' Random effects can also be specified for each model parameter.
#' @param type Type of missingness mechanism assumed. Choices are Missing At Random (MAR) and Missing Not At Random (MNAR).
#' @param center Logical. If \code{center} is \code{TRUE} all the covariates in the model are centered.
#' @keywords read data
#' @importFrom stats na.omit sd as.formula model.matrix model.frame model.response terms aggregate reshape
#' @export
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

data_read_selection_long <- function(data, model.eff, model.cost, model.mu, model.mc, type, center) {
  if(is.data.frame(data) == FALSE) {
    stop("object data must be provided as data frame")
  }
  if(any(names(data) == "u") == TRUE & any(names(data) == "c") == TRUE) {
    u <- as.name("u")
    c <- as.name("c")
  }
  data.wide <- reshape(data, v.names = c("u", "c"), timevar = c("time"),
                       direction = "wide", idvar = "id", sep = ".")
  data.wide_u <- data.wide[, grepl('u\\.|id', colnames(data.wide))]
  wide_u_index <- grepl('u\\.|id', colnames(data.wide_u))
  data.wide_u <- data.wide_u[, wide_u_index]
  data.wide_u_obs_index <- as.data.frame(which(is.na(data.wide_u) == FALSE, arr.ind = TRUE))
  data.wide_u_obs_index <- merge(data.wide_u_obs_index, aggregate(col ~ row, data.wide_u_obs_index, max))
  data.wide_u_obs_index <- unique(data.wide_u_obs_index)
  names(data.wide_u_obs_index) <- c("id", "drop_u")
  data.wide_u_drop <- merge(data.wide, data.wide_u_obs_index, by.x = "id", all.x = TRUE)
  data.wide_u_drop$drop_u <- data.wide_u_drop$drop_u - 1
  data.wide_c <- data.wide[, grepl('c\\.|id', colnames(data.wide))]
  wide_c_index <- grepl('c\\.|id', colnames(data.wide_c))
  data.wide_c <- data.wide_c[, wide_c_index]
  data.wide_c_obs_index <- as.data.frame(which(is.na(data.wide_c) == FALSE, arr.ind = TRUE))
  data.wide_c_obs_index <- merge(data.wide_c_obs_index, aggregate(col ~ row, data.wide_c_obs_index, max))
  data.wide_c_obs_index <- unique(data.wide_c_obs_index)
  names(data.wide_c_obs_index) <- c("id", "drop_c")
  data.wide_uc_drop <- merge(data.wide_u_drop, data.wide_c_obs_index, by.x = "id", all.x = TRUE)
  data.wide_uc_drop$drop_c <- data.wide_uc_drop$drop_c - 1
  data.wide_uc_drop <- data.wide_uc_drop[, which(names(data.wide_uc_drop) %in% c("id","drop_u","drop_c"))]
  data.wide_drop <- merge(data.wide, data.wide_uc_drop, by = "id", all = TRUE)
  u_index_names <- colnames(data.wide_drop)[grepl('u\\.', colnames(data.wide_drop))]
  c_index_names <- colnames(data.wide_drop)[grepl('c\\.', colnames(data.wide_drop))]
  data.long <- reshape(data.wide_drop, varying = c(u_index_names,c_index_names),
                     direction = "long", idvar = "id", sep = ".")
  time_max <- max(data.long$time)
  count_mis_u <- rep(NA, length(data.long$drop_u))
  count_mis_u <- ifelse(is.na(data.long$u) == FALSE, 1, count_mis_u)
  count_mis_u <- ifelse(data.long$drop_u == 1 & is.na(data.long$u) == TRUE, 3, count_mis_u)
  count_mis_u <- ifelse(data.long$drop_u <= data.long$time & is.na(data.long$u) == TRUE, 3, count_mis_u)
  count_mis_u <- ifelse(data.long$drop_u > data.long$time & is.na(data.long$u) == TRUE, 2, count_mis_u)
  count_mis_c <- rep(NA, length(data.long$drop_c))
  count_mis_c <- ifelse(is.na(data.long$c) == FALSE, 1, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c == 1 & is.na(data.long$c) == TRUE, 3, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c <= data.long$time & is.na(data.long$c) == TRUE, 3, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c > data.long$time & is.na(data.long$c) == TRUE, 2, count_mis_c)
  data <- data.long
  my_subset <- function() {
    drop_u <- drop_c <- NULL
    cov_matrix <- subset(data, select = -c(u, c, drop_u, drop_c, time, t))[data$time == 1,]
    cov_matrix
  }
  cov_matrix <- my_subset()
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  is.formula <- function (x) { inherits(x, "formula") }
  if(is.formula(model.eff) == FALSE | is.formula(model.cost) == FALSE) {
    stop("model.eff and/or model.cost must be formula objects")
  }
  if(is.logical(center) == FALSE) { stop("center must be either TRUE or FALSE") }
  fixed_u <- nobars_(model.eff)
  fixed_c <- nobars_(model.cost)
  random_u <- fb(model.eff)
  random_c <- fb(model.cost)
  fname_re_u_coeff <- as.formula(paste("u", "0", sep=" ~ "))
  fname_re_c_coeff <- as.formula(paste("c", "0", sep=" ~ "))
  fname_re_mu_coeff <- as.formula(paste("mu", "0", sep=" ~ "))
  fname_re_mc_coeff <- as.formula(paste("mc", "0", sep=" ~ "))
  clusn_u <- clusn_c <- NULL
  clusn_mu <- clusn_mc <- NULL
  if(!is.null(random_u) & length(random_u) > 1 | !is.null(random_c) & length(random_c) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(all(names(model.frame(fixed_u, data = data)) %in% c("u", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_c, data = data)) %in% c("c", "u", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the fixed effects model")
  }
  if(all(names(model.frame(fixed_u, data = data)) %in% names(data)) == FALSE | 
     all(names(model.frame(fixed_c, data = data)) %in% names(data)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if("u" %in% labels(terms(fixed_u)) | "c" %in% labels(terms(fixed_c))) {
    stop("please remove 'u' from the right hand side of model.eff and/or 'c' from the right hand side of model.cost")
  }
  if(names(model.frame(fixed_u, data = data)[1]) != "u") {
    stop("you must set 'u' as the response in the formula model.eff")
  }
  if("c" %in% names(model.frame(fixed_u, data = data))) {
    stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
  }
  if(names(model.frame(fixed_c, data = data)[1]) != "c") {
    stop("you must set 'c' as the response in the formula model.cost")
  }
  if("u" %in% labels(terms(fixed_c))) {
    if(length(grep(":u", labels(terms(fixed_c)))) != 0 | length(grep("u:", labels(terms(fixed_c)))) != 0) {
      stop("no interaction effects for 'u' is allowed")
    } 
  }
  if("t" %in% names(model.frame(fixed_c, data = data)) | "t" %in% names(model.frame(fixed_u, data = data))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  index_mis_u <- which(is.na(data$u[data$time == 1]))
  index_mis_c <- which(is.na(data$c[data$time == 1]))
  data2 <- data
  data2$mu <- data2$drop_u
  data2$mc <- data2$drop_c
  data.t1 <- data[data$time == 1, ]
  data.t1$u[is.na(data.t1$u)] <- -999999
  data.t1$c[is.na(data.t1$c)] <- -999999
  mf_u_fixed <- model.frame(formula = fixed_u, data = data.t1)
  mf_c_fixed <- model.frame(formula = fixed_c, data = data.t1)
  terms <- NULL
  x_u_fixed <- model.matrix(attr(mf_u_fixed, "terms"), data = mf_u_fixed)
  x_c_fixed <- model.matrix(attr(mf_c_fixed, "terms"), data = mf_c_fixed)
  if("u" %in% names(mf_c_fixed)){
    mf_c_fixed$u[index_mis_u] <- NA
  }
  name_re_u_coeff <- NULL
  name_re_c_coeff <- NULL
  if(!is.null(random_u)){
    name_re_u_coeff <- sub("\\|.*", "", random_u)
    if(grepl("0 + 1", name_re_u_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_u <- sub('.*\\|', '', random_u)
    if(lengths(strsplit(name_clus_u, " ")) > 2) { stop("a single clustering variable must selected for each formula") }
    name_clus_u <- gsub(" ", "", name_clus_u, fixed = TRUE)
    if(!name_clus_u %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_u_coeff, "")[[1]][1] == 0) {
      no_random_int_u <- TRUE} else {no_random_int_u <- FALSE }
    if(no_random_int_u == TRUE) { 
      name_re_u_coeff <- sub("[0]", "", name_re_u_coeff) 
      name_re_u_coeff <- sub("[+]", "", name_re_u_coeff) 
    }
    if(name_re_u_coeff == "" | name_re_u_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    fname_re_u_coeff <- as.formula(paste("u", name_re_u_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_u_coeff, data = data.t1)) %in% c("0", "1", names(model.frame(fixed_u, data = data.t1)))) == FALSE) {
      stop("only covariates defined as fixed effects can be included in the random effects model")
    }
    if("u" %in% labels(terms(fname_re_u_coeff))) {
      stop("please remove 'u' from the random effects expression of model.eff")
    }
    if("c" %in% labels(terms(fname_re_u_coeff))) {
      stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
    }
    mf_u_random <- model.frame(formula = fname_re_u_coeff, data = data.t1)
    x_u_random <- model.matrix(attr(mf_u_random, "terms"), data = mf_u_random)
    if(no_random_int_u == TRUE) {
      x_u_random <- as.matrix(x_u_random[, !colnames(x_u_random) == "(Intercept)"])
      if(is.null(colnames(x_u_random)) == TRUE & dim(x_u_random)[2] == 1) {
        colnames(x_u_random) <- gsub(" ", "", name_re_u_coeff)
      }
    }
    clus_u <- data[, name_clus_u]
    if(!is.factor(clus_u)) { stop("clustering variables must be defined as factors") }
    clusn_u <- as.numeric(clus_u)
    if(!all(diff(sort(unique(clusn_u))) == 1) | !min(clusn_u) == 1) {
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
    if(gsub(" ", "", name_re_c_coeff) == "u" & no_random_int_c == FALSE) {name_re_c_coeff <- "1 + u" }
    fname_re_c_coeff <- as.formula(paste("c", name_re_c_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_c_coeff, data = data.t1)) %in% c("0", "1", names(model.frame(fixed_c, data = data.t1)))) == FALSE) {
      stop("only covariates defined as fixed effects can be included in the random effects model")
    }
    if("c" %in% labels(terms(fname_re_c_coeff))) {
      stop("please remove 'c' from the random effects expression of model.cost")
    }
    if("u" %in% labels(terms(fname_re_c_coeff))) {
      if(length(grep(":u", labels(terms(fname_re_c_coeff)))) != 0 | length(grep("u:", labels(terms(fname_re_c_coeff)))) != 0) {
        stop("no interaction effects for 'u' is allowed")
      } 
    }
    mf_c_random <- model.frame(formula = fname_re_c_coeff, data = data.t1)
    x_c_random <- model.matrix(attr(mf_c_random, "terms"), data = mf_c_random)
    if("u" %in% labels(terms(fname_re_c_coeff)) & length(labels(terms(fname_re_c_coeff))) == 1) {
      x_c_random <- subset(x_c_random, select = -c(u))
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
  N1 <- N2 <- c() 
  N1 <- length(data.t1$u[data.t1$t == 1])
  N2 <- length(data.t1$u[data.t1$t == 2])
  N <- c(N1, N2)
  n.time <- length(unique(data$time))
  eff1.mat <- cost1.mat <- matrix(NA, nrow = N1, ncol = n.time)
  eff2.mat <- cost2.mat <- matrix(NA, nrow = N2, ncol = n.time)
  m_eff <- count_mis_u
  m_cost <- count_mis_c
  time1.mat <- matrix(NA, nrow = N1, ncol = n.time)
  time2.mat <- matrix(NA, nrow = N2, ncol = n.time)
  m_eff1.mat <- m_cost1.mat <- matrix(NA, nrow = N1, ncol = n.time)
  m_eff2.mat <- m_cost2.mat <- matrix(NA, nrow = N2, ncol = n.time)
  for(i in 1:n.time) {
    eff1.mat[, i] <- data$u[data$time == i & data$t == 1] 
    eff2.mat[, i] <- data$u[data$time == i & data$t == 2] 
    cost1.mat[, i] <- data$c[data$time == i & data$t == 1] 
    cost2.mat[, i] <- data$c[data$time == i & data$t == 2] 
    time1.mat[, i] <- data$time[data$time == i & data$t == 1] 
    time2.mat[, i] <- data$time[data$time == i & data$t == 2] 
    m_eff1.mat[, i] <- count_mis_u[data$time == i & data$t == 1] 
    m_eff2.mat[, i] <- count_mis_u[data$time == i & data$t == 2] 
    m_cost1.mat[, i] <- count_mis_c[data$time == i & data$t == 1] 
    m_cost2.mat[, i] <- count_mis_c[data$time == i & data$t == 2] 
  }
  effects <- list(eff1.mat, eff2.mat)
  costs <- list(cost1.mat, cost2.mat)
  time <- list(time1.mat, time2.mat)
  m_eff <- list(m_eff1.mat, m_eff2.mat) 
  m_cost <- list(m_cost1.mat, m_cost2.mat) 
  N1_cc <- N2_cc <- N1_mis <- N2_mis <- matrix(NA, nrow = 2, ncol = n.time) 
  eff1_cc <- eff2_cc <- cost1_cc <- cost2_cc <- vector(mode = "list", length = n.time)
  for(i in 1:n.time){
  N1_cc[1, i] <- length(as.vector(apply(eff1.mat, 2, na.omit)[[i]]))
  N1_cc[2, i] <- length(as.vector(apply(cost1.mat, 2, na.omit)[[i]]))
  N2_cc[1, i] <- length(as.vector(apply(eff2.mat, 2, na.omit)[[i]]))
  N2_cc[2, i] <- length(as.vector(apply(cost2.mat, 2, na.omit)[[i]]))
  eff1_cc[[i]] <- as.vector(apply(eff1.mat, 2, na.omit)[[i]]) 
  cost1_cc[[i]] <- as.vector(apply(cost1.mat, 2, na.omit)[[i]]) 
  eff2_cc[[i]] <- as.vector(apply(eff2.mat, 2, na.omit)[[i]]) 
  cost2_cc[[i]] <- as.vector(apply(cost2.mat, 2, na.omit)[[i]]) 
  }
  N_cc <- list(N1_cc, N2_cc) 
  N1_mis <- N1 - N1_cc 
  N2_mis <- N2 - N2_cc 
  N_mis <- list(N1_mis, N2_mis) 
  eff_cc <- list(eff1_cc, eff2_cc) 
  cost_cc <- list(cost1_cc, cost2_cc)
  cov1_u_fixed <- as.data.frame(x_u_fixed[data.t1$t == 1, ])
  names(cov1_u_fixed) <- colnames(x_u_fixed)
  cov2_u_fixed <- as.data.frame(x_u_fixed[data.t1$t == 2, ])
  names(cov2_u_fixed) <- colnames(x_u_fixed)
  cov_u_fixed <- list(cov1_u_fixed, cov2_u_fixed)
  x_c_hold_fixed <- x_c_fixed
  if("u" %in% colnames(x_c_hold_fixed)) {
    x_c_fixed <- subset(x_c_hold_fixed, select = -c(u))
  }
  cov1_c_fixed <- as.data.frame(x_c_fixed[data.t1$t == 1, ])
  names(cov1_c_fixed) <- colnames(x_c_fixed)
  cov2_c_fixed <- as.data.frame(x_c_fixed[data.t1$t == 2, ])
  names(cov2_c_fixed) <- colnames(x_c_fixed)
  cov_c_fixed <- list(cov1_c_fixed, cov2_c_fixed)
  covu_fixed <- list(cov1_u_fixed, cov2_u_fixed) 
  mean_cov_u_fixed <- list(apply(as.matrix(cov1_u_fixed), 2, mean), apply(as.matrix(cov2_u_fixed), 2, mean))
  covc_fixed <- list(cov1_c_fixed, cov2_c_fixed) 
  mean_cov_c_fixed <- list(apply(as.matrix(cov1_c_fixed), 2, mean), apply(as.matrix(cov2_c_fixed), 2, mean))
  cov1_u_center_fixed <- as.data.frame(scale(cov1_u_fixed, scale = FALSE))
  cov2_u_center_fixed <- as.data.frame(scale(cov2_u_fixed, scale = FALSE))
  cov1_u_center_fixed[, 1] <- rep(1, nrow(cov1_u_fixed))
  cov2_u_center_fixed[, 1] <- rep(1, nrow(cov2_u_fixed))
  cov_u_center_fixed <- list(cov1_u_center_fixed, cov2_u_center_fixed)
  mean_cov_u_center_fixed <- list(apply(as.matrix(cov1_u_center_fixed), 2, mean), apply(as.matrix(cov2_u_center_fixed), 2, mean))
  cov1_c_center_fixed <- as.data.frame(scale(cov1_c_fixed, scale = FALSE))
  cov2_c_center_fixed <- as.data.frame(scale(cov2_c_fixed, scale = FALSE))
  cov1_c_center_fixed[, 1] <- rep(1, nrow(cov1_c_fixed))
  cov2_c_center_fixed[, 1] <- rep(1, nrow(cov2_c_fixed))
  cov_c_center_fixed <- list(cov1_c_center_fixed, cov2_c_center_fixed)
  mean_cov_c_center_fixed <- list(apply(as.matrix(cov1_c_center_fixed), 2, mean), apply(as.matrix(cov2_c_center_fixed), 2, mean))
  if(center == TRUE) {
    cov_u_fixed <- cov_u_center_fixed
    cov_c_fixed <- cov_c_center_fixed
    mean_cov_u_fixed <- mean_cov_u_center_fixed
    mean_cov_c_fixed <- mean_cov_c_center_fixed
  }
  if(!is.null(random_u)){
    cov1_u_random <- as.data.frame(x_u_random[data.t1$t == 1, ])
    names(cov1_u_random) <- colnames(x_u_random)
    cov2_u_random <- as.data.frame(x_u_random[data.t1$t == 2, ])
    names(cov2_u_random) <- colnames(x_u_random)
    cov_u_random <- list(cov1_u_random, cov2_u_random)
    covu_random <- list(cov1_u_random, cov2_u_random) 
    mean_cov_u_random <- list(apply(as.matrix(cov1_u_random), 2, mean), apply(as.matrix(cov2_u_random), 2, mean))
    cov1_u_center_random <- as.data.frame(scale(cov1_u_random, scale = FALSE))
    cov2_u_center_random <- as.data.frame(scale(cov2_u_random, scale = FALSE))
    if(no_random_int_u == FALSE) {
      cov1_u_center_random[, 1] <- rep(1, nrow(cov1_u_random))
      cov2_u_center_random[, 1] <- rep(1, nrow(cov2_u_random))
    }
    cov_u_center_random <- list(cov1_u_center_random, cov2_u_center_random)
    mean_cov_u_center_random <- list(apply(as.matrix(cov1_u_center_random), 2, mean), apply(as.matrix(cov2_u_center_random), 2, mean))
    if(center == TRUE) {
      cov_u_random <- cov_u_center_random
      mean_cov_u_random <- mean_cov_u_center_random
    }
    clusn_u1 <- clusn_u[data.t1$t == 1]
    clusn_u1 <- factor(clusn_u1, levels = unique(clusn_u1))
    clusn_u2 <- clusn_u[data.t1$t == 2]
    clusn_u2 <- factor(clusn_u2, levels = unique(clusn_u2))
  }
  if(!is.null(random_c)){
    x_c_hold_random <- x_c_random
    if("u" %in% colnames(x_c_hold_random)) {
      x_c_random <- subset(x_c_hold_random, select = -c(u))
    }
    cov1_c_random <- as.data.frame(x_c_random[data.t1$t == 1, ])
    names(cov1_c_random) <- colnames(x_c_random)
    cov2_c_random <- as.data.frame(x_c_random[data.t1$t == 2, ])
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
    clusn_c1 <- clusn_c[data.t1$t == 1]
    clusn_c1 <- factor(clusn_c1, levels = unique(clusn_c1))
    clusn_c2 <- clusn_c[data.t1$t == 2]
    clusn_c2 <- factor(clusn_c2, levels = unique(clusn_c2))
  }
  data2.t1 <- data[data$time == 1, ]
  data2.t1$u[is.na(data2.t1$u) == TRUE] <- -999999
  data2.t1$c[is.na(data2.t1$c) == TRUE] <- -999999
  data2.t1$mu <- c(m_eff1.mat[,1], m_eff2.mat[,1])
  data2.t1$mc <- c(m_cost1.mat[,1], m_cost2.mat[,1])
  if(!is.formula(model.mu) | !is.formula(model.mc)) {
    stop("model.mu and/or model.mc must be formula objects")
  }
  fixed_mu <- nobars_(model.mu)
  fixed_mc <- nobars_(model.mc)
  random_mu <- fb(model.mu)
  random_mc <- fb(model.mc)
  if(!is.null(random_mu) & length(random_mu) > 1 | !is.null(random_mc) & length(random_mc) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(all(names(model.frame(fixed_mu, data = data2.t1)) %in% c("mu", "u", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_mc, data = data2.t1)) %in% c("mc", "c", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the model")
  }
  if(all(names(model.frame(fixed_mu, data = data2.t1)) %in% names(data2.t1)) == FALSE | 
     all(names(model.frame(fixed_mc, data = data2.t1)) %in% names(data2.t1)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if(names(model.frame(fixed_mu, data = data2.t1)[1]) != "mu") {
    stop("you must set 'mu' as the response in the formula model.mu")
  }
  if(names(model.frame(fixed_mc, data = data2.t1)[1]) != "mc") {
    stop("you must set 'mc' as the response in the formula model.mc")
  }
  if("t" %in% names(model.frame(fixed_mc, data = data2.t1)) | "t" %in% names(model.frame(fixed_mu, data = data2.t1))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.mu' and/or 'model.mc'")
  }
  if("c" %in% names(model.frame(fixed_mu, data = data2.t1)) | "u" %in% names(model.frame(fixed_mc, data = data2.t1))) {
    stop("please remove 'u' from model.mc and/or remove 'c' from model.mu")
  }
  if("u" %in% labels(terms(fixed_mu))) {
    if(length(grep(":u", labels(terms(fixed_mu)))) != 0 | length(grep("u:", labels(terms(fixed_mu)))) != 0) {
      stop("no interaction effects for 'u' is allowed")
    } 
  }
  if("c" %in% labels(terms(fixed_mc))) {
    if(length(grep(":c", labels(terms(fixed_mc)))) != 0 | length(grep("c:", labels(terms(fixed_mc)))) != 0) {
      stop("no interaction effects for 'c' is allowed")
    } 
  }
  name_re_mu_coeff <- NULL
  name_re_mc_coeff <- NULL
  if(!is.null(random_mu)){
    name_re_mu_coeff <- sub("\\|.*", "", random_mu)
    if(grepl("0 + 1", name_re_mu_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_mu <- sub('.*\\|', '', random_mu)
    if(lengths(strsplit(name_clus_mu, " ")) > 2) {stop("a single clustering variable must be selected for each formula") }
    name_clus_mu <- gsub(" ", "", name_clus_mu, fixed = TRUE)
    if(!name_clus_mu %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_mu_coeff, "")[[1]][1] == 0) {
      no_random_int_mu <- TRUE} else {no_random_int_mu <- FALSE }
    if(no_random_int_mu == TRUE) { 
      name_re_mu_coeff <- sub("[0]", "", name_re_mu_coeff) 
      name_re_mu_coeff <- sub("[+]", "", name_re_mu_coeff) 
      }
    if(name_re_mu_coeff == "" | name_re_mu_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    if(gsub(" ", "", name_re_mu_coeff) == "u" & no_random_int_mu == FALSE) {name_re_mu_coeff <- "1 + u" }
    fname_re_mu_coeff <- as.formula(paste("mu", name_re_mu_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_mu_coeff, data = data2)) %in% c("0","1", names(model.frame(fixed_mu, data = data2)))) == FALSE) {
      stop("only covariates inlcued as fixed effects can be included in the random effects model")
    }
    if("mu" %in% labels(terms(fname_re_mu_coeff))) {
      stop("please remove 'mu' from the right hand side of model.mu")
    }
    if("u" %in% labels(terms(fname_re_mu_coeff))) {
      if(length(grep(":u", labels(terms(fname_re_mu_coeff)))) != 0 | length(grep("u:", labels(terms(fname_re_mu_coeff)))) != 0) {
        stop("no interaction effects for 'u' are allowed")
      } 
    }
    clus_mu <- data[, name_clus_mu]
    if(!is.factor(clus_mu)) { stop("clustering variables must be defined as factors") }
    clusn_mu <- as.numeric(clus_mu)
    if(!all(diff(sort(unique(clusn_mu))) == 1) | !min(clusn_mu) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  if(!is.null(random_mc)){
    name_re_mc_coeff <- sub("\\|.*", "", random_mc)
    if(grepl("0 + 1", name_re_mc_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_mc <- sub('.*\\|', '', random_mc)
    if(lengths(strsplit(name_clus_mc, " ")) > 2) {stop("a single clustering variable must be selected for each formula") }
    name_clus_mc <- gsub(" ", "", name_clus_mc, fixed = TRUE)
    if(!name_clus_mc %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_mc_coeff, "")[[1]][1] == 0) {
      no_random_int_mc <- TRUE} else {no_random_int_mc <- FALSE }
    if(no_random_int_mc == TRUE) { 
      name_re_mc_coeff <- sub("[0]", "", name_re_mc_coeff) 
      name_re_mc_coeff <- sub("[+]", "", name_re_mc_coeff) 
      }
    if(name_re_mc_coeff == "" | name_re_mc_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    if(gsub(" ", "", name_re_mc_coeff) == "c" & no_random_int_mc == FALSE) {name_re_mc_coeff <- "1 + c" }
    fname_re_mc_coeff <- as.formula(paste("mc", name_re_mc_coeff, sep=" ~ "))
    if(all(names(model.frame(fname_re_mc_coeff, data = data2)) %in% c("0","1", names(model.frame(fixed_mc, data = data2)))) == FALSE) {
      stop("only covariates inlcued as fixed effects can be included in the random effects model")
    }
    if("mc" %in% labels(terms(fname_re_mc_coeff))) {
      stop("please remove 'mc' from the right hand side of model.mc")
    }
    if("c" %in% labels(terms(fname_re_mc_coeff))) {
      if(length(grep(":c", labels(terms(fname_re_mc_coeff)))) != 0 | length(grep("c:", labels(terms(fname_re_mc_coeff)))) != 0) {
        stop("no interaction effects for 'c' is allowed")
      } 
    }
    clus_mc <- data[, name_clus_mc]
    if(!is.factor(clus_mc)) { stop("clustering variables must be defined as factors") }
    clusn_mc <- as.numeric(clus_mc)
    if(!all(diff(sort(unique(clusn_mc))) == 1) | !min(clusn_mc) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  mf_mu_fixed <- model.frame(formula = fixed_mu, data = data2.t1)
  mf_mc_fixed <- model.frame(formula = fixed_mc, data = data2.t1)
  z_u_fixed <- model.matrix(attr(mf_mu_fixed, "terms"), data = mf_mu_fixed)
  z_c_fixed <- model.matrix(attr(mf_mc_fixed, "terms"), data = mf_mc_fixed)
  z_u_hold_fixed <- z_u_fixed
  if("u" %in% colnames(z_u_hold_fixed)) {
    z_u_fixed <- subset(z_u_hold_fixed, select = -c(u))
  }
  z_c_hold_fixed <- z_c_fixed
  if("c" %in% colnames(z_c_hold_fixed)) {
    z_c_fixed <- subset(z_c_hold_fixed, select = -c(c))
  }
  covz1_u_fixed <- as.data.frame(z_u_fixed[data2.t1$t == 1, ])
  names(covz1_u_fixed) <- colnames(z_u_fixed)
  covz2_u_fixed <- as.data.frame(z_u_fixed[data2.t1$t == 2, ])
  names(covz2_u_fixed) <- colnames(z_u_fixed)
  covz_u_fixed <- list(covz1_u_fixed, covz2_u_fixed)
  covzu_fixed <- list(covz1_u_fixed, covz2_u_fixed) 
  mean_covz_u_fixed <- list(apply(as.matrix(covz1_u_fixed), 2, mean), apply(as.matrix(covz2_u_fixed), 2, mean))
  names(covz_u_fixed) <- names(mean_covz_u_fixed) <- c("Control", "Intervention")
  covz1_c_fixed <- as.data.frame(z_c_fixed[data2.t1$t == 1, ])
  names(covz1_c_fixed) <- colnames(z_c_fixed)
  covz2_c_fixed <- as.data.frame(z_c_fixed[data2.t1$t == 2, ])
  names(covz2_c_fixed) <- colnames(z_c_fixed)
  covz_c_fixed <- list(covz1_c_fixed, covz2_c_fixed)
  covzc_fixed <- list(covz1_c_fixed, covz2_c_fixed) 
  mean_covz_c_fixed <- list(apply(as.matrix(covz1_c_fixed), 2, mean), apply(as.matrix(covz2_c_fixed), 2, mean))
  covz1_u_center_fixed <- as.data.frame(scale(covz1_u_fixed, scale = FALSE))
  covz2_u_center_fixed <- as.data.frame(scale(covz2_u_fixed, scale = FALSE))
  covz1_u_center_fixed[, 1] <- rep(1, nrow(covz1_u_fixed))
  covz2_u_center_fixed[, 1] <- rep(1, nrow(covz2_u_fixed))
  covz_u_center_fixed <- list(covz1_u_center_fixed, covz2_u_center_fixed)
  mean_covz_u_center_fixed <- list(apply(as.matrix(covz1_u_center_fixed), 2, mean), apply(as.matrix(covz2_u_center_fixed), 2, mean))
  covz1_c_center_fixed <- as.data.frame(scale(covz1_c_fixed, scale = FALSE))
  covz2_c_center_fixed <- as.data.frame(scale(covz2_c_fixed, scale = FALSE))
  covz1_c_center_fixed[, 1] <- rep(1, nrow(covz1_c_fixed))
  covz2_c_center_fixed[, 1] <- rep(1, nrow(covz2_c_fixed))
  covz_c_center_fixed <- list(covz1_c_center_fixed, covz2_c_center_fixed)
  mean_covz_c_center_fixed <- list(apply(as.matrix(covz1_c_center_fixed), 2, mean), apply(as.matrix(covz2_c_center_fixed), 2, mean))
  if(center == TRUE) {
    covz_u_fixed <- covz_u_center_fixed
    covz_c_fixed <- covz_c_center_fixed
    mean_covz_u_fixed <- mean_covz_u_center_fixed
    mean_covz_c_fixed <- mean_covz_c_center_fixed
  }
  if(!is.null(random_mu)){
    mf_mu_random <- model.frame(formula = fname_re_mu_coeff, data = data2.t1)
    z_u_random <- model.matrix(attr(mf_mu_random, "terms"), data = mf_mu_random)
    if(no_random_int_mu == TRUE) {
      z_u_random <- as.matrix(z_u_random[, !colnames(z_u_random) == "(Intercept)"])
      if(dim(z_u_random)[2] == 1) { colnames(z_u_random) <- colnames(model.matrix(attr(mf_mu_random, "terms"), data = mf_mu_random))[2] }
    }
    z_u_hold_random <- z_u_random
    if("u" %in% colnames(z_u_hold_random)) {
      z_u_random <- subset(z_u_hold_random, select = -c(u))
    }
    covz1_u_random <- as.data.frame(z_u_random[data2.t1$t == 1, ])
    names(covz1_u_random) <- colnames(z_u_random)
    covz2_u_random <- as.data.frame(z_u_random[data2.t1$t == 2, ])
    names(covz2_u_random) <- colnames(z_u_random)
    covz_u_random <- list(covz1_u_random, covz2_u_random)
    covzu_random <- list(covz1_u_random, covz2_u_random) 
    mean_covz_u_random <- list(apply(as.matrix(covz1_u_random), 2, mean), apply(as.matrix(covz2_u_random), 2, mean))
    names(covz_u_random) <- names(mean_covz_u_random) <- c("Control", "Intervention")
    covz1_u_center_random <- as.data.frame(scale(covz1_u_random, scale = FALSE))
    covz2_u_center_random <- as.data.frame(scale(covz2_u_random, scale = FALSE))
    if(no_random_int_mu == FALSE) {
      covz1_u_center_random[, 1] <- rep(1, nrow(covz1_u_random))
      covz2_u_center_random[, 1] <- rep(1, nrow(covz2_u_random))
    }
    covz_u_center_random <- list(covz1_u_center_random, covz2_u_center_random)
    mean_covz_u_center_random <- list(apply(as.matrix(covz1_u_center_random), 2, mean), apply(as.matrix(covz2_u_center_random), 2, mean))
    if(center == TRUE) {
      covz_u_random <- covz_u_center_random
      mean_covz_u_random <- mean_covz_u_center_random
    }
    clusn_mu1 <- clusn_mu[data2.t1$t == 1]
    clusn_mu1 <- factor(clusn_mu1, levels = unique(clusn_mu1))
    clusn_mu2 <- clusn_mu[data2.t1$t == 2]
    clusn_mu2 <- factor(clusn_mu2, levels = unique(clusn_mu2))
  }
  if(!is.null(random_mc)){
    mf_mc_random <- model.frame(formula = fname_re_mc_coeff, data = data2.t1)
    z_c_random <- model.matrix(attr(mf_mc_random, "terms"), data = mf_mc_random)
    if(no_random_int_mc == TRUE) {
      z_c_random <- as.matrix(z_c_random[, !colnames(z_c_random) == "(Intercept)"])
      if(dim(z_c_random)[2] == 1) { colnames(z_c_random) <- colnames(model.matrix(attr(mf_mc_random, "terms"), data = mf_mc_random))[2] }
    }
    z_c_hold_random <- z_c_random
    if("c" %in% colnames(z_c_hold_random)) {
      z_c_random <- subset(z_c_hold_random, select = -c(c))
    }
    covz1_c_random <- as.data.frame(z_c_random[data2.t1$t == 1, ])
    names(covz1_c_random) <- colnames(z_c_random)
    covz2_c_random <- as.data.frame(z_c_random[data2.t1$t == 2, ])
    names(covz2_c_random) <- colnames(z_c_random)
    covz_c_random <- list(covz1_c_random, covz2_c_random)
    covzc_random <- list(covz1_c_random, covz2_c_random) 
    mean_covz_c_random <- list(apply(as.matrix(covz1_c_random), 2, mean), apply(as.matrix(covz2_c_random), 2, mean))
    names(covz_c_random) <- names(mean_covz_c_random) <- c("Control", "Intervention")
    covz1_c_center_random <- as.data.frame(scale(covz1_c_random, scale = FALSE))
    covz2_c_center_random <- as.data.frame(scale(covz2_c_random, scale = FALSE))
    if(no_random_int_mc == FALSE) {
      covz1_c_center_random[, 1] <- rep(1, nrow(covz1_c_random))
      covz2_c_center_random[, 1] <- rep(1, nrow(covz2_c_random))
    }
    covz_c_center_random <- list(covz1_c_center_random, covz2_c_center_random)
    mean_covz_c_center_random <- list(apply(as.matrix(covz1_c_center_random), 2, mean), apply(as.matrix(covz2_c_center_random), 2, mean))
    if(center == TRUE) {
      covz_c_random <- covz_c_center_random
      mean_covz_c_random <- mean_covz_c_center_random
    }
    clusn_mc1 <- clusn_mc[data2.t1$t == 1]
    clusn_mc1 <- factor(clusn_mc1, levels = unique(clusn_mc1))
    clusn_mc2 <- clusn_mc[data2.t1$t == 2]
    clusn_mc2 <- factor(clusn_mc2, levels = unique(clusn_mc2))
  }
  names(covz_u_fixed) <- names(mean_covz_u_fixed) <- names(covz_c_fixed) <- names(mean_covz_c_fixed) <- c("Control", "Intervention")
  names(cov_u_fixed) <- names(cov_c_fixed) <- names(mean_cov_u_fixed) <- names(mean_cov_c_fixed) <- c("Control", "Intervention")
  if(!is.null(random_c)) {
    names(cov_c_random) <- names(mean_cov_c_random) <- c("Control", "Intervention")
    clusn_c <- list("Control" = clusn_c1, "Intervention" = clusn_c2)
  } else {cov_c_random <- mean_cov_c_random <- NULL}
  if(!is.null(random_mc)) {
    names(covz_c_random) <- names(mean_covz_c_random) <- c("Control", "Intervention")
    clusn_mc <- list("Control" = clusn_mc1, "Intervention" = clusn_mc2)
  } else {covz_c_random <- mean_covz_c_random <- NULL}
  if(!is.null(random_u)) {
    names(cov_u_random) <- names(mean_cov_u_random) <- c("Control", "Intervention")
    clusn_u <- list("Control" = clusn_u1, "Intervention" = clusn_u2)
  } else {cov_u_random <- mean_cov_u_random <- NULL}
  if(!is.null(random_mu)) {
    names(covz_u_random) <- names(mean_covz_u_random) <- c("Control", "Intervention")
    clusn_mu <- list("Control" = clusn_mu1, "Intervention" = clusn_mu2)
  } else {covz_u_random <- mean_covz_u_random <- NULL}
  names(m_eff) <- names(m_cost) <- names(time) <- c("Control", "Intervention")
  names(effects) <- names(costs) <- names(eff_cc) <- names(cost_cc) <- c("Control", "Intervention")
  data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                   "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "missing_effects" = m_eff, "missing_costs" = m_cost, "time" = time,
                   "covariates_effects_fixed" = cov_u_fixed, "covariates_costs_fixed" = cov_c_fixed, "mean_cov_effects_fixed" = mean_cov_u_fixed, "mean_cov_costs_fixed" = mean_cov_c_fixed, 
                   "covariates_missing_effects_fixed" = covz_u_fixed, "mean_cov_missing_effects_fixed" = mean_covz_u_fixed, "covariates_missing_costs_fixed" = covz_c_fixed, 
                   "mean_cov_missing_costs_fixed" = mean_covz_c_fixed, "covariates_effects_random" = cov_u_random, "covariates_costs_random" = cov_c_random, "mean_cov_effects_random" = mean_cov_u_random, "mean_cov_costs_random" = mean_cov_c_random, 
                   "covariates_missing_effects_random" = covz_u_random, "mean_cov_missing_effects_random" = mean_covz_u_random, "covariates_missing_costs_random" = covz_c_random, 
                   "mean_cov_missing_costs_random" = mean_covz_c_random, "clus_u" = clusn_u, "clus_c" = clusn_c, "clus_mu" = clusn_mu, "clus_mc" = clusn_mc, "data_ind" = data2.t1, "data_long" = data) 
  model_formula <- list("mf_model.u_fixed" = fixed_u, "mf_model.c_fixed" = fixed_c, "mf_model.mu_fixed" = fixed_mu, "mf_model.mc_fixed" = fixed_mc,
                        "mf_model.u_random" = fname_re_u_coeff, "mf_model.c_random" = fname_re_c_coeff, "mf_model.mu_random" = fname_re_mu_coeff, "mf_model.mc_random" = fname_re_mc_coeff)
  data_list <- list("data_raw" = data_raw, "model_formula" = model_formula)
  return(data_list) 
}