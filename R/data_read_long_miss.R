#' A function to read and re-arrange the data in different ways
#'
#' This internal function imports the data and outputs only those variables that are needed to run the model
#' according to the information provided by the user.
#' @param data A data frame in which to find variables supplied in \code{model.eff} and \code{model.cost}. Among these,
#' effectiveness, cost and treatment indicator (only two arms) variables must always be provided and named 'e', 'c' and 't' respectively. 
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

data_read_long_miss <- function(data, model.eff, model.cost, model.me, model.mc, type, center) {
  if(is.data.frame(data) == FALSE) {
    stop("object data must be provided as data frame")
  }
  if(any(names(data) == "e") == TRUE & any(names(data) == "c") == TRUE) {
    e <- as.name("e")
    c <- as.name("c")
  }
  data.wide <- reshape(data, v.names = c("e", "c"), timevar = c("time"),
                       direction = "wide", idvar = "id", sep = ".")
  data.wide_e <- data.wide[, grepl('e\\.|id', colnames(data.wide))]
  wide_e_index <- grepl('e\\.|id', colnames(data.wide_e))
  data.wide_e <- data.wide_e[, wide_e_index]
  data.wide_e_obs_index <- as.data.frame(which(is.na(data.wide_e) == FALSE, arr.ind = TRUE))
  data.wide_e_obs_index <- merge(data.wide_e_obs_index, aggregate(col ~ row, data.wide_e_obs_index, max))
  data.wide_e_obs_index <- unique(data.wide_e_obs_index)
  names(data.wide_e_obs_index) <- c("id", "drop_e")
  data.wide_e_drop <- merge(data.wide, data.wide_e_obs_index, by.x = "id", all.x = TRUE)
  data.wide_e_drop$drop_e <- data.wide_e_drop$drop_e - 1
  data.wide_c <- data.wide[, grepl('c\\.|id', colnames(data.wide))]
  wide_c_index <- grepl('c\\.|id', colnames(data.wide_c))
  data.wide_c <- data.wide_c[, wide_c_index]
  data.wide_c_obs_index <- as.data.frame(which(is.na(data.wide_c) == FALSE, arr.ind = TRUE))
  data.wide_c_obs_index <- merge(data.wide_c_obs_index, aggregate(col ~ row, data.wide_c_obs_index, max))
  data.wide_c_obs_index <- unique(data.wide_c_obs_index)
  names(data.wide_c_obs_index) <- c("id", "drop_c")
  data.wide_cc_drop <- merge(data.wide_e_drop, data.wide_c_obs_index, by.x = "id", all.x = TRUE)
  data.wide_cc_drop$drop_c <- data.wide_cc_drop$drop_c - 1
  data.wide_cc_drop <- data.wide_cc_drop[, which(names(data.wide_cc_drop) %in% c("id","drop_e","drop_c"))]
  data.wide_drop <- merge(data.wide, data.wide_cc_drop, by = "id", all = TRUE)
  e_index_names <- colnames(data.wide_drop)[grepl('e\\.', colnames(data.wide_drop))]
  c_index_names <- colnames(data.wide_drop)[grepl('c\\.', colnames(data.wide_drop))]
  data.long <- reshape(data.wide_drop, varying = c(e_index_names,c_index_names),
                     direction = "long", idvar = "id", sep = ".")
  time_max <- max(data.long$time)
  count_mis_e <- rep(NA, length(data.long$drop_e))
  count_mis_e <- ifelse(is.na(data.long$e) == FALSE, 1, count_mis_e)
  count_mis_e <- ifelse(data.long$drop_e == 1 & is.na(data.long$e) == TRUE, 3, count_mis_e)
  count_mis_e <- ifelse(data.long$drop_e <= data.long$time & is.na(data.long$e) == TRUE, 3, count_mis_e)
  count_mis_e <- ifelse(data.long$drop_e > data.long$time & is.na(data.long$e) == TRUE, 2, count_mis_e)
  count_mis_c <- rep(NA, length(data.long$drop_c))
  count_mis_c <- ifelse(is.na(data.long$c) == FALSE, 1, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c == 1 & is.na(data.long$c) == TRUE, 3, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c <= data.long$time & is.na(data.long$c) == TRUE, 3, count_mis_c)
  count_mis_c <- ifelse(data.long$drop_c > data.long$time & is.na(data.long$c) == TRUE, 2, count_mis_c)
  data <- data.long
  my_subset <- function() {
    drop_e <- drop_c <- NULL
    cov_matrix <- subset(data, select = -c(e, c, drop_e, drop_c, time, t))[data$time == 1, ]
    cov_matrix
  }
  cov_matrix <- my_subset()
  cov_matrix <- cov_matrix[!unlist(vapply(cov_matrix, anyNA, logical(1)))]
  is.formula <- function (x) { inherits(x, "formula") }
  if(is.formula(model.eff) == FALSE | is.formula(model.cost) == FALSE) {
    stop("model.eff and/or model.cost must be formula objects")
  }
  if(is.logical(center) == FALSE) { stop("center must be either TRUE or FALSE") }
  fixed_e <- nobars_(model.eff)
  fixed_c <- nobars_(model.cost)
  random_e <- fb(model.eff)
  random_c <- fb(model.cost)
  fname_re_e_coeff <- as.formula(paste("e", "0", sep=" ~ "))
  fname_re_c_coeff <- as.formula(paste("c", "0", sep=" ~ "))
  fname_re_me_coeff <- as.formula(paste("me", "0", sep=" ~ "))
  fname_re_mc_coeff <- as.formula(paste("mc", "0", sep=" ~ "))
  clusn_e <- clusn_c <- NULL
  clusn_me <- clusn_mc <- NULL
  if(!is.null(random_e) & length(random_c) > 1 | !is.null(random_e) & length(random_c) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(!is.null(random_e)){
    if(!grepl(paste0("\\+\\s+\\("), model.eff[3])) {
      stop("random effects for model.eff are not correctly specified")
    }
  }
  if(!is.null(random_c)){
    if(!grepl(paste0("\\+\\s+\\("), model.cost[3])) {
      stop("random effects for model.eff are not correctly specified")
    }
  }
  if(all(names(model.frame(fixed_e, data = data)) %in% c("e", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_c, data = data)) %in% c("c", "e", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the fixed effects model")
  }
  if(all(names(model.frame(fixed_e, data = data)) %in% names(data)) == FALSE | 
     all(names(model.frame(fixed_c, data = data)) %in% names(data)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
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
  if("t" %in% names(model.frame(fixed_e, data = data)) | "t" %in% names(model.frame(fixed_c, data = data))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.eff' and/or 'model.cost'")
  }
  index_mis_e <- which(is.na(data$e[data$time == 1]))
  index_mis_c <- which(is.na(data$c[data$time == 1]))
  data2 <- data
  data2$me <- data2$drop_e
  data2$mc <- data2$drop_c
  data.t1 <- data[data$time == 1, ]
  data.t1$e[is.na(data.t1$e)] <- -999999
  data.t1$c[is.na(data.t1$c)] <- -999999
  mf_e_fixed <- model.frame(formula = fixed_e, data = data.t1)
  mf_c_fixed <- model.frame(formula = fixed_c, data = data.t1)
  terms <- NULL
  x_e_fixed <- model.matrix(attr(mf_e_fixed, "terms"), data = mf_e_fixed)
  x_c_fixed <- model.matrix(attr(mf_c_fixed, "terms"), data = mf_c_fixed)
  if("e" %in% names(mf_c_fixed)){
    mf_e_fixed$e[index_mis_e] <- NA
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
    if(all(names(model.frame(fname_re_e_coeff, data = data.t1)) %in% c("0", "1", names(model.frame(fixed_e, data = data.t1)))) == FALSE) {
      stop("only covariates defined as fixed effects can be included in the random effects model")
    }
    if("e" %in% labels(terms(fname_re_e_coeff))) {
      stop("please remove 'e' from the random effects expression of model.eff")
    }
    if("c" %in% labels(terms(fname_re_e_coeff))) {
      stop("dependence allowed only through the cost model; please remove 'c' from model.eff")
    }
    mf_e_random <- model.frame(formula = fname_re_e_coeff, data = data.t1)
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
    if(all(names(model.frame(fname_re_c_coeff, data = data.t1)) %in% c("0", "1", names(model.frame(fixed_c, data = data.t1)))) == FALSE) {
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
    mf_c_random <- model.frame(formula = fname_re_c_coeff, data = data.t1)
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
  N1 <- N2 <- c() 
  N1 <- length(data.t1$e[data.t1$t == 1])
  N2 <- length(data.t1$e[data.t1$t == 2])
  N <- c(N1, N2)
  n.time <- length(unique(data$time))
  eff1.mat <- cost1.mat <- matrix(NA, nrow = N1, ncol = n.time)
  eff2.mat <- cost2.mat <- matrix(NA, nrow = N2, ncol = n.time)
  mpatt_e1 <- unique(count_mis_e[data$t == 1]) 
  mpatt_e2 <- unique(count_mis_e[data$t == 2]) 
  mpatt_c1 <- unique(count_mis_c[data$t == 1]) 
  mpatt_c2 <- unique(count_mis_c[data$t == 2]) 
  n_patt_e1 <- max(mpatt_e1)
  n_patt_e2 <- max(mpatt_e2)
  n_patt_c1 <- max(mpatt_c1)
  n_patt_c2 <- max(mpatt_c2)
  if(!any(mpatt_e1 %in% c(2, 3)) | any(mpatt_e1 < 1) | any(mpatt_e1 > 3) | 
     !any(mpatt_e2 %in% c(2, 3)) | any(mpatt_e2 < 1) | any(mpatt_e2 > 3)) { 
    stop("A minimum of two and a maximum of 3 different missingness patterns can be modelled for each outcome and arm")
  }
  if(!any(mpatt_c1 %in% c(2, 3)) | any(mpatt_c1 < 1) | any(mpatt_c1 > 3) | 
     !any(mpatt_c2 %in% c(2, 3)) | any(mpatt_c2 < 1) | any(mpatt_c2 > 3)) { 
    stop("A minimum of two and a maximum of 3 different missingness patterns can be modelled for each outcome and arm")
  }
  if(length(mpatt_e1) == 2){
    if(all(mpatt_e1 == c(1, 2))) {
      n_patt_e1 <- 2
    }
    if(all(mpatt_e1 == c(1, 3))) {
      n_patt_e1 <- 2
      count_mis_e[data$t == 1] <- ifelse(count_mis_e[data$t == 1] == 3, 2, count_mis_e[data$t == 1])
    }
    if(all(mpatt_e1 == c(2, 3))) {
      n_patt_e1 <- 2
      count_mis_e[data$t == 1] <- ifelse(count_mis_e[data$t == 1] == 2, 1, count_mis_e[data$t == 1])
      count_mis_e[data$t == 1] <- ifelse(count_mis_e[data$t == 1] == 3, 2, count_mis_e[data$t == 1])
    }
  }
  if(length(mpatt_e2) == 2){
    if(all(mpatt_e2 == c(1, 2))) {
      n_patt_e2 <- 2
    }
    if(all(mpatt_e2 == c(1, 3))) {
      n_patt_e2 <- 2
      count_mis_e[data$t == 2] <- ifelse(count_mis_e[data$t == 2] == 3, 2, count_mis_e[data$t == 2])
    }
    if(all(mpatt_e2 == c(2, 3))) {
      n_patt_e2 <- 2
      count_mis_e[data$t == 2] <- ifelse(count_mis_e[data$t == 2] == 2, 1, count_mis_e[data$t == 2])
      count_mis_e[data$t == 2] <- ifelse(count_mis_e[data$t == 2] == 3, 2, count_mis_e[data$t == 2])
    }
  }
  if(length(mpatt_c1) == 2){
    if(all(mpatt_c1 == c(1, 2))) {
      n_patt_c1 <- 2
    }
    if(all(mpatt_c1 == c(1, 3))) {
      n_patt_c1 <- 2
      count_mis_c[data$t == 1] <- ifelse(count_mis_c[data$t == 1] == 3, 2, count_mis_c[data$t == 1])
    }
    if(all(mpatt_c1 == c(2, 3))) {
      n_patt_c1 <- 2
      count_mis_c[data$t == 1] <- ifelse(count_mis_c[data$t == 1] == 2, 1, count_mis_c[data$t == 1])
      count_mis_c[data$t == 1] <- ifelse(count_mis_c[data$t == 1] == 3, 2, count_mis_c[data$t == 1])
    }
  }
  if(length(mpatt_c2) == 2){
    if(all(mpatt_c2 == c(1, 2))) {
      n_patt_c2 <- 2
    }
    if(all(mpatt_c2 == c(1, 3))) {
      n_patt_c2 <- 2
      count_mis_c[data$t == 2] <- ifelse(count_mis_c[data$t == 2] == 3, 2, count_mis_c[data$t == 2])
    }
    if(all(mpatt_c2 == c(2, 3))) {
      n_patt_c2 <- 2
      count_mis_c[data$t == 2] <- ifelse(count_mis_c[data$t == 2] == 2, 1, count_mis_c[data$t == 2])
      count_mis_c[data$t == 2] <- ifelse(count_mis_c[data$t == 2] == 3, 2, count_mis_c[data$t == 2])
    }
  }
  count_me1 <- array(NA, dim = c(N1, n.time, n_patt_e1)) 
  for(i in 1:n.time) {
    count_mis_e1_patt1 <- ifelse(count_mis_e[data$time == i & data$t == 1] == 1, 1, 0)
    count_mis_e1_patt2 <- ifelse(count_mis_e[data$time == i & data$t == 1] == 2, 1, 0)
    count_mis_e1_patt <- cbind(count_mis_e1_patt1, count_mis_e1_patt2)
    if(n_patt_e1 == 3){
    count_mis_e1_patt3 <- ifelse(count_mis_e[data$time == i & data$t == 1] == 3, 1, 0)
    count_mis_e1_patt <- cbind(count_mis_e1_patt, count_mis_e1_patt3)
    }
    count_me1[, i, ] <- count_mis_e1_patt
  }
  count_me2 <- array(NA, dim = c(N2, n.time, n_patt_e2)) 
  for(i in 1:n.time) {
    count_mis_e2_patt1 <- ifelse(count_mis_e[data$time == i & data$t == 2] == 1, 1, 0)
    count_mis_e2_patt2 <- ifelse(count_mis_e[data$time == i & data$t == 2] == 2, 1, 0)
    count_mis_e2_patt <- cbind(count_mis_e2_patt1, count_mis_e2_patt2)
    if(n_patt_e2 == 3){
      count_mis_e2_patt3 <- ifelse(count_mis_e[data$time == i & data$t == 2] == 3, 1, 0)
      count_mis_e2_patt <- cbind(count_mis_e2_patt, count_mis_e2_patt3)
    }
    count_me2[, i, ] <- count_mis_e2_patt
  }
  count_mc1 <- array(NA, dim = c(N1, n.time, n_patt_c1)) 
  for(i in 1:n.time) {
    count_mis_c1_patt1 <- ifelse(count_mis_c[data$time == i & data$t == 1] == 1, 1, 0)
    count_mis_c1_patt2 <- ifelse(count_mis_c[data$time == i & data$t == 1] == 2, 1, 0)
    count_mis_c1_patt <- cbind(count_mis_c1_patt1, count_mis_c1_patt2)
    if(n_patt_c1 == 3){
      count_mis_c1_patt3 <- ifelse(count_mis_c[data$time == i & data$t == 1] == 3, 1, 0)
      count_mis_c1_patt <- cbind(count_mis_c1_patt, count_mis_c1_patt3)
    }
    count_mc1[, i, ] <- count_mis_c1_patt
  }
  count_mc2 <- array(NA, dim = c(N2, n.time, n_patt_c2)) 
  for(i in 1:n.time) {
    count_mis_c2_patt1 <- ifelse(count_mis_c[data$time == i & data$t == 2] == 1, 1, 0)
    count_mis_c2_patt2 <- ifelse(count_mis_c[data$time == i & data$t == 2] == 2, 1, 0)
    count_mis_c2_patt <- cbind(count_mis_c2_patt1, count_mis_c2_patt2)
    if(n_patt_c2 == 3){
      count_mis_c2_patt3 <- ifelse(count_mis_c[data$time == i & data$t == 2] == 3, 1, 0)
      count_mis_c2_patt <- cbind(count_mis_c2_patt, count_mis_c2_patt3)
    }
    count_mc2[, i, ] <- count_mis_c2_patt
  }
  m_eff <- count_mis_e
  m_cost <- count_mis_c
  time1.mat <- matrix(NA, nrow = N1, ncol = n.time)
  time2.mat <- matrix(NA, nrow = N2, ncol = n.time)
  for(i in 1:n.time) {
    eff1.mat[, i] <- data$e[data$time == i & data$t == 1] 
    eff2.mat[, i] <- data$e[data$time == i & data$t == 2] 
    cost1.mat[, i] <- data$c[data$time == i & data$t == 1] 
    cost2.mat[, i] <- data$c[data$time == i & data$t == 2] 
    time1.mat[, i] <- data$time[data$time == i & data$t == 1] 
    time2.mat[, i] <- data$time[data$time == i & data$t == 2] 
  }
  effects <- list(eff1.mat, eff2.mat)
  costs <- list(cost1.mat, cost2.mat)
  time <- list(time1.mat, time2.mat)
  m_eff <- list(count_me1, count_me2) 
  m_cost <- list(count_mc1, count_mc2) 
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
  n_patt_e <- c(n_patt_e1, n_patt_e2)
  n_patt_c <- c(n_patt_c1, n_patt_c2)
  eff_cc <- list(eff1_cc, eff2_cc) 
  cost_cc <- list(cost1_cc, cost2_cc)
  cov1_e_fixed <- as.data.frame(x_e_fixed[data.t1$t == 1, ])
  names(cov1_e_fixed) <- colnames(x_e_fixed)
  cov2_e_fixed <- as.data.frame(x_e_fixed[data.t1$t == 2, ])
  names(cov2_e_fixed) <- colnames(x_e_fixed)
  cov_e_fixed <- list(cov1_e_fixed, cov2_e_fixed)
  x_c_hold_fixed <- x_c_fixed
  if("e" %in% colnames(x_c_hold_fixed)) {
    x_c_fixed <- subset(x_c_hold_fixed, select = -c(e))
  }
  cov1_c_fixed <- as.data.frame(x_c_fixed[data.t1$t == 1, ])
  names(cov1_c_fixed) <- colnames(x_c_fixed)
  cov2_c_fixed <- as.data.frame(x_c_fixed[data.t1$t == 2, ])
  names(cov2_c_fixed) <- colnames(x_c_fixed)
  cov_c_fixed <- list(cov1_c_fixed, cov2_c_fixed)
  cov_e_fixed <- list(cov1_e_fixed, cov2_e_fixed) 
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
    cov1_e_random <- as.data.frame(x_e_random[data.t1$t == 1, ])
    names(cov1_e_random) <- colnames(x_e_random)
    cov2_e_random <- as.data.frame(x_e_random[data.t1$t == 2, ])
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
    clusn_e1 <- clusn_e[data.t1$t == 1]
    clusn_e1 <- factor(clusn_e1, levels = unique(clusn_e1))
    clusn_e2 <- clusn_e[data.t1$t == 2]
    clusn_e2 <- factor(clusn_e2, levels = unique(clusn_e2))
  }
  if(!is.null(random_c)){
    x_c_hold_random <- x_c_random
    if("e" %in% colnames(x_c_hold_random)) {
      x_c_random <- subset(x_c_hold_random, select = -c(e))
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
  data2.t1$e[is.na(data2.t1$e) == TRUE] <- -999999
  data2.t1$c[is.na(data2.t1$c) == TRUE] <- -999999
  data2.t1$me <- data$drop_e[data$time == 1]
  data2.t1$mc <- data$drop_c[data$time == 1]
  if(!is.formula(model.me) | !is.formula(model.mc)) {
    stop("model.me and/or model.mc must be formula objects")
  }
  fixed_me <- nobars_(model.me)
  fixed_mc <- nobars_(model.mc)
  random_me <- fb(model.me)
  random_mc <- fb(model.mc)
  if(!is.null(random_me) & length(random_me) > 1 | !is.null(random_mc) & length(random_mc) > 1) {
    stop("random effects can be included in the formula only through a single expression within brackets")
  }
  if(all(names(model.frame(fixed_me, data = data2.t1)) %in% c("me", "e", names(cov_matrix))) == FALSE | 
     all(names(model.frame(fixed_mc, data = data2.t1)) %in% c("mc", "c", names(cov_matrix))) == FALSE) {
    stop("partially-observed covariates cannot be included in the model")
  }
  if(all(names(model.frame(fixed_me, data = data2.t1)) %in% names(data2.t1)) == FALSE | 
     all(names(model.frame(fixed_mc, data = data2.t1)) %in% names(data2.t1)) == FALSE) {
    stop("you must provide names in the formula that correspond to those in the data")
  }
  if(names(model.frame(fixed_me, data = data2.t1)[1]) != "me") {
    stop("you must set 'me' as the response in the formula model.me")
  }
  if(names(model.frame(fixed_mc, data = data2.t1)[1]) != "mc") {
    stop("you must set 'mc' as the response in the formula model.mc")
  }
  if("t" %in% names(model.frame(fixed_mc, data = data2.t1)) | "t" %in% names(model.frame(fixed_me, data = data2.t1))) {
    stop("treatment indicator must be provided only in the data. Please remove 't' from 'model.me' and/or 'model.mc'")
  }
  if("c" %in% names(model.frame(fixed_me, data = data2.t1)) | "e" %in% names(model.frame(fixed_mc, data = data2.t1))) {
    stop("please remove 'e' from model.mc and/or remove 'c' from model.me")
  }
  if("e" %in% labels(terms(fixed_me))) {
    if(length(grep(":e", labels(terms(fixed_me)))) != 0 | length(grep("e:", labels(terms(fixed_me)))) != 0) {
      stop("no interaction effects for 'e' is allowed")
    } 
  }
  if("c" %in% labels(terms(fixed_mc))) {
    if(length(grep(":c", labels(terms(fixed_mc)))) != 0 | length(grep("c:", labels(terms(fixed_mc)))) != 0) {
      stop("no interaction effects for 'c' is allowed")
    } 
  }
  name_re_me_coeff <- NULL
  name_re_mc_coeff <- NULL
  if(!is.null(random_me)){
    if(!grepl(paste0("\\+\\s+\\("), model.me[3])) {
        stop("random effects for model.me are not correctly specified")
    }
    name_re_me_coeff <- sub("\\|.*", "", random_me)
    if(grepl("0 + 1", name_re_me_coeff, fixed = TRUE) == TRUE) { stop("Either remove or add the random intercept")}
    name_clus_me <- sub('.*\\|', '', random_me)
    if(lengths(strsplit(name_clus_me, " ")) > 2) {stop("a single clustering variable must be selected for each formula") }
    name_clus_me <- gsub(" ", "", name_clus_me, fixed = TRUE)
    if(!name_clus_me %in% names(cov_matrix)) { stop("the clustering variable must be among the variables in the dataset") }
    if(strsplit(name_re_me_coeff, "")[[1]][1] == 0) {
      no_random_int_me <- TRUE} else {no_random_int_me <- FALSE }
    if(no_random_int_me == TRUE) { 
      name_re_me_coeff <- sub("[0]", "", name_re_me_coeff) 
      name_re_me_coeff <- sub("[+]", "", name_re_me_coeff) 
      }
    if(name_re_me_coeff == "" | name_re_me_coeff == " ") { stop("please state for which variables the random effects are assumed") }
    if(gsub(" ", "", name_re_me_coeff) == "e" & no_random_int_me == FALSE) {name_re_me_coeff <- "1 + e" }
    fname_re_me_coeff <- as.formula(paste("me", name_re_me_coeff, sep = " ~ "))
    if(all(names(model.frame(fname_re_me_coeff, data = data2)) %in% c("0","1", names(model.frame(fixed_me, data = data2)))) == FALSE) {
      stop("only covariates inlcued as fixed effects can be included in the random effects model")
    }
    if("me" %in% labels(terms(fname_re_me_coeff))) {
      stop("please remove 'me' from the right hand side of model.me")
    }
    if("e" %in% labels(terms(fname_re_me_coeff))) {
      if(length(grep(":e", labels(terms(fname_re_me_coeff)))) != 0 | length(grep("e:", labels(terms(fname_re_me_coeff)))) != 0) {
        stop("no interaction effects for 'e' are allowed")
      } 
    }
    clus_me <- data[, name_clus_me]
    if(!is.factor(clus_me)) { stop("clustering variables must be defined as factors") }
    clusn_me <- as.numeric(clus_me)
    if(!all(diff(sort(unique(clusn_me))) == 1) | !min(clusn_me) == 1) {
      stop("ordered levels of clustering variables must not have gaps and must start from 1")
    }
  }
  if(!is.null(random_mc)){
    if(!grepl(paste0("\\+\\s+\\("), model.mc[3])) {
      stop("random effects for model.mc are not correctly specified")
    }
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
  mf_me_fixed <- model.frame(formula = fixed_me, data = data2.t1)
  mf_mc_fixed <- model.frame(formula = fixed_mc, data = data2.t1)
  z_e_fixed <- model.matrix(attr(mf_me_fixed, "terms"), data = mf_me_fixed)
  z_c_fixed <- model.matrix(attr(mf_mc_fixed, "terms"), data = mf_mc_fixed)
  z_e_hold_fixed <- z_e_fixed
  if("e" %in% colnames(z_e_hold_fixed)) {
    z_e_fixed <- subset(z_e_hold_fixed, select = -c(e))
  }
  z_c_hold_fixed <- z_c_fixed
  if("c" %in% colnames(z_c_hold_fixed)) {
    z_c_fixed <- subset(z_c_hold_fixed, select = -c(c))
  }
  covz1_e_fixed <- as.data.frame(z_e_fixed[data2.t1$t == 1, ])
  names(covz1_e_fixed) <- colnames(z_e_fixed)
  covz2_e_fixed <- as.data.frame(z_e_fixed[data2.t1$t == 2, ])
  names(covz2_e_fixed) <- colnames(z_e_fixed)
  covz_e_fixed <- list(covz1_e_fixed, covz2_e_fixed)
  covze_fixed <- list(covz1_e_fixed, covz2_e_fixed) 
  mean_covz_e_fixed <- list(apply(as.matrix(covz1_e_fixed), 2, mean), apply(as.matrix(covz2_e_fixed), 2, mean))
  names(covz_e_fixed) <- names(mean_covz_e_fixed) <- c("Control", "Intervention")
  covz1_c_fixed <- as.data.frame(z_c_fixed[data2.t1$t == 1, ])
  names(covz1_c_fixed) <- colnames(z_c_fixed)
  covz2_c_fixed <- as.data.frame(z_c_fixed[data2.t1$t == 2, ])
  names(covz2_c_fixed) <- colnames(z_c_fixed)
  covz_c_fixed <- list(covz1_c_fixed, covz2_c_fixed)
  covzc_fixed <- list(covz1_c_fixed, covz2_c_fixed) 
  mean_covz_c_fixed <- list(apply(as.matrix(covz1_c_fixed), 2, mean), apply(as.matrix(covz2_c_fixed), 2, mean))
  covz1_e_center_fixed <- as.data.frame(scale(covz1_e_fixed, scale = FALSE))
  covz2_e_center_fixed <- as.data.frame(scale(covz2_e_fixed, scale = FALSE))
  covz1_e_center_fixed[, 1] <- rep(1, nrow(covz1_e_fixed))
  covz2_e_center_fixed[, 1] <- rep(1, nrow(covz2_e_fixed))
  covz_e_center_fixed <- list(covz1_e_center_fixed, covz2_e_center_fixed)
  mean_covz_e_center_fixed <- list(apply(as.matrix(covz1_e_center_fixed), 2, mean), apply(as.matrix(covz2_e_center_fixed), 2, mean))
  covz1_c_center_fixed <- as.data.frame(scale(covz1_c_fixed, scale = FALSE))
  covz2_c_center_fixed <- as.data.frame(scale(covz2_c_fixed, scale = FALSE))
  covz1_c_center_fixed[, 1] <- rep(1, nrow(covz1_c_fixed))
  covz2_c_center_fixed[, 1] <- rep(1, nrow(covz2_c_fixed))
  covz_c_center_fixed <- list(covz1_c_center_fixed, covz2_c_center_fixed)
  mean_covz_c_center_fixed <- list(apply(as.matrix(covz1_c_center_fixed), 2, mean), apply(as.matrix(covz2_c_center_fixed), 2, mean))
  if(center == TRUE) {
    covz_e_fixed <- covz_e_center_fixed
    covz_c_fixed <- covz_c_center_fixed
    mean_covz_e_fixed <- mean_covz_e_center_fixed
    mean_covz_c_fixed <- mean_covz_c_center_fixed
  }
  if(!is.null(random_me)){
    mf_me_random <- model.frame(formula = fname_re_me_coeff, data = data2.t1)
    z_e_random <- model.matrix(attr(mf_me_random, "terms"), data = mf_me_random)
    if(no_random_int_me == TRUE) {
      z_e_random <- as.matrix(z_e_random[, !colnames(z_e_random) == "(Intercept)"])
      if(dim(z_e_random)[2] == 1) { colnames(z_e_random) <- colnames(model.matrix(attr(mf_me_random, "terms"), data = mf_me_random))[2] }
    }
    z_e_hold_random <- z_e_random
    if("e" %in% colnames(z_e_hold_random)) {
      z_e_random <- subset(z_e_hold_random, select = -c(e))
    }
    covz1_e_random <- as.data.frame(z_e_random[data2.t1$t == 1, ])
    names(covz1_e_random) <- colnames(z_e_random)
    covz2_e_random <- as.data.frame(z_e_random[data2.t1$t == 2, ])
    names(covz2_e_random) <- colnames(z_e_random)
    covz_e_random <- list(covz1_e_random, covz2_e_random)
    covze_random <- list(covz1_e_random, covz2_e_random) 
    mean_covz_e_random <- list(apply(as.matrix(covz1_e_random), 2, mean), apply(as.matrix(covz2_e_random), 2, mean))
    names(covz_e_random) <- names(mean_covz_e_random) <- c("Control", "Intervention")
    covz1_e_center_random <- as.data.frame(scale(covz1_e_random, scale = FALSE))
    covz2_e_center_random <- as.data.frame(scale(covz2_e_random, scale = FALSE))
    if(no_random_int_me == FALSE) {
      covz1_e_center_random[, 1] <- rep(1, nrow(covz1_e_random))
      covz2_e_center_random[, 1] <- rep(1, nrow(covz2_e_random))
    }
    covz_e_center_random <- list(covz1_e_center_random, covz2_e_center_random)
    mean_covz_e_center_random <- list(apply(as.matrix(covz1_e_center_random), 2, mean), apply(as.matrix(covz2_e_center_random), 2, mean))
    if(center == TRUE) {
      covz_e_random <- covz_e_center_random
      mean_covz_e_random <- mean_covz_e_center_random
    }
    clusn_me1 <- clusn_me[data2.t1$t == 1]
    clusn_me1 <- factor(clusn_me1, levels = unique(clusn_me1))
    clusn_me2 <- clusn_me[data2.t1$t == 2]
    clusn_me2 <- factor(clusn_me2, levels = unique(clusn_me2))
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
  names(covz_e_fixed) <- names(mean_covz_e_fixed) <- names(covz_c_fixed) <- names(mean_covz_c_fixed) <- c("Control", "Intervention")
  names(cov_e_fixed) <- names(cov_c_fixed) <- names(mean_cov_e_fixed) <- names(mean_cov_c_fixed) <- c("Control", "Intervention")
  if(!is.null(random_c)) {
    names(cov_c_random) <- names(mean_cov_c_random) <- c("Control", "Intervention")
    clusn_c <- list("Control" = clusn_c1, "Intervention" = clusn_c2)
  } else {cov_c_random <- mean_cov_c_random <- NULL}
  if(!is.null(random_mc)) {
    names(covz_c_random) <- names(mean_covz_c_random) <- c("Control", "Intervention")
    clusn_mc <- list("Control" = clusn_mc1, "Intervention" = clusn_mc2)
  } else {covz_c_random <- mean_covz_c_random <- NULL}
  if(!is.null(random_e)) {
    names(cov_e_random) <- names(mean_cov_e_random) <- c("Control", "Intervention")
    clusn_e <- list("Control" = clusn_e1, "Intervention" = clusn_e2)
  } else {cov_e_random <- mean_cov_e_random <- NULL}
  if(!is.null(random_me)) {
    names(covz_e_random) <- names(mean_covz_e_random) <- c("Control", "Intervention")
    clusn_me <- list("Control" = clusn_me1, "Intervention" = clusn_me2)
  } else {covz_e_random <- mean_covz_e_random <- NULL}
  names(m_eff) <- names(m_cost) <- names(time) <- c("Control", "Intervention")
  names(effects) <- names(costs) <- names(eff_cc) <- names(cost_cc) <- c("Control", "Intervention")
  data_raw <- list("raw_effects" = effects, "raw_costs" = costs, "raw_effects_cc" = eff_cc, "raw_costs_cc" = cost_cc, "arm_lengths" = N, 
                   "arm_lengths_cc" = N_cc, "arm_missing_data" = N_mis, "arm_missing_effects_pattern" = n_patt_e, "arm_missing_costs_pattern" = n_patt_c,
                   "missing_effects" = m_eff, "missing_costs" = m_cost, "time" = time, "covariates_effects_fixed" = cov_e_fixed, "covariates_costs_fixed" = cov_c_fixed, 
                   "mean_cov_effects_fixed" = mean_cov_e_fixed, "mean_cov_costs_fixed" = mean_cov_c_fixed, "covariates_missing_effects_fixed" = covz_e_fixed, 
                   "mean_cov_missing_effects_fixed" = mean_covz_e_fixed, "covariates_missing_costs_fixed" = covz_c_fixed, 
                   "mean_cov_missing_costs_fixed" = mean_covz_c_fixed, "covariates_effects_random" = cov_e_random, "covariates_costs_random" = cov_c_random, 
                   "mean_cov_effects_random" = mean_cov_e_random, "mean_cov_costs_random" = mean_cov_c_random, 
                   "covariates_missing_effects_random" = covz_e_random, "mean_cov_missing_effects_random" = mean_covz_e_random, "covariates_missing_costs_random" = covz_c_random, 
                   "mean_cov_missing_costs_random" = mean_covz_c_random, "clus_e" = clusn_e, "clus_c" = clusn_c, "clus_me" = clusn_me, "clus_mc" = clusn_mc, "data_ind" = data2.t1, "data_long" = data) 
  model_formula <- list("mf_model.e_fixed" = fixed_e, "mf_model.c_fixed" = fixed_c, "mf_model.me_fixed" = fixed_me, "mf_model.mc_fixed" = fixed_mc,
                        "mf_model.e_random" = fname_re_e_coeff, "mf_model.c_random" = fname_re_c_coeff, "mf_model.me_random" = fname_re_me_coeff, "mf_model.mc_random" = fname_re_mc_coeff)
  data_list <- list("data_raw" = data_raw, "model_formula" = model_formula)
  return(data_list) 
}