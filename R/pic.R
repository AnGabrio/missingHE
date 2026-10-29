#' Predictive information criteria for a Bayesian model
#'
#' Efficient approximate leave-one-out cross validation (`LOO`), deviance information criterion (`DIC`) and widely applicable information criterion (`WAIC`) for Bayesian models,
#' calculated on the observed data.
#' @keywords loo waic dic JAGS
#' @param x A `missingHE` object containing the results of a Bayesian model fitted in cost-effectiveness analysis using the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#' @param criterion type of information criteria to be produced. Available choices are `dic` for the Deviance Information Criterion,
#' `waic` for the Widely Applicable Information Criterion, and `looic` for the Leave-One-Out Information Criterion.
#' @param cases group of cases for which information criteria should be computed for: either `all` for all cases, `cc` for complete cases,
#' `ac_e` and `ac_c` for only the observed effect and cost cases, respectively.
#' @param \dots Additional parameters that can be provided to manage the output of `pic` when `looic` is selected. For more details see [bayesplot].
#' @return A named list containing different predictive information criteria results and quantities according to the value of `criterion`. In all cases, the measures are
#' computed on the observed data for the specific value selected in `cases`. The elements contained include:
#'
#'
#' - **d_bar**: Posterior mean deviance (only if `criterion = "dic"`).
#' - **pD**: Effective number of parameters calculated with the formula used by `JAGS` (only if `criterion = "dic"`).
#' - **dic**: Deviance Information Criterion calculated with the formula used by `JAGS` (only if `criterion = "dic"`).
#' - **d_hat**: Deviance evaluated at the posterior mean of the parameters and calculated with the formula used by `JAGS` (only if `criterion = "dic"`).
#' - **elpd**: Expected log pointwise predictive density calculated on the selected cases (only if `criterion = "waic"` or `criterion = "looic"`).
#' - **elpd_se**: `elpd` standard error calculated on the selected cases (only if `criterion = "waic"` or `criterion = "looic"`).
#' - **p**: Effective number of parameters calculated on the selected cases (only if `criterion = "waic"` or `criterion = "looic"`).
#' - **p_se**: `p` standard error calculated on the selected cases (only if `criterion = "waic"` or `criterion = "looic"`).
#' - **looic**: The leave-one-out information criterion calculated on the selected cases (only if `criterion = "looic"`).
#' - **looic_se**: The `looic` standard error calculated on the selected cases (only if `criterion = "looic"`).
#' - **waic**: The widely applicable information criterion calculated on the selected cases (only if `criterion = "waic"`).
#' - **waic_se**: The `waic` standard error calculated on the selected cases (only if `criterion = "waic"`).
#' - **pointwise**: A matrix containing the pointwise contributions of each of the above measures calculated on the selected cases (only if `criterion = "waic"` or `criterion = "looic"`).
#' - **diagnostics**: A named list containing additional diagnostic measures (only if `criterion = "looic"`). See [loo::loo()] for details about interpreting the list elements.
#' - **psis_object**: A named list containing the matrix of (smoothed) log weights (only if `criterion = "looic"` and `save_psis = TRUE`).
#'
#' @seealso [R2jags::jags()], [loo::loo()], [loo::waic()]
#' @author Andrea Gabrio
#' @importFrom stats var
#' @details The Deviance Information Criterion (`DIC`), Leave-One-Out Information Criterion (`LOOIC`) and the Widely Applicable Information Criterion (`WAIC`) are methods for estimating
#' out-of-sample predictive accuracy from a Bayesian model using the log-likelihood evaluated at the posterior simulations of the parameters.
#' `DIC` is computationally simple to calculate but it is known to have some problems, arising in part from it not being fully Bayesian in that it is based on a point estimate.
#' `LOOIC` can be computationally expensive but can be easily approximated using importance weights that are smoothed by fitting a generalised Pareto distribution to the upper tail
#' of the distribution of the importance weights.
#' `WAIC` is fully Bayesian and closely approximates Bayesian cross-validation. Unlike `DIC`, `WAIC` is invariant to parameterisation and also works for singular models.
#' In finite cases, `WAIC` and `LOOIC` give similar estimates, but for influential observations `WAIC` underestimates the effect of leaving out one observation.
#'
#' @references
#' Plummer, M. *JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling*. (2003).
#'
#' Vehtari, A. Gelman, A. Gabry, J. (2016a) *Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC*.
#' Statistics and Computing. Advance online publication.
#'
#' Vehtari, A. Gelman, A. Gabry, J. (2016b) *Pareto smoothed importance sampling*. ArXiv preprint.
#'
#' Gelman, A. Hwang, J. Vehtari, A. (2014) *Understanding predictive information criteria for Bayesian models*.
#' Statistics and Computing 24, 997-1016.
#'
#' Watanable, S. (2010). *Asymptotic equivalence of Bayes cross validation and widely application information criterion in singular learning theory*.
#' Journal of Machine Learning Research 11, 3571-3594.
#'
#' @export
#' @examples
#' # For examples see the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#' #
#' #

pic <- function(x, criterion = "dic", cases = "cc", ...) {
  if(!isTRUE(requireNamespace("loo", quietly = TRUE))) {
    stop("You need to install the R package 'loo'. Please run in your R terminal:\n install.packages('loo')")}
  if(!is.character(criterion)) {
    stop("Criterion must be a character name")}
  criterion <- tolower(criterion)
  if(!criterion %in% c("dic", "waic", "looic")) {
    stop("Please, select a character name for 'criterion' among those available. For details type 'help(pic)'")}
  if(!cases %in% c("all", "cc", "ac_e", "ac_c")) {
    stop("Please, select a character name for 'cases' among those available. For details type 'help(pic)'")}
  if(!is.list(x)) { stop("Only an object or list of objects of class 'missingHE' can be used")}
  if(!inherits(x, "missingHE")) {
    if(!all(unlist(lapply(x, function(m) { inherits(m, "missingHE")})))) {
      stop("Only a list of objects of class 'missingHE' can be used with same 'method' and 'data_format' argument values")}
      if(length(unique(lapply(x, function(x) x$model_output$method))) != 1) {
        stop("Only a list of objects of class 'missingHE' can be used with same 'method' and 'data_format' argument values")}
    if(length(unique(lapply(x, function(x) x$data_format))) != 1) {
      stop("Only a list of objects of class 'missingHE' can be used with same 'method' and 'data_format' argument values")}
  }
  exArgs <- list(...)
  if(exists("r_eff", where = exArgs)) { r_eff = exArgs$r_eff } else { r_eff = 1}
  if(exists("save_psis", where = exArgs)) { save_psis = exArgs$save_psis } else { save_psis = FALSE}
  if(exists("cores", where = exArgs)) { cores = exArgs$cores } else { cores = getOption("mc.cores", 1)}
  if(exists("is_method", where = exArgs)) { is_method = exArgs$is_method } else { is_method = c("psis", "tis", "sis")}
  if(inherits(x, "missingHE")) {
    if(x$model_output$model$BUGSoutput$n.chains == 1){
      stop("Information criteria cannot be computed if 'n.chain' = 1")}
    if(x$data_format == "wide") {
      loglik_e_all <- x$model_output$loglik$effects
      loglik_c_all <- x$model_output$loglik$costs
      }
    if(x$data_format == "long") {
      loglik_e_all <- as.matrix(as.data.frame(x$model_output$loglik$effects))
      loglik_c_all <- as.matrix(as.data.frame(x$model_output$loglik$costs))
      }
      loglik_ec_all <- loglik_e_all + loglik_c_all
      m_e <- x$data_set$data_raw$me
      m_c <- x$data_set$data_raw$mc
      loglik_e_ac <- loglik_e_all[, m_e == 0]
      loglik_c_ac <- loglik_c_all[, m_c == 0]
      loglik_e_cc <- loglik_e_all[, m_e == 0 & m_c == 0]
      loglik_c_cc <- loglik_c_all[, m_e == 0 & m_c == 0]
      loglik_ec_cc <- loglik_e_cc + loglik_c_cc
      if(x$model_output$method %in% c("SELECTION", "LMDM")) {
        if(x$model_output$method == "SELECTION") {
          loglik_me_all <- x$model_output$loglik$`missing indicators effects`
          loglik_mc_all <- x$model_output$loglik$`missing indicators costs`
        }
        if(x$model_output$method == "LMDM") {
          loglik_me_all <- as.matrix(as.data.frame(x$model_output$loglik$`missing indicators effects`))
          loglik_mc_all <- as.matrix(as.data.frame(x$model_output$loglik$`missing indicators costs`))
        }
      loglik_me_ac <- loglik_me_all[, m_e == 0]
      loglik_mc_ac <- loglik_mc_all[, m_c == 0]
      loglik_me_cc <- loglik_me_all[, m_e == 0 & m_c == 0]
      loglik_mc_cc <- loglik_mc_all[, m_e == 0 & m_c == 0]
      if(cases == "all") { loglik_total <- loglik_ec_all + loglik_me_all + loglik_mc_all}
      if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_me_ac}
      if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_mc_ac}
      if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_me_cc + loglik_mc_cc}
      }
      if(x$model_output$method == "HURDLE") {
        if(!is.null(x$model_output$mean_nostr_effects)) {
        loglik_se_all <- x$model_output$loglik$`structural indicators effects`
        loglik_se_ac <- loglik_se_all[, m_e == 0]
        loglik_se_cc <- loglik_se_all[, m_e == 0 & m_c == 0]
        } else { loglik_se_all <- loglik_se_ac <- loglik_se_cc <- 0}
        if(!is.null(x$model_output$mean_nostr_costs)) {
          loglik_sc_all <- x$model_output$loglik$`structural indicators costs`
          loglik_sc_ac <- loglik_sc_all[, m_c == 0]
          loglik_sc_cc <- loglik_sc_all[, m_e == 0 & m_c == 0]
        } else { loglik_sc_all <- loglik_sc_ac <- loglik_sc_cc <- 0}
        if(cases == "all") { loglik_total <- loglik_ec_all + loglik_se_all + loglik_sc_all}
        if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_se_ac}
        if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_sc_ac}
        if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_se_cc + loglik_sc_cc}
      }
      if(x$model_output$method == "PATTERN") {
        loglik_d_all <- x$model_output$loglik$`pattern indicators`
        loglik_d_e_ac <- loglik_d_all[, m_e == 0]
        loglik_d_c_ac <- loglik_d_all[, m_c == 0]
        loglik_d_cc <- loglik_d_all[, m_e == 0 & m_c == 0]
        if(cases == "all") { loglik_total <- loglik_ec_all + loglik_d_all}
        if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_d_e_ac}
        if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_d_c_ac}
        if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_d_cc}
      }
      if(criterion == "dic") {
        d_bar <- mean(-2 * rowSums(loglik_total))
        pD <- var(-2* c(loglik_total)) / 2
        dic <- d_bar + pD
        d_hat <- dic - 2 * pD
        ic <- list("d_bar" = d_bar, "pD" = pD, "dic" = dic, "d_hat" = d_hat)
      }
      if(criterion == "waic") {
        waic_l <- suppressWarnings(loo::waic(loglik_total))
        elpd <- waic_l[[1]]["elpd_waic", "Estimate"]
        elpd_se <- waic_l[[1]]["elpd_waic", "SE"]
        p <- waic_l[[1]]["p_waic", "Estimate"]
        p_se <- waic_l[[1]]["p_waic", "SE"]
        waic <- waic_l[[1]]["waic", "Estimate"]
        waic_se <- waic_l[[1]]["waic", "SE"]
        pointwise <- waic_l[[2]]
        ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se,
                   "waic" = waic, "waic_se" = waic_se, "pointwise" = pointwise)
      }
      if(criterion == "looic") {
        loo_l <- suppressWarnings(loo::loo(loglik_total, r_eff = r_eff,
                                           save_psis = save_psis, cores = cores,
                                           is_method = is_method))
        elpd <- loo_l[[1]]["elpd_loo", "Estimate"]
        elpd_se <- loo_l[[1]]["elpd_loo", "Estimate"]
        p <- loo_l[[1]]["p_loo", "Estimate"]
        p_se <- loo_l[[1]]["p_loo", "Estimate"]
        looic <- loo_l[[1]]["looic", "Estimate"]
        looic_se <- loo_l[[1]]["looic", "Estimate"]
        pointwise <- loo_l[[2]]
        diagnostics <- loo_l[[3]]
        psis_object <- loo_l[[4]]
        ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se,
                   "looic" = looic, "looic_se" = looic_se, "pointwise" = pointwise,
                   "diagnostics" = diagnostics, "psis_object" = psis_object)
      }
  }
  if(!inherits(x, "missingHE")) {
    if(any(lapply(x, function(x) x$model_output$model$BUGSoutput$n.chains) == 1)){
      stop("Information criteria cannot be computed if 'n.chain' = 1")}
    if(all(lapply(x, function(x) x$data_format) %in% c("wide", "long"))) {
        ic_list <- lapply(x, function(x, criterion, cases) {
          if(x$data_format == "wide") {
            loglik_e_all <- x$model_output$loglik$effects
            loglik_c_all <- x$model_output$loglik$costs
          }
          if(x$data_format == "long") {
            loglik_e_all <- as.matrix(as.data.frame(x$model_output$loglik$effects))
            loglik_c_all <- as.matrix(as.data.frame(x$model_output$loglik$costs))
          }
          loglik_ec_all <- loglik_e_all + loglik_c_all
          m_e <- x$data_set$data_raw$me
          m_c <- x$data_set$data_raw$mc
          loglik_e_ac <- loglik_e_all[, m_e == 0]
          loglik_c_ac <- loglik_c_all[, m_c == 0]
          loglik_e_cc <- loglik_e_all[, m_e == 0 & m_c == 0]
          loglik_c_cc <- loglik_c_all[, m_e == 0 & m_c == 0]
          loglik_ec_cc <- loglik_e_cc + loglik_c_cc
          if(x$model_output$method %in% c("SELECTION", "LMDM")) {
            if(x$model_output$method == "SELECTION") {
              loglik_me_all <- x$model_output$loglik$`missing indicators effects`
              loglik_mc_all <- x$model_output$loglik$`missing indicators costs`
            }
            if(x$model_output$method == "LMDM") {
              loglik_me_all <- as.matrix(as.data.frame(x$model_output$loglik$`missing indicators effects`))
              loglik_mc_all <- as.matrix(as.data.frame(x$model_output$loglik$`missing indicators costs`))
            }
          loglik_me_ac <- loglik_me_all[, m_e == 0]
          loglik_mc_ac <- loglik_mc_all[, m_c == 0]
          loglik_me_cc <- loglik_me_all[, m_e == 0 & m_c == 0]
          loglik_mc_cc <- loglik_mc_all[, m_e == 0 & m_c == 0]
          if(cases == "all") { loglik_total <- loglik_ec_all + loglik_me_all + loglik_mc_all}
          if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_me_ac}
          if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_mc_ac}
          if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_me_cc + loglik_mc_cc}
          }
          if(x$model_output$method == "HURDLE") {
            if(!is.null(x$model_output$mean_nostr_effects)) {
              loglik_se_all <- x$model_output$loglik$`structural indicators effects`
              loglik_se_ac <- loglik_se_all[, m_e == 0]
              loglik_se_cc <- loglik_se_all[, m_e == 0 & m_c == 0]
            } else { loglik_se_all <- loglik_se_ac <- loglik_se_cc <- 0}
            if(!is.null(x$model_output$mean_nostr_costs)) {
              loglik_sc_all <- x$model_output$loglik$`structural indicators costs`
              loglik_sc_ac <- loglik_sc_all[, m_c == 0]
              loglik_sc_cc <- loglik_sc_all[, m_e == 0 & m_c == 0]
            } else { loglik_sc_all <- loglik_sc_ac <- loglik_sc_cc <- 0}
            if(cases == "all") { loglik_total <- loglik_ec_all + loglik_se_all + loglik_sc_all}
            if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_se_ac}
            if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_sc_ac}
            if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_se_cc + loglik_sc_cc}
          }
          if(x$model_output$method == "PATTERN") {
            loglik_d_all <- x$model_output$loglik$`pattern indicators`
            loglik_d_e_ac <- loglik_d_all[, m_e == 0]
            loglik_d_c_ac <- loglik_d_all[, m_c == 0]
            loglik_d_cc <- loglik_d_all[, m_e == 0 & m_c == 0]
            if(cases == "all") { loglik_total <- loglik_ec_all + loglik_d_all}
            if(cases == "ac_e") { loglik_total <- loglik_e_ac + loglik_d_e_ac}
            if(cases == "ac_c") { loglik_total <- loglik_c_ac + loglik_d_c_ac}
            if(cases == "cc") { loglik_total <- loglik_ec_cc + loglik_d_cc}
          }
          if(criterion == "dic") {
            d_bar <- mean(-2 * rowSums(loglik_total))
            pD <- var(-2* c(loglik_total)) / 2
            dic <- d_bar + pD
            d_hat <- dic - 2 * pD
            ic <- list("d_bar" = d_bar, "pD" = pD, "dic" = dic, "d_hat" = d_hat)
          }
          if(criterion == "waic") {
            waic_l <- suppressWarnings(loo::waic(loglik_total))
            elpd <- waic_l[[1]]["elpd_waic", "Estimate"]
            elpd_se <- waic_l[[1]]["elpd_waic", "SE"]
            p <- waic_l[[1]]["p_waic", "Estimate"]
            p_se <- waic_l[[1]]["p_waic", "SE"]
            waic <- waic_l[[1]]["waic", "Estimate"]
            waic_se <- waic_l[[1]]["waic", "SE"]
            pointwise <- waic_l[[2]]
            ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se,
                       "waic" = waic, "waic_se" = waic_se, "pointwise" = pointwise)
          }
          if(criterion == "looic") {
            loo_l <- suppressWarnings(loo::loo(loglik_total, r_eff = r_eff,
                                               save_psis = save_psis, cores = cores,
                                               is_method = is_method))
            elpd <- loo_l[[1]]["elpd_loo", "Estimate"]
            elpd_se <- loo_l[[1]]["elpd_loo", "Estimate"]
            p <- loo_l[[1]]["p_loo", "Estimate"]
            p_se <- loo_l[[1]]["p_loo", "Estimate"]
            looic <- loo_l[[1]]["looic", "Estimate"]
            looic_se <- loo_l[[1]]["looic", "Estimate"]
            pointwise <- loo_l[[2]]
            diagnostics <- loo_l[[3]]
            psis_object <- loo_l[[4]]
            ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se,
                       "looic" = looic, "looic_se" = looic_se, "pointwise" = pointwise,
                       "diagnostics" = diagnostics, "psis_object" = psis_object)
          }
        return(ic)
      }, criterion = criterion, cases = cases)
        if(criterion == "dic") {
         pic_compare <-  unlist(lapply(ic_list, function(x) x$dic))
         pic_se_compare <-  NULL
         peff_compare <-  unlist(lapply(ic_list, function(x) x$pD))
         peff_se_compare <-  NULL
         names(pic_compare) <- paste("DIC model", seq(1:length(x)))
         names(peff_compare) <- paste("peff model", seq(1:length(x)))
        }
        if(criterion == "waic") {
          pic_compare <-  unlist(lapply(ic_list, function(x) x$waic))
          pic_se_compare <-  unlist(lapply(ic_list, function(x) x$waic_se))
          peff_compare <-  unlist(lapply(ic_list, function(x) x$p))
          peff_se_compare <-  unlist(lapply(ic_list, function(x) x$p_se))
          names(pic_compare) <- paste("WAIC model", seq(1:length(x)))
          names(pic_se_compare) <- paste("SE(WAIC) model", seq(1:length(x)))
          names(peff_compare) <- paste("peff model", seq(1:length(x)))
          names(peff_se_compare) <- paste("SE(peff) model", seq(1:length(x)))
        }
        if(criterion == "looic") {
          pic_compare <-  unlist(lapply(ic_list, function(x) x$looic))
          pic_se_compare <-  unlist(lapply(ic_list, function(x) x$looic_se))
          peff_compare <-  unlist(lapply(ic_list, function(x) x$p))
          peff_se_compare <-  unlist(lapply(ic_list, function(x) x$p_se))
          names(pic_compare) <- paste("LOOIC model", seq(1:length(x)))
          names(pic_se_compare) <- paste("SE(LOOIC) model", seq(1:length(x)))
          names(peff_compare) <- paste("peff model", seq(1:length(x)))
          names(peff_se_compare) <- paste("SE(peff) model", seq(1:length(x)))
        }
    }
    ic <- list("pic" = pic_compare, "pic_se" = pic_se_compare,
               "peff" = peff_compare, "peff_se" = peff_se_compare, "ic_res" = ic_list)
  }
  return(ic)
}
