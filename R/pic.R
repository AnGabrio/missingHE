#' Predictive information criteria for Bayesian models fitted in \code{JAGS} using the function \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}
#' 
#' Efficient approximate leave-one-out cross validation (LOO), deviance information criterion (DIC) and widely applicable information criterion (WAIC) for Bayesian models, 
#' calculated on the observed data.
#' @keywords loo waic dic JAGS   
#' @param x A \code{missingHE} object containing the results of a Bayesian model fitted in cost-effectiveness analysis using the function \code{\link{selection}}, \code{\link{pattern}} or \code{\link{hurdle}}.
#' @param criterion type of information criteria to be produced. Available choices are \code{'dic'} for the Deviance Information Criterion, 
#' \code{'waic'} for the Widely Applicable Information Criterion, and \code{'looic'} for the Leave-One-Out Information Criterion.
#' @param module The modules with respect to which the information criteria should be computed. Available choices are \code{'total'} for the whole model, 
#' \code{'e'} for the effectiveness variables only, \code{'c'} for the cost variables only, and \code{'both'} for both outcome variables.
#' @return A named list containing different predictive information criteria results and quantities according to the value of \code{criterion}. In all cases, the measures are 
#' computed on the observed data for the specific modules of the model selected in \code{module}.
#' \describe{
#'   \item{d_bar}{Posterior mean deviance (only if \code{criterion} is \code{'dic'}).}
#'   \item{pD}{Effective number of parameters calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}.
#'   \item{dic}{Deviance Information Criterion calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}. 
#'   \item{d_hat}{Deviance evaluated at the posterior mean of the parameters and calculated with the formula used by \code{JAGS} (only if \code{criterion} is \code{'dic'})}
#'   \item{elpd, elpd_se}{Expected log pointwise predictive density and standard error calculated on the observed data for the model nodes indicated in \code{module}
#'    (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
#'   \item{p, p_se}{Effective number of parameters and standard error calculated on the observed data for the model nodes indicated in \code{module}
#'    (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
#'   \item{looic, looic_se}{The leave-one-out information criterion and standard error calculated on the observed data for the model nodes indicated in \code{module}
#'    (only if \code{criterion} is \code{'loo'}).}
#'   \item{waic, waic_se}{The widely applicable information criterion and standard error calculated on the observed data for the model nodes indicated in \code{module}
#'    (only if \code{criterion} is \code{'waic'}).}
#'   \item{pointwise}{A matrix containing the pointwise contributions of each of the above measures calculated on the observed data for the model nodes indicated in \code{module}
#'    (only if \code{criterion} is \code{'waic'} or \code{'loo'}).}
#'   \item{pareto_k}{A vector containing the estimates of the shape parameter \eqn{k} for the generalised Pareto fit to the importance ratios for each leave-one-out distribution 
#'    calculated on the observed data for the model nodes indicated in \code{module} (only if \code{criterion} is \code{'loo'}). 
#'    See \code{\link[loo]{loo}} for details about interpreting \eqn{k}.}
#' }
#' @seealso \code{\link[R2jags]{jags}}, \code{\link[loo]{loo}}, \code{\link[loo]{waic}}
#' @author Andrea Gabrio
#' @importFrom stats var
#' @details The Deviance Information Criterion (DIC), Leave-One-Out Information Criterion (LOOIC) and the Widely Applicable Information Criterion (WAIC) are methods for estimating 
#' out-of-sample predictive accuracy from a Bayesian model using the log-likelihood evaluated at the posterior simulations of the parameters. 
#' DIC is computationally simple to calculate but it is known to have some problems, arising in part from it not being fully Bayesian in that it is based on a point estimate.
#' LOOIC can be computationally expensive but can be easily approximated using importance weights that are smoothed by fitting a generalised Pareto distribution to the upper tail 
#' of the distribution of the importance weights. For more details about the methods used to compute LOOIC see the PSIS-LOO section in \code{\link{loo-package}}.
#' WAIC is fully Bayesian and closely approximates Bayesian cross-validation. Unlike DIC, WAIC is invariant to parameterisation and also works for singular models. 
#' In finite cases, WAIC and LOO give similar estimates, but for influential observations WAIC underestimates the effect of leaving out one observation.
#' 
#' @references  
#' Plummer, M. \emph{JAGS: A program for analysis of Bayesian graphical models using Gibbs sampling.} (2003).
#' 
#' Vehtari, A. Gelman, A. Gabry, J. (2016a) Practical Bayesian model evaluation using leave-one-out cross-validation and WAIC. 
#' \emph{Statistics and Computing}. Advance online publication.
#' 
#' Vehtari, A. Gelman, A. Gabry, J. (2016b) Pareto smoothed importance sampling. \emph{ArXiv} preprint.
#' 
#' Gelman, A. Hwang, J. Vehtari, A. (2014) Understanding predictive information criteria for Bayesian models. 
#' \emph{Statistics and Computing} 24, 997-1016. 
#' 
#' Watanable, S. (2010). Asymptotic equivalence of Bayes cross validation and widely application information 
#' criterion in singular learning theory. \emph{Journal of Machine Learning Research} 11, 3571-3594.
#' 
#' @export
#' @examples  
#' # For examples see the function \code{\link{selection}}, 
#' # \code{\link{pattern}} or \code{\link{hurdle}}
#' # 
#' # 

pic <- function(x, criterion = "dic", module = "total") {
  if(!isTRUE(requireNamespace("loo", quietly = TRUE))) {
    stop("You need to install the R package 'loo'. Please run in your R terminal:\n install.packages('loo')")
  }
  if(class(x) != "missingHE") { 
    stop("Only objects of class 'missingHE' can be used") 
  }
  if(x$model_output$`model summary`$BUGSoutput$n.chains == 1){
    stop("information criteria cannot be computed if n.chain=1")
  }
  if(is.character(criterion) == FALSE | is.character(module) == FALSE) {
    stop("criterion and module must be character names")
  }
  criterion <- tolower(criterion)
  module <- tolower(module)
  if(!criterion %in% c("dic", "waic", "looic")) {
    stop("you must select a character name among those available. For details type 'help(pic)'")
  }
  if(!module %in% c("total", "e", "c", "both")) {
    stop("you must select a character name among those available. For details type 'help(pic)'")
  }
  m_e1 <- x$data_set$missing_effects$Control
  loglik_e1 <- x$model_output$loglik$effects$control[, m_e1 == 0]
  m_e2 <- x$data_set$missing_effects$Intervention
  loglik_e2 <- x$model_output$loglik$effects$intervention[, m_e2 == 0]
  m_c1 <- x$data_set$missing_costs$Control
  loglik_c1 <- x$model_output$loglik$costs$control[, m_c1 == 0]
  m_c2 <- x$data_set$missing_costs$Intervention
  loglik_c2 <- x$model_output$loglik$costs$intervention[, m_c2 == 0]
  if(x$model_output$type == "SELECTION" | x$model_output$type == "SELECTION_e" | 
     x$model_output$type == "SELECTION_c" | x$model_output$type == "SELECTION_ec" ) {
  loglik_me1 <- x$model_output$loglik$`missing indicators effects`$control
  loglik_me2 <- x$model_output$loglik$`missing indicators effects`$intervention
  loglik_mc1 <- x$model_output$loglik$`missing indicators costs`$control
  loglik_mc2 <- x$model_output$loglik$`missing indicators costs`$intervention
  if(criterion == "dic") {
  loglik_e <- rowSums(cbind(loglik_e1, loglik_e2))
  loglik_c <- rowSums(cbind(loglik_c1, loglik_c2))
  loglik_me <- rowSums(cbind(loglik_me1, loglik_me2))
  loglik_mc <- rowSums(cbind(loglik_mc1, loglik_mc2))
  loglik_e_c <- rowSums(cbind(loglik_e, loglik_c))
  loglik_total <- rowSums(cbind(loglik_e_c, loglik_me, loglik_mc))
  } else if(criterion == "looic" | criterion == "waic"){
    loglik_e <- cbind(loglik_e1, loglik_e2)
    loglik_c <- cbind(loglik_c1, loglik_c2)
    loglik_me <- cbind(loglik_me1, loglik_me2)
    loglik_mc <- cbind(loglik_mc1, loglik_mc2)
    loglik_e_c <- cbind(loglik_e, loglik_c)
    loglik_total <- cbind(loglik_e_c, loglik_me, loglik_mc)
   }
  }
  if(x$model_output$type == "PATTERN" | x$model_output$type == "PATTERN_e" | 
     x$model_output$type == "PATTERN_c" | x$model_output$type == "PATTERN_ec" ) {
  loglik_d1 <- x$model_output$loglik$`pattern indicators`$control
  loglik_d2 <- x$model_output$loglik$`pattern indicators`$intervention
  if(criterion == "dic") {
    loglik_e <- rowSums(cbind(loglik_e1, loglik_e2))
    loglik_c <- rowSums(cbind(loglik_c1, loglik_c2))
    loglik_d <- rowSums(cbind(loglik_d1, loglik_d2))
    loglik_e_c <- rowSums(cbind(loglik_e, loglik_c))
    loglik_total <- rowSums(cbind(loglik_e_c, loglik_d))
  } else if(criterion == "looic" | criterion == "waic"){
    loglik_e <- cbind(loglik_e1, loglik_e2)
    loglik_c <- cbind(loglik_c1, loglik_c2)
    loglik_d <- cbind(loglik_d1, loglik_d2)
    loglik_e_c <- cbind(loglik_e, loglik_c)
    loglik_total <- cbind(loglik_e_c, loglik_d)
   }
  }
  if(x$model_output$type == "HURDLE" | x$model_output$type == "HURDLE_e" | 
     x$model_output$type == "HURDLE_c" | x$model_output$type == "HURDLE_ec") {
    if(x$model_output$type == "HURDLE_ec" | x$model_output$type == "HURDLE_e") {
    m_de1 <- ifelse(is.na(x$data_set$structural_effects$Control) == TRUE, 1, 0)
    loglik_e1 <- x$model_output$loglik$effects$control[, m_e1 == 0]
    loglik_de1 <- x$model_output$loglik$`structural indicators effects`$control[, m_e1 == 0]
    m_de2 <- ifelse(is.na(x$data_set$structural_effects$Intervention) == TRUE, 1, 0)
    loglik_e2 <- x$model_output$loglik$effects$intervention[, m_e2 == 0]
    loglik_de2 <- x$model_output$loglik$`structural indicators effects`$intervention[, m_e2 == 0]
    } 
    if(x$model_output$type == "HURDLE_ec" | x$model_output$type == "HURDLE_c") {
    m_dc1 <- ifelse(is.na(x$data_set$structural_costs$Control) == TRUE, 1, 0)
    loglik_c1 <- x$model_output$loglik$costs$control[, m_c1 == 0]
    loglik_dc1 <- x$model_output$loglik$`structural indicators costs`$control[, m_c1 == 0]
    m_dc2 <- ifelse(is.na(x$data_set$structural_costs$Intervention) == TRUE, 1, 0)
    loglik_c2 <- x$model_output$loglik$costs$intervention[, m_c2 == 0]
    loglik_dc2 <- x$model_output$loglik$`structural indicators costs`$intervention[, m_c2 == 0]
    }
    if(criterion == "dic") {
      loglik_e <- rowSums(cbind(loglik_e1, loglik_e2))
      loglik_c <- rowSums(cbind(loglik_c1, loglik_c2))
      loglik_e_c <- rowSums(cbind(loglik_e, loglik_c))
      if(x$model_output$type == "HURDLE_ec") {
        loglik_de <- rowSums(cbind(loglik_de1, loglik_de2))
        loglik_dc <- rowSums(cbind(loglik_dc1, loglik_dc2))
        loglik_total <- rowSums(cbind(loglik_e_c, loglik_de, loglik_dc))
      } else if(x$model_output$type == "HURDLE_e") {
        loglik_de <- rowSums(cbind(loglik_de1, loglik_de2))
        loglik_total <- rowSums(cbind(loglik_e_c, loglik_de))
      } else if(x$model_output$type == "HURDLE_c") {
        loglik_dc <- rowSums(cbind(loglik_dc1, loglik_dc2))
        loglik_total <- rowSums(cbind(loglik_e_c, loglik_dc))
      }
    } else if(criterion == "looic" | criterion == "waic"){
      loglik_e <- cbind(loglik_e1, loglik_e2)
      loglik_c <- cbind(loglik_c1, loglik_c2)
      loglik_e_c <- cbind(loglik_e, loglik_c)
      if(x$model_output$type == "HURDLE_ec") {
        loglik_de <- cbind(loglik_de1, loglik_de2)
        loglik_dc <- cbind(loglik_dc1, loglik_dc2)
        loglik_total <- cbind(loglik_e_c, loglik_de, loglik_dc)
      } else if(x$model_output$type == "HURDLE_e") {
        loglik_de <- cbind(loglik_de1, loglik_de2)
        loglik_total <- cbind(loglik_e_c, loglik_de)
      } else if(x$model_output$type == "HURDLE_c") {
        loglik_dc <- cbind(loglik_dc1, loglik_dc2)
        loglik_total <- cbind(loglik_e_c, loglik_dc)
      }
    }
  }
  if(module == "total") { 
    loglik <- loglik_total
  } else if(module == "e") {
    loglik <- loglik_e
  } else if(module == "c") {
    loglik <- loglik_c
  } else if(module == "both") {
    loglik <- loglik_e_c
  }
  if(criterion == "dic") {
      d_bar <- mean(-2 * loglik)
      pD <- var(-2 * loglik) / 2
      dic <- d_bar + pD
      d_hat <- 2 * d_bar - dic 
      ic <- list("d_bar" = d_bar, "pD" = pD, "dic" = dic, "d_hat" = d_hat)
  }
  if(criterion == "waic") {
    waic_l <- suppressWarnings(loo::waic(loglik))
    elpd <- waic_l$estimates[1, 1]
    elpd_se <- waic_l$estimates[1, 2]
    p <- waic_l$estimates[2, 1]
    p_se <- waic_l$estimates[2, 2]
    waic <- waic_l$estimates[3, 1]
    waic_se <- waic_l$estimates[3, 2]
    pointwise <- waic_l$pointwise
    ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se, 
               "waic" = waic, "waic_se" = waic_se, "pointwise" = pointwise)
  }
  if(criterion == "looic") {
    loo_l <- suppressWarnings(loo::loo(loglik))
    elpd <- loo_l$estimates[1, 1]
    elpd_se <- loo_l$estimates[1, 2]
    p <- loo_l$estimates[2, 1]
    p_se <- loo_l$estimates[2, 2]
    looic <- loo_l$estimates[3, 1]
    looic_se <- loo_l$estimates[3, 2]
    pointwise <- loo_l$pointwise
    pareto_k <- loo_l$pareto_k
    ic <- list("elpd" = elpd, "elpd_se" = elpd_se, "p" = p, "p_se" = p_se, 
               "looic" = looic, "looic_se" = looic_se, "pointwise" = pointwise, "pareto_k" = pareto_k)
  }
  return(ic)
}















