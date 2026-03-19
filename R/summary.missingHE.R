#' Summary method for objects in the class \code{missingHE}
#'
#' Produces a table printout with some summary results of the health economic evaluation probabilistic model
#' run using the function \code{\link{selection}}, \code{\link{pattern}}, \code{\link{hurdle}} or \code{\link{lmdm}}.
#' @param object A \code{missingHE} object containing the results of the Bayesian modelling and the economic evaluation
#' @param incremental Logical. If \code{incremental} is \code{TRUE}, incremental CE results are printed.
#' @param prob A numeric vector of probabilities within the range (0, 1), representing the upper and lower
#'  CI quantiles to be calculated and returned for the posterior estimates.
#' @param digits Integer indicating the number of decimal places to be used for rounding (default = 3).
#' @param ... Additional arguments affecting the summary produced.
#' @return Prints tables with information on the CE results based on a model fitted using the function \code{selection}, \code{pattern} or \code{hurdle}. 
#' Summary information on the main parameters of interests is provided.
#' @seealso \code{\link{selection}} \code{\link{pattern}} \code{\link{hurdle}} \code{\link{lmdm}} \code{\link{diagnostic}} \code{\link{plot.missingHE}}
#' @author Andrea Gabrio
#' @references 
#' Baio, G.(2012). \emph{Bayesian Methods in Health Economcis}. CRC/Chapman Hall, London.
#' @importFrom stats quantile
#' @export
#' @examples 
#' # For examples see the function selection, pattern or hurdle
#' #
#' #

summary.missingHE <- function(object, incremental = FALSE, prob = c(0.025, 0.975),
                              digits = 3, ...) {
  exArgs <- list(...)
  if(!inherits(object, "missingHE")) {
    stop("Only objects of class 'missingHE' can be used")}
  if(!is.logical(incremental)) {
    stop("Please provide 'incremental' as a logical value")}
  if(!is.vector(prob) | length(prob) != 2 | !is.numeric(prob) | any(prob <= 0 | prob >= 1)) {
    stop("Please provide a lower and an upper quantile for the posterior results")}
  if(exists("ref", where = exArgs)) {
    if(length(exArgs$ref) != 1 | !is.numeric(exArgs$ref)) { stop("Please provide a single numeric indicator for the reference treatment")}
    if(exArgs$ref <= 0 | !exArgs$ref %% 1 == 0) { stop("Please provide a valid indicator value for the reference treatment")}
    ref = exArgs$ref } else { ref = object$cea$ref}
  if(exists("wtp", where = exArgs)) {
    if(length(exArgs$wtp) != 1 | !is.numeric(exArgs$wtp)) { stop("Please provide a single numeric value for the acceptance threshold")}
    if(exArgs$wtp <= 0 | !exArgs$wtp %% 1 == 0) { stop("Please provide a valid value for the acceptance threshold")}
    wtp = exArgs$wtp } else { wtp = object$cea$Kmax}
  n_trt <- length(object$data_set$data_raw$n)
  trt_lev <- names(object$data_set$data_raw$n)
  mu_e <- as.matrix(object$cea$e)
  mu_c <- as.matrix(object$cea$c)
  nmb <- mu_e * wtp - mu_c
  colnames(mu_e) <- colnames(mu_c) <- colnames(nmb) <- trt_lev
  delta_e <- as.matrix(mu_e[, ref] - mu_e[, -ref])
  delta_c <- as.matrix(mu_c[, ref] - mu_c[, -ref])
  inmb <- delta_e * wtp - delta_c
  incr_lev <- paste(trt_lev[ref], "vs", trt_lev[-ref])
  colnames(delta_e) <- colnames(delta_c) <- colnames(inmb) <- incr_lev
  mean_mu_e <- round(apply(mu_e, 2, mean, na.rm = TRUE), digits = digits)
  sd_mu_e <- round(apply(mu_e, 2, sd, na.rm = TRUE), digits = digits)
  ql_mu_e <- round(apply(mu_e, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_mu_e <- round(apply(mu_e, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  mean_mu_c <- round(apply(mu_c, 2, mean, na.rm = TRUE), digits = digits)
  sd_mu_c <- round(apply(mu_c, 2, sd, na.rm = TRUE), digits = digits)
  ql_mu_c <- round(apply(mu_c, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_mu_c <- round(apply(mu_c, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  mean_nmb <- round(apply(nmb, 2, mean, na.rm = TRUE), digits = digits)
  sd_nmb <- round(apply(nmb, 2, sd, na.rm = TRUE), digits = digits)
  ql_nmb <- round(apply(nmb, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_nmb <- round(apply(nmb, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  mean_delta_e <- round(apply(delta_e, 2, mean, na.rm = TRUE), digits = digits)
  sd_delta_e <- round(apply(delta_e, 2, sd, na.rm = TRUE), digits = digits)
  ql_delta_e <- round(apply(delta_e, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_delta_e <- round(apply(delta_e, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  mean_delta_c <- round(apply(delta_c, 2, mean, na.rm = TRUE), digits = digits)
  sd_delta_c <- round(apply(delta_c, 2, sd, na.rm = TRUE), digits = digits)
  ql_delta_c <- round(apply(delta_c, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_delta_c <- round(apply(delta_c, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  mean_inmb <- round(apply(inmb, 2, mean, na.rm = TRUE), digits = digits)
  sd_inmb <- round(apply(inmb, 2, sd, na.rm = TRUE), digits = digits)
  ql_inmb <- round(apply(inmb, 2, quantile, probs = prob[1], na.rm = TRUE), digits = digits)
  qu_inmb <- round(apply(inmb, 2, quantile, probs = prob[2], na.rm = TRUE), digits = digits)
  icer <- round(mean_delta_c / mean_delta_e, digits = digits)
  tbl_abs_e <- t(rbind(mean_mu_e, sd_mu_e, ql_mu_e, qu_mu_e))
  tbl_abs_c <- t(rbind(mean_mu_c, sd_mu_c, ql_mu_c, qu_mu_c))
  tbl_abs_nmb <- t(rbind(mean_nmb, sd_nmb, ql_nmb, qu_nmb))
  colnames(tbl_abs_e) <- colnames(tbl_abs_c) <- colnames(tbl_abs_nmb) <- c("Mean", "SD", "QL", "QU")
  tbl_inc_e <- t(rbind(mean_delta_e, sd_delta_e, ql_delta_e, qu_delta_e))
  tbl_inc_c <- t(rbind(mean_delta_c, sd_delta_c, ql_delta_c, qu_delta_c))
  tbl_inc_nmb <- t(rbind(mean_inmb, sd_inmb, ql_inmb, qu_inmb, icer))
  colnames(tbl_inc_e) <- colnames(tbl_inc_c) <- c("Mean", "SD", "QL", "QU")
  colnames(tbl_inc_nmb) <- c("Mean", "SD", "QL", "QU", "ICER")
  if(!incremental) {
    tbl_res1 <- tbl_abs_e; tbl_res2 <- tbl_abs_c; tbl_res3 <- tbl_abs_nmb
    text1 <- paste("\n \n Mean effects by intervention \n")
    text2 <- paste("\n \n Mean costs by intervention \n")
    text3 <- paste("\n \n Mean net monetary benefit by intervention and wtp =", wtp, "\n")
  }
  if(incremental) {
    tbl_res1 <- tbl_inc_e; tbl_res2 <- tbl_inc_c; tbl_res3 <- tbl_inc_nmb
    text1 <- paste("\n \n Mean incremental effects \n")
    text2 <- paste("\n \n Mean incremental costs \n")
    text3 <- paste("\n \n Mean incremental net monetary benefit and wtp =", wtp, "\n")
  }
  cat(paste("\n Cost-effectiveness analysis summary \n \n", "CE parameter estimates under", object$type, "assumption"))
  cat(text1)
  print(tbl_res1, quote = FALSE, digits = digits, justify = "center")
  cat(text2)
  print(tbl_res2, quote = FALSE, digits = digits, justify = "center")
  cat(text3)
  print(tbl_res3, quote = FALSE, digits = digits, justify = "center")
 }