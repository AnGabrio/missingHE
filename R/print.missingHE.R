#' Print method for objects in the class `missingHE`
#'
#' Produces a printout of the (posterior mean) estimate for each regression coefficient included in the outcome and missingness/structural value
#' models based on the health economic evaluation probabilistic model run using the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#' @param x A `missingHE` object containing the results of the Bayesian modelling and the economic evaluation.
#' @param digits Integer indicating the number of decimal places to be used for rounding (default = 3).
#' @param \dots Additional arguments affecting the summary produced.
#' @return Prints all fixed regression coefficients based on the type of model fitted using the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#' Summary estimates based on posterior means of each fixed regression coefficient included into the fitted models are provided.
#' @seealso [selection()] [pattern()] [hurdle()] [lmdm()] [diagnostic()] [print.missingHE()]
#' @author Andrea Gabrio
#' @references
#' Baio, G.(2012). *Bayesian Methods in Health Economcis*. CRC/Chapman Hall, London.
#' @importFrom stats quantile
#' @export
#' @examples
#' # For examples see the function [selection()], [pattern()], [hurdle()] or [lmdm()].
#' #
#' #

print.missingHE <- function(x, digits = 3, ...) {
  exArgs <- list(...)
  if(!inherits(x, "missingHE")) {
    stop("Only objects of class 'missingHE' can be used")}
  mod_est <- lapply(x$mod_est, round, digits = digits)
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(mod_est[names(mod_est)[1]])
  print(mod_est[names(mod_est)[2]])
  print(mod_est[names(mod_est)[3]])
  if(x$model_output$method != "PATTERN") {
    print(mod_est[names(mod_est)[4]])}
  invisible(mod_est)
 }
