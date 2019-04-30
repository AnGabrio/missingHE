#' An internal function to summarise results from BUGS model
#'
#' This function hides missing data distribution from summary results of BUGS models
#' @keywords summary JAGS model
#' @param x The \code{rjags}, \code{rjags.parallel}, or \code{mcmc.list} object
#'   for which results will be printed.
#' @param params Character vector or a regular expression pattern. The 
#'   parameters for which results will be printed (unless \code{invert} is
#'   \code{FALSE}, in which case results for all parameters other than those
#'   given in \code{params} will be returned). If \code{regex} is \code{FALSE},
#'   only those parameters that match \code{params} exactly will be returned. If 
#'   \code{regex} is \code{TRUE}, \code{param} should be a character string 
#'   giving the regular expression pattern to be matched.
#' @param regex If \code{regex} is \code{TRUE}, then \code{param} is 
#'   expected to be a single string giving a text pattern to be matched. 
#'   Parameters with names matching the pattern will be returned (unless 
#'   \code{invert} is \code{TRUE}, which results in all parameters that do not 
#'   match the pattern being returned). Text pattern matching uses regular 
#'   expressions (\code{\link{regex}}). 
#' @param invert Logical. If \code{invert} is \code{TRUE}, only those parameters
#'   that do not match elements of \code{params} will be returned.
#' @param probs A numeric vector of probabilities within range [0, 1],
#'   representing the sample quantiles to be calculated and returned.
#' @param signif If supplied, all columns other than \code{n.eff} will have 
#'   their values rounded such that the most extreme value has the specified
#'   number of significant digits.
#' @param ... Additional arguments accepted by \code{\link{grep}}, e.g. 
#'   \code{perl=TRUE}, to allow look-around pattern matching.
#' @importFrom stats quantile
#' @importFrom methods is
#' @export 
#' @examples
#' # Internal function only
#' # no examples

jagsresults <-
  function(x, params, regex=FALSE, invert=FALSE, 
           probs=c(0.025, 0.25, 0.5, 0.75, 0.975), signif, ...) {
    if(!isTRUE(requireNamespace("methods"))) {
      stop("You need to install the R package 'methods'. 
           Please run in your R terminal:\n install.packages('methods')")
    }
    if(!regex) {
      params <- paste0(gsub('(?=\\.|\\[|\\])', '\\1\\\\', params, perl=TRUE),
                       '(\\[.*\\])?', collapse='|')
      params <- paste("^", gsub("\\|", "$|^", params), '$', sep = "")
    } else if(length(params) > 1) {
      stop("If 'regex' is TRUE, 'params' must be a single regex string.",
           call.=FALSE)
    }
    if(any(methods::is(x) %in% c('rjags.parallel', 'rjags'))) {
      nm <- dimnames(x$BUGSoutput$sims.array)[[3]]
      i <- grep(params, nm, invert=invert, ...)
      if(length(i) == 0) stop("No parameters match 'params'", call.=FALSE)
      samp <- x$BUGSoutput$sims.array[, , i, drop=FALSE] 
      rhat_neff <- x$BUGSoutput$summary[i, c('Rhat', 'n.eff'), drop=FALSE]
      out <- cbind(t(apply(
        samp, 3, function(x) 
          c(mean = mean(x), sd = sd(x), quantile(x, probs = probs)))), rhat_neff)
    } else if(any(is(x) == 'mcmc.list')) {
      nm <- colnames(x[[1]])
      i <- grep(params, nm, invert = invert, ...)
      if(length(i) == 0) stop("No parameters match 'params'", call. = FALSE)
      out <- t(apply(do.call(rbind, x), 2, function(z) {
        c(mean = mean(z), sd = sd(z), quantile(z, probs))
      }))[i, , drop = FALSE]
    } else {
      stop("x must be an 'mcmc.list' or 'rjags'  object.")
    }
    if(!missing(signif)) {
      out[, colnames(out) != 'n.eff'] <- 
        signif(out[, colnames(out) != 'n.eff'], signif)  
      out
    } else out
  }
