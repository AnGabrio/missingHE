#' An internal function to detect the random effects component from an object of class formula
#' 
#' @keywords internal
#' @param term formula to be processed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

anyBars <- function(term) {
  any(c('|','||') %in% all.names(term))
}

