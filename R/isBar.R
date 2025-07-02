#' An internal function to detect the random effects component from an object of class formula
#' 
#' @keywords random effects models
#' @param term formula to be processed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

isBar <- function(term) {
  if(is.call(term)) {
    if((term[[1]] == as.name("|")) || (term[[1]] == as.name("||"))) {
      return(TRUE)
    }
  }
  FALSE
}