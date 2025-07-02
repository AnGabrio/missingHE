#' An internal function to detect the random effects component from an object of class formula
#' 
#' @keywords random effects models
#' @param term formula to be processed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

isAnyArgBar <- function(term) {
  if ((term[[1]] != as.name("~")) && (term[[1]] != as.name("("))) {
    for(i in seq_along(term)) {
      if(isBar(term[[i]])) return(TRUE)
    }
  }
  FALSE
}
