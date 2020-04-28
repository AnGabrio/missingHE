#' An internal function to extract the random effects component from an object of class formula
#' 
#' @keywords random effects models
#' @param term formula to be processed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

fb <- function(term) {
  if (is.name(term) || !is.language(term)) return(NULL)
  if (term[[1]] == as.name("(")) return(fb(term[[2]]))
  stopifnot(is.call(term))
  if (term[[1]] == as.name('|')) return(term)
  if (length(term) == 2) return(fb(term[[2]]))
  c(fb(term[[2]]), fb(term[[3]]))
}





