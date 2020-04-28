#' An internal function to separate the fixed and random effects components from an object of class formula
#' 
#' @keywords random effects models
#' @param term formula to be processed
#' @examples
#' #Internal function only
#' #no examples
#' #
#' #

nobars_ <- function(term)
{
  if (!anyBars(term)) return(term)
  if (isBar(term)) return(NULL)
  if (isAnyArgBar(term)) return(NULL)
  if (length(term) == 2) {
    nb <- nobars_(term[[2]])
    if(is.null(nb)) return(NULL)
    term[[2]] <- nb
    return(term)
  }
  nb2 <- nobars_(term[[2]])
  nb3 <- nobars_(term[[3]])
  if (is.null(nb2)) return(nb3)
  if (is.null(nb3)) return(nb2)
  term[[2]] <- nb2
  term[[3]] <- nb3
  term
}
