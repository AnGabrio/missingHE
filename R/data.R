#' MenSS economic data on STIs
#'
#' Data from a pilot RCT trial (The MenSS trial) on youn men at risk of Sexually Trasmitted Infections (STIs).
#' A total of 159 individuals were enrolled in trial: 75 in the control (t=1) and 84 in the active intervention (t=2).
#' Health related quality of life and cost data were collected via self-reported questionnaires at 4 time points throughout the study: at baseline,
#' 3 months, 6 months and 12 months follow-up. QALYs and total costs were then computed by combining the utility scores using the area under the curve method
#' and by summing up the cost components at each time point. 
#'
#' \describe{
#'   \item{e}{id}{id number}
#'   \item{u.0}{baseline utilities}
#'   \item{u.1}{utilities at 3-months follow-up}
#'   \item{u.2}{utilities at 6-months follow-up}
#'   \item{u.3}{utilities at 12-months follow-up}
#'   \item{e}{Quality Adjusted Life Years (QALYs)}
#'   \item{c.1}{costs in pounds at 3-months follow-up}
#'   \item{c.2}{costs in pounds at 6-months follow-up}
#'   \item{c.3}{costs in pounds at 12-months follow-up}
#'   \item{c}{Total costs in pounds}
#'   \item{age}{Age in years}
#'   \item{ethnicity}{binary: white (1) and other (0)}
#'   \item{employment}{binary: working (1) and other (0)}
#'   \item{t}{Treatment arm indicator for the control (t=1) and the active intervention (t=2)}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name MenSS
#' @usage data(MenSS)
#' @format A data frame with 159 rows and 14 variables
#' @references Bailey et al. (2016) Health Technology Assessment 20
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204131/}{PubMed})
#'
#'
#' @examples
#' data(MenSS)
#' summary(MenSS)
#' str(MenSS)
"MenSS"