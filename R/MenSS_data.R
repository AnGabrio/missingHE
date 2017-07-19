#' MenSS economic data on STIs
#'
#' Data from a pilot RCT trial (The MenSS trial) on Sexually Trasmitted Infections (STIs)
#' QALY and cost data on 159 individuals (75 in the control and 84 in the active intervention) of whom
#' 113 are missing (48 in the control and 65 in the active intervention).
#' Data were collected via self-reported utility/cost questionnaires collected at 4 time points throughout the study: at baseline,
#' 3 months, 6 months and 12 months. QALYs were computed by combining the utility scores at each time point using the AUC method while
#' Total Costs were obtained by summing up all the cost components. 
#'
#' @docType data
#'
#' @usage data(MenSS_data)
#'
#' @format A data frame with 159 rows and 4 variables:
#' \describe{
#'   \item{e}{QALYs}
#'   \item{c}{Total Costs, in pounds}
#'   \item{t}{Treatment Arm indicator, 1 for the control and 2 for the active intervention}
#'   \item{X1}{Age in years}
#' }
#'
#' @keywords datasets
#'
#' @references Bailey et al. (2016) Health Technology Assessment 20
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204131/}{PubMed})
#'
#'
#' @examples
#' data(MenSS_data)
#' summary(MenSS_data)
#' str(MenSS_data)
"MenSS_data"
