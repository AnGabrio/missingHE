#' PBS economic data on intellectual disability and challenging behaviour
#'
#' Longitudinal data from a cluster RCT trial (The PBS trial) on people suffering from intellectual disability and challenging behaviour.
#' A total of 244 individuals across 23 sites were enrolled in the trial: 136 in the control (t=1) and 108 in the active intervention (t=2).
#' Health economic outcome data were collected via self-reported questionnaires at three time points throughout the study: baseline (time=1),
#' 6 months (time=2) and 12 months (time=3) follow-up, and included utility scores related to quality of life and costs. 
#' Baseline data are available for age, gender, ethnicity, living status, type of carer, marital status, and disability level variables. 
#'
#'
#' \describe{
#'   \item{id}{id number}
#'   \item{time}{time indicator}
#'   \item{u}{utilities}
#'   \item{c}{costs (in pounds)}
#'   \item{age}{Age in years}
#'   \item{gender}{binary: male (1) and female (0)}
#'   \item{ethnicity}{binary: white (1) and other (0)}
#'   \item{carer}{binary: paid carer (1) and family carer (0)}
#'   \item{marital}{binary: single (1) and married (0)}
#'   \item{living}{categorical: alone (1), with partner (2) and with parents (3)}
#'   \item{disability}{categorical: mild (1), moderate (2) and severe (3)}
#'   \item{site}{site number}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name PBS
#' @usage data(PBS)
#' @format A data frame with 732 rows and 16 variables
#' @references Hassiotis et al. (2014) BMC Psychiatry 14
#' (\href{https://pubmed.ncbi.nlm.nih.gov/25927187/}{PubMed})
#'
#'
#' @examples
#' PBS <- data(PBS)
#' summary(PBS)
#' str(PBS)
"PBS"