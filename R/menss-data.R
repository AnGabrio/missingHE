#' MenSS economic data on STIs
#'
#' Data from a pilot RCT trial (The MenSS trial) on youn men at risk of Sexually Trasmitted Infections (STIs).
#' A total of `159` individuals were enrolled in trial: `75` in the control (`t = 1`) and `84` in the active intervention (`t = 2`).
#' Clinical and health economic outcome data were collected via self-reported questionnaires at four time points throughout the study: baseline,
#' `3` months, `6` months and `12` months follow-up. Health economic data include utility scores related to quality of life and costs, from which QALYs and total costs were then computed using the area under the curve method
#' and by summing up the cost components at each time point. Clinical data include the total number of instances of unprotected sex and whether the individual was associated with an STI diagnosis or not.
#' Baseline data are available for the utilities (no baseline costs collected), instances of unprotected sex, sti diagnosis, age, ethnicity and employment variables.
#'
#' - **id**: id number.
#' - **e**: Quality Adjusted Life Years (QALYs).
#' - **c**: Total costs in pounds.
#' - **u.0**: baseline utilities.
#' - **age**: Age in years.
#' - **ethnicity**: binary: white (`1`) and other (`0`).
#' - **employment**: binary: working (`1`) and other (`0`).
#' - **trt**: Treatment arm indicator for the control (`t = 1`) and the active intervention (`t = 2`).
#' - **sex_inst.0**: baseline number of instances of unprotected sex.
#' - **sex_inst**: number of instances of unprotected sex at `12` months follow-up.
#' - **sti.0**: binary : baseline sti diagnosis (`1`) and no baseline sti diagnosis (`0`).
#' - **sti**: binary : sti diagnosis (`1`) and no sti diagnosis (`0`) at `12` months follow-up.
#' - **site**: site number.
#'
#' @docType data
#' @keywords datasets
#' @name MenSS
#' @usage data(MenSS)
#' @format A data frame with `159` rows and `12` variables
#' @references Bailey et al. (2016) Health Technology Assessment 20
#' (\href{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5204131/}{PubMed})
#'
#'
#' @examples
#' MenSS <- data(MenSS)
#' summary(MenSS)
#' str(MenSS)
"MenSS"
