# PBS economic data on intellectual disability and challenging behaviour

Longitudinal data from a cluster RCT trial (The PBS trial) on people
suffering from intellectual disability and challenging behaviour. A
total of `244` individuals across `23` sites were enrolled in the trial:
`136` in the control (`t = 1`) and `108` in the active intervention
(`t = 2`). Health economic outcome data were collected via self-reported
questionnaires at three time points throughout the study: baseline
(`time = 1`), `6` months (`time = 2`) and `12` months (`time = 3`)
follow-up, and included utility scores related to quality of life and
costs. Baseline data are available for age, gender, ethnicity, living
status, type of carer, marital status, and disability level variables.

## Usage

``` r
data(PBS)
```

## Format

A data frame with `732` rows and `16` variables

## Details

- **id**: id number.

- **time**: time indicator.

- **e**: Health-Related Quality of Life (HRQoL) utilities.

- **c**: Costs in pounds.

- **age**: Age in years.

- **gender**: binary: male (`1`) and female (`0`).

- **ethnicity**: binary: white (`1`) and other (`0`).

- **trt**: Treatment arm indicator for the control (`t = 1`) and the
  active intervention (`t = 2`).

- **carer**: binary: paid carer (`1`) and family carer (`0`).

- **marital**: binary: single (`1`) and married (`0`).

- **living**: categorical: alone (`1`), with partner (`2`) and with
  parents (`3`).

- **disability**: categorical: mild (`1`), moderate (`2`) and severe
  (`3`).

- **site**: site number.

## References

Hassiotis et al. (2014) BMC Psychiatry 14
([PubMed](https://pubmed.ncbi.nlm.nih.gov/25927187/))

## Examples

``` r
PBS <- data(PBS)
summary(PBS)
#>    Length  N.unique   N.blank Min.nchar Max.nchar 
#>         1         1         0         3         3 
str(PBS)
#>  chr "PBS"
```
