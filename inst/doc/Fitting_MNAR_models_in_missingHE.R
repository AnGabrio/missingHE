## ----echo = FALSE, message = FALSE--------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
set.seed(1014)

## ----selection1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.sel1=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1, 
  model.me = me ~  e, model.mc = mc ~ 1, type = "MNAR", 
  n.iter = 1000, dist_e = "norm", dist_c = "norm")

## ----selection1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.sel1=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   model.me = me ~  e, model.mc = mc ~ 1, type = "MNAR",
#   n.iter = 1000, dist_e = "norm", dist_c = "norm")

## ----print_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
print(NN.sel1)

## ----prior_selection2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
my.prior <- list(
  "delta.prior.e" = c(10, 1)
)

## ----selection2_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.sel2=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1, 
  model.me = me ~  e, model.mc = mc ~ 1, type = "MNAR", 
  n.iter = 1000, dist_e = "norm", dist_c = "norm", prior = my.prior)

## ----selection2, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.sel2=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   model.me = me ~  e, model.mc = mc ~ 1, type = "MNAR",
#   n.iter = 1000, dist_e = "norm", dist_c = "norm", prior = my.prior)

## ----print_selection2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
print(NN.sel2)

## ----BCEA_selection, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(1,2))
BCEA::ceac.plot(NN.sel1$cea)
BCEA::ceac.plot(NN.sel2$cea)

## ----sp_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
Delta_e <- matrix(NA, 2, 2)
Delta_e[1, ] <- c(- 0.3, - 0.2) 
Delta_e[2, ] <- c(-0.1, 0)

## ----pattern2_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.pat2=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1, 
  type = "MNAR", restriction = "CC", n.iter = 1000, Delta_e = Delta_e, Delta_c = 0, 
  dist_e = "norm", dist_c = "norm")

## ----pattern2, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.pat2=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   type = "MNAR", restriction = "CC", n.iter = 1000, Delta_e = Delta_e, Delta_c = 0,
#   dist_e = "norm", dist_c = "norm")

## ----print_pattern2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
print(NN.pat2)

## ----pattern1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.pat1=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1, 
  type = "MAR", restriction = "CC", n.iter = 1000, Delta_e = 0, Delta_c = 0, 
  dist_e = "norm", dist_c = "norm")

## ----pattern1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.pat1=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   type = "MAR", restriction = "CC", n.iter = 1000, Delta_e = 0, Delta_c = 0,
#   dist_e = "norm", dist_c = "norm")

## ----BCEA_pattern, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(1,2))
BCEA::ceac.plot(NN.pat1$cea)
BCEA::ceac.plot(NN.pat2$cea)

## ----hurdle1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.hur1=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
  model.se = se ~ 1, model.sc = sc ~ 1, type = "SCAR", se = 1, sc = 0,
  n.iter = 1000, dist_e = "norm", dist_c = "norm")

## ----hurdle1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.hur1=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   model.se = se ~ 1, model.sc = sc ~ 1, type = "SCAR", se = 1, sc = 0,
#   n.iter = 1000, dist_e = "norm", dist_c = "norm")

## ----d_hur2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
d_e <- ifelse(MenSS$e == 1, 1, 0)

#number of ones
sum(d_e == 1, na.rm = T)

## ----d_age_hur2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
myd_e <- ifelse(is.na(d_e) & MenSS$age < 22, 1, d_e)

#number of ones
sum(myd_e == 1, na.rm = T)

## ----hurdle2_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.hur2=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
  model.se = se ~ 1, model.sc = sc ~ 1, type = "SCAR", se = 1, sc = 0,
  n.iter = 1000, dist_e = "norm", dist_c = "norm", d_e = myd_e)

## ----hurdle2, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.hur2=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ 1,
#   model.se = se ~ 1, model.sc = sc ~ 1, type = "SCAR", se = 1, sc = 0,
#   n.iter = 1000, dist_e = "norm", dist_c = "norm", d_e = myd_e)

## ----print_hurdle2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
print(NN.hur2)

## ----plot_hurdle2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.hur2, outcome = "effects")

## ----BCEA_hurdle, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(1,2))
BCEA::ceac.plot(NN.hur1$cea)
BCEA::ceac.plot(NN.hur2$cea)

