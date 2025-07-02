## ----echo = FALSE, message = FALSE--------------------------------------------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
set.seed(1234)

## ----menss_data2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE------------------------
MenSS2 <- MenSS
MenSS2$e <- MenSS$sex_inst

#first 10 entries of e
head(MenSS2$e, n = 10)

## ----hist_sex, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(1,2))
hist(MenSS2$e[MenSS2$t==1], main = "N sex instances - Control")
hist(MenSS2$e[MenSS2$t==2], main = "N sex instances - Intervention")

## ----mv, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE---------------------------------
#proportions of missing values in the control group
sum(is.na(MenSS2$e[MenSS$t==1])) / length(MenSS2$e[MenSS$t==1])  

#proportions of missing values in the intervention group
sum(is.na(MenSS2$e[MenSS$t==2])) / length(MenSS2$e[MenSS$t==2])  

## ----costs_const, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE------------------------
 MenSS2$c <- MenSS2$c + 0.01

## ----selection1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE------
PG.sel=selection(data = MenSS2, model.eff = e ~ sex_inst.0, model.cost = c ~ 1, 
  model.me = me ~ age + ethnicity + employment, 
  model.mc = mc ~ age + ethnicity + employment, type = "MAR", 
  n.iter = 1000, dist_e = "pois", dist_c = "gamma")

## ----selection1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE------------------------
# PG.sel=selection(data = MenSS2, model.eff = e ~ sex_inst.0, model.cost = c ~ 1,
#   model.me = me ~ age + ethnicity + employment,
#   model.mc = mc ~ age + ethnicity + employment, type = "MAR",
#   n.iter = 1000, dist_e = "pois", dist_c = "gamma")

## ----plot_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(PG.sel, outcome = "effects")

## ----coef_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE--------------------
coef(PG.sel, prob = c(0.05, 0.95))

## ----summary_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE-----------------
summary(PG.sel)

## ----pattern1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE--------
PG.pat=pattern(data = MenSS2, model.eff = e ~ sex_inst.0 + (1 | site), 
               model.cost = c ~ 1 + (1 | site), type = "MAR", restriction = "AC", 
               n.iter = 1000, Delta_e = 0, Delta_c = 0, dist_e = "pois", dist_c = "gamma")

## ----pattern1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE--------------------------
# PG.pat=pattern(data = MenSS2, model.eff = e ~ sex_inst.0 + (1 | site),
#                model.cost = c ~ 1 + (1 | site), type = "MAR", restriction = "AC",
#                n.iter = 1000, Delta_e = 0, Delta_c = 0, dist_e = "pois", dist_c = "gamma")

## ----coef_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----------------------
coef(PG.pat, random = TRUE)

## ----prior_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----------------------
my.prior <- list(
  "alpha0.prior" = c(0 , 0.0000001),
  "alpha.prior" = c(0, 0.0000001),
  "beta0.prior" = c(0, 0.0000001),
  "gamma0.prior.c"= c(0, 1),
  "gamma.prior.c" = c(0, 0.01),
  "mu.b0.prior" = c(0, 0.001),
  "mu.g0.prior.c"= c(0, 0.001),
  "s.b0.prior" = c(0, 100),
  "s.g0.prior.c"= c(0, 100),
  "sigma.prior.c" = c(0, 10000)
)

## ----hurdle1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE---------
#remove added constant from costs
#MenSS2$c <- MenSS2$c - 0.01

PG.hur=hurdle(data = MenSS2, model.eff = e ~ sex_inst.0, model.cost = c ~ 1 + (1 | site),
  model.se = se ~ 1, model.sc = sc ~ age + (1 | site), type = "SAR", se = 1, sc = 0,
  n.iter = 1000, dist_e = "pois", dist_c = "gamma", prior = my.prior)

## ----hurdle1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE---------------------------
# #remove added constant from costs
# #MenSS2$c <- MenSS2$c - 0.01
# 
# PG.hur=hurdle(data = MenSS2, model.eff = e ~ sex_inst.0, model.cost = c ~ 1 + (1 | site),
#   model.se = se ~ 1, model.sc = sc ~ age + (1 | site), type = "SAR", se = NULL, sc = 0,
#   n.iter = 1000, dist_e = "pois", dist_c = "gamma", prior = my.prior)

## ----coef_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE-----------------------
coef(PG.hur, random = FALSE)

## ----summary_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE--------------------
coef(PG.hur, random = TRUE)

