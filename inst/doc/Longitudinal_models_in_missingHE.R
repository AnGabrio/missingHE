## ----echo = FALSE, message = FALSE--------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
set.seed(1234)

## ----menss--------------------------------------------------------------------
str(PBS)

## ----hist, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(2,2))
hist(PBS$e[PBS$time==2], main = "utilities - 6 months")
hist(PBS$e[PBS$time==2], main = "utilities - 6 months")
hist(PBS$c[PBS$time==3], main = "costs - 12 months")
hist(PBS$c[PBS$time==3], main = "costs - 12 months")

## ----mv, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#proportions of missing values in the control group at 6 months
c(sum(is.na(PBS$e[PBS$time==2 & PBS$t==1])) / length(PBS$e[PBS$time==2 & PBS$t==1]),  
sum(is.na(PBS$c[PBS$time==2 & PBS$t==1])) / length(PBS$c[PBS$time==2 & PBS$t==1]))

#proportions of missing values in the intervention group at 12 months
c(sum(is.na(PBS$e[PBS$time==2 & PBS$t==2])) / length(PBS$e[PBS$time==2 & PBS$t==2]),  
sum(is.na(PBS$c[PBS$time==2 & PBS$t==2])) / length(PBS$c[PBS$time==2 & PBS$t==2]))

## ----costs_const, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
PBS$c <- PBS$c + 0.01

## ----long1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.long=long_miss(data = PBS, model.eff = e ~ 1, model.cost = c ~ e, 
    model.me = me ~ gender, model.mc = mc ~ gender, type = "MAR",
    time_dep = "AR1", n.iter = 500, 
    dist_e = "normal", dist_c = "normal", ppc = TRUE)

## ----long1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# NN.long=long_miss(data = PBS, model.eff = e ~ 1, model.cost = c ~ e,
#     model.me = me ~ gender, model.mc = mc ~ gender, type = "MAR",
#     time_dep = "AR1", n.iter = 500,
#     dist_e = "normal", dist_c = "normal", ppc = TRUE)

## ----plot_long1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.long, outcome = "effects", time_plot = 3)

## ----coef_long1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.long, time = 2)

## ----summary_long1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.long)

## ----diag_long1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
diagnostic(NN.long, type = "traceplot", param = "mu.c")

## ----diag_long2, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
diagnostic(NN.long, type = "Rhat", param = "mu.c")

## ----ppc_long, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
ppc(NN.long, type = "dens", outcome = "effects", time_plot = 3, ndisplay = 3)

## ----long2_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
my.prior <- list(
  "delta.prior.e" = c(0, 1),
  "sigma.prior.c" = c(0, 20000)
)
PBS$site <- factor(PBS$site)

LG.long=long_miss(data = PBS, model.eff = e ~ 1 + (1|site), model.cost = c ~ 1 + (1|site), 
    model.me = me ~ e, model.mc = mc ~ 1, type = "MNAR",
    time_dep = "none", n.iter = 500, prior = my.prior,
    dist_e = "logistic", dist_c = "gamma", ppc = FALSE)

## ----long2, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
# my.prior <- list(
#   "delta.prior.e" = c(0, 1),
#   "sigma.prior.c" = c(0, 20000)
# )
# PBS$site <- factor(PBS$site)
# 
# LG.long=long_miss(data = PBS, model.eff = e ~ 1 + (1|site), model.cost = c ~ 1 + (1|site),
#     model.me = me ~ e, model.mc = mc ~ 1, type = "MNAR",
#     time_dep = "none", n.iter = 500, prior = my.prior,
#     dist_e = "logistic", dist_c = "gamma", ppc = FALSE)

## ----coef_sl1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(LG.long, time = 3, random = TRUE)

## ----summ_sl1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.long)

## ----bcea_sl1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
library(BCEA)
ceplane.plot(NN.long$cea)
ceac.plot(NN.long$cea)

