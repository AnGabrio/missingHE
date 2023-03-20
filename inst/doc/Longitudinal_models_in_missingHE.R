## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
set.seed(1234)

## ----pbs_data, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#first 10 entries of u and c at time 1
head(PBS$u[PBS$time == 1], n = 10)
head(PBS$c[PBS$time == 1], n = 10)

## ----hist_u, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(2, 2))
hist(PBS$u[PBS$t==1 & PBS$time == 1], main = "Utility at time 1 - Control")
hist(PBS$u[PBS$t==2 & PBS$time == 1], main = "Utility at time 1 - Intervention")
hist(PBS$u[PBS$t==1 & PBS$time == 2], main = "Utility at time 2 - Control")
hist(PBS$u[PBS$t==2 & PBS$time == 2], main = "Utility at time 2 - Intervention")

## ----mv, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#proportions of missing values in the control group 
sum(is.na(PBS$u[PBS$t==1 & PBS$time == 1])) / length(PBS$u[PBS$t==1 & PBS$time == 1])  
sum(is.na(PBS$u[PBS$t==1 & PBS$time == 2])) / length(PBS$u[PBS$t==1 & PBS$time == 2])
sum(is.na(PBS$u[PBS$t==1 & PBS$time == 3])) / length(PBS$u[PBS$t==1 & PBS$time == 3])

#proportions of missing values in the intervention group
sum(is.na(PBS$u[PBS$t==2 & PBS$time == 1])) / length(PBS$u[PBS$t==2 & PBS$time == 1])
sum(is.na(PBS$u[PBS$t==2 & PBS$time == 2])) / length(PBS$u[PBS$t==2 & PBS$time == 2])  
sum(is.na(PBS$u[PBS$t==2 & PBS$time == 3])) / length(PBS$u[PBS$t==2 & PBS$time == 3])  

## ----costs_const, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
PBS$c <- PBS$c + 0.01

## ----selection1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.sel=selection_long(data = PBS, model.eff = u ~ 1, model.cost = c ~ u, 
  model.mu = mu ~ gender, 
  model.mc = mc ~ gender, type = "MAR", 
  n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none")

## ----selection1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  NN.sel=selection_long(data = PBS, model.eff = u ~ 1, model.cost = c ~ u,
#    model.mu = mu ~ gender,
#    model.mc = mc ~ gender, type = "MAR",
#    n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none")

## ----plot_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.sel, outcome = "effects", time_plot = 3)

## ----coef_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.sel, time = 2)

## ----summary_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.sel)

## ----sl1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
PBS$site <- factor(PBS$site)
NN.re=selection_long(data = PBS, model.eff = u ~ 1 + (1 | site), 
               model.cost = c ~ u + (1 | site), model.mu = mu ~ gender, 
               model.mc = mc ~ gender, type = "MAR", 
               n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none")

## ----pattern1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  PBS$site <- factor(PBS$site)
#  NN.re=selection_long(data = PBS, model.eff = u ~ 1 + (1 | site),
#                 model.cost = c ~ u + (1 | site), model.mu = mu ~ gender,
#                 model.mc = mc ~ gender, type = "MAR",
#                 n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none")

## ----coef_sl1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.re, time = 3, random = TRUE)

## ----prior_sl1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
my.prior <- list(
  "alpha0.prior" = c(0 , 0.0000001),
  "beta0.prior" = c(0, 0.0000001),
  "gamma0.prior.c"= c(0, 1),
  "gamma.prior.c" = c(0, 0.01),
  "sigma.prior.c" = c(0, 100)
)

## ----sl12_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.prior=selection_long(data = PBS, model.eff = u ~ 1 + (1 | site), 
               model.cost = c ~ u + (1 | site), model.mu = mu ~ gender, 
               model.mc = mc ~ gender, type = "MAR", 
               n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none",
               prior = my.prior)

## ----sl12, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  NN.prior=selection_long(data = PBS, model.eff = u ~ 1 + (1 | site),
#                 model.cost = c ~ u + (1 | site), model.mu = mu ~ gender,
#                 model.mc = mc ~ gender, type = "MAR",
#                 n.iter = 500, dist_u = "norm", dist_c = "norm", time_dep = "none",
#                 prior = my.prior)

## ----coef_sl12, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.prior, random = FALSE, time = 3)

## ----summary_sl12, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.prior, random = TRUE, time = 3)

