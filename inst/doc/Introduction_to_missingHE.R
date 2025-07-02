## ---- echo = FALSE, message = FALSE-------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
set.seed(1014)

## ----menss--------------------------------------------------------------------
str(MenSS)

## ----hist, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(2,2))
hist(MenSS$e[MenSS$t==1], main = "QALYs - Control")
hist(MenSS$e[MenSS$t==2], main = "QALYs - Intervention")
hist(MenSS$c[MenSS$t==1], main = "Costs - Control")
hist(MenSS$c[MenSS$t==2], main = "Costs - Intervention")

## ----sv, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#proportions of ones and zeros in the control group
c(sum(MenSS$e[MenSS$t==1]==1, na.rm = TRUE) / length(MenSS$e[MenSS$t==1]),  
sum(MenSS$c[MenSS$t==1]==0, na.rm = TRUE) / length(MenSS$e[MenSS$t==1]))

#proportions of ones and zeros in the intervention group
c(sum(MenSS$e[MenSS$t==2]==1, na.rm = TRUE) / length(MenSS$e[MenSS$t==2]),  
sum(MenSS$c[MenSS$t==2]==0, na.rm = TRUE) / length(MenSS$e[MenSS$t==2]))

## ----mv, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#proportions of missing values in the control group
c(sum(is.na(MenSS$e[MenSS$t==1])) / length(MenSS$e[MenSS$t==1]),  
sum(is.na(MenSS$c[MenSS$t==1])) / length(MenSS$e[MenSS$t==1]))

#proportions of missing values in the intervention group
c(sum(is.na(MenSS$e[MenSS$t==2])) / length(MenSS$e[MenSS$t==2]),  
sum(is.na(MenSS$c[MenSS$t==2])) / length(MenSS$e[MenSS$t==2]))

## ----selection1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.sel=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e, 
  model.me = me ~ age + ethnicity + employment, 
  model.mc = mc ~ age + ethnicity + employment, type = "MAR", 
  n.iter = 1000, dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----selection1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  NN.sel=selection(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e,
#    model.me = me ~ age + ethnicity + employment,
#    model.mc = mc ~ age + ethnicity + employment, type = "MAR",
#    n.iter = 1000, dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----print_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
print(NN.sel, value.mis = FALSE, only.means = TRUE)

## ----coef_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.sel, random = FALSE)

## ----summary_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.sel)

## ----BCEA_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
par(mfrow=c(1,2))
BCEA::ceplane.plot(NN.sel$cea)
BCEA::ceac.plot(NN.sel$cea)

## ----pattern1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.pat=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e, 
  type = "MAR", restriction = "CC", n.iter = 1000, Delta_e = 0, Delta_c = 0, 
  dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----pattern1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  NN.pat=pattern(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e,
#    type = "MAR", restriction = "CC", n.iter = 1000, Delta_e = 0, Delta_c = 0,
#    dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----coef_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.pat, random = FALSE)

## ----summary_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.pat)

## ----hurdle1_no, eval=TRUE, echo=FALSE, include=FALSE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
NN.hur=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e,
  model.se = se ~ 1, model.sc = sc ~ age, type = "SAR", se = 1, sc = 0,
  n.iter = 1000, dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----hurdle1, eval=FALSE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
#  NN.hur=hurdle(data = MenSS, model.eff = e ~ u.0, model.cost = c ~ e,
#    model.se = se ~ 1, model.sc = sc ~ age, type = "SAR", se = 1, sc = 0,
#    n.iter = 1000, dist_e = "norm", dist_c = "norm", ppc = TRUE)

## ----coef_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
coef(NN.hur, random = FALSE)

## ----summary_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE----
summary(NN.hur)

## ----diag_sel1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
diagnostic(NN.sel, type = "denplot", param = "mu.e", theme = NULL)

## ----diag_pat1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
diagnostic(NN.pat, type = "traceplot", param = "mu.c", theme = NULL)

## ----diag_hur1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
diagnostic(NN.hur, type = "acf", param = "p.c", theme = "base")

## ----plot_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.sel, class = "scatter", outcome = "all")

## ----plot_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.pat, class = "histogram", outcome = "all")

## ----plot_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
plot(NN.hur, class = "scatter", outcome = "costs_arm1")

## ----ppc_selection1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
ppc(NN.sel, type = "histogram", outcome = "effects_arm1", ndisplay = 8)

## ----ppc_pattern1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
ppc(NN.pat, type = "dens", outcome = "effects_arm1", ndisplay = 8)

## ----ppc_hurdle1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
ppc(NN.hur, type = "dens_overlay", outcome = "all", ndisplay = 25)

## ----pic_model1, eval=TRUE, echo=TRUE, comment=NA,warning=FALSE,error=FALSE,message=FALSE, fig.width=15,fig.height=9,out.width='65%',fig.align='center'----
pic_sel <- pic(NN.sel, criterion = "waic", module = "both")
pic_pat <- pic(NN.pat, criterion = "waic", module = "both")
pic_hur <- pic(NN.hur, criterion = "waic", module = "both")

#print results
c(pic_sel$waic, pic_pat$waic, pic_hur$waic)

