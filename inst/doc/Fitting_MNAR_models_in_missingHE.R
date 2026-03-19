## ----echo = FALSE, message = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
library(bookdown)
library(ggplot2)

MenSS$trt <- factor(MenSS$trt)
levels(MenSS$trt) <- c("SoC", "MenSS")

options(width = 300)
set.seed(1014)

## ----mnar-sel1-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.sel1 = selection(data = MenSS, model.eff = e ~ trt + u.0, 
                    model.cost = c ~ trt, model.me = me ~ e, 
                    model.mc = mc ~ 1, type = "MNAR", 
                    n.iter = 1000, ref = 2, 
                    dist_e = "norm", dist_c = "norm")

## ----mnar-sel1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.sel1 = selection(data = MenSS, model.eff = e ~ trt + u.0,
#                     model.cost = c ~ trt, model.me = me ~ e,
#                     model.mc = mc ~ 1, type = "MNAR",
#                     n.iter = 1000, ref = 2,
#                     dist_e = "norm", dist_c = "norm")

## ----print-sel1, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(NN.sel1)

## ----prior-sel2, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my.prior <- list(
  "delta.prior.e" = c("norm", 10, 1)
)

## ----mnar-sel2-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.sel2 = selection(data = MenSS, model.eff = e ~ trt + u.0, 
                    model.cost = c ~ trt, model.me = me ~ e, 
                    model.mc = mc ~ 1, type = "MNAR", 
                    n.iter = 1000, ref = 2, prior = my.prior,
                    dist_e = "norm", dist_c = "norm")

## ----mnar-sel2, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.sel2 = selection(data = MenSS, model.eff = e ~ trt + u.0,
#                     model.cost = c ~ trt, model.me = me ~ e,
#                     model.mc = mc ~ 1, type = "MNAR",
#                     n.iter = 1000, ref = 2, prior = my.prior,
#                     dist_e = "norm", dist_c = "norm")

## ----print-sel2, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(NN.sel2)

## ----figplotceamnarsel1, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from model 1 fitted under default MNAR prior specification for the missing effectiveness mechanism.",out.width='100%', fig.pos='h', out.extra=''----
ceac_sm1_mnar <- BCEA::ceac.plot(NN.sel1$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_sm2_mnar <- BCEA::ceac.plot(NN.sel2$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_sm1_mnar

## ----figplotceamnarsel2, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from model 2 fitted under custom MNAR prior specification for the missing effectiveness mechanism.",out.width='100%', fig.pos='h', out.extra=''----
ceac_sm2_mnar

## ----prior-pm2, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my.prior <- list(
  "delta.prior.e" = c("unif", -0.3, -0.2)
)

## ----mnar-pat2-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.pat2 = pattern(data = MenSS, model.eff = e ~ trt + u.0, 
                model.cost = c ~ trt, type = "MNAR", 
                restriction = "CC", n.iter = 1000, 
                ref = 2, prior = my.prior,
                dist_e = "norm", dist_c = "norm")

## ----mnar-pat2, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.pat2 = pattern(data = MenSS, model.eff = e ~ trt + u.0,
#                 model.cost = c ~ trt, type = "MNAR",
#                 restriction = "CC", n.iter = 1000,
#                 ref = 2, prior = my.prior,
#                 dist_e = "norm", dist_c = "norm")

## ----print-pm2, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(NN.pat2)

## ----prior-pm3, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
my.prior2 <- list(
  "delta.prior.e" = c("norm", -0.25, 1/(0.05^2))
)

## ----mnar-pat3-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.pat3 = pattern(data = MenSS, model.eff = e ~ trt + u.0, 
                model.cost = c ~ trt, type = "MNAR", 
                restriction = "CC", n.iter = 1000, 
                ref = 2, prior = my.prior2,
                dist_e = "norm", dist_c = "norm")

## ----mnar-pat3, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.pat3 = pattern(data = MenSS, model.eff = e ~ trt + u.0,
#                 model.cost = c ~ trt, type = "MNAR",
#                 restriction = "CC", n.iter = 1000,
#                 ref = 2, prior = my.prior2,
#                 dist_e = "norm", dist_c = "norm")

## ----figplotceamnarpat1, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from model 1 fitted under uniform MNAR prior specification for the effectiveness sensitivity parameter.",out.width='100%', fig.pos='h', out.extra=''----
ceac_pm2_mnar <- BCEA::ceac.plot(NN.pat2$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_pm3_mnar <- BCEA::ceac.plot(NN.pat3$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_pm2_mnar

## ----figplotceamnarpat2, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from model 2 fitted under normal MNAR prior specification for the effectiveness sensitivity parameter.",out.width='100%', fig.pos='h', out.extra=''----
ceac_pm3_mnar

## ----mar-hm1-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.hur1 = hurdle(data = MenSS, model.eff = e ~ trt + u.0, 
                 model.cost = c ~ trt, model.se = se ~ 1, 
                 model.sc = sc ~ 1, type = "SCAR", 
                 se = 1, sc = 0, n.iter = 1000, ref = 2, 
                 dist_e = "beta", dist_c = "gamma")

## ----mar-hm1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.hur1 = hurdle(data = MenSS, model.eff = e ~ trt + u.0,
#                  model.cost = c ~ trt, model.se = se ~ 1,
#                  model.sc = sc ~ 1, type = "SCAR",
#                  se = 1, sc = 0, n.iter = 1000, ref = 2,
#                  dist_e = "beta", dist_c = "gamma")

## ----d-mar-hm, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
s_e <- ifelse(MenSS$e == 1, 1, 0)

#number of ones
sum(s_e == 1, na.rm = T)

## ----d-mnar-hm, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
mys_e <- ifelse(is.na(s_e) & MenSS$age < 22, 1, s_e)

#number of ones
sum(mys_e == 1, na.rm = T)

## ----mnar-hm2-hide, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
NN.hur2 = hurdle(data = MenSS, model.eff = e ~ trt + u.0, 
                 model.cost = c ~ trt, model.se = se ~ 1, 
                 model.sc = sc ~ 1, type = "SCAR", s_e = mys_e,
                 se = 1, sc = 0, n.iter = 1000, ref = 2, 
                 dist_e = "beta", dist_c = "gamma")

## ----mnar-hm2, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# NN.hur2 = hurdle(data = MenSS, model.eff = e ~ trt + u.0,
#                  model.cost = c ~ trt, model.se = se ~ 1,
#                  model.sc = sc ~ 1, type = "SCAR", s_e = mys_e,
#                  se = 1, sc = 0, n.iter = 1000, ref = 2,
#                  dist_e = "beta", dist_c = "gamma")

## ----print-hm2, echo=TRUE, eval=TRUE, message=FALSE, error=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(NN.hur2)

## ----figplothm2, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Scatter plots of observed and imputed effectiveness data in the SoC (trt=1) arm based on an hurdle model assuming all individuals younger than 22 years old are associated with a perfect healf status.",out.width='100%', fig.pos='h', out.extra=''----
plot(NN.hur2, class = "scatter", outcome = "effects", trt = "SoC")

## ----figplotceamnarhm1, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from the hurdle model fitted under a MAR assumption about the number of structural ones in SoC.",out.width='100%', fig.pos='h', out.extra=''----
ceac_hm1_mnar <- BCEA::ceac.plot(NN.hur1$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_hm2_mnar <- BCEA::ceac.plot(NN.hur2$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_hm1_mnar

## ----figplotceamnarhm2, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Cost-effectiveness acceptability curves (CEAC) derived from the hurdle models fitted under a MNAR assumption about the number of structural ones in SoC.",out.width='100%', fig.pos='h', out.extra=''----
ceac_hm2_mnar

