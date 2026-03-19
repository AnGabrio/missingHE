## ----echo = FALSE, message = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
library(bookdown)
library(ggplot2)

options(width = 300)
set.seed(1014)

## ----seldist, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#rename trt levels
MenSS$trt <- factor(MenSS$trt)
levels(MenSS$trt) <- c("SoC", "MenSS")
MenSS$e <- MenSS$e - 0.01 #ensure no ones QALYs occur
MenSS$c <- MenSS$c + 0.01 #ensure no zero costs

#fit models with different distributions for outcomes
#1=Normal-Normal
sm1_nn <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                    model.eff = e ~ trt + u.0, model.cost = c ~ trt + e,
                    model.me = me ~ 1, model.mc = mc ~ 1,
                    type = "MAR", n.iter = 1000, ref = 2)
#2=Normal-Gamma
sm1_ng <- selection(data = MenSS, dist_e = "norm", dist_c = "gamma",
                    model.eff = e ~ trt + u.0, model.cost = c ~ trt + e,
                    model.me = me ~ 1, model.mc = mc ~ 1,
                    type = "MAR", n.iter = 1000, ref = 2)
#3=Beta-Gamma
sm1_bg <- selection(data = MenSS, dist_e = "beta", dist_c = "gamma",
                    model.eff = e ~ trt + u.0, model.cost = c ~ trt + e,
                    model.me = me ~ 1, model.mc = mc ~ 1,
                    type = "MAR", n.iter = 1000, ref = 2)

## ----selpic, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#estimate looic for all models based on complete cases
looic_m123 <- pic(x = list(sm1_nn, sm1_ng, sm1_bg),
                  criterion = "looic", cases = "cc")
#print criteria for each model
looic_m123$pic

## ----selpeff, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# #print peff values
# looic_m123$peff

## ----selppc, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#compare densities of observed vs replicated outcome data under each model
ppc_dens_m1 <- ppc(x = sm1_nn, type = "dens_overlay", outcome = "both",
                   ndisplay = 20, trt = "none")
ppc_dens_m2 <- ppc(x = sm1_ng, type = "dens_overlay", outcome = "costs",
                   ndisplay = 20, trt = "none")
ppc_dens_m3 <- ppc(x = sm1_bg, type = "dens_overlay", outcome = "effects",
                   ndisplay = 20, trt = "none")

## ----figplotppc1, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Posterior densities of 20 replicated data under model 1 compared to the empirical density of the original data.", out.width='100%', fig.pos='h', out.extra=''---------------------------
ppc_dens_m1

## ----figplotppc2, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Posterior densities of 20 replicated data under model 2 compared to the empirical density of the original data.", out.width='100%', fig.pos='h', out.extra=''---------------------------
ppc_dens_m2

## ----figplotppc3, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Posterior densities of 20 replicated data under model 3 compared to the empirical density of the original data.", out.width='100%', fig.pos='h', out.extra=''---------------------------
ppc_dens_m3

## ----selmar, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#fit selection model 1 (Normal-Normal) under MAR conditional on u.0 and age
sm1_nn_mar <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                    model.eff = e ~ trt + u.0, model.cost = c ~ trt + e,
                    model.me = me ~ u.0 + age, model.mc = mc ~ u.0 + age,
                    type = "MAR", n.iter = 1000, ref = 2)

## ----selprint, echo=TRUE, eval=TRUE, tidy=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#print posterior summaries for fixed effects
print(x = sm1_nn_mar, display = "fixed")

## ----selpe, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#store posterior samples of p_e and compute some custom summaries
p_e <- sm1_nn_mar$model_output$model$BUGSoutput$sims.list$p_e
summary(c(p_e))

## ----selcoef, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#print coefficient estimates from model.eff and model.cost
coef(sm1_nn_mar, random = FALSE, digits = 2)

## ----selcov, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#fit selection model 1 (Normal-Normal) with covariates into both outcome and missingness models
sm1_nn_cov <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                    model.eff = e ~ trt + u.0 + age, 
                    model.cost = c ~ trt + age + employment + e,
                    model.me = me ~ u.0 + age, model.mc = mc ~ age + employment,
                    type = "MAR", n.iter = 1000, ref = 2)

## ----selcovre, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#fit selection model 1 (Normal-Normal) with covariates and random intercepts
#into the model
sm1_nn_cov_re <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                    model.eff = e ~ trt + u.0 + age + (1 | site), 
                    model.cost = c ~ trt + age + employment + e + (1 | site),
                    model.me = me ~ u.0 + age, model.mc = mc ~ age + employment,
                    type = "MAR", n.iter = 1000, ref = 2)

## ----figplotdiag, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Checking convergence using the diagnostic function for a family of model parameters estimated in missingHE, for example through inspection of the autocorrelation plots.", out.width='100%', fig.pos='h', out.extra=''----
diagnostic(x = sm1_nn_cov, type = "acf", param = "alpha")

## ----selprior, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#update priors on outcome regression coefficients and standard deviations
myprior <- list("beta.prior" = c("norm", 0, 0.01), 
                "beta_f.prior" = c("norm", 0, 0.001),
                "alpha.prior" = c("norm", 0, 0.01),
                "sigma.prior.e" = c("unif", 0, 5),
                "sigma.prior.c" = c("unif", 0, 1000))
#update model with new priors
sm2_nn_cov <- selection(data = MenSS, dist_e = "norm", dist_c = "norm",
                    model.eff = e ~ trt + u.0 + age, 
                    model.cost = c ~ trt + age + employment + e,
                    model.me = me ~ u.0 + age, model.mc = mc ~ age + employment,
                    type = "MAR", n.iter = 1000, ref = 2,
                    prior = myprior)

