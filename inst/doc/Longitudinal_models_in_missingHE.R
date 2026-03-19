## ----echo = FALSE, message = FALSE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(prompt = TRUE, highlight = F, background = '#FFFFFF',
                      collapse = T, comment = "#>")
library(missingHE)
library(bookdown)
library(ggplot2)

options(width = 300)
set.seed(1234)

## ----pbs, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
n_long <- nrow(PBS)
n <- nrow(PBS[PBS$time == 1, ])
n1 <- nrow(PBS[PBS$trt=="1" & PBS$time == 1,])
n2 <- nrow(PBS[PBS$trt=="2" & PBS$time == 1,])
PBS_outcomes <- PBS[,c("e","c","trt","time")]
PBS_outcomes$time <- factor(PBS_outcomes$time)
PBS_outcomes$trt <- factor(PBS_outcomes$trt)
levels(PBS_outcomes$trt) <- c("TAU", "PBS")
PBS_ec <- PBS[,c("id","e","c","trt","time")]
PBS_ec_cc <- PBS_ec[complete.cases(PBS_ec),]
PBS_ec_cc_unique <- PBS_ec_cc[unique(PBS_ec_cc$id),]
n_ec_cc <- nrow(PBS_ec_cc_unique)
n1_ec_cc <- nrow(PBS_ec_cc_unique[PBS_ec_cc_unique$trt==1,])
n2_ec_cc <- nrow(PBS_ec_cc_unique[PBS_ec_cc_unique$trt==2,])
n_long_ec_cc <- nrow(PBS_ec_cc)
n1_long_ec_cc <- nrow(PBS_ec_cc[PBS_ec_cc$trt==1,])
n2_long_ec_cc <- nrow(PBS_ec_cc[PBS_ec_cc$trt==2,])
n1_ec_cc_t1 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="1" & PBS_ec_cc$time==1,])
n2_ec_cc_t1 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="2" & PBS_ec_cc$time==1,])
n1_ec_cc_t2 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="1" & PBS_ec_cc$time==2,])
n2_ec_cc_t2 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="2" & PBS_ec_cc$time==2,])
n1_ec_cc_t3 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="1" & PBS_ec_cc$time==3,])
n2_ec_cc_t3 <- nrow(PBS_ec_cc[PBS_ec_cc$trt=="2" & PBS_ec_cc$time==3,])
n_ec_cc_t1 <- n1_ec_cc_t1 + n2_ec_cc_t1
n_ec_cc_t2 <- n1_ec_cc_t2 + n2_ec_cc_t2
n_ec_cc_t3 <- n1_ec_cc_t3 + n2_ec_cc_t3

hist_e <- ggplot(PBS_outcomes, aes(x=e, fill=trt)) + geom_histogram(color="black", binwidth = 0.05) + facet_wrap(~trt+time, labeller=label_parsed) + scale_fill_manual(values=c("grey","grey")) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(panel.grid.major = element_blank(), legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position="none",
axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+ ylab("Frequency") + xlab("utilities") + scale_x_continuous(breaks=seq(-0.5,1,0.5))
hist_c <- ggplot(PBS_outcomes, aes(x=c, fill=trt)) + geom_histogram(color="black", binwidth = 2000) + facet_wrap(~trt+time) + scale_fill_manual(values=c("grey","grey")) + scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + theme(panel.grid.major = element_blank(), legend.key = element_rect(fill = "white"), panel.grid.minor = element_blank(), panel.background = element_blank(), legend.position="none",
axis.line = element_line(colour = "black"), plot.title = element_text(hjust = 0.5))+ ylab("Frequency") + xlab("costs") + scale_x_continuous(breaks=seq(0,40000,40000))

## ----tbl1-hide, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# rbind(head(PBS), tail(PBS))

## ----tbl1, echo=FALSE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(rbind(head(PBS), tail(PBS)), row.names = FALSE)

## ----fig1e, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Histograms of the observed data distributions for the utilities in the TAU and PBS arm by time point.", out.width='100%', fig.pos='h', out.extra=''-------------------------------------------
hist_e

## ----fig1c, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Histograms of the observed data distributions for the costs in the TAU and PBS arm by time point.", out.width='100%', fig.pos='h', out.extra=''-----------------------------------------------
hist_c

## ----lmdm1, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
lmdm1_mar <- lmdm(data = PBS, dist_e = "norm", dist_c = "norm",
             model.eff = e ~ trt, model.cost = c ~ trt + e,
             model.me = me ~ 1, model.mc = mc ~ 1, time_dep = "biv",
             type = "MAR", n.iter = 1000, ref = 2)

## ----diag1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, tidy=TRUE-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# diagnostic(x = lmdm1_mar, type = "denplot", param = "beta.f")

## ----fig2, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Checking convergence using the diagnostic function for a model fitted in missingHE, for example through inspection of the density plots.", out.width='100%', fig.pos='h', out.extra=''---------
diagnostic(x = lmdm1_mar, type = "denplot", param = "beta.f")

## ----ppc1, echo=TRUE, eval=FALSE, message=FALSE, warning=FALSE, error=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# ppc(x = lmdm1_mar, type = "histogram", outcome = "effects",
#     ndisplay = 2, time.plot = 3, trt = "2")

## ----fig3, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Posterior predictive graphical checks using the ppc function for a model fitted in missingHE, for example through inspection of histograms.", out.width='100%', fig.pos='h', out.extra=''------
ppc(x = lmdm1_mar, type = "histogram", outcome = "effects", 
    ndisplay = 2, time.plot = 3, trt = "2")

## ----pic1, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE, error=FALSE, results='hide'----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pic_lmdm1_mar <- pic(x = lmdm1_mar, criterion = "looic", cases = "cc")

## ----pic1show, echo=TRUE, eval=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
pic_lmdm1_mar$looic

## ----pic2, echo=TRUE, eval=FALSE, message=FALSE, error=FALSE, warning=FALSE, results='hide'---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# pic_compare_mar <- pic(x = list(lmdm1_mar, lmdm2_mar),
#                        criterion = "looic", cases = "cc")

## ----print, echo=TRUE, eval=FALSE, tidy=TRUE--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# print(lmdm1_mar, display = "fixed")

## ----printshow, echo=FALSE, eval=TRUE, tidy=TRUE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
print(lmdm1_mar, display = "fixed")

## ----coef, echo=TRUE, eval=FALSE, tidy=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# coef(lmdm1_mar)

## ----coefshow, echo=FALSE, eval=TRUE, tidy=TRUE-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
coef(lmdm1_mar)

## ----plot1, echo=TRUE, eval=FALSE, tidy=TRUE, error=FALSE, warning=FALSE----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# plot(lmdm1_mar, class = "boxplot", outcome = "costs",
#      trt = "1", time.plot = 2)

## ----figplot, echo=FALSE, eval=TRUE, tidy=TRUE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Boxplots of observed and imputed effectiveness data in the TAU arm at time 2 under normal distributions using a longitudinal model.", out.width='100%', fig.pos='h', out.extra=''-----------
plot(lmdm1_mar, class = "boxplot", outcome = "costs", 
     trt = "1", time.plot = 2)

## ----summary, echo=TRUE, eval=FALSE, tidy=TRUE------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# summary(lmdm1_mar, incremental = TRUE)

## ----bcea, echo=TRUE, eval=FALSE, tidy=TRUE---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# BCEA::ceplane.plot(lmdm1_mar$cea, graph = "ggplot2")
# BCEA::ceac.plot(lmdm1_mar$cea, graph = "ggplot2")

## ----figplotcea1, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Inspecting probabilistic sensitivity analysis using BCEA built-in functions, for example in terms of cost-effectiveness plane.",out.width='100%', fig.pos='h', out.extra=''------------
cep_sm1_mar <- BCEA::ceplane.plot(lmdm1_mar$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
ceac_sm1_mar <- BCEA::ceac.plot(lmdm1_mar$cea, 
            graph = "ggplot2") + ggplot2::ggtitle("")
cep_sm1_mar

## ----figplotcea2, echo=FALSE, eval=TRUE, tidy=FALSE, error=FALSE, warning=FALSE, dpi=300, fig.show='hold',fig.cap="Inspecting probabilistic sensitivity analysis using BCEA built-in functions, for example in terms of cost-effectiveness acceptability curve.",out.width='100%', fig.pos='h', out.extra=''----
ceac_sm1_mar

