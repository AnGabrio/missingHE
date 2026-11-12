# test diagnostic/pic/ppc/plot functions

test_that("output of diagnostic/pic/ppc/plot functions", {

  testthat::skip_on_cran()

  # test for input errors

  load(test_path("testdata/MenSS.RData"))

  model.eff <- e ~ trt
  model.cost <- c ~ trt
  type <- "MAR"
  dist_e <- "norm"
  dist_c <- "norm"
  n.iter <- 10
  n.trt <- length(unique(MenSS$trt))
  n.trt.names <- unique(MenSS$trt)
  n.trt.names.char <- as.character(n.trt.names)
  n.ind <- nrow(MenSS)
  n.cc <- nrow(MenSS[complete.cases(MenSS[,c("e","c")]),])

  sm <- selection(data = MenSS, model.eff = model.eff, dist_e = dist_e,
                  model.cost = model.cost, dist_c = dist_c, type = type, n.iter = n.iter)

  # test for error/check output of diagnostic

  expect_error(
    diagnostic(sm$data_set, type = "denplot", param = "all", theme = NULL),
    regexp = "Only objects of class 'missingHE' can be used")
  expect_error(
    diagnostic(sm, type = "denplot", param = "all", theme = "gg"),
    regexp = "Please provide a valid theme type among those available")
  expect_error(
    diagnostic(sm, type = "denplot", param = c("alpha", "beta"), theme = NULL),
    regexp = "You can only visualise diagnostic checks for one family of parameters at a time
         or for all parameters together by setting param = 'all'")
  expect_error(
    diagnostic(sm, type = "denplots", param = "all", theme = NULL),
    regexp = "Types of diagnostics available for use are: 'histogram', 'running', 'denplot', 'compare',
         'traceplot', 'acf', 'cross', 'Rhat', 'geweke', 'caterpillar', 'pairs'")
  expect_error(
    diagnostic(sm, type = "denplot", param = "feta", theme = NULL),
    regexp = "Please provide valid parameter names based on the fitted model")
  expect_error(
    diagnostic(sm, type = "denplot", param = "random.alpha", theme = NULL),
    regexp = "no random effects for effects found")
  expect_error(
    diagnostic(sm, type = "denplot", param = "random.beta", theme = NULL),
    regexp = "no random effects for costs found")
  expect_error(
    diagnostic(sm, type = "denplot", param = "random.beta.f", theme = NULL),
    regexp = "no random effects for costs found related effects predictor")
  expect_error(
    diagnostic(sm, type = "denplot", param = "random.gamma.e", theme = NULL),
    regexp = "no random effects for missing effects found")
  expect_error(
    diagnostic(sm, type = "denplot", param = "random.gamma.c", theme = NULL),
    regexp = "no random effects for missing costs found")
  expect_error(
    diagnostic(sm, type = "denplot", param = "beta.f", theme = NULL),
    regexp = "no fixed effects for costs found related effects predictor")

  diag_sm <- diagnostic(sm, type = "denplot", param = "alpha", theme = NULL)

  expect_equal(is_ggplot(diag_sm), TRUE)
  expect_length(diag_sm, 1)

  # test for error/check output of pic

  sm2 <- selection(data = MenSS, model.eff = model.eff, dist_e = "logis",
                  model.cost = model.cost, dist_c = dist_c, type = type, n.iter = n.iter)

  expect_error(
    pic(sm, criterion = 2, cases = "cc"),
    regexp = "Criterion must be a character name")
  expect_error(
    pic(sm, criterion = "cc", cases = "cc"),
    regexp = "Please, select a character name for 'criterion' among those available")
  expect_error(
    pic(sm, criterion = "waic", cases = "cas"),
    regexp = "Please, select a character name for 'cases' among those available")
  expect_error(
    pic(c(1,2), criterion = "waic", cases = "cc"),
    regexp = "Only an object or list of objects of class 'missingHE' can be used")
  expect_error(
    pic(sm$data_set, criterion = "waic", cases = "cc"),
    regexp = "Only a list of objects of class 'missingHE' can be used with same 'method' and 'data_format' argument values")

  dic_sm <- pic(sm, criterion = "dic", cases = "cc")
  looic_sm <- pic(sm, criterion = "looic", cases = "cc", save_psis = TRUE)
  sm_list <- list(sm, sm2)
  n.sm <- length(sm_list)
  waic_sm_list <- pic(sm_list, criterion = "waic", cases = "cc")

  expect_type(dic_sm, "list")
  expect_type(looic_sm, "list")
  expect_type(waic_sm_list, "list")
  expect_length(dic_sm, 4)
  expect_length(looic_sm, 9)
  expect_length(waic_sm_list, 5)
  expect_named(dic_sm, c("d_bar", "pD", "dic", "d_hat"))
  expect_named(looic_sm, c("elpd", "elpd_se", "p", "p_se", "looic", "looic_se", "pointwise", "diagnostics", "psis_object"))
  expect_named(waic_sm_list, c("pic", "pic_se", "peff", "peff_se", "ic_res"))
  expect_vector(dic_sm$d_bar)
  expect_length(dic_sm$d_bar, 1)
  expect_vector(dic_sm$pD)
  expect_length(dic_sm$pD, 1)
  expect_vector(dic_sm$dic)
  expect_length(dic_sm$dic, 1)
  expect_vector(dic_sm$d_hat)
  expect_length(dic_sm$d_hat, 1)

  expect_vector(looic_sm$elpd)
  expect_length(looic_sm$elpd, 1)
  expect_vector(looic_sm$elpd_se)
  expect_length(looic_sm$elpd_se, 1)
  expect_vector(looic_sm$p)
  expect_length(looic_sm$p, 1)
  expect_vector(looic_sm$p_se)
  expect_length(looic_sm$p_se, 1)
  expect_vector(looic_sm$looic)
  expect_length(looic_sm$looic, 1)
  expect_vector(looic_sm$looic_se)
  expect_length(looic_sm$looic_se, 1)
  expect_shape(looic_sm$pointwise, nrow = n.cc)
  expect_equal(colnames(looic_sm$pointwise), c("elpd_loo", "mcse_elpd_loo",
                                               "p_loo", "looic", "influence_pareto_k"))
  expect_type(looic_sm$pointwise, "double")
  expect_type(looic_sm$diagnostics, "list")
  expect_named(looic_sm$diagnostics, c("pareto_k", "r_eff", "n_eff"), ignore.order = TRUE)
  expect_length(looic_sm$diagnostics, 3)
  expect_vector(looic_sm$diagnostics$r_eff)
  expect_length(looic_sm$diagnostics$r_eff, n.cc)
  expect_vector(looic_sm$diagnostics$n_eff)
  expect_length(looic_sm$diagnostics$n_eff, n.cc)
  expect_vector(looic_sm$diagnostics$pareto_k)
  expect_length(looic_sm$diagnostics$pareto_k, n.cc)
  expect_type(looic_sm$psis_object, "list")
  expect_named(looic_sm$psis_object, c("log_weights", "diagnostics"))
  expect_length(looic_sm$psis_object, 2)
  expect_shape(looic_sm$psis_object$log_weights, dim = c(n.iter, n.cc))
  expect_type(looic_sm$psis_object$log_weights, "double")
  expect_length(looic_sm$psis_object$diagnostics, 3)
  expect_named(looic_sm$psis_object$diagnostics, c("pareto_k", "r_eff", "n_eff"), ignore.order = TRUE)

  expect_type(waic_sm_list, "list")
  expect_length(waic_sm_list, 5)
  expect_named(waic_sm_list, c("pic", "pic_se", "peff", "peff_se", "ic_res"))
  expect_length(waic_sm_list$pic, n.sm)
  expect_vector(waic_sm_list$pic)
  expect_length(waic_sm_list$pic_se, n.sm)
  expect_vector(waic_sm_list$pic_se)
  expect_length(waic_sm_list$peff, n.sm)
  expect_vector(waic_sm_list$peff)
  expect_length(waic_sm_list$peff_se, n.sm)
  expect_vector(waic_sm_list$peff_se)
  expect_type(waic_sm_list$ic_res, "list")
  expect_length(waic_sm_list$ic_res, 2)


  # test for error/check output of plot

  expect_error(
    plot(sm, prob = c(0.025, 0.975),
        class = "scatter", outcome = "both", trt = "all",
        time_plot = NULL, theme = "ff", only.plot = TRUE),
    regexp = "Please provide a valid theme type among those available")
  expect_error(
    plot(sm, prob = c(0.025, 0.975),
         class = "scatter", outcome = "both", trt = "all",
         time_plot = NULL, theme = NULL, only.plot = 1),
    regexp = "Please provide a logical value to the 'only.plot' argument")
  expect_error(
    plot(sm, prob = c(0.025, 0.975),
         class = "dr", outcome = "both", trt = "all",
         time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Please provide a plot class name among those available")
  expect_error(
    plot(sm, prob = c(0.025, 0.975),
         class = "scatter", outcome = "fg", trt = "all",
         time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Please provide an outcome name among those available")
  expect_error(
    plot.missingHE(sm$data_set, prob = c(0.025, 0.975),
         class = "scatter", outcome = "both", trt = "all",
         time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Only objects of class 'missingHE' can be used")
  expect_error(
    plot.missingHE(sm, prob = 0.025,
                   class = "scatter", outcome = "both", trt = "all",
                   time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Please provide a lower and an upper quantile for the imputed data distribution")
  expect_error(
    plot.missingHE(sm, prob = c(0.025, 0.975),
                   class = "scatter", outcome = "both", trt = "gg",
                   time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Please provide a trt name among those available")
  expect_error(
    plot.missingHE(sm, prob = c(0.025, 0.975),
                   class = "scatter", outcome = "both", trt = 4,
                   time_plot = NULL, theme = NULL, only.plot = TRUE),
    regexp = "Please provide a trt value among those available")


  plot.sm <- plot(sm, prob = c(0.025, 0.975),
                      class = "scatter", outcome = "both", trt = "all",
                      time_plot = NULL, theme = NULL, only.plot = TRUE)
  noplot.sm <- plot(sm, prob = c(0.025, 0.975),
                  class = "scatter", outcome = "both", trt = "all",
                  time_plot = NULL, theme = NULL, only.plot = FALSE)

  expect_equal(is_ggplot(plot.sm), TRUE)
  expect_length(plot.sm, 1)
  expect_type(noplot.sm, "list")
  expect_length(noplot.sm, 2)
  expect_named(noplot.sm, c("plot", "list.data"))
  expect_equal(is_ggplot(noplot.sm$plot), TRUE)
  expect_length(noplot.sm$plot, 1)
  expect_type(noplot.sm$list.data, "list")
  expect_length(noplot.sm$list.data, 3)
  expect_named(noplot.sm$list.data, c("complete", "observed", "imputed"))
  expect_type(noplot.sm$list.data$complete, "list")
  expect_type(noplot.sm$list.data$observed, "list")
  expect_type(noplot.sm$list.data$imputed, "list")
  expect_length(noplot.sm$list.data$complete, 2)
  expect_length(noplot.sm$list.data$observed, 2)
  expect_length(noplot.sm$list.data$imputed, 2)
  expect_named(noplot.sm$list.data$complete, c("effects", "costs"))
  expect_named(noplot.sm$list.data$observed, c("effects", "costs"))
  expect_named(noplot.sm$list.data$imputed, c("effects", "costs"))
  expect_type(noplot.sm$list.data$complete$effects, "list")
  expect_type(noplot.sm$list.data$observed$effects, "list")
  expect_type(noplot.sm$list.data$imputed$effects, "list")
  expect_type(noplot.sm$list.data$complete$costs, "list")
  expect_type(noplot.sm$list.data$observed$costs, "list")
  expect_type(noplot.sm$list.data$imputed$costs, "list")
  expect_length(noplot.sm$list.data$complete$effects, 2)
  expect_length(noplot.sm$list.data$observed$effects, 2)
  expect_length(noplot.sm$list.data$imputed$effects, 2)
  expect_named(noplot.sm$list.data$complete$effects, n.trt.names.char)
  expect_named(noplot.sm$list.data$observed$effects, n.trt.names.char)
  expect_named(noplot.sm$list.data$imputed$effects, n.trt.names.char)
  expect_length(noplot.sm$list.data$complete$costs, 2)
  expect_length(noplot.sm$list.data$observed$costs, 2)
  expect_length(noplot.sm$list.data$imputed$costs, 2)
  expect_named(noplot.sm$list.data$complete$costs, n.trt.names.char)
  expect_named(noplot.sm$list.data$observed$costs, n.trt.names.char)
  expect_named(noplot.sm$list.data$imputed$costs, n.trt.names.char)
  expect_type(noplot.sm$list.data$complete$effects$`1`, "list")
  expect_type(noplot.sm$list.data$complete$effects$`2`, "list")
  expect_shape(noplot.sm$list.data$complete$effects$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$complete$effects$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$complete$effects$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$complete$effects$`2`, c("Mean", "QL", "QU"))
  expect_type(noplot.sm$list.data$observed$effects$`1`, "list")
  expect_type(noplot.sm$list.data$observed$effects$`2`, "list")
  expect_shape(noplot.sm$list.data$observed$effects$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$observed$effects$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$observed$effects$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$observed$effects$`2`, c("Mean", "QL", "QU"))
  expect_type(noplot.sm$list.data$imputed$effects$`1`, "list")
  expect_type(noplot.sm$list.data$imputed$effects$`2`, "list")
  expect_shape(noplot.sm$list.data$imputed$effects$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$imputed$effects$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$imputed$effects$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$imputed$effects$`2`, c("Mean", "QL", "QU"))
  expect_type(noplot.sm$list.data$complete$costs$`1`, "list")
  expect_type(noplot.sm$list.data$complete$costs$`2`, "list")
  expect_shape(noplot.sm$list.data$complete$costs$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$complete$costs$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$complete$costs$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$complete$costs$`2`, c("Mean", "QL", "QU"))
  expect_type(noplot.sm$list.data$observed$costs$`1`, "list")
  expect_type(noplot.sm$list.data$observed$costs$`2`, "list")
  expect_shape(noplot.sm$list.data$observed$costs$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$observed$costs$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$observed$costs$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$observed$costs$`2`, c("Mean", "QL", "QU"))
  expect_type(noplot.sm$list.data$imputed$costs$`1`, "list")
  expect_type(noplot.sm$list.data$imputed$costs$`2`, "list")
  expect_shape(noplot.sm$list.data$imputed$costs$`1`, ncol = 3)
  expect_shape(noplot.sm$list.data$imputed$costs$`2`, ncol = 3)
  expect_named(noplot.sm$list.data$imputed$costs$`1`, c("Mean", "QL", "QU"))
  expect_named(noplot.sm$list.data$imputed$costs$`2`, c("Mean", "QL", "QU"))


  # test for error/check output of ppc

  expect_error(
    ppc(sm$data_set, type = "histogram", outcome = "both", ndisplay = 5,
         trt = "all", theme = NULL, scheme_set = NULL),
    regexp = "Only objects of class 'missingHE' can be used")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = 5,
        trt = "all", theme = "fg", scheme_set = NULL),
    regexp = "Please provide a valid theme type among those available")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = "sd"),
    regexp = "Please provide a valid scheme set among those available")
  expect_error(
    ppc(sm, type = "sd", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = "Please provide a plot type among those available")
  expect_error(
    ppc(sm, type = "histogram", outcome = "sd", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = "Please provide an outcome name among those available")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL, corr = 2),
    regexp = "The argument 'corr' must be a function and is available only when outcome is set to 'both'")
  expect_error(
    ppc(sm, type = "stat_2d", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL, bpp = TRUE),
    regexp = "Bayesian p-values are not available for 'stat_2d' type")
  expect_error(
    ppc(sm, type = "stat", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL, bpp = TRUE),
    regexp = "Please use the argument 'stat' to provide a single function for computing statistics")
  expect_error(
    ppc(sm, type = "stat_2d", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL, bpp = FALSE),
    regexp = "Please use the argument 'stat' to provide a vector of two functions for computing statistics")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = 5,
        trt = 5, theme = NULL, scheme_set = NULL),
    regexp = "Please provide a trt value among those available")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = 5,
        trt = "sd", theme = NULL, scheme_set = NULL),
    regexp = "Please provide a trt name among those available")
  expect_error(
    ppc(sm, type = "histogram", outcome = "both", ndisplay = -2,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = "Number of replications to display is not valid")
  expect_error(
    ppc(sm, type = "loo_pit_overlay", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = NULL)
  expect_error(
    ppc(sm, type = "loo_pit_ecdf", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = NULL)
  expect_error(
    ppc(sm, type = "loo_pit_ecdf", outcome = "both", ndisplay = 5,
        trt = "all", theme = NULL, scheme_set = NULL),
    regexp = NULL)

  ppc_sm <- ppc(sm, type = "histogram", outcome = "both", ndisplay = 3,
                    trt = "all", theme = NULL, scheme_set = NULL)

  expect_equal(is_ggplot(ppc_sm), TRUE)
  expect_length(ppc_sm, 1)

})


