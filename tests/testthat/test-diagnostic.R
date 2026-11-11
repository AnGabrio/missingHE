# test diagnostic/plot functions

test_that("output of diagnostic function", {

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
  n.ind <- nrow(MenSS)

  sm <- selection(data = MenSS, model.eff = model.eff, dist_e = dist_e,
                  model.cost = model.cost, dist_c = dist_c, type = type, n.iter = n.iter)

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

  # test for output

  diag_sm <- diagnostic(sm, type = "denplot", param = "alpha", theme = NULL)

  expect_equal(is_ggplot(diag_sm), TRUE)
  expect_length(diag_sm, 1)


})


