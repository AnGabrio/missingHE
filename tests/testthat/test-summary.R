# test summary function

test_that("output of summary function", {

  load(test_path("testdata/sm.RData"))

  expect_output(summary(sm, incremental = FALSE, prob = c(0.025, 0.975),
                        digits = 3),
                regexp = "Cost-effectiveness analysis summary")
  expect_output(summary(sm, incremental = TRUE, prob = c(0.025, 0.975),
                        digits = 3),
                regexp = "Cost-effectiveness analysis summary")
  expect_output(summary(sm, incremental = FALSE, prob = c(0.1, 0.9),
                        digits = 3),
                regexp = "Cost-effectiveness analysis summary")
  expect_output(summary(sm, incremental = FALSE, prob = c(0.025, 0.975),
                        digits = 1),
                regexp = "Cost-effectiveness analysis summary")

  expect_error(
    summary.missingHE(sm$data_set),
    regexp = "Only objects of class 'missingHE' can be used")

  expect_error(
    summary(sm, incremental = 1),
    regexp = "Please provide 'incremental' as a logical value")

  expect_error(
    summary(sm, prob = c(-0.2, 0.6)),
    regexp = "Please provide a lower and an upper quantile for the posterior results")

  expect_error(
    summary(sm, ref = "trt"),
    regexp = "Please provide a single numeric indicator for the reference treatment")

  expect_error(
    summary(sm, ref = 0),
    regexp = "Please provide a valid indicator value for the reference treatment")

  expect_error(
    summary(sm, wtp = "wtp"),
    regexp = "Please provide a single numeric value for the acceptance threshold")

  expect_error(
    summary(sm, wtp = -2),
    regexp = "Please provide a valid value for the acceptance threshold")

})

# test postsum function

test_that("output of postsum function", {

  load(test_path("testdata/pmm.RData"))

  expect_output(postsum(pmm, display = "fixed", digits = 3),
                regexp = "mean        sd     2.5%       50%    97.5%   Rhat n.eff")
  expect_output(postsum(pmm, display = "loglik", digits = 3),
                regexp = "mean      sd     2.5%     50%   97.5%  Rhat n.eff")

  expect_error(
    postsum.missingHE(pmm$data_set),
    regexp = "Only objects of class 'missingHE' can be used")

  expect_error(
    postsum(pmm, display = "no"),
    regexp = "Please select a valid input for 'display' argument")

})


# test print function

test_that("output of print function", {

  load(test_path("testdata/hm.RData"))

  expect_output(print(hm, digits = 3),
                regexp = "Call:
hurdle")
  expect_output(print(hm, digits = 1),
                regexp = "Call:
hurdle")

  expect_error(
    print.missingHE(hm$data_set),
    regexp = "Only objects of class 'missingHE' can be used")

})


# test coef function

test_that("output of coef function", {

  load(test_path("testdata/lmm.RData"))

  expect_output(coef(lmm, prob = c(0.025, 0.975),
                     random = FALSE, digits = 3),
                regexp = "Mean     SD      QL     QU")
  expect_output(coef(lmm, prob = c(0.1, 0.9),
                     random = FALSE, digits = 3),
                regexp = "Mean     SD      QL     QU")
  expect_output(coef(lmm, prob = c(0.025, 0.975),
                     random = FALSE, digits = 1),
                regexp = "Mean   SD    QL   QU")

  expect_error(
    coef.missingHE(lmm$data_set),
    regexp = "Only objects of class 'missingHE' can be used")
  expect_error(
    coef(lmm, prob = c(0.025, 1.2)),
    regexp = "You must provide valid lower/upper quantiles for the coefficients distributions")
  expect_error(
    coef(lmm, random = TRUE),
    regexp = "No random effects estimates found")

})
