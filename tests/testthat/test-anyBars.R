# test internal functions

test_that("anyBars() divide input formula into a fixed and random part", {
  expect_equal(anyBars(y ~ 1 | z) , TRUE)
})

test_that("fb() extracts the random effects component from an object of class formula", {
  expect_type(fb(y ~ 1 + (1 | x)), "list")
})

test_that("nobars_() separates the fixed and random effects components from an object of class formula", {
  expect_length(nobars_(y ~ 1), 3)
})

test_that("isAnyArgBar() detect the random effects component from an object of class formula", {
  test_term1 <- quote(foo(x, 1 | group))
  expect_equal(isAnyArgBar(test_term1), TRUE)
})

test_that("isBar() detect the random effects component from an object of class formula", {
  test_term1 <- quote(1 | group)
  expect_equal(isBar(test_term1) , TRUE)
})


test_that("jagsresults() summarise results from a fitted model in JAGS", {

  load(test_path("testdata/sm.RData"))

  jags_sm_obj <- sm$model_output$model
  params_sm <- c("s_e", "s_c", "beta", "alpha", "p_e", "p_c")
  jags_sm <- jagsresults(jags_sm_obj, params = params_sm, regex = FALSE, invert = FALSE,
                         probs = c(0.025, 0.25, 0.5, 0.75, 0.975), signif = 2)
  expect_type(jags_sm, "double")
  expect_shape(jags_sm, ncol = 9)

  load(test_path("testdata/pmm.RData"))

  jags_pmm_obj <- pmm$model_output$model
  params_pmm <- c("s_e_p", "s_c_p", "beta_p", "alpha_p", "p_prob")
  jags_pmm <- jagsresults(jags_pmm_obj, params = params_pmm, regex = FALSE, invert = FALSE,
                         probs = c(0.025, 0.5, 0.975), signif = 1)
  expect_type(jags_pmm, "double")
  expect_shape(jags_pmm, ncol = 7)

})
