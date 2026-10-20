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



