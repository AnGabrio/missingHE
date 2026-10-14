test_that("anyBars() divide input formula into a fixed and random part", {
  expect_equal(anyBars(y ~ 1 | z) , TRUE)
})
