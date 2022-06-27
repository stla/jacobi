test_that("wsigma values", {
  expect_equal(
    wsigma(1, g = c(12, -8)), 
    sin(1i*sqrt(3))/(1i*sqrt(3)*sqrt(exp(1)))
  )
  expect_equal(
    wsigma(2, g = c(1, 2i)), 
    1.8646253716-0.3066001355i
  )
})
