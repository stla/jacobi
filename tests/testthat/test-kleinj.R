test_that("kleinj at 2i", {
  expect_equal(kleinj(2i), as.complex(66^3))
})

test_that("kleinj alternate expression", {
  tau <- 0.5 + 0.3i
  lambda <- (jtheta2(0, tau = tau) / jtheta3(0, tau = tau))^4 
  x <- lambda * (1- lambda)
  expect_equal(kleinj(tau), 256*(1-x)^3 / x^2)
})
