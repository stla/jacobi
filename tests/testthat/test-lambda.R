test_that("lambda in terms of the 'Weierstrass roots'", {
  omega1 <- 1.4 - 1i
  omega2 <- 1.6 + 0.5i
  omega <- c(omega1, omega2)
  e1 <- pweierstrass(omega1, omega = omega)
  e2 <- pweierstrass(omega2, omega = omega)
  e3 <- pweierstrass(omega1+omega2, omega = omega)
  tau <- omega2 / omega1
  expect_equal(lambda(tau), (e3-e2)/(e1-e2))
})
