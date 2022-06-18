test_that("Sum of the e_i is zero.", {
  omega1 <- 1.4 - 1i
  omega2 <- 1.6 + 0.5i
  omega <- c(omega1, omega2)
  e1 <- wp(omega1, omega = omega)
  e2 <- wp(omega2, omega = omega)
  e3 <- wp(-omega1-omega2, omega = omega)
  expect_equal(e1 + e2 + e3, 0i)
})

test_that("Differential equation.", {
  z <- 1 + 1i
  g2 <- 5 + 3i
  g3 <- 2 + 7i
  g <- c(g2, g3)
  pw <- wp(z, g)
  pwprimesquared <- wp(z, g, derivative = 1)^2
  expect_equal(pwprimesquared, 4*pw**3 - g2*pw - g3)
})


