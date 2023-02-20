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
  expect_equal(
    pwprimesquared, 
    4*pw**3 - g2*pw - g3
  )
})

test_that("wpprime value.", {
  z <- 0.1 + 0.1i
  g2 <- 2 + 1i
  g3 <- 2 - 1i
  g <- c(g2, g3)
  wpprime <- wp(z, g, derivative = 1)
  expect_equal(
    wpprime,
    500.009714501399 + 500.03085572727i,
    tolerance = 1e-4
  )
})

test_that("Equianharmonic case.", {
  omega2 <- gamma(1/3)^3 / 4 / pi
  z0 <- omega2 * (1 + 1i/sqrt(3))
  expect_equal(
    wp(z0, g = c(0, 1)), 0i
  )
  expect_equal(
    wp(Conj(z0), g = c(0, 1)), 0i
  )
})