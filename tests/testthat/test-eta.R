test_that("eta values", {
  expect_equal(eta(1i/2), as.complex(gamma(1/4) / 2^(7/8) / pi^(3/4)))
  expect_equal(eta(2i), as.complex(gamma(1/4) / 2^(11/8) / pi^(3/4)))
})

test_that("eta is correctly vectorized", {
  tau1 <- 1 + 1i
  tau2 <- -2 + 3i
  tau <- c(tau1, tau2)
  expect_equal(eta(tau), c(eta(tau1), eta(tau2)))
})

test_that("Alternative expression of Dedekind eta", {
  tau <- 0.2 + 0.2i
  expect_equal(
    eta(tau),
    -1i * exp(1i*pi*tau/3) * jtheta1(tau*pi, 3*tau)
  )
})

test_that("Other alternative expression of Dedekind eta", {
  tau <- 9.2 + 9.2i
  qbar <- exp(2i*pi*tau)
  expect_equal(
    eta(tau),
    jtheta2(pi/6, tau = tau/3) / sqrt(3)
  )
})
