test_that("kleinj at 2i", {
  expect_equal(
    kleinj(2i), 
    as.complex(66^3)
  )
})

test_that("kleinj alternate expression", {
  tau <- 0.5 + 0.3i
  j2 <- jtheta2_cpp(0, tau)
  j3 <- jtheta3_cpp(0, tau)
  g2 <- 4/3 * (pi/2)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  g3 <- 8/27 * (pi/2)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  ) + j3**12)
  g2cube <- g2*g2*g2
  expect_equal(
    kleinj(tau), 
    1728 * g2cube / (g2cube - 27*g3*g3)
  )
})

test_that("kleinj inverse", {
  tau <- 0.5 + 0.3i
  j <- kleinj(tau)
  x <- polyroot(c(256, -3*256, 3*256-j, -256))[1L]
  lbd <- polyroot(c(-x, 1, -1))[1L]
  tau <- 1i*agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd))
  expect_equal(
    j, 
    kleinj(tau)
  )
})

test_that("kleinj and Dedekind eta", {
  tau <- 0.1 + 2i
  x <- (eta(tau)/eta(5*tau))^6 / 125
  expected <- (5*x*x + 10*x + 1)^3 / x^5
  obtained <- kleinj(tau)
  expect_equal(obtained, expected)
})