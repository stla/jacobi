test_that("kleinj at 2i", {
  expect_equal(
    kleinj(2i), 
    as.complex(66^3)
  )
})

test_that("kleinj alternate expression", {
  tau <- 0.5 + 0.3i
  lbd <- lambda(tau) 
  x <- lbd * (1- lbd)
  expect_equal(
    kleinj(tau), 
    256*(1-x)^3 / x^2
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