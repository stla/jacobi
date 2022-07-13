test_that("wpinv works.", {
  omega <- c(1.4 - 1i, 1.6 + 0.5i)
  w <- 1 + 1i
  z <- wpinv(w, omega = omega)
  expect_equal(w, wp(z, omega = omega))
})
