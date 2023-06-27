test_that("agm value", {
  expect_equal(
    agm(1, sqrt(2)), 
    as.complex(2 * pi^(3/2) * sqrt(2) / gamma(1/4)^2)
  )
})

test_that("agm functional equation", {
  a <- 1 + 2i
  b <- 2 + 3i
  expect_equal(
    agm(a, b), 
    agm((a+b)/2, sqrt(a*b))
  )
})

test_that("agm - Gauss theorem", {
  a <- 2
  b <- 3
  f <- function(phi) {
    1 / sqrt(a^2*cos(phi)^2 + b^2*sin(phi)^2)
  }
  I <- integrate(f, lower = 0, upper = pi/2)
  expect_equal(
    1 / agm(a, b), 
    as.complex(2/pi * I$value)
  )
})