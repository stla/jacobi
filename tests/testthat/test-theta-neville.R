test_that("Neville theta values", {
  expect_equal(theta.c(2.5, m = 0.3), as.complex(-0.65900466676738154967))
  expect_equal(theta.d(2.5, m = 0.3), as.complex(0.95182196661267561994))
  expect_equal(theta.n(2.5, m = 0.3), as.complex(1.0526693354651613637))
  expect_equal(theta.s(2.5, m = 0.3), as.complex(0.82086879524530400536))
})

test_that("Neville theta values involving the nome", {
  m <- -2 + 3i
  actual <- theta.c(1i*Carlson::elliptic_F(pi/2, 1-m), m = m)
  expected <- 1 / (m^0.25 * nome(m)^0.25)
  expect_equal(actual, expected)
  #
  actual <- theta.s(
    Carlson::elliptic_F(pi/2, m) + 1i*Carlson::elliptic_F(pi/2, 1-m), m = m
  )
  expected <- 1 / ((1-m)^0.25 * m^0.25 * nome(m)^0.25)
  expect_equal(actual, expected)
  
})
