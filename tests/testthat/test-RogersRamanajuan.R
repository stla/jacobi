test_that("A value of the Rogers-Ramanajuan function", {
  expected <- tan(atan(2)/4)
  expect_equal(
    RR(exp(-2*pi)), as.complex(expected)
  )
})

test_that("A value of the alternating Rogers-Ramanajuan function", {
  expected <- tan(pi/4 - atan(2)/4)
  expect_equal(
    RRa(exp(-pi)), as.complex(expected)
  )
})

test_that("Value of `q` for which Rogers-Ramanajuan is `i`", {
  tau <- (7 + 1i)/10
  q <- exp(1i * pi * tau)
  expect_equal(
    RR(q^2)^5, 1i, tolerance = 1e-7
  )
})

test_that("Relation to Dedekind eta 1", {
  tau <- 0.1 + 1i
  q <- exp(1i * pi * tau)
  expect_equal(
    1/RR(q^2) - RR(q^2), 1 + eta(tau/5) / eta(5*tau)
  )
})

test_that("Relation to Dedekind eta 2", {
  tau <- 0.1 + 1i
  q <- exp(1i * pi * tau)
  expect_equal(
    1/RR(q^2)^5 - RR(q^2)^5, 11 + (eta(tau) / eta(5*tau))^6
  )
})