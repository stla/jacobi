test_that("eta values", {
  expect_equal(eta(1i/2), as.complex(gamma(1/4) / 2^(7/8) / pi^(3/4)))
  expect_equal(eta(2i), as.complex(gamma(1/4) / 2^(11/8) / pi^(3/4)))
})
