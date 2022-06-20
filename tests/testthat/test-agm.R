test_that("agm value", {
  expect_equal(
    agm(1, sqrt(2)), 
    as.complex(2 * pi^(3/2) * sqrt(2) / gamma(1/4)^2)
  )
})
