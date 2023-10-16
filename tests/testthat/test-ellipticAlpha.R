test_that("A value of the elliptic alpha function", {
  expected <- (sqrt(7) - 2) / 2
  expect_equal(
    ellipticAlpha(7), as.complex(expected)
  )
})
