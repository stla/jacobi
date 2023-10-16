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