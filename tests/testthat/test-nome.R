test_that("Some value of the nome", {
  expect_equal(nome(-2), as.complex(-0.06823783), tolerance = 1e-6)
  expect_equal(nome(0.5), as.complex(exp(-pi)))
  expect_equal(nome(5/16), as.complex(0.02340125), tolerance = 1e-6)
})