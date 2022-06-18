test_that("jtheta3 value.", {
  expect_equal(jtheta3(0, q = exp(-pi)), as.complex(pi^(1/4) / gamma(3/4)))
})
