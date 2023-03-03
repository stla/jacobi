test_that("jtheta3 value.", {
  expect_equal(jtheta3(0, q = exp(-pi)), as.complex(pi^(1/4) / gamma(3/4)))
})

test_that("Some values of jtheta functions.", {
  q <- 0.556 + 0.283i
  expect_equal(jtheta2(0, q = q), 2.0062976673+0.8296971428i)
  expect_equal(jtheta3(0, q = q), 2.0061026281+0.8298431382i)
  expect_equal(jtheta4(0, q = q), -0.19172221286-0.25120246965i)
})

test_that("jtheta1prime0 value.", {
  q <- 0.556 + 0.283i
  expect_equal(jtheta1prime0(q = q), 0.1966992019-1.4764061381i)
})
