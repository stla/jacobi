test_that("sigmaw values", {
  expect_equal(
    sigmaw(1, g = c(12, -8)), 
    sin(1i*sqrt(3))/(1i*sqrt(3)*sqrt(exp(1)))
  )
  expect_equal(
    sigmaw(2, g = c(1, 2i)), 
    1.864625371572 - 0.306600135476i
  )
})
