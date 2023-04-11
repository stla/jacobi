test_that("Dixon functions relation", {
  z <- 1 + 2i
  expect_equal(
    cm(z)^3 + sm(z)^3, as.complex(1)
  )
})
