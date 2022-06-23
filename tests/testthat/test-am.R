test_that("am works", {
  phi <- 1 + 1i
  k <- 2
  u <- Carlson::elliptic_F(phi, k^2)
  expect_equal(
    am(u, k), 
    phi
  )
})
