test_that("am works", {
  phi <- 1 + 1i
  m <- 2
  u <- Carlson::elliptic_F(phi, m)
  expect_equal(
    am(u, m), 
    phi
  )
})
