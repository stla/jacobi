test_that("G4 is modular", {
  z <- 1 + 1i
  expect_equal(
    G4(z), G4(z + 1)
  )
  expect_equal(
    G4(-1/z), z^4 * G4(z)
  )
})

test_that("G6 is modular", {
  z <- 1 + 1i
  expect_equal(
    G6(z), G6(z + 1)
  )
  expect_equal(
    G6(-1/z), z^6 * G6(z)
  )
})

test_that("Klein-j in terms of Eisenstein", {
  z <- 1 + 1i
  expect_equal(
    kleinj(z),
    1728 * E4(z)^3 / (E4(z)^3 - E6(z)^2)
  )
})


# test_that("Relation G2(tau/2), G2(2tau), theta3(tau)", {
#   tau <- 0.5 + 0.7i
#   expect_equal(
#     1/3 * jtheta3(0, tau)^4,
#     4*E2(2*tau) - E2(tau/2)
#   )
# })

test_that("Value of E2(2i)", {
  expect_equal(
    E2(2i), 3/8*jtheta3(0,1i)^4 + 3/(2*pi)
  )
})