test_that("sine lemniscate is the inverse of arcsine lemniscate", {
  expect_equal(
    asl(sl(1+1i)), 1+1i
  )
})

test_that("cosine lemniscate is the inverse of arccosine lemniscate", {
  expect_equal(
    acl(cl(1+1i)), 1+1i
  )
})

test_that("Relation hyperbolic lemniscate sine and Jacobi elliptic functions", {
  z <- 1 + 1i
  expect_equal(
    slh(z), 
    jellip("sn", z, m = 0.5) / jellip("cd", z, m = 0.5)
  )
})

test_that("Relation hyperbolic lemniscate cosine and Jacobi elliptic functions", {
  z <- 0.5 + 1i
  expect_equal(
    clh(z), 
    jellip("cd", z, m = 0.5) / jellip("sn", z, m = 0.5)
  )
})