test_that("jellip relations", {
  u <- 0.3 + 0.7i
  m <- 0.4
  expect_equal(jellip("cn", u, m = m), jellip("nc", 1i*u, m = 1-m))
  expect_equal(jellip("sn", u, m = m), -1i*jellip("sc", 1i*u, m = 1-m))
  expect_equal(jellip("dn", u, m = m), jellip("dc", 1i*u, m = 1-m))
})

test_that("jellip identities", {
  u <- 0.3 + 0.7i
  m <- 0.4
  cn <- jellip("cn", u, m = m)
  sn <- jellip("sn", u, m = m)
  dn <- jellip("dn", u, m = m)
  expect_equal(sn^2 + cn^2, as.complex(1))
  expect_equal(dn^2 + m*sn^2, as.complex(1))
})

test_that("relation sn & ellipticK", {
  m <- 0.4
  K <- Carlson::elliptic_F(pi/2, m)
  sn <- jellip("sn", K, m = m)
  expect_equal(sn, as.complex(1))
  sn <- jellip("sn", K/2, m = m)
  expect_equal(sn^2, as.complex(1/(1+sqrt(1-m))))
})