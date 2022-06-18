test_that("jellip relations", {
  u <- 0.3 + 0.7i
  m <- 0.4
  expect_equal(jellip("cn", u, m = m), jellip("nc", 1i*u, m = 1-m))
  expect_equal(jellip("sn", u, m = m), -1i*jellip("sc", 1i*u, m = 1-m))
  expect_equal(jellip("dn", u, m = m), jellip("dc", 1i*u, m = 1-m))
})
