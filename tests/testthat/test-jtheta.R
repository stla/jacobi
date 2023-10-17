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

test_that("jtheta1prime0 triple product.", {
  tau <- 0.556 + 0.283i
  # !! true only if tau is derived from a nome
  expect_equal(
    jtheta1prime0(tau), 
    jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau)
  )
})

test_that("jtheta1prime0 relation with jtheta_ab.", {
  tau <- 3.556 + 0.1283i
  q <- exp(1i*pi*tau)
  tau <- check_and_get_tau(NULL, q)
  expect_equal(
    jtheta1prime0(tau)^(1/3), 
    (-2i)^(1/3) * jtheta_ab(1/6, 1/2, 0, 3*tau)
  )
  expect_equal(
    jtheta1prime0(tau), 
    -2i * jtheta_ab(1/6, 1/2, 0, 3*tau)^3
  )
})

test_that("jtheta1prime0 relation with Dedekind eta.", {
  tau <- 5.556 + 0.1283i
  q <- exp(1i*pi*tau)
  tau <- check_and_get_tau(NULL, q)
  expect_equal(
    jtheta1prime0(tau), 
    2 * eta(tau)^3
  )
})

test_that("AGM relations for theta(0,q).", {
  z <- 0
  q <- 0.556 + 0.283i
  theta3q <- jtheta3(z, q = q) 
  theta4q <- jtheta4(z, q = q) 
  theta3q2 <- jtheta3(z, q = q^2) 
  theta4q2 <- jtheta4(z, q = q^2) 
  expect_equal((theta3q^2+theta4q^2)/2, theta3q2^2)
  expect_equal(theta3q^2*theta4q^2, theta4q2^4)
})

test_that("Jacobi identity.", {
  z <- 0
  q <- 0.556 + 0.283i
  theta2q <- jtheta2(z, q = q) 
  theta3q <- jtheta3(z, q = q) 
  theta4q <- jtheta4(z, q = q) 
  expect_equal(theta2q^4 + theta4q^4, theta3q^4)
})

test_that("An edge case for jtheta2.", {
  tau <- 0.7792256 + 1e-7i
  expected <- 27.7468161 + 31.2412167i
  expect_equal(jtheta2(0, tau), expected)
})

test_that("Periodicity-like properties of jtheta_ab.", {
  a   <- 2 + 0.3i
  b   <- 1 - 0.6i
  z   <- 0.1 + 0.4i
  tau <- 0.2 + 0.3i
  jab <- jtheta_ab(a, b, z, tau)
  # first property
  expect_equal(
    jtheta_ab(a, b, z + pi, tau), 
    jab * exp(2i*pi*a)
  )
  # second property
  expect_equal(
    jtheta_ab(a, b, z + pi*tau, tau), 
    jab * exp(-1i*(pi*tau + 2*z + 2*pi*b))
  )
})

test_that("A formula involving some jtheta_ab at z=0", {
  tau <- 1 + 0.1i
  expect_equal(
    jtheta_ab(1/6, 1/2, 0, tau)^3,
    jtheta_ab(1/6, 1/6, 0, tau)^3 + jtheta_ab(1/6, 5/6, 0, tau)^3
  )
})