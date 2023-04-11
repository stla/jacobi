library(jacobi)

sm <- function(z) {
  stopifnot(jacobi:::isComplex(z))
  logm <- complex(real = 0, imaginary = 5 * pi / 3)
  m <- exp(logm)
  theta <- 3^0.25 * exp(logm / 4)
  xi <- function(u) {
    s <- jellip("sn", u, m = m)
    c <- jellip("cn", u, m = m)
    d <- jellip("dn", u, m = m)
    x <- theta * s * c * d
    (x - 1) / (x + 1)
  }
  pi3 <- beta(1/3, 1/3)
  xi((z + pi3/6) / (2^(1/3) * theta))
}

cm <- function(z) {
  stopifnot(jacobi:::isComplex(z))
  logm <- complex(real = 0, imaginary = 5 * pi / 3)
  m <- exp(logm)
  theta <- 3^0.25 * exp(logm / 4)
  xi <- function(u) {
    s <- jellip("sn", u, m = m)
    c <- jellip("cn", u, m = m)
    d <- jellip("dn", u, m = m)
    2^(1/3) * (1 + theta*theta*s*s) / (1 + theta*s*c*d)
  }
  pi3 <- beta(1/3, 1/3)
  xi((z + pi3/6) / (2^(1/3) * theta))
}


p <- function(z) {
  wp(z, g = c(0, 1/27))
}

pdash <- function(z) {
  wp(z, g = c(0, 1/27), derivative = 1L)
}

cm <- function(z) {
  x <- 3 * pdash(z)
  (x + 1) / (x - 1)
}

sm <- function(z) {
  -6 * p(z) / (3 * pdash(z) - 1)
}


# cubic Fermat curve x^3+y^3=1
pi3 <- beta(1/3, 1/3)
epsilon <- 0.7
t_ <- seq(-pi3/3 + epsilon, 2*pi3/3 - epsilon, length.out = 100)
pts <- t(vapply(t_, function(t) {
  c(Re(cm(t)), Re(sm(t)))
}, FUN.VALUE = numeric(2L)))
plot(pts, type = "l", asp = 1)

z <- 5 - 5i
cm(z)^3 + sm(z)^3