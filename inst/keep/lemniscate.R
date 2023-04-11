library(jacobi)
library(Carlson)

# lemniscate constant
ombar <- gamma(1/4)^2 / (2 * sqrt(2*pi))

sl <- function(z) {
  w <- pi / ombar * z
  q <- exp(-pi)
  jtheta1(w, q = q) / jtheta3(w, q = q)
}

cl <- function(z) {
  w <- pi / ombar * z
  q <- exp(-pi)
  jtheta2(w, q = q) / jtheta4(w, q = q)
}

arcsl <- function(x) {
  x <- as.complex(x)
  elliptic_F(asin(sqrt(2)*x / sqrt(1+x*x)), 1/2)   / sqrt(2)
}
arcsl(sl(1+1i))

arccl <- function(x) {
  x <- as.complex(x)
  elliptic_F(acos(x), 0.5) / sqrt(2)
}
arccl(cl(1+1i))


# hyperbolic

slh <- function(z) {
  z <- as.complex(z / sqrt(2))
  cosl <- cl(z)
  (1 + cosl*cosl) * sl(z) / (sqrt(2) * cosl)
}

slh(1+1i)
jellip("sn", 1+1i, m = 0.5) / jellip("cd", 1+1i, m = 0.5)

clh <- function(z) {
  z <- as.complex(z / sqrt(2))
  sinl <- sl(z)
  (1 + sinl*sinl) * cl(z) / (sqrt(2) * sinl)
}

clh(1+1i)
jellip("cd", 1+1i, m = 0.5) / jellip("sn", 1+1i, m = 0.5)
