library(jacobi)
library(Carlson)

al <- function(z) {
  z <- as.complex(z)
  kz <- lambda(1i*sqrt(z))
  K <- elliptic_F(pi/2, kz)
  pi / (4 * K*K) + sqrt(z)*(1 - elliptic_E(pi/2, kz)/K)
}

al(2)
