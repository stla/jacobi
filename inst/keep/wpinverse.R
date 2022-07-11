library(jacobi)

omega1 <- 1.4 - 1i
omega2 <- 1.6 + 0.5i
omega <- c(omega1, omega2)
e1 <- wp(omega1, omega = omega)
e2 <- wp(omega2, omega = omega)
e3 <- wp(-omega1-omega2, omega = omega)

z <- 1+1i

f <- Carlson::Carlson_RF(z-e1, z-e2, z-e3)
wp(f, omega = c(omega1, omega2))

