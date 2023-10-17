library(jacobi)

Theta <- function(a, b, z, tau) {
  alpha <- a * tau * pi
  beta  <- z + b*pi
  exp(1i*a*(alpha + 2*beta)) * jtheta3(alpha + beta, tau) 
}

z <- 0.1 + 0.4i
tau <- 0.2 + 0.3i

jtheta1(z, tau)
-Theta(0.5, 0.5, z, tau)
#
jtheta2(z, tau)
Theta(0.5, 0, z, tau)
#
jtheta3(z, tau)
Theta(0, 0, z, tau)
#
jtheta4(z, tau)
Theta(0, 0.5, z, tau)



a   <- 2 + 0.3i
b   <- 1 - 0.6i
z   <- 0.1 + 0.4i
tau <- 0.2 + 0.3i
# first property
Theta(a, b, z + pi, tau) # is equal to:
Theta(a, b, z, tau) * exp(2i*a*pi)
# second property
Theta(a, b, z + pi*tau, tau) # is equal to:
Theta(a, b, z, tau) * exp(-1i*(pi*tau + 2*z + 2*pi*b))

a   <- 2 + 0.3i
b   <- 1 - 0.6i
z   <- 0.1 + 0.4i
tau <- 0.2 + 0.3i
# first property
jtheta_ab(a, b, z + pi, tau) # is equal to:
jtheta_ab(a, b, z, tau) * exp(2i*a*pi)
# second property
jtheta_ab(a, b, z + pi*tau, tau) # is equal to:
jtheta_ab(a, b, z, tau) * exp(-1i*(pi*tau + 2*z + 2*pi*b))


jtheta_ab(1/6, 1/2, 0, tau)^3
jtheta_ab(1/6, 1/6, 0, tau)^3 + jtheta_ab(1/6, 5/6, 0, tau)^3

#
tau <- 2.5 + .1i
q <- exp(1i*pi*tau)
( tau <- jacobi:::check_and_get_tau(NULL, q) )
jacobi:::jtheta1prime0(tau)^(1/3) / jtheta_ab(1/6, 1/2, 0, 3*tau)
(-2i)^(1/3)
