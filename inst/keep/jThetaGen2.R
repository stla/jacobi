library(jacobi)

Theta <- function(a, b, z, tau) {
  alpha <- a * tau * pi
  beta  <- z + b*pi
  exp(1i*a*(alpha + 2*beta)) * 
    jtheta3(alpha + beta, tau) 
  # z <- z + 2*pi*b
  # exp(0.5*1i*a*(2*z+a*tau)) * jacobi:::jtheta3(z+a*tau, tau)
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