library(jacobi)

modulo <- function(a, p) {
  i <- if(a > 0) floor(a/p) else ceiling(a/p)
  a - i * p;
}

jthetagen_cpp <- function(a, z, tau) {
  exp(1i*pi*a*(2*z+a*tau)) * jacobi:::jtheta3_cpp(z+a*tau, tau)
}

jthetagen <- function(a, z, tau) {
  jthetagen_cpp(a, z/pi, tau)
}

jthetaGen <- function(a, b, z, tau) {
  #a <- modulo(a, 1)
  jthetagen(a, z+b, tau)
}

z <- .2 - .1i
tau <- .2 + .3i
m <- 2
j <- 4
T(m, j, z, tau) # ok pour m = 3 et m = 1
T(m, -j, -z, tau)
q <- exp(1i*pi*tau)
q^((j*m+1)/3)*T(m, 2*m+j, z, tau) # ok quand m = 2

# m = 2 => m + j
# m = 3 

(log(T(m, 2*m+j, z, tau)) - log(T(m,j,z,tau))) /
  log(q)

T <- function(m, j, z, tau) {
  q <- exp(1i * pi *tau)
  exp(-1i*2*pi*j/m) / q^((j/m)^2/2) * jthetaGen(j/m/2, 1/2, m*z, m*tau)
}

1i* q^((m-2*j)/4) * exp(-1i*z*(j-m)) * 
  jacobi:::jtheta3_cpp(m*z + (j-m)*pi*tau/4, m*tau)


jthetaGen(1/2, 1/2, z, tau)
jacobi:::jtheta1_cpp(z, tau)

a <- 0.25
b <- 1/2
lq <- -1i*tau
exp(-0.5*a^2*lq + 1i*a*(z + 2*pi*b)) * 
  jacobi:::jtheta3(z + 2*pi*b + 1i*a*lq, tau) 
jthetaGen(a, b, z, tau)

test <- function(z) { # c'est ça
  exp(0.5*1i*a*(2*z+a*tau)) * jacobi:::jtheta3(z+a*tau, tau)
}
test(z+2*pi*b)
jtheta1(z, tau)


Test <- function(z) {
  test(z + 2*pi*b)
}
Test(z) * exp(2*1i*pi *a)
Test(z+2*pi)

Theta <- function(a, b, z, tau) {
  # lq <- -1i*tau
  # exp(-0.5*a^2*lq + 1i*a*(z + 2*pi*b)) * 
  #   jacobi:::jtheta3(z + 2*pi*b + 1i*a*lq, tau) 
  z <- z + 2*pi*b
  exp(0.5*1i*a*(2*z+a*tau)) * jacobi:::jtheta3(z+a*tau, tau)
}

Theta <- function(a, b, z, R) {
  tau <- R * 1i / pi
  exp(-0.5*a^2*R + 1i*a*(z + 2*pi*b)) * 
     jtheta3(z/2 + pi*b + 1i*a*R/2, tau/2) 
  # z <- z + 2*pi*b
  # exp(0.5*1i*a*(2*z+a*tau)) * jacobi:::jtheta3(z+a*tau, tau)
}

R <- 0.2 + 0.3i
tau <- 1i * R / pi

jtheta1(z, tau)
jtheta2(z, tau)
jtheta3(z, tau)
jtheta4(z, tau)
Theta(0, 0.5, 2*z, 2*R)

a <- .2 + .3i
Theta(a, b, z, R) * exp(2*1i*pi * a)
Theta(a, b, z+2*pi, R)
