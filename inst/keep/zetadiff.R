E2 <- function (tau) 
{
  q <- exp(1i*pi*tau)
  q3 <- jacobi:::jtheta3_cpp(0, tau)^2
  q4 <- jacobi:::jtheta4_cpp(0, tau)^4
  6/pi * Carlson::elliptic_E(pi/2, 1-q4/q3^2, minerror = 1000 * .Machine$double.eps) * 
    q3 - q3^2 - q4
}

G2 <- function(tau){
  E2(tau) * pi^2 / 3
}

G2 <- function(tau){
  pi^2/3 - 8*pi^2*(E2(tau)-1)/(-24)
}

n <- 1L:60L
tau <- 1i/2
q <- exp(1i*pi*tau)
1 - 24*sum(n*q^n/(1-q^n))
E2(tau/2)



library(jacobi)
library(jacobi)

d <- 1
z <- 0.7+0.6i
om <- c(1/2, 1i/2)
elliptic::zeta(z+d, Om=om) - 2*elliptic::zeta(d/2, Om=om)
elliptic::zeta(z, Om=om)

library(elliptic)

z <- 0.27+0.16i
omega1 <- 1i
omega2 <- 6/7
Omega <- c(omega1, omega2)/2
Omega <- as.primitive(Omega)
tau <- Omega[2L] / Omega[1L]
Omega <- c(1, tau)/2
zeta(z+1, Omega=Omega) - zeta(z, Omega=Omega)
# 3.289593+0i
G2(tau)#/omega1^2
# 3.289593+0i

a <- 6
zetaw(z, om=c(1, tau)/a) 
# w1 - w3 = 
1/(a*pi) # = 
# tau de info = 1 + 1/tau

z <- 0.17+0.26i
omega1 <- 1i/3
omega2 <- 1
Omega <- c(omega1, omega2)/2
tau <- omega1/omega2
elliptic::zeta(z+(tau), Om=om) - elliptic::zeta(z, Om=om) 
1/tau*G2(tau) #+ 2i*pi

z <- 0.7+0.6i
omega1 <- 1i
omega2 <- 6
Omega <- c(omega1, omega2)/2
Omega <- as.primitive(Omega)
tau <- Omega[2L] / Omega[1L]
Omega <- c(1, tau)
eta1 <- zeta(z+2*Omega[1L], Omega=Omega) - zeta(z, Omega=Omega)
eta2 <- zeta(z+2*tau, Omega=Omega) - zeta(z, Omega=Omega)
omega2/2*eta1 - omega1/2*eta2 # 1i*2*pi
eta2
-2*tau*G2(tau)
eta1
3*tau*G2(tau)

zetaw(1i+0.5, omega = c(1/2,1)) - zetaw(1i-0.5, omega = c(1i/2,-1)) 
1i/2*G2(1i/2) -2i*pi
kappa <- 1/0.1591549
zetaw(kappa*(1+1i/2), omega = kappa*c(1, 1i/2)) - zetaw(kappa, omega = kappa*c(1, 1i/2))

2 * zetaw(2 * z, omega=c(2,4i)) 
zetaw(z, omega=c(1,2i))
k <- 3.5 # real only
k * zetaw(k*z, omega=k*c(1,2i))


zetaw(0.5+1i, omega = c(1/2, 1i/2))
k <- 3.5
k * zetaw(k*z, omega=k*c(1,2i))




tau <- .2i
m <- n <- c((-25):(-1), 1:25)
m <- 1:85
2*sum(pi^2/sin(pi*m*tau)^2)
sum(1/(4*tau+m))
pi/tan(pi*4*tau)
sum(1/(4*tau+m)^2)
pi^2/sin(pi*4*tau)^2

m <- 1:65
E2(tau) / (sum(pi^2/sin(pi*m*tau)^2) + sum(pi^2/sin(-pi*m*tau)^2) + pi^2/3)

s = 0
for(k in 1:185){
  s <- s + sum(1/(k*tau+m)^2) + sum(1/(-k*tau+m)^2)
}


s = 0
for(k in 1:185){
  s <- s + sum(1/(k*tau+m)^2) + sum(1/(-k*tau+m)^2)
}
s + pi^2/3
E2(tau)












