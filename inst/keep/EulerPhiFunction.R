library(jacobi)
# https://mathworld.wolfram.com/q-PochhammerSymbol.html
# https://dlmf.nist.gov/20.4
tau0 <- 0.15+0.15i
q <- exp(1i * pi * tau0)
( tau0 <- -1i * log(q) / pi )
q2 <- q*q
tau <- -1i * log(q2) / pi


(eta(tau0))^3 * 2 
(1/sqrt(3)  * (jtheta2(pi/6, q = q^(1/3))))^3 * 2 # pb choix de q^(1/3)

jacobi:::jtheta1prime0(q = q)
jacobi:::theta1dash(0, tau0) # le bon selon mpmath

1/2 * q^(-1/4) * (theta1dash(0, tau))

q <- 0.556 + 0.283i
jacobi:::jtheta1prime0(q = q)
tau <- -1i * log(q) / pi
jacobi:::theta1dash(0, tau) # le bon selon mpmath

2*eta(tau)^3 # should be the same, and is not! (works for other values of q)


q <- 0.556 + 0.283i
jtheta2(0, q=q)
jtheta3(0, q=q)
jtheta4(0, q=q)
tau <- -1i * log(q) / pi
jacobi:::jtheta3_fremling(-0.0000001, tau)
