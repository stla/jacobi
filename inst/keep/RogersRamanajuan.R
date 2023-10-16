library(jacobi)

RR <- function(q) {
  tau <- jacobi:::check_and_get_tau(NULL, q)
  x <- atan(0.5 - 0.5*jacobi:::jtheta4_cpp(0, tau)^2/jacobi:::jtheta4_cpp(0, 5*tau)^2)
  tan(0.5*x)^0.2 * tan(0.5 * (pi/2 - x))^0.4
}

RRa <- function(q) {
  tau <- jacobi:::check_and_get_tau(NULL, q)
  x <- atan(0.5*jacobi:::jtheta3_cpp(0, tau)^2/jacobi:::jtheta3_cpp(0, 5*tau)^2 - 0.5)
  tan(0.5*x)^0.2 / tan(0.5 * (pi/2 - x))^0.4
}

RR(exp(-2*pi))
tan(atan(2)/4)
RRa(exp(-pi))
tan(pi/4 - atan(2)/4)
