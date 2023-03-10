G4 <- function(tau){
  pi^4/45 * E4(tau)
}

G6 <- function(tau){
  2*pi^6/945 * E6(tau)
}

omega1_and_tau <- function(g) {
  g2 <- g[1L]
  g3 <- g[2L]
  if(g2 == 0) {
    omega1 <- gamma(1/3)^3 / 4 / pi / g3^(1/6)
    tau <- complex(real = 0.5, imaginary = sqrt(3)/2)
  } else {
    g2cube <- g2*g2*g2
    j <- 1728 * g2cube / (g2cube - 27*g3*g3)
    if(is.infinite(j)){
      return(c(-1i*pi/2/sqrt(3), complex(real = Inf, imaginary = Inf)))
    } 
    tau <- kleinjinv(j) 
    if(g3 == 0) {
      omega1 <- 1i * sqrt(sqrt(15 / 4 / g2 * G4(tau))) 
    } else {
      omega1 <- sqrt(as.complex(7 * G6(tau) * g2 / (12 * G4(tau) * g3)))
    }
  }
  c(omega1, tau)
}

#' @title Half-periods
#' @description Half-periods from elliptic invariants.
#'
#' @param g2g3 the elliptic invariants, a vector of two complex numbers
#'
#' @return The half-periods, a vector of two complex numbers.
#' @export
halfPeriods <- function(g2g3) {
  stopifnot(isComplexPair(g2g3))
  omega1_tau <- omega1_and_tau(g2g3)
  omega1 <- omega1_tau[1L]
  c(omega1, omega1 * omega1_tau[2L])
}


g2_from_omega1_and_tau <- function(omega1, tau){
  # if(Im(w2)*Re(w1) <= Im(w1)*Re(w2)){
  #   stop(
  #     "Invalid `omega` values. Do you want to exchange `omega1` and `omega2`?"
  #   )
  # }
  j2 <- jtheta2_cpp(0, tau)
  j3 <- jtheta3_cpp(0, tau)
  4/3 * (pi/2/omega1)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
}

g_from_omega1_and_tau <- function(omega1, tau){ # used in zetaw
  # if(Im(w2)*Re(w1) <= Im(w1)*Re(w2)){
  #   stop(
  #     "Invalid `omega` values. Do you want to exchange `omega1` and `omega2`?"
  #   )
  # }
  # tau <- w2 / w1
  # if(Im(tau) <= 0){
  #   stop("The ratio `omega2/omega1` must have a positive imaginary part.")
  # }
  j2 <- jtheta2_cpp(0, tau)
  j3 <- jtheta3_cpp(0, tau)
  g2 <- 4/3 * (pi/2/omega1)**4 * (j2**8 - (j2*j3)**4 + j3**8)
  g3 <- 8/27 * (pi/2/omega1)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8)
  ) + j3**12)
  c(g2, g3)
}

#' @title Elliptic invariants
#' @description Elliptic invariants from half-periods
#'
#' @param omega1omega2 the half-periods, a vector of two complex numbers
#'
#' @return The elliptic invariants, a vector of two complex numbers.
#' @export
ellipticInvariants <- function(omega1omega2) {
  stopifnot(isComplexPair(omega1omega2))
  omega1 <- omega1omega2[1L]
  tau <- omega1omega2[2L] / omega1
  if(Im(tau) <= 0) {
    stop("The ratio `omega2/omega1` must have a positive imaginary part.")
  }
  g_from_omega1_and_tau(omega1, tau)
}