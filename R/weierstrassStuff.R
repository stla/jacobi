G4 <- function(tau){
  pi^4/45 * E4(tau)
}

halfPeriods <- function(g){
  g2 <- g[1L]
  g2cube <- g2*g2*g2
  g3 <- g[2L]
  j <- 1728 * g2cube / (g2cube - 27*g3*g3)
  if(is.infinite(j)){
    return(c(-1i*pi/2/sqrt(3), complex(real = Inf, imaginary = Inf)))
  } 
  tau <- kleinjinv(j) 
  omega1 <- 1i * sqrt(sqrt(15 / 4 / g2 * G4(tau)))
  c(omega1, tau)
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

