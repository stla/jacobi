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

wpinv <- function(w, g = NULL, omega = NULL, tau = NULL){
  stopifnot(isComplex(w))
  if((is.null(g) + is.null(omega) + is.null(tau)) != 2L){
    stop("You must supply exactly one of `g`, `omega` or `tau`.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
    om1_tau <- halfPeriods(g)
    omega1 <- om1_tau[1L]
    omega2 <- omega1 * om1_tau[2L]
    omega <- c(omega1, omega2)
  }else if(!is.null(tau)){
    stopifnot(isComplexNumber(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be nonnegative.")
    }
    omega1 <- 1/2
    omega2 <- tau/2
    omega <- c(omega1, omega2)
  }else{ # omega is given
    stopifnot(isComplexPair(omega))
    omega1 <- omega[1L]
    omega2 <- omega[2L]
    if(Im(omega2/omega1) <= 0){
      stop(
        "The imaginary part of the `omega[2]/omega[1]` must be nonnegative."
      )
    }
  }
  e1 <- wp(omega1, omega = omega)
  e2 <- wp(omega2, omega = omega)
  e3 <- wp(-omega1-omega2, omega = omega)
  Carlson_RF(z-e1, z-e2, z-e3)
}