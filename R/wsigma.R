#' @title Weierstrass sigma function
#' @description Evaluation of the Weierstrass sigma function.
#'
#' @param z a complex number
#' @param g the elliptic invariants, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param omega the half-periods, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param tau the half-periods ratio; supplying \code{tau} is equivalent to 
#'   supply \code{omega = c(1/2, tau/2)}
#'
#' @return A complex number.
#' @export
#' 
#' @examples
#' wsigma(1, g = c(12, -8))
#' # should be equal to:
#' sin(1i*sqrt(3))/(1i*sqrt(3)) / sqrt(exp(1))
wsigma <- function(z, g = NULL, omega = NULL, tau = NULL){
  stopifnot(isComplex(z))
  if((is.null(g) + is.null(omega) + is.null(tau)) != 2L){
    stop("You must supply exactly one of `g`, `omega` or `tau`.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
    om1_tau <- halfPeriods(g)
    omega1 <- om1_tau[1L]
    tau <- om1_tau[2L]
    if(is.infinite(tau)){
      return(2*omega1/pi * exp(1/6*(pi*z/2/omega1)^2) * sin(pi*z/2/omega1))
    }
  }else if(!is.null(tau)){
    stopifnot(isComplex(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be nonnegative.")
    }
    omega1 <- 1/2
    g <- g_from_omega1_and_tau(omega1, tau)
  }else{ # omega is given
    stopifnot(isComplexPair(omega))
    omega1 <- omega[1L]
    tau <- omega[2L]/omega1
    if(Im(tau) <= 0){
      stop(
        "The imaginary part of the `omega[2]/omega[1]` must be nonnegative."
      )
    }
    g <- g_from_omega1_and_tau(omega1, tau)
  }
  if(g[1L] == 0 && g[2L] == 0){
    return(1/z)
  }
  if(g[1L] == 3 && g[2L] == 1){
    return(z/2 + sqrt(3/2)/tan(sqrt(3/2)*z))
  }
  w1 <- -2 * omega1 / pi
  f <- jtheta1prime0(tau = tau)
  h <- -pi/6/w1 * jtheta1primeprimeprime0(tau) / f
  w1 * exp(h*z*z/w1/pi) * jtheta1_cpp(z/w1/pi, tau) / f
}