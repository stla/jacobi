g_from_omega1_and_tau <- function(omega1, tau){
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


#' @title Weierstrass zeta function
#' @description Evaluation of the Weierstrass zeta function.
#'
#' @param z complex number
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
#' # Mirror symmetry property:
#' z <- 1 + 1i
#' g <- c(1i, 1+2i)
#' zetaw(Conj(z), Conj(g))
#' Conj(zetaw(z, g))
zetaw <- function(z, g = NULL, omega = NULL, tau = NULL){
  stopifnot(isComplex(z))
  if((is.null(g) + is.null(omega) + is.null(tau)) != 2L){
    stop("You must supply exactly one of `g`, `omega` or `tau`.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
    om1_tau <- halfPeriods(g)
    omega1 <- om1_tau[1L]
    tau <- om1_tau[2L]
  }
  if(!is.null(tau)){
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
  w1 <- - omega1 / pi
  q <- exp(1i * pi * tau)
  p <- 1 / w1 / 2
  eta1 <- p / 6 / w1 * jtheta1primeprimeprime0(tau) / jtheta1prime0(tau)
  - eta1 * z + p * dljtheta1(p * z, tau, q)
  # if(fix && (is.nan(out))){
  #   out <- zetaw(z-1, g = g, fix = FALSE) + 2*zetaw(1/2, g)
  #   if(is.nan(out)){
  #     out <- zetaw(z+1, g = g, fix = FALSE) - 2*zetaw(1/2, g)
  #   }
  # }
  # attr(out, "info") <- c(tau = tau, eta1 = eta1, p = p, w1 = w1, w3 = w3)
}
