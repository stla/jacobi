#' @title Weierstrass zeta function
#' @description Evaluation of the Weierstrass zeta function.
#'
#' @param z complex number, vector or matrix
#' @param g the elliptic invariants, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param omega the half-periods, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param tau the half-periods ratio; supplying \code{tau} is equivalent to 
#'   supply \code{omega = c(1/2, tau/2)}
#'
#' @return A complex number, vector or matrix.
#' @export
#' 
#' @examples
#' # Mirror symmetry property:
#' z <- 1 + 1i
#' g <- c(1i, 1+2i)
#' wzeta(Conj(z), Conj(g))
#' Conj(wzeta(z, g))
wzeta <- function(z, g = NULL, omega = NULL, tau = NULL){
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
    stopifnot(isComplexNumber(tau))
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
        "The imaginary part of the ratio `omega[2]/omega[1]` must be nonnegative."
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
  # q <- exp(1i * pi * tau)
  p <- 1 / w1 / 2
  eta1 <- p / 6 / w1 * jtheta1primeprimeprime0(tau) / jtheta1prime0(tau)
  - eta1 * z + p * dljtheta1(p * z, tau)
  # if(fix && (is.nan(out))){
  #   out <- zetaw(z-1, g = g, fix = FALSE) + 2*zetaw(1/2, g)
  #   if(is.nan(out)){
  #     out <- zetaw(z+1, g = g, fix = FALSE) - 2*zetaw(1/2, g)
  #   }
  # }
  # attr(out, "info") <- c(tau = tau, eta1 = eta1, p = p, w1 = w1, w3 = w3)
}
