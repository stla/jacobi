#' @title Inverse of Weierstrass elliptic function
#' @description Evaluation of the inverse of the Weierstrass elliptic function.
#'
#' @param w complex number
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
#' @importFrom Carlson Carlson_RF
#'
#' @examples
#' library(jacobi)
#' omega <- c(1.4 - 1i, 1.6 + 0.5i)
#' w <- 1 + 1i
#' z <- wpinv(w, omega = omega)
#' wp(z, omega = omega) # should be w
wpinv <- function(w, g = NULL, omega = NULL, tau = NULL){
  stopifnot(isComplex(w))
  if((is.null(g) + is.null(omega) + is.null(tau)) != 2L){
    stop("You must supply exactly one of `g`, `omega` or `tau`.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
    om1_tau <- omega1_and_tau(g)
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
  Carlson_RF(w-e1, w-e2, w-e3)
}
