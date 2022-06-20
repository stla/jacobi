#' @title Neville theta functions
#' @description Evaluation of the Neville theta functions.
#'
#' @param z a complex number
#' @param tau complex number with strictly positive imaginary part; it is 
#'   related to \code{m} and only one of them must be supplied
#' @param m the "parameter", square of the elliptic modulus; it is related to 
#'   \code{tau} and only one of them must be supplied
#'
#' @return A complex number.
#' @export
#' @rdname neville
theta.s <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi
  jtheta3_cpp(0, tau)^2 * jtheta1_cpp(zprime, tau) / jtheta1prime0(tau)
  # jtheta3_cpp(0, tau) * jtheta1_cpp(zprime, tau) / 
  #   jtheta1_cpp(0, tau) / jtheta2_cpp(0, tau)
}

#' @rdname neville
#' @export
theta.c <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi
  jtheta2_cpp(zprime, tau) / jtheta2_cpp(0, tau)
}

#' @rdname neville
#' @export
theta.n <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi 
  jtheta4_cpp(zprime, tau) / jtheta4_cpp(0, tau)
}

#' @rdname neville
#' @export
theta.d <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi
  jtheta3_cpp(zprime, tau) / jtheta3_cpp(0, tau)
}
