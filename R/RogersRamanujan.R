#' @title Rogers-Ramanujan continued fraction
#' @description Evaluates the Rogers-Ramanujan continued fraction.
#'
#' @param q the nome, a complex number whose modulus is strictly less than one, 
#'   and which is not zero nor a negative real number
#'
#' @return A complex number
#' @export
#' @note This function is sometimes denoted by \eqn{R}.
RR <- function(q) {
  tau <- check_and_get_tau(NULL, q)
  x <- atan(0.5 - 0.5 * jtheta4_cpp(0, tau)^2 / jtheta4_cpp(0, 5*tau)^2)
  tan(0.5 * x)^0.2 * tan(0.5 * (pi/2 - x))^0.4
}

#' @title Alternating Rogers-Ramanujan continued fraction
#' @description Evaluates the alternating Rogers-Ramanujan continued fraction.
#'
#' @param q the nome, a complex number whose modulus is strictly less than one, 
#'   and which is not zero nor a negative real number
#'
#' @return A complex number
#' @export
#' @note This function is sometimes denoted by \eqn{S}.
RRa <- function(q) {
  tau <- check_and_get_tau(NULL, q)
  x <- atan(0.5 * jtheta3_cpp(0, tau)^2 / jtheta3_cpp(0, 5*tau)^2 - 0.5)
  tan(0.5 * x)^0.2 / tan(0.5 * (pi/2 - x))^0.4
}
