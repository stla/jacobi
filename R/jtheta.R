#' @title Jacobi theta function one
#' @description Evaluates the first Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number; \code{jtheta1} evaluates the first Jacobi theta 
#'   function and \code{ljtheta1} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta1(1 + 1i, q = exp(-pi/2))
jtheta1 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta1_cpp(z/pi, tau)
}

#' @rdname jtheta1
#' @export
ljtheta1 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  ljtheta1_cpp(z/pi, tau)
}


#' @title Jacobi theta function two
#' @description Evaluates the second Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number; \code{jtheta2} evaluates the second Jacobi theta 
#'   function and \code{ljtheta2} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta2(1 + 1i, q = exp(-pi/2))
jtheta2 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta2_cpp(z/pi, tau)
}

#' @rdname jtheta2
#' @export
ljtheta2 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  ljtheta2_cpp(z/pi, tau)
}


#' @title Jacobi theta function three
#' @description Evaluates the third Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number; \code{jtheta3} evaluates the third Jacobi theta 
#'   function and \code{ljtheta3} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta3(1 + 1i, q = exp(-pi/2))
jtheta3 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta3_cpp(z/pi, tau)
}

#' @rdname jtheta3
#' @export
ljtheta3 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  ljtheta3_cpp(z/pi, tau)
}


#' @title Jacobi theta function four
#' @description Evaluates the fourth Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number; \code{jtheta4} evaluates the fourth Jacobi theta 
#'   function and \code{ljtheta4} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta4(1 + 1i, q = exp(-pi/2))
jtheta4 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta4_cpp(z/pi, tau)
}

#' @rdname jtheta4
#' @export
ljtheta4 <- function(z, tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  ljtheta4_cpp(z/pi, tau)
}

jtheta1prime0 <- function(tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau)
}