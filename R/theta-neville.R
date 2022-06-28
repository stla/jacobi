#' @title Neville theta functions
#' @description Evaluation of the Neville theta functions.
#'
#' @param z a complex number, vector, or matrix
#' @param tau complex number with strictly positive imaginary part; it is 
#'   related to \code{m} and only one of them must be supplied
#' @param m the "parameter", square of the elliptic modulus; it is related to 
#'   \code{tau} and only one of them must be supplied
#'
#' @return A complex number, vector or matrix.
#' @export
#' @rdname neville
theta.s <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  j3sq <- jtheta3_cpp(0, tau)^2
  zprime <- z / j3sq / pi
  if(length(z) == 1L){
    j3sq * jtheta1_cpp(zprime, tau) / jtheta1prime0(tau)
  }else{
    if(!is.matrix(z)){
      j3sq * JTheta1(cbind(zprime), tau)[, 1L] / jtheta1prime0(tau)
    }else{
      j3sq * JTheta1(zprime, tau) / jtheta1prime0(tau)
    }
  }
  # jtheta3_cpp(0, tau) * jtheta1_cpp(zprime, tau) / 
  #   jtheta1_cpp(0, tau) / jtheta2_cpp(0, tau)
}

#' @rdname neville
#' @export
theta.c <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplex(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi
  if(length(z) == 1L){
    jtheta2_cpp(zprime, tau) / jtheta2_cpp(0, tau)
  }else{
    if(!is.matrix(z)){
      JTheta2(cbind(zprime), tau)[, 1L] / jtheta2_cpp(0, tau)
    }else{
      JTheta2(zprime, tau) / jtheta2_cpp(0, tau)
    }
  }
}

#' @rdname neville
#' @export
theta.n <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplexNumber(z))
  tau <- check_and_get_tau_from_m(tau, m)
  zprime <- z / jtheta3_cpp(0, tau)^2 / pi 
  if(length(z) == 1L){
    jtheta4_cpp(zprime, tau) / jtheta4_cpp(0, tau)
  }else{
    if(!is.matrix(z)){
      JTheta4(cbind(zprime), tau)[, 1L] / jtheta4_cpp(0, tau)
    }else{
      JTheta4(zprime, tau) / jtheta4_cpp(0, tau)
    }
  }
}

#' @rdname neville
#' @export
theta.d <- function(z, tau = NULL, m = NULL){
  stopifnot(isComplexNumber(z))
  tau <- check_and_get_tau_from_m(tau, m)
  j3 <- jtheta3_cpp(0, tau)
  zprime <- z / j3 / j3 / pi
  if(length(z) == 1L){
    jtheta3_cpp(zprime, tau) / jtheta3_cpp(0, tau)
  }else{
    if(!is.matrix(z)){
      JTheta3(cbind(zprime), tau)[, 1L] / j3
    }else{
      JTheta3(zprime, tau) / j3
    }
  }
}
