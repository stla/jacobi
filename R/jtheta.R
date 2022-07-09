#' @title Jacobi theta function one
#' @description Evaluates the first Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number, vector or matrix; \code{jtheta1} evaluates the 
#'   first Jacobi theta function and \code{ljtheta1} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta1(1 + 1i, q = exp(-pi/2))
jtheta1 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    jtheta1_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      JTheta1(cbind(z/pi), tau)[, 1L]
    }else{
      JTheta1(z/pi, tau)
    }
  }
}

#' @rdname jtheta1
#' @export
ljtheta1 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    ljtheta1_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      LJTheta1(cbind(z/pi), tau)[, 1L]
    }else{
      LJTheta1(z/pi, tau)
    }
  }
}


#' @title Jacobi theta function two
#' @description Evaluates the second Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number, vector or matrix; \code{jtheta2} evaluates the 
#'   second Jacobi theta function and \code{ljtheta2} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta2(1 + 1i, q = exp(-pi/2))
jtheta2 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    jtheta2_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      LJTheta2(cbind(z/pi), tau)[, 1L]
    }else{
      LJTheta2(z/pi, tau)
    }
  }
}

#' @rdname jtheta2
#' @export
ljtheta2 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    ljtheta2_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      LJTheta2(cbind(z/pi), tau)[, 1L]
    }else{
      LJTheta2(z/pi, tau)
    }
  }
}

#' @title Jacobi theta function three
#' @description Evaluates the third Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number, vector or matrix; \code{jtheta3} evaluates the 
#'   third Jacobi theta function and \code{ljtheta3} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta3(1 + 1i, q = exp(-pi/2))
jtheta3 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    jtheta3_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      JTheta3(cbind(z/pi), tau)[, 1L]
    }else{
      JTheta3(z/pi, tau)
    }
  }
}

#' @rdname jtheta3
#' @export
ljtheta3 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    ljtheta3_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      LJTheta3(cbind(z/pi), tau)[, 1L]
    }else{
      LJTheta3(z/pi, tau)
    }
  }
}

#' @title Jacobi theta function four
#' @description Evaluates the fourth Jacobi theta function.
#'
#' @template jthetaTemplate
#'
#' @return A complex number, vector or matrix; \code{jtheta4} evaluates the 
#'   fourth Jacobi theta function and \code{ljtheta4} evaluates its logarithm.
#' @export
#'
#' @examples
#' jtheta4(1 + 1i, q = exp(-pi/2))
jtheta4 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    jtheta4_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      JTheta4(cbind(z/pi), tau)[, 1L]
    }else{
      JTheta4(z/pi, tau)
    }
  }
}

#' @rdname jtheta4
#' @export
ljtheta4 <- function(z, tau = NULL, q = NULL){
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  if(length(z) == 1L){
    ljtheta4_cpp(z[1L]/pi, tau)
  }else{
    if(!is.matrix(z)){
      LJTheta4(cbind(z/pi), tau)[, 1L]
    }else{
      LJTheta4(z/pi, tau)
    }
  }
}
