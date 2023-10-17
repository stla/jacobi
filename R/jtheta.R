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
      JTheta2(cbind(z/pi), tau)[, 1L]
    }else{
      JTheta2(z/pi, tau)
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


#' @title Jacobi theta function with characteristics
#' @description Evaluates the Jacobi theta function with characteristics.
#'
#' @param a,b the characteristics, two complex numbers
#' @template jthetaTemplate
#'
#' @return A complex number, vector or matrix, like \code{z}.
#' @export
#'
#' @details
#' The Jacobi theta function with characteristics generalizes the four Jacobi 
#' theta functions. It is denoted by 
#' \ifelse{html}{\out{&#120579;[a,b](z|&tau;)}}{\eqn{\theta[a,b](z|\tau)}{theta(z|tau)}}.
#' One gets the four Jacobi theta functions when \code{a} and \code{b} take the 
#' values \code{0} or \code{0.5}:
#' \describe{
#' \item{if \code{a=b=0.5}}{then one gets
#'   \ifelse{html}{\out{-&#120599;<sub>1</sub>(z|&tau;)}}{\eqn{\vartheta_1(z|\tau)}{theta_1(z|tau)}}
#' }
#' \item{if \code{a=0.5} and \code{b=0}}{then one gets
#'   \ifelse{html}{\out{&#120599;<sub>2</sub>(z|&tau;)}}{\eqn{\vartheta_2(z|\tau)}{theta_2(z|tau)}}
#' }
#' \item{if \code{a=b=0}}{then one gets
#'   \ifelse{html}{\out{&#120599;<sub>3</sub>(z|&tau;)}}{\eqn{\vartheta_3(z|\tau)}{theta_3(z|tau)}}
#' }
#' \item{if \code{a=0} and \code{b=0.5}}{then one gets
#'   \ifelse{html}{\out{&#120599;<sub>4</sub>(z|&tau;)}}{\eqn{\vartheta_4(z|\tau)}{theta_4(z|tau)}}
#' }
#' }
#' Both 
#' \ifelse{html}{\out{&#120579;[a,b](z+&pi;|&tau;)}}{\eqn{\theta[a,b](z+\pi|\tau)}{theta(z+pi|tau)}}
#' and
#' \ifelse{html}{\out{&#120579;[a,b](z+&pi;&times;&tau;|&tau;)}}{\eqn{\theta[a,b](z+\pi\tau|\tau)}{theta(z+pi*tau|tau)}}
#' are equal to 
#' \ifelse{html}{\out{&#120579;[a,b](z|&tau;)}}{\eqn{\theta[a,b](z|\tau)}{theta(z|tau)}}
#' up to a factor - see the examples for the details.
#' 
#' @note
#' Different conventions are used in the book cited as reference.
#' 
#' @references 
#' Hershel M. Farkas, Irwin Kra. 
#' \emph{Theta Constants, Riemann Surfaces and the Modular Group: An Introduction with Applications to Uniformization Theorems, Partition Identities and Combinatorial Number Theory}.
#' Graduate Studies in Mathematics, volume 37, 2001.
#' 
#' @examples
#' a   <- 2 + 0.3i
#' b   <- 1 - 0.6i
#' z   <- 0.1 + 0.4i
#' tau <- 0.2 + 0.3i
#' jab <- jtheta_ab(a, b, z, tau) 
#' # first property ####
#' jtheta_ab(a, b, z + pi, tau) # is equal to:
#' jab * exp(2i*pi*a)
#' # second property ####
#' jtheta_ab(a, b, z + pi*tau, tau) # is equal to:
#' jab * exp(-1i*(pi*tau + 2*z + 2*pi*b))
jtheta_ab <- function(a, b, z, tau = NULL, q = NULL) {
  stopifnot(isComplexNumber(a))
  stopifnot(isComplexNumber(b))
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  tau <- check_and_get_tau(tau, q)
  alpha <- a * tau
  beta  <- z/pi + b
  C <- exp(1i*pi*a*(alpha + 2*beta)) 
  if(length(z) == 1L) {
    C * jtheta3_cpp(alpha + beta, tau)
  } else {
    if(!is.matrix(z)){
      C * JTheta3(cbind(alpha + beta), tau)[, 1L]
    } else {
      C * JTheta3(alpha + beta, tau)
    }
  }
}
