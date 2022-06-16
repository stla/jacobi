isComplex <- function(z){
  (is.complex(z) || is.numeric(z)) && length(z) == 1L && !is.na(z)
}

check_and_get_tau <- function(tau, q){
  if(is.null(tau) && is.null(q)){
    stop("You must supply either `tau` or `q`.")
  }
  if(!is.null(tau) && !is.null(q)){
    stop("You must supply either `tau` or `q`, not both.")
  }
  if(!is.null(tau)){
    stopifnot(isComplex(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be strictly positive.")
    }
  }
  if(!is.null(q)){
    stopifnot(isComplex(q))
    if(Mod(q) >= 1){
      stop("The modulus of `q` must be strictly less than one.")
    }
    if(Im(q) == 0 && Re(q) <= 0){
      stop("If `q` is real, it must be strictly positive.")
    }
    tau <- -1i * log(q) / pi
  }
  tau
}

#' @title Jacobi theta function one
#' @description Evaluates the first Jacobi theta function.
#'
#' @param z complex number
#' @param tau complex number with strictly positive imaginary part; 
#'   the two complex numbers \code{tau} and \code{q} are related by 
#'   \code{q = exp(1i*pi*tau)}, and only one of them must be supplied
#' @param q the nome, a complex number whose modulus is strictly less than one, 
#'   and which is not zero nor a negative real number
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