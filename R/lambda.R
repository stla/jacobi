#' @title Lambda modular function
#' @description Evaluation of the lambda modular function.
#'
#' @param tau a complex number with strictly positive imaginary part, or 
#'   a vector or matrix of such complex numbers
#'
#' @return A complex number, vector or matrix.
#' @export
#' 
#' @note The lambda function is the square of the elliptic modulus.
#'
#' @examples
#' x <- 2
#' lambda(1i*sqrt(x)) + lambda(1i*sqrt(1/x)) # should be one
lambda <- function(tau){
  stopifnot(isComplex(tau))
  if(any(Im(tau) <= 0)){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  if(length(tau) == 1L){
    (jtheta2_cpp(0, tau) / jtheta3_cpp(0, tau))^4
  }else{
    if(is.matrix(tau)){
      lambda_cpp(tau)
    }else{
      lambda_cpp(cbind(tau))[, 1L]
    }
  }
}
