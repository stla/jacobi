#' @title Lambda modular function
#' @description Evaluation of the lambda modular function.
#'
#' @param tau a complex number with strictly positive imaginary part, or 
#'   a vector or matrix of such complex numbers; missing values allowed
#' @param transfo Boolean, whether to use a transformation of the values  
#'   of \code{tau} close to the real line; using this option can fix some 
#'   failures of the computation (at the cost of speed), e.g. when the 
#'   algorithm reaches the maximal number of iterations
#'
#' @return A complex number, vector or matrix.
#' @export
#' 
#' @note The lambda function is the square of the elliptic modulus.
#'
#' @examples
#' x <- 2
#' lambda(1i*sqrt(x)) + lambda(1i*sqrt(1/x)) # should be one
lambda <- function(tau, transfo = FALSE){
  stopifnot(isBoolean(transfo))
  stopifnot(isComplexObject(tau))
  storage.mode(tau) <- "complex"
  nonmissing <- !is.na(tau)
  if(any(Im(tau[nonmissing]) <= 0)){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  if(length(tau) == 1L){
    if(is.na(tau)){
      NA_complex_
    }else{
      (jtheta2_cpp(0, tau) / jtheta3_cpp(0, tau))^4
    }
  }else{
    if(is.matrix(tau)){
      if(transfo){
        lambda_transfo(tau)
      }else{
        lambda_cpp(tau)
      }
    }else{
      if(transfo){
        lambda_transfo(cbind(tau))[, 1L]
      }else{
        lambda_cpp(cbind(tau))[, 1L]
      }
    }
  }
}
