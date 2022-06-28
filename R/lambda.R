#' @title Lambda modular function
#' @description Evaluation of the lambda modular function.
#'
#' @param tau a complex number with strictly positive imaginary part
#'
#' @return A complex number.
#' @export
#' 
#' @note The lambda function is the square of the elliptic modulus.
#'
#' @examples
#' x <- 2
#' lambda(1i*sqrt(x)) + lambda(1i*sqrt(1/x)) # should be one
lambda <- function(tau){
  stopifnot(isComplexNumber(tau))
  if(Im(tau) <= 0){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  (jtheta2_cpp(0, tau) / jtheta3_cpp(0, tau))^4 
}
