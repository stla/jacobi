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
#' eta(2i)
#' gamma(1/4) / 2^(11/8) / pi^(3/4)
lambda <- function(tau){
  stopifnot(isComplex(tau))
  if(Im(tau) <= 0){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  (jtheta2_cpp(0, tau) / jtheta3_cpp(0, tau))^4 
}
