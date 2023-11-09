#' @title Dedekind eta function
#' @description Evaluation of the Dedekind eta function.
#'
#' @param tau a vector of complex numbers with strictly positive 
#'   imaginary parts
#'
#' @return A vector of complex numbers.
#' @export
#'
#' @examples
#' eta(2i)
#' gamma(1/4) / 2^(11/8) / pi^(3/4)
eta <- function(tau) {
  stopifnot(isComplex(tau))
  storage.mode(tau) <- "complex"
  if(any(Im(tau) <= 0)){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  if(length(tau) == 1L){
    exp(1i * pi * tau/12) * jtheta3_cpp((tau+1)/2, 3*tau)
  } else if(!is.matrix(tau)) {
    JTheta2_tau(1/6, cbind(tau/3))[, 1L] / sqrt(3)
  } else {
    JTheta2_tau(1/6, tau/3) / sqrt(3)
  }
}
