#' @title Dedekind eta function
#' @description Evaluation of the Dedekind eta function.
#'
#' @param tau a complex number with strictly positive imaginary part
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' eta(2i)
#' gamma(1/4) / 2^(11/8) / pi^(3/4)
eta <- function(tau){
  stopifnot(isComplexNumber(tau))
  if(Im(tau) <= 0){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  exp(1i * pi * tau/12) * jtheta3_fremling((tau+1)/2, 3*tau)
}
