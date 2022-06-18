#' @title Klein j-function
#' @description Evaluation of the Klein j-invariant function.
#'
#' @param tau a complex number with strictly positive imaginary part
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' kleinj(2i)
#' 66^3
kleinj <- function(tau){
  stopifnot(isComplex(tau))
  if(Im(tau) <= 0){
    stop("The complex number `tau` must have a positive imaginary part.")
  }
  j2 <- jtheta2(0, tau = tau)
  j3 <- jtheta3(0, tau = tau)
  g2 <- 4/3 * (pi/2)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  g3 <- 8/27 * (pi/2)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  ) + j3**12)
  g2cube <- g2*g2*g2
  1728 * g2cube / (g2cube - 27*g3*g3)
}
