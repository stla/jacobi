#' @title Amplitude function
#' @description Evaluation of the amplitude function.
#'
#' @param u complex number
#' @param k elliptic modulus, a complex number
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' library(Carlson)
#' phi <- 1 + 1i
#' k <- 2
#' u <- elliptic_F(phi, k^2)
#' am(u, k) # should be phi
am <- function(u, k){
  asin(jellip("sn", u, m = k*k))
}