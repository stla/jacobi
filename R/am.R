#' @title Amplitude function
#' @description Evaluation of the amplitude function.
#'
#' @param u complex number
#' @param m square of elliptic modulus, a complex number
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' library(Carlson)
#' phi <- 1 + 1i
#' m <- 2
#' u <- elliptic_F(phi, m)
#' am(u, m) # should be phi
am <- function(u, m){
  w <- asin(jellip("sn", u, m=m))
  k <- round(Re(u)/pi) + round(Re(w)/pi)
  (-1)^k * w + k * pi
}