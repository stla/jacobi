#' @title Arithmetic-geometric mean
#' @description Evaluation of the arithmetic-geometric mean of two 
#'   complex numbers.
#'
#' @param x,y complex numbers
#'
#' @return A complex number, the arithmetic-geometric mean of \code{x} 
#'   and \code{y}.
#' @export
#' 
#' @importFrom Carlson elliptic_F
#'
#' @examples
#' agm(1, sqrt(2))
#' 2*pi^(3/2)*sqrt(2) / gamma(1/4)^2
agm <- function(x, y){
  if(x + y == 0 || x == 0 || y == 0){
    return(0)
  }
  pi/4 * (x + y) / elliptic_F(pi/2, ((x-y)/(x+y))^2)
}