#' @title Elliptic alpha function
#' @description Evaluates the elliptic alpha function.
#'
#' @param z a complex number
#'
#' @return A complex number.
#' @export
#' @importFrom Carlson elliptic_F elliptic_E
#'
#' @references Weisstein, Eric W. 
#' \href{https://mathworld.wolfram.com/EllipticAlphaFunction.html}{"Elliptic Alpha Function"}.
ellipticAlpha <- function(z) {
  stopifnot(isComplexNumber(z))
  z <- as.complex(z)
  sqrtz <- sqrt(z)
  kz <- lambda(1i*sqrtz)
  K <- elliptic_F(pi/2, kz)
  pi / (4 * K*K) + sqrtz*(1 - elliptic_E(pi/2, kz)/K)
}

