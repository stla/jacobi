#' @title Nome
#' @description The nome in function of the parameter \eqn{m}.
#'
#' @param m the parameter, square of elliptic modulus
#'
#' @return A complex number.
#' @export
#' @importFrom Carlson elliptic_F
#'
#' @examples
#' nome(-2)
nome <- function(m) {
  if(m == 0) {
    return(as.complex(0))
  }
  if(m == 1) {
    return(as.complex(1))
  }
  exp(-pi * elliptic_F(pi/2, 1-m) / elliptic_F(pi/2, m))
}