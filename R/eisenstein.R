#' @importFrom Carlson elliptic_E
#' @noRd
E2 <- function(tau){
  q3 <- jtheta3_cpp(0, tau)^2
  6/pi * elliptic_E(pi/2, lambda(tau), minerror = 100 * .Machine$double.eps) * 
    q3 - q3^2 - jtheta4_cpp(0, tau)^4
}

E4 <- function(tau){
  (exp(8*ljtheta2_cpp(0, tau)) + 
     exp(8*ljtheta3_cpp(0, tau)) + 
     exp(8*ljtheta4_cpp(0, tau))) / 2
}

E6 <- function(tau){
  (jtheta3_cpp(0, tau)^12 + jtheta4_cpp(0, tau)^12 - 3*jtheta2_cpp(0, tau)^8 * 
      (jtheta3_cpp(0, tau)^4 + jtheta4_cpp(0, tau)^4)) / 2
}

#' @title Eisenstein series
#' @description Evaluation of Eisenstein series with weight 2, 4 or 6.
#'
#' @param n the weight, can be 2, 4 or 6
#' @param q nome, complex number with modulus smaller than one, but not zero
#'
#' @return A complex number, the value of the Eisenstein series.
#' @export
EisensteinE <- function(n, q){
  stopifnot(n %in% c(2, 4, 6))
  stopifnot(isComplexNumber(q))
  tau <- check_and_get_tau(NULL, q) / 2
  switch(
    as.character(n),
    "2" = E2(tau),
    "4" = E4(tau),
    "6" = E6(tau)
  )
}
