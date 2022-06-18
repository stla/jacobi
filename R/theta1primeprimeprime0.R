#' @importFrom Carlson elliptic_E
#' @noRd
E2 <- function(tau){
  q3 <- jtheta3_cpp(0, tau)^2
  6/pi * elliptic_E(pi/2, lambda(tau)) * q3 - q3^2 - jtheta4_cpp(0, tau)^4
}

theta1primeprimeprime0 <- function(tau){
  -2 * eta(tau)^3 * E2(tau)
}