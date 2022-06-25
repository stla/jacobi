jtheta1prime0 <- function(tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau)
}

jtheta1primeprimeprime0 <- function(tau){
  -2 * eta(tau)^3 * E2(tau)
}

dljtheta1 <- function(z, tau, q){
  if(z == 0){
    return(jtheta1prime0(tau) / jtheta1_cpp(0, tau))
  }
  out <- dlogjtheta1(z, q)
  if(is.nan(out) || is.infinite(out)){
    out <- theta1dash(z, q=q) / jtheta1_cpp(z/pi, tau)
  }
  out
}
