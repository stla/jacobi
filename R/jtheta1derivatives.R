jtheta1prime0 <- function(tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  exp(ljtheta2_cpp(0, tau) + ljtheta3_cpp(0, tau) + ljtheta4_cpp(0, tau))
}

jtheta1primeprimeprime0 <- function(tau){
  -2 * eta(tau)^3 * E2(tau)
}

dljtheta1 <- function(z, tau){
  # out <- z
  # if(z == 0){
  #   return(jtheta1prime0(tau) / jtheta1_cpp(0, tau))
  # }
  if(length(z) == 1L){
    theta1dash(z, tau) / jtheta1_cpp(z/pi, tau)
  }else{
    if(!is.matrix(z)){
      dLTheta1(cbind(z), tau)[, 1L]
    }else{
      dLTheta1(z, tau)
    }
  }
  # out <- dlogjtheta1(z, q)
  # if(is.nan(out) || is.infinite(out)){
  #   out <- theta1dash(z, q=q) / jtheta1_cpp(z/pi, tau)
  # }
  # out
}
