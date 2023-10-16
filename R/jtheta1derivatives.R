jtheta1prime0 <- function(tau = NULL, q = NULL){
  tau <- check_and_get_tau(tau, q)
  return(theta1dash(0, tau))
  #jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau)
}

jtheta1primeprimeprime0 <- function(tau){
  -theta1dash(0, tau) * E2(tau)
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

# jtheta4prime0 <- function(tau){
#   tau <- check_and_get_tau(tau, NULL)
#   r <- -tau*tau
#   # 1i * tau = - sqrt(r)
#   # - tau^2 = r
#   q <- exp(-pi*sqrt(r))
#   tau <- check_and_get_tau(NULL, q)
#   z <- ellipticAlpha(r)
#   - pi * z * jtheta3_cpp(0, tau)^4 / (4 * sqrt(r) / jtheta4_cpp(0, tau))
# }
# Wolfram: it's derivative with respect to q, not z!