jtheta1prime <- function(z, tau){
  jtheta1_cpp(z/pi, tau) * dlogjtheta1(z, exp(1i * pi * tau))
}

#' @title Weierstrass zeta function
#' @description Evaluation of the Weierstrass zeta function.
#'
#' @param z complex number
#' @param g the elliptic invariants, a vector of two complex numbers; they are 
#'   related to the half-periods (\code{omega}) and only one of \code{g} 
#'   and \code{omega} must be given
#' @param omega the half-periods, a vector of two complex numbers; they are 
#'   related to the elliptic invariants (\code{g}) and only one of \code{g} 
#'   and \code{omega} must be given
#'
#' @return A complex number.
#' @export
#' 
#' @examples
#' 2+2
zetaw <- function(z, g = NULL, omega = NULL){
  stopifnot(isComplex(z))
  if(is.null(g) && is.null(omega)){
    stop("You must supply either `g` or `omega`.")
  }
  if(!is.null(g) && !is.null(omega)){
    stop("You must supply either `g` or `omega`, not both.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
  }
  if(!is.null(omega)){
    stopifnot(isComplexPair(omega))
    g <- g_from_omega(omega[1L], omega[2L])
  }
  # g2 <- g[1L]
  # g3 <- g[2L]
  # r <- sort(polyroot(c(-g3, -g2, 0, 4)))
  r <- e3e2e1(g)
  r1 <- r[1L]
  r2 <- r[2L]
  r3 <- r[3L]
  a <- sqrt(r1 - r3)
  b <- sqrt(r1 - r2)
  c <- sqrt(r2 - r3)
  if(abs(a + b) < abs(a - b)) b <- -b
  if(abs(a + c) < abs(a - c)) c <- -c
  if(abs(c + 1i*b) < abs(c - 1i*b)){
    a <- sqrt(r3 - r1)
    b <- sqrt(r3 - r2)
    c <- sqrt(r2 - r1)
    w1 <- 1 / agm(1i*b, c) 
  }else{
    w1 <- 1 / agm(a, b)
  }
  w3 <- 1i / agm(a, c)
  tau <- w3 / w1
  if(Im(tau) <= 0){
    stop("Invalid values of the parameters.")
  }
  q <- exp(1i * pi * tau)
  p <- 1 / w1 / 2
  eta1 <- p / 6 / w1 * jtheta1primeprimeprime0(tau) / jtheta1prime0(tau)
  - eta1 * z + p * jtheta1prime(p*z, tau) / jtheta1_cpp(p*z/pi, tau)
}
