dljtheta1 <- function(z, tau, q){
  if(z == 0){
    return(jtheta1prime0(tau) / jtheta1_cpp(0, tau))
  }
  dlogjtheta1(z, q)
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
#' # Mirror symmetry property:
#' z <- 1 + 1i
#' g <- c(1i, 1+2i)
#' zetaw(Conj(z), Conj(g))
#' Conj(zetaw(z, g))
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
  if(g[1L] == 0 && g[2L] == 0){
    return(1/z)
  }
  if(g[1L] == 3 && g[2L] == 1){
    return(z/2 + sqrt(3/2)/tan(sqrt(3/2)*z))
  }
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
  p <- 1 / w1 
  eta1 <- p / 3 / w1 * jtheta1primeprimeprime0(tau) / jtheta1prime0(tau)
  - eta1 * z + p * dljtheta1(p*z, tau, q)
}
