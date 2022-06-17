isComplexPair <- function(x){
  (is.complex(x) || is.numeric(x)) && length(x) == 2L && !anyNA(x)
}

g_from_omega <- function(w1, w2){
  if(Im(w2)*Re(w1) < Im(w1)*Re(w2)){
    stop(
      "Invalid `omega` values. Do you want to exchange `omega1` and `omega2`?"
    )
  }
  ratio <- w2 / w1
  if(Im(ratio) <= 0){
    stop("The ratio `omega2/omega1` must have a positive imaginary part.")
  }
  q <- exp(1i * pi * ratio)
  j2 <- jtheta2(0, q = q)
  j3 <- jtheta3(0, q = q)
  g2 <- 4/3 * (pi/2/w1)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  g3 <- 8/27 * (pi/2/w1)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  ) + j3**12)
  c(g2, g3)
}

#' @title Weierstrass elliptic function
#' @description Evaluation of the Weierstrass elliptic function and its 
#'   derivatives.
#'
#' @param z complex number
#' @param g the elliptic invariants, a vector of two complex numbers; they are 
#'   related to the half-periods (\code{omega}) and only one of \code{g} 
#'   and \code{omega} must be given
#' @param omega the half-periods, a vector of two complex numbers; they are 
#'   related to the elliptic invariants (\code{g}) and only one of \code{g} 
#'   and \code{omega} must be given
#' @param derivative differentiation order, an integer between 0 and 3
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' omega1 <- 1.4 - 1i
#' omega2 <- 1.6 + 0.5i
#' omega <- c(omega1, omega2)
#' e1 <- pweierstrass(omega1, omega = omega)
#' e2 <- pweierstrass(omega2, omega = omega)
#' e3 <- pweierstrass(-omega1-omega2, omega = omega)
#' e1 + e2 + e3 # should be 0
pweierstrass <- function(z, g = NULL, omega = NULL, derivative = 0L){
  if(!is.element(derivative, 0L:3L)){
    stop("`derivative` must be an integer between 0 and 3.") 
  }
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
  g2 <- g[1L]
  g3 <- g[2L]
  r <- polyroot(c(-g3, -g2, 0, 4))
  r1 <- r[1L]
  r2 <- r[2L]
  r3 <- r[3L]
  e3 <- r3
  a <- sqrt(r1 - r3)
  b <- sqrt(r1 - r2)
  c <- sqrt(r2 - r3)
  if(abs(a + b) < abs(a - b)) b <- -b
  if(abs(a + c) < abs(a - c)) c <- -c
  if(abs(c + 1i*b) < abs(c - 1i*b)){
    e3 <- r1
    a <- sqrt(r3 - r1)
    b <- sqrt(r3 - r2)
    c <- sqrt(r2 - r1)
    w1 <- 1 / agm(1i*b, c) 
  }else{
    w1 <- 1 / agm(a, b)
  }
  w3 <- 1i / agm(a, c)
  q <- exp(1i * pi * w3/w1)
  if(derivative != 1){
    pw0 <- e3 + 
      (jtheta2(0, q = q) * jtheta3(0, q = q) * jtheta4(z/w1, q = q) /
        (w1 * jtheta1(z/w1, q = q)))**2
    if(derivative == 0) return(pw0)
    if(derivative == 2) return(6*pw0**2 - g2/2)
  }
  f <- jtheta1prime0(q = q)**3 /
    (jtheta2(0, q = q) * jtheta3(0, q = q) * jtheta4(0, q = q))
  pw1 <- -2*(1/w1)**3 * jtheta2(z/w1, q = q) * jtheta3(z/w1, q = q) *
    jtheta4(z/w1, q = q) * f / jtheta1(z/w1, q = q)**3
  if(derivative == 1) return(pw1)
  12 * pw0 * pw1 
}
