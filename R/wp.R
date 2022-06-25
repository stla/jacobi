wp_from_tau <- function(z, tau){ # wp(z, omega1 = 1/2, omega2 = tau/2)
  j2 <- jtheta2_cpp(0, tau)
  j3 <- jtheta3_cpp(0, tau)
  (pi * j2 * j3 * jtheta4_cpp(z, tau) / jtheta1_cpp(z, tau))^2 -
    pi^2/3 * (j2^4 + j3^4)
}

wp_from_omega <- function(z, omega1, omega2){
  wp_from_tau(z/omega1/2, omega2/omega1) / omega1 / omega1 / 4
}

wp_from_omega1_and_tau <- function(z, omega1, tau){
  wp_from_tau(z/omega1/2, tau) / omega1 / omega1 / 4
}

wp_from_g <- function(z, g){
  om1_tau <- halfPeriods(g)
  wp_from_omega1_and_tau(z, om1_tau[1L], om1_tau[2L])
}

weierDerivative <- function(z, omega1, tau){
  f <- jtheta1prime0(tau = tau)**3 /
    (jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau))
  w1 <- 2 * omega1 / pi
  z1 <- -z / 2 / omega1
  2*(1/w1)**3 * jtheta2_cpp(z1, tau) * jtheta3_cpp(z1, tau) *
    jtheta4_cpp(z1, tau) * f / jtheta1_cpp(z1, tau)**3
}

#' @title Weierstrass elliptic function
#' @description Evaluation of the Weierstrass elliptic function and its 
#'   derivatives.
#'
#' @param z complex number
#' @param g the elliptic invariants, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param omega the half-periods, a vector of two complex numbers; only 
#'   one of \code{g}, \code{omega} and \code{tau} must be given
#' @param tau the half-periods ratio; supplying \code{tau} is equivalent to 
#'   supply \code{omega = c(1/2, tau/2)}
#' @param derivative differentiation order, an integer between 0 and 3
#'
#' @return A complex number.
#' @export
#'
#' @examples
#' omega1 <- 1.4 - 1i
#' omega2 <- 1.6 + 0.5i
#' omega <- c(omega1, omega2)
#' e1 <- wp(omega1, omega = omega)
#' e2 <- wp(omega2, omega = omega)
#' e3 <- wp(-omega1-omega2, omega = omega)
#' e1 + e2 + e3 # should be 0
wp <- function(z, g = NULL, omega = NULL, tau = NULL, derivative = 0L){
  stopifnot(isComplex(z))
  if(!is.element(derivative, 0L:3L)){
    stop("`derivative` must be an integer between 0 and 3.") 
  }
  if((is.null(g) + is.null(omega) + is.null(tau)) != 2L){
    stop("You must supply exactly one of `g`, `omega` or `tau`.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
    if(derivative != 1){
      weier <- wp_from_g(z, g)
      if(derivative == 0){
        return(weier)
      } 
      if(derivative == 2){
        return(6*weier*weier - g[1L]/2)
      }
    }
    om1_tau <- halfPeriods(g)
    omega1 <- om1_tau[1L]
    tau <- om1_tau[2L]
    weierPrime <- weierDerivative(z, omega1, tau)
    if(derivative == 1){
      return(weierPrime)
    } 
    return(12 * weier * weierPrime) # derivative = 3
  }
  if(!is.null(tau)){
    stopifnot(isComplex(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be nonnegative.")
    }
    omega1 <- 1/2
    if(derivative != 1){
      weier <- wp_from_tau(z, tau)
    }
  }else{ # omega is given
    stopifnot(isComplexPair(omega))
    omega1 <- omega[1L]
    tau <- omega[2L]/omega1
    if(Im(tau) <= 0){
      stop(
        "The imaginary part of the `omega[2]/omega[1]` must be nonnegative."
      )
    }
    if(derivative != 1){
      weier <- wp_from_omega1_and_tau(z, omega1, tau)
    }
  }
  if(derivative == 0){
    return(weier)
  } 
  if(derivative == 2){
    g2 <- g2_from_omega1_and_tau(omega1, tau)
    return(6*weier*weier - g2/2)
  }
  weierPrime <- weierDerivative(z, omega1, tau)
  if(derivative == 1){
    return(weierPrime)
  } 
  12 * weier * weierPrime # derivative = 3
}

# wp <- function(z, g = NULL, omega = NULL, derivative = 0L){
#   stopifnot(isComplex(z))
#   if(!is.element(derivative, 0L:3L)){
#     stop("`derivative` must be an integer between 0 and 3.") 
#   }
#   if(is.null(g) && is.null(omega)){
#     stop("You must supply either `g` or `omega`.")
#   }
#   if(!is.null(g) && !is.null(omega)){
#     stop("You must supply either `g` or `omega`, not both.")
#   }
#   if(!is.null(g)){
#     stopifnot(isComplexPair(g))
#   }
#   if(!is.null(omega)){
#     stopifnot(isComplexPair(omega))
#     g <- g_from_omega(omega[1L], omega[2L])
#   }
#   # g2 <- g[1L]
#   # g3 <- g[2L]
#   # r <- sort(polyroot(c(-g3, -g2, 0, 4)))
#   r <- e3e2e1(g)
#   if(isReal(g)) r <- r[c(1L, 3L, 2L)]# unname(elliptic::e1e2e3(g))
#   # Delta <- g2^3 - 27*g3^2
#   # if(FALSE && Im(Delta) == 0 && Re(Delta) > 0){
#   #   r <- sort(Re(r), decreasing = TRUE)
#   #   r1 <- r[1L]
#   #   r2 <- r[2L]
#   #   r3 <- r[3L]
#   #   tau <- elliptic_F(pi/2, (r2-r3) / (r1-r2))
#   #   e3 <- r3
#   #   w1 <- 1/sqrt(r1 - r3) * elliptic_F(pi/2, (r2 - r3)/ (r1 - r3))
#   #   w3 <- 1i/sqrt(r1 - r3) * elliptic_F(pi/2, 1 - (r2 - r3)/ (r1 - r3))
#   #   tau <- w3/w1
#   # }else{
#     r1 <- r[1L]
#     r2 <- r[2L]
#     r3 <- r[3L]
#     e3 <- r3
#     a <- sqrt(r1 - r3)
#     b <- sqrt(r1 - r2)
#     c <- sqrt(r2 - r3)
#     if(abs(a + b) < abs(a - b)) b <- -b
#     if(abs(a + c) < abs(a - c)) c <- -c
#     if(abs(c + 1i*b) < abs(c - 1i*b)){
#       e3 <- r1
#       a <- sqrt(r3 - r1)
#       b <- sqrt(r3 - r2)
#       c <- sqrt(r2 - r1)
#       w1 <- 1 / agm(1i*b, c) 
#     }else{
#       w1 <- 1 / agm(a, b)
#     }
#     w3 <- 1i / agm(a, c)
#     tau <- w3 / w1
#   # }
#   if(Im(tau) <= 0){
#     stop("Invalid values of the parameters.")
#   }
#   z1 <- z/w1/pi #  w1 = -omega1/pi*2
#   if(derivative != 1){
#     pw0 <- e3 + 
#       (jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(z1, tau) /
#         (w1 * jtheta1_cpp(z1, tau)))**2
#     if(derivative == 0) return(pw0)
#     if(derivative == 2) return(6*pw0**2 - g2/2)
#   }
#   f <- jtheta1prime0(tau = tau)**3 /
#     (jtheta2_cpp(0, tau) * jtheta3_cpp(0, tau) * jtheta4_cpp(0, tau))
#   pw1 <- -2*(1/w1)**3 * jtheta2_cpp(z1, tau) * jtheta3_cpp(z1, tau) *
#     jtheta4_cpp(z1, tau) * f / jtheta1_cpp(z1, tau)**3
#   if(derivative == 1) return(pw1)
#   12 * pw0 * pw1 
# }
