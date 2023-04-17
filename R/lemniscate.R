#' @title Lemniscate functions
#' @description Lemniscate sine, cosine, arcsine, arccosine, hyperbolic sine, 
#'   and hyperbolic cosine functions.
#' 
#' @param z a real number or a complex number
#'
#' @return A complex number. 
#' @importFrom Carlson elliptic_F
#' @export
#' @name lemniscate
#' @rdname lemniscate
#'
#' @examples
#' sl(1+1i) * cl(1+1i) # should be 1
#' ## | the lemniscate ####
#' # lemniscate parameterization
#' p <- Vectorize(function(s) {
#'   a <- Re(cl(s))
#'   b <- Re(sl(s))
#'   c(a, a * b) / sqrt(1 + b*b)
#' })
#' # lemnniscate constant
#' ombar <- 2.622 # gamma(1/4)^2 / (2 * sqrt(2*pi))
#' # plot
#' s_ <- seq(0, ombar, length.out = 100)
#' lemniscate <- t(p(s_))
#' plot(lemniscate, type = "l", col = "blue", lwd = 3)
#' lines(cbind(lemniscate[, 1L], -lemniscate[, 2L]), col="red", lwd = 3)
sl <- function(z) {
  stopifnot(isComplex(z))
  ombar <- gamma(1/4)^2 / (2 * sqrt(2*pi)) # lemniscate constant
  w <- as.complex(z / ombar)
  jtheta1_cpp(w, 1i) / jtheta3_cpp(w, 1i)
}

#' @export
#' @rdname lemniscate
cl <- function(z) {
  stopifnot(isComplex(z))
  ombar <- gamma(1/4)^2 / (2 * sqrt(2*pi)) # lemniscate constant
  w <- as.complex(z / ombar)
  jtheta2_cpp(w, 1i) / jtheta4_cpp(w, 1i)
}

#' @export
#' @rdname lemniscate
asl <- function(z) {
  z <- as.complex(z)
  elliptic_F(asin(sqrt(2)*z / sqrt(1+z*z)), 0.5) / sqrt(2)
}

#' @export
#' @rdname lemniscate
acl <- function(z) {
  z <- as.complex(z)
  elliptic_F(acos(z), 0.5) / sqrt(2)
}

#' @export
#' @rdname lemniscate
slh <- function(z) {
  z <- z / sqrt(2)
  cosl <- cl(z)
  (1 + cosl*cosl) * sl(z) / (sqrt(2) * cosl)
}

#' @export
#' @rdname lemniscate
clh <- function(z) {
  z <- z / sqrt(2)
  sinl <- sl(z)
  (1 + sinl*sinl) * cl(z) / (sqrt(2) * sinl)
}
