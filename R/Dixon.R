#' @title Dixon elliptic functions
#' @description The Dixon elliptic functions.
#' 
#' @param z a real or complex number 
#'
#' @return A complex number.
#' @export
#' @name Dixon
#' @rdname Dixon
#'
#' @examples
#' # cubic Fermat curve x^3+y^3=1
#' pi3 <- beta(1/3, 1/3)
#' epsilon <- 0.7
#' t_ <- seq(-pi3/3 + epsilon, 2*pi3/3 - epsilon, length.out = 100)
#' pts <- t(vapply(t_, function(t) {
#'   c(Re(cm(t)), Re(sm(t)))
#' }, FUN.VALUE = numeric(2L)))
#' plot(pts, type = "l", asp = 1)
sm <- function(z) {
  stopifnot(isComplex(z))
  -6 * wp(z, g = c(0, 1/27)) / 
    (3 * wp(z, g = c(0, 1/27), derivative = 1L) - 1)
}

#' @rdname Dixon
#' @export
cm <- function(z) {
  stopifnot(isComplex(z))
  x <- 3 * wp(z, g = c(0, 1/27), derivative = 1L)
  (x + 1) / (x - 1)
}