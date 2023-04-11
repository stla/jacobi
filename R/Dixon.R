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
  logm <- complex(real = 0, imaginary = 5 * pi / 3)
  m <- exp(logm)
  theta <- 3^0.25 * exp(logm / 4)
  xi <- function(u) {
    s <- jellip("sn", u, m = m)
    c <- jellip("cn", u, m = m)
    d <- jellip("dn", u, m = m)
    x <- theta * s * c * d
    (x - 1) / (x + 1)
  }
  pi3 <- beta(1/3, 1/3)
  xi((z + pi3/6) / (2^(1/3) * theta))
}

#' @rdname Dixon
#' @export
cm <- function(z) {
  stopifnot(isComplex(z))
  logm <- complex(real = 0, imaginary = 5 * pi / 3)
  m <- exp(logm)
  theta <- 3^0.25 * exp(logm / 4)
  xi <- function(u) {
    s <- jellip("sn", u, m = m)
    c <- jellip("cn", u, m = m)
    d <- jellip("dn", u, m = m)
    2^(1/3) * (1 + theta*theta*s*s) / (1 + theta*s*c*d)
  }
  pi3 <- beta(1/3, 1/3)
  xi((z + pi3/6) / (2^(1/3) * theta))
}