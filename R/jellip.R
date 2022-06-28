isString <- function(x){
  is.character(x) && length(x) == 1L && !is.na(x)
}

#' @title Jacobi elliptic functions
#' @description Evaluation of the Jacobi elliptic functions.
#'
#' @param kind a string with two characters among \code{"s"}, \code{"c"}, 
#'   \code{"d"} and \code{"n"}; this string specifies the function: the two 
#'   letters respectively denote the basic functions \eqn{sn}, \eqn{cn}, 
#'   \eqn{dn} and \eqn{1}, and the string specifies the ratio of two such 
#'   functions, e.g. \eqn{ns = 1/sn} and \eqn{cd = cn/dn}   
#' @param u a complex number, vector or matrix
#' @param tau complex number with strictly positive imaginary part; it is 
#'   related to \code{m} and only one of them must be supplied
#' @param m the "parameter", square of the elliptic modulus; it is related to 
#'   \code{tau} and only one of them must be supplied
#'
#' @return A complex number, vector or matrix.
#' @export
#' 
#' @examples 
#' u <- 2 + 2i
#' tau <- 1i
#' jellip("cn", u, tau)^2 + jellip("sn", u, tau)^2 # should be 1
jellip <- function(kind, u, tau = NULL, m = NULL){
  stopifnot(isString(kind))
  stopifnot(nchar(kind) == 2L)
  f1 <- substr(kind, 1L, 1L)
  f2 <- substr(kind, 2L, 2L)
  stopifnot(f1 %in% c("s", "c", "d", "n"), f2 %in% c("s", "c", "d", "n"))
  ntheta1 <- switch(
    f1,
    s = theta.s(u, tau = tau, m = m),
    c = theta.c(u, tau = tau, m = m),
    d = theta.d(u, tau = tau, m = m),
    n = theta.n(u, tau = tau, m = m)
  )
  ntheta2 <- switch(
    f2,
    s = theta.s(u, tau = tau, m = m),
    c = theta.c(u, tau = tau, m = m),
    d = theta.d(u, tau = tau, m = m),
    n = theta.n(u, tau = tau, m = m)
  )
  ntheta1 / ntheta2
}