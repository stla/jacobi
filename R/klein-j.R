#' @title Klein j-function and its inverse
#' @description Evaluation of the Klein j-invariant function and its inverse.
#'
#' @param tau a complex number with strictly positive imaginary part, or a 
#'   vector or matrix of such complex numbers; missing values allowed
#' @param j a complex number
#' @param transfo Boolean, whether to use a transformation of the values  
#'   of \code{tau} close to the real line; using this option can fix some 
#'   failures of the computation (at the cost of speed), e.g. when the 
#'   algorithm reaches the maximal number of iterations
#'
#' @return A complex number, vector or matrix.
#' @export
#' 
#' @note The Klein-j function is the one with the factor 1728.
#'
#' @examples
#' ( j <- kleinj(2i) )
#' 66^3
#' kleinjinv(j)
kleinj <- function(tau, transfo = FALSE){
  stopifnot(isBoolean(transfo))
  lbd <- lambda(tau, transfo) 
  x <- lbd * (1 - lbd)
  256 * (1-x)^3 / x^2
  # stopifnot(isComplexNumber(tau))
  # if(Im(tau) <= 0){
  #   stop("The complex number `tau` must have a positive imaginary part.")
  # }
  # j2 <- jtheta2_cpp(0, tau)
  # j3 <- jtheta3_cpp(0, tau)
  # g2 <- 4/3 * (pi/2)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  # g3 <- 8/27 * (pi/2)**6 * (j2**12 - (
  #   (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  # ) + j3**12)
  # g2cube <- g2*g2*g2
  # 1728 * g2cube / (g2cube - 27*g3*g3)
}

#' @rdname kleinj
#' @export
kleinjinv <- function(j) {
  stopifnot(isComplexNumber(j))
  if(is.infinite(j)) {
    x <- 0
  } else if(j == 0) {
    return(0.5 + sqrt(3)/2*1i)
  } else {
    # x <- polyroot(c(1, -3, 3-j/256, -1))[1L]
    t <- -j*j*j + 2304 * j*j + 12288 * 
      sqrt(3 * (1728 * j*j - j*j*j)) - 884736 * j
    u <- t^(1/3)
    x <- 1/768 * u - (1536 * j - j*j) / (768 * u) + (1 - j/768)
  }
  # lbd <- polyroot(c(-x, 1, -1))[1L]
  lbd <- -(-1 - sqrt(1-4*x)) / 2
  1i * agm(1, sqrt(1-lbd)) / agm(1, sqrt(lbd))
}