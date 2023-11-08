#' @title Square to disk
#' @description Conformal map from the unit square to the unit disk. The 
#'   function is vectorized.
#' 
#' @param z a complex number in the unit square \eqn{[0,1] \times [0,1]}
#'
#' @return A complex number in the unit disk.
#' @export
#'
#' @examples
#' x <- y <- seq(0, 1, length.out = 25L)
#' Grid <- transform(
#'   expand.grid(X = x, Y = y), 
#'   Z = complex(real = X, imaginary = Y)
#' )
#' u <- square2disk(Grid$Z)
#' plot(u, pch = 19, asp = 1)
square2disk <- function(z) {
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  if(length(z) == 1L) {
    w <- -jtheta2_cpp(z, 2i) / jtheta3_cpp(z, 2i)
  } else {
    if(!is.matrix(z)) {
      w <- -JTheta2(cbind(z), 2i)[, 1L] / JTheta3(cbind(z), 2i)[, 1L]
    } else {
      w <- -JTheta2(z, 2i) / JTheta3(z, 2i)
    }
  }
  (1i - w)/(1i + w)
}

#' @title Square to upper half-plane
#' @description Conformal map from the unit square to the upper half-plane. 
#'   The function is vectorized.
#' 
#' @param z a complex number in the unit square \eqn{[0,1] \times [0,1]}
#'
#' @return A complex number in the upper half-plane.
#' @export
#' 
#' @examples
#' n <- 1024L
#' x <- y <- seq(0.0001, 0.9999, length.out = n)
#' Grid <- transform(
#'   expand.grid(X = x, Y = y), 
#'   Z = complex(real = X, imaginary = Y)
#' )
#' \donttest{K <- kleinj(square2H(Grid$Z))
#' dim(K) <- c(n, n)
#' # plot
#' if(require("RcppColors")) {
#'   img <- colorMap5(K)
#' } else {
#'   img <- as.raster((Arg(K) + pi)/(2*pi))
#' }
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
#'      axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
#' rasterImage(img, 0, 0, 1, 1)
#' par(opar)}
square2H <- function(z) {
  stopifnot(isComplex(z))
  storage.mode(z) <- "complex"
  if(length(z) == 1L) {
    -jtheta2_cpp(z, 2i) / jtheta3_cpp(z, 2i)
  } else {
    if(!is.matrix(z)) {
      -JTheta2(cbind(z), 2i)[, 1L] / JTheta3(cbind(z), 2i)[, 1L]
    } else {
      -JTheta2(z, 2i) / JTheta3(z, 2i)
    }
  }
}

#' @title Disk to square
#' @description Conformal map from the unit disk to the square
#'   \eqn{[-1,1] \times [-1,1]}. The function is not vectorized.
#'   
#' @param z a complex number in the unit disk
#'
#' @return A complex number in the square \eqn{[-1,1] \times [-1,1]}.
#' @export
#' @importFrom Carlson elliptic_F
#'
#' @examples
#' n <- 70L
#' r <- seq(0, 1, length.out = n)
#' theta <- seq(0, 2*pi, length.out = n+1L)[-1L]
#' Grid <- transform(
#'   expand.grid(R = r, Theta = theta),
#'   Z = R*exp(1i*Theta)
#' )
#' s <- vapply(Grid$Z, disk2square, complex(1L))
#' plot(Re(s), Im(s), pch = ".", asp = 1, cex = 2)
#' #
#' # a more insightful plot ####
#' r_ <- seq(0, 1, length.out = 10L)
#' theta_ <- seq(0, 2*pi, length.out = 33)[-1L]
#' plot(
#'   NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, xlab = "x", ylab = "y"
#' )
#' for(r in r_) {
#'   theta <- sort(
#'     c(seq(0, 2, length.out = 200L), c(1/4, 3/4, 5/4, 7/4))
#'   )
#'   z <- r*(cospi(theta) + 1i*sinpi(theta))
#'   s <- vapply(z, disk2square, complex(1L))
#'   lines(Re(s), Im(s), col = "blue", lwd = 2)
#' }
#' for(theta in theta_) {
#'   r <- seq(0, 1, length.out = 30L)
#'   z <- r*exp(1i*theta)
#'   s <- vapply(z, disk2square, complex(1L))
#'   lines(Re(s), Im(s), col = "green", lwd = 2)
#' }
disk2square <- function(z) {
  C <- 2*gamma(5/4)^2 / sqrt(pi)
  sqrt2half <- sqrt(2) / 2
  w <- complex(real = sqrt2half, imaginary = sqrt2half)
  -w * elliptic_F(1i * asinh(w * z), -1) / C  
}


#' @title Disk to upper half-plane
#' @description Conformal map from the unit disk to the upper half-plane. 
#'   The function is vectorized.
#'   
#' @param z a complex number in the unit disk
#'
#' @return A complex number in the upper half-plane.
#' @export
#'
#' @examples
#' # map the disk to H and calculate kleinj
#' f <- function(x, y) {
#'   z <- complex(real = x, imaginary = y)
#'   K <- rep(NA_complex_, length(x))
#'   inDisk <- Mod(z) < 1
#'   K[inDisk] <- kleinj(disk2H(z[inDisk]))
#'   K
#' }
#' n <- 1024L
#' x <- y <- seq(-1, 1, length.out = n)
#' Grid <- expand.grid(X = x, Y = y)
#' \donttest{K <- f(Grid$X, Grid$Y)
#' dim(K) <- c(n, n)
#' # plot
#' if(require("RcppColors")) {
#'   img <- colorMap5(K)
#' } else {
#'   img <- as.raster(1 - abs(Im(K))/Mod(K))
#' }
#' opar <- par(mar = c(0, 0, 0, 0))
#' plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
#'      axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
#' rasterImage(img, 0, 0, 1, 1)
#' par(opar)}
disk2H <- function(z) {
  1i + (2i*z) / (1i - z)
}