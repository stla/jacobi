library(jacobi)
library(RcppColors)


square2H <- function(z) {
  z <- z*pi
  tau <- -jtheta2(z, 2i) / jtheta3(z, 2i) # or jtheta_cpp
  kleinj(tau)
}

n <- 1024L
x <- y <- seq(0.0001, 0.9999, length.out = n)
Grid <- transform(
  expand.grid(X = x, Y = y), 
  Z = complex(real = X, imaginary = Y)
)
K <- square2H(Grid$Z)
dim(K) <- c(n, n)
# plot
if(require("RcppColors")) {
  img <- colorMap5(K)
} else {
  img <- as.raster((Arg(K) + pi)/(2*pi))
}
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)

