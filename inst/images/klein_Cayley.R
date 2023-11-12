library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
disk2H <- function(z) {
  1i + (2i*z) / (1i - z)
}
disk2H <- function(z) {
  1i * (1-z)/(1+z)
}

# background color
bkgcol <- "#ffffff"

# map the disk to H and calculate kleinj
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  K <- rep(NA_complex_, length(x))
  inDisk <- Mod(z) < 1
  K[inDisk] <- kleinj(disk2H(z[inDisk]))
  K
}
n <- 1024L
x <- y <- seq(-1, 1, length.out = n)
Grid <- expand.grid(X = x, Y = y)
K <- f(Grid$X, Grid$Y)
dim(K) <- c(n, n)
# plot
if(require("RcppColors")) {
  img <- colorMap5(K)
} else {
  img <- as.raster(1 - abs(Im(K))/Mod(K))
}
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)

