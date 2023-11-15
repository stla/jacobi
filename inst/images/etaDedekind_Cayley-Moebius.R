library(jacobi)
library(RcppColors)

# background color
bkgcol <- "#ffffff"

# Möbius transformation of order 3
Mob <- function(z, t) {
  a <- pi*t/3
  ((sqrt(3)*cos(a) - sin(a)) * z - 2*sin(a))/
    (2*sin(a) * z + sqrt(3)*cos(a) + sin(a))
}

# map the disk to H and calculate eta
f <- function(x, y, t) {
  z <- complex(real = x, imaginary = y)
  K <- rep(NA_complex_, length(z))
  inDisk <- Mod(z) < 1
  K[inDisk] <- eta(disk2H(z[inDisk]))
  K
}
n <- 512L
x <- y <- seq(-1, 1, length.out = n)
Grid <- expand.grid(X = x, Y = y)
K <- f(Grid$X, Grid$Y, 0.05)^24
dim(K) <- c(n, n)

# plot
img <- colorMap5(K)
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)


svg("x.svg", width = 15, height = 15)
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "etaDedekindOnCircle.png", 512, 512)
