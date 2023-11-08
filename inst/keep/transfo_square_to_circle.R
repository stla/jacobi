library(jacobi)
library(RcppColors)


square2circle <- function(z) {
  z <- z*pi
  w <- -jtheta2(z, 2i) / jtheta3(z, 2i) # or jtheta_cpp
  (1i - w)/(1i + w)
}

x <- y <- seq(0, 1, length.out = 25L)
Grid <- transform(
  expand.grid(X = x, Y = y), 
  Z = complex(real = X, imaginary = Y)
)
u <- square2circle(Grid$Z)
plot(u, pch = 19, asp = 1)


# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  ifelse(
    x %in% c(0, 1) || y %in% c(0,1), 
    NA_complex_,
    square2circle(z)
  )
}
x <- y <- seq(0, 1, length.out = 56L)
L <- outer(x, y, Vectorize(f))

img <- colorMap5(L)
img <- permuteRGB(img, "brg")

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1, interpolate = FALSE)
par(opar)

