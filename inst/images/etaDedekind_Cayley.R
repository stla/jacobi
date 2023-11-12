library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

# background color
bkgcol <- "#ffffff"


# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  tau <- PhiInv(z)
  ifelse(
    Mod(z) >= 1, 
    NA_complex_,
    eta(tau)
  )
}
x <- y <- seq(-1, 1, length.out = 512L)
S <- outer(x, y, Vectorize(f))

img <- colorMap5(S^24)
img <- permuteRGB(img, "gbr")

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)

