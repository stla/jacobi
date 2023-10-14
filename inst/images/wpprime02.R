library(jacobi)
library(RcppColors)

x <- seq(-1, 1, len = 1024)
y <- seq(-1, 1, len = 1024)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, g = c(189, 0), derivative = 0L)
}

Z <- outer(x, y, f)
img <- colorMap1(Z)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)


svg("x.svg", width = 16, height = 8)
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(
  c(-200, 200), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(img, -200, -100, 200, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "wp2x1.png", width = 1024, height = 512)

