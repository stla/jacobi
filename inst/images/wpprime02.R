library(jacobi)
library(RcppColors)

x <- seq(0, 1, len = 1024)
y <- seq(0, 1, len = 1024)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(1/2, 1i/2), derivative = 0L)
}

Z <- outer(x, y, f)
img <- colorMap1(Z)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)


svg("x.svg", width = 12, height = 12)
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(
  c(-100, 100), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(img, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "wp.png", width = 1024, height = 1024)

