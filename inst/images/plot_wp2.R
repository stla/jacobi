library(jacobi)

# background color
bkgcol <- rgb(21, 25, 30, maxColorValue = 255)

x <- seq(-4, 4, len = 1500)
y <- seq(-4, 4, len = 1500)
W <- outer(y, x, function(x, y) complex(real = x, imaginary = y))

# Z <- wzeta(W, omega = c(-1, -1i))

Z <- matrix(Bessel::BesselY(W, nu = 3), nrow = nrow(W), ncol = ncol(W))
image <- jacobi:::ColorMap1(Z)


opar <- par(mar = c(0,0,0,0), bg = "#15191E")
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)

svg("BesselY.svg")
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png(
  "BesselY.svg", "BesselY.png", width = 512, height = 512
)

