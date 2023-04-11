library(jacobi)
library(RcppColors)

f <- Vectorize(function(x, y){
  z <- x + 1i*y
  cm(z)
})

h <- beta(1/3, 1/3)
x <- seq(-h, h, len = 1024)
y <- seq(-h, h, len = 1024)
Z <- outer(x, y, f)

img <- colorMap1(Z)

opar <- par(mar = c(0,0,0,0))
plot(c(-100, 100), c(-100, 100), type = "n", asp = 1, 
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
rasterImage(img, -100, -100, 100, 100)
par(opar)


svg("cm.svg")
opar <- par(mar = c(0,0,0,0))
plot(c(-100, 100), c(-100, 100), type = "n", asp = 1, 
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
rasterImage(img, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("cm.svg", "cm.png", width = 512, height = 512)
