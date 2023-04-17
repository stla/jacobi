library(jacobi)
library(RcppColors)

f <- Vectorize(function(x, y){
  w <- x + 1i*y
  sl(w)
})

x <- seq(-3, 3, len = 1024)
y <- seq(-3, 3, len = 1024)
Z <- outer(x, y, f)

img <- colorMap1(Z)

opar <- par(mar = c(0,0,0,0))
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)


svg("sl.svg")
opar <- par(mar = c(0,0,0,0))
plot(c(-100, 100), c(-100, 100), type = "n", asp = 1, 
     xlab = "", ylab = "", axes = FALSE, xaxs = "i", yaxs = "i")
rasterImage(img, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("sl.svg", "sl.png", width = 512, height = 512)
