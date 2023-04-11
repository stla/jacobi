library(jacobi)
library(RcppColors)

f <- Vectorize(function(x, y){
  w <- x + 1i*y
  slh(w)
})

x <- seq(-1, 1, len = 200)
y <- seq(-1, 1, len = 200)
Z <- outer(x, y, f)

img <- colorMap1(Z)

opar <- par(mar = c(0,0,0,0))
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)
