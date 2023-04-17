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


# lemniscate parameterization
p <- Vectorize(function(s) {
  a <- Re(cl(s))
  b <- Re(sl(s))
  c(a, a * b) / sqrt(1 + b*b)
})
# lemnniscate constant
ombar <- 2.622 # gamma(1/4)^2 / (2 * sqrt(2*pi))
# plot
s_ <- seq(0, ombar, length.out = 100)
lemniscate <- t(p(s_))
plot(lemniscate, type = "l", col = "blue", lwd = 3)
lines(cbind(lemniscate[, 1L], -lemniscate[, 2L]), col="red", lwd = 3)

