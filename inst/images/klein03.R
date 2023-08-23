library(jacobi)
library(RcppColors)

# modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i) / (-1 + 1i/z)
}

kleinj3 <- function(tau) {
  j2 <- jacobi:::jtheta2_cpp(0, tau)
  j3 <- jacobi:::jtheta3_cpp(0, tau)
  g2 <- 4/3 * (pi/2)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  g3 <- 8/27 * (pi/2)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  ) + j3**12)
  g2cube <- g2*g2*g2
  1728 * g2cube / (g2cube - 27*g3*g3)
}

f <- Vectorize(function(x, y){
  q <- x + 1i*y
  if(Mod(q) >= 0.96){
    return(bkgcol) # || (Im(q) == 0 && Re(q) <= 0)) return(bkgcol)
#  tau <- -1i * log(q) / pi
  } else {
    z <- if(y < 0.8) kleinj(-1/PhiInv(q)) else kleinj(PhiInv(q))
    colorMap1(1/z)
  }
})

x <- seq(-0.99, 0.99, len = 1024)
y <- seq(-0.99, 0.99, len = 1024)
image <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)


svg("x.svg", width = 10, height = 10)
x <- seq(-0.95, 0.95, len = 2048)
y <- seq(-0.95, 0.95, len = 2048)
image <- outer(x, y, f)
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "Klein01.png", width = 512, height = 512)
