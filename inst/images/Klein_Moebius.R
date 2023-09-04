library(jacobi)
library(RcppColors)

# modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i) / (-1 + 1i/z)
}

bkgcol <- "#15191e"

M_power_t <- function(t, gamma = 0.2 + 0.2i){
  h <- sqrt(1-Mod(gamma)^2)
  d2 <- h^t * (cos(t*pi/2) + 1i*sin(t*pi/2))
  d1 <- Conj(d2)
  a <- Re(d1) - 1i*Im(d1)/h
  b <- gamma * Im(d2)/h
  c <- Conj(b)
  d <- Conj(a)
  c(a = a, b = b, c = c, d = d)
}

t_ <- seq(0, 2, length.out = 121L)[-1L]
for(i in seq_along(t_)) {
  Mobius <- M_power_t(t_[i])
  a <- Mobius["a"]
  b <- Mobius["b"]
  c <- Mobius["c"]
  d <- Mobius["d"]
  f <- Vectorize(function(x, y){
    z <- complex(real = x, imaginary = y)
    z <- (a*z + b) / (c*z + d)
    ifelse(
      Mod(z) > 0.96, 
      NA_complex_,
      ifelse(
        y < 0, 
        kleinj(-1/PhiInv(z)), kleinj(PhiInv(z))))
  })
  
  x <- seq(-0.99, 0.99, len = 2048)
  y <- seq(-0.99, 0.99, len = 2048)
  Z <- outer(x, y, f)
  image <- colorMap3(Z, s = 95, n = 1)
  #
  svg("x.svg", width = 15, height = 15)
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(
    c(-100, 100), c(-100, 100), type = "n", 
    xlab = NA, ylab = NA, axes = FALSE, asp = 1
  )
  rasterImage(image, -100, -100, 100, 100)
  par(opar)
  dev.off()
  pngf <- sprintf("zzpic%03d.png", i)
  rsvg::rsvg_png(
    "x.svg", pngf, width = 512, height = 512
  )
}



