library(jacobi)
library(RcppColors)

bkgcol <- "white"

M_power_t <- function(t, gamma = -0.7 - 0.7i){
  h <- sqrt(1-Mod(gamma)^2)
  d2 <- (exp(pi*1i/2)*h)^t#h^t * (cospi(t/2) + 1i*sinpi(t/2))
  d1 <- Conj(d2)
  a <- Re(d1) - 1i*Im(d1)/h
  b <- gamma * Im(d2)/h
  c <- Conj(b)
  d <- Conj(a)
  c(a = a, b = b, c = c, d = d)
}



t_ <- seq(0, 2, length.out = 121L)[-1L]
for(i in seq_along(t_)) {
  Mt <- M_power_t(t_[i])
  a <- Mt["a"]
  b <- Mt["b"]
  c <- Mt["c"]
  d <- Mt["d"]
  f <- function(x, y) {
    z <- complex(real = x, imaginary = y)
    z <- (a*z + b) / (c*z + d)
    ifelse(
      Mod(z) >= 0.9999 | (Im(z) == 0 & Re(z) <= 0), 
      NA_complex_,
      -1i * log(z) / pi)
  }
  x <- seq(-1, 1, len = 1024)
  y <- seq(-1, 1, len = 1024)
  Z <- outer(x, y, f)
  K <- kleinj(Z) / 1728
  G <- 1 / (1/K - 1 - K)
  image <- colorMap4(G, bkgcolor = bkgcol)
  #
  svg("x.svg", width = 15, height = 15)
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(
    c(-1, 1), c(-1, 1), type = "n", 
    xlab = NA, ylab = NA, axes = FALSE, asp = 1
  )
  rasterImage(image, -1, -1, 1, 1)
  par(opar)
  dev.off()
  pngf <- sprintf("zzpic%03d.png", i)
  rsvg::rsvg_png(
    "x.svg", pngf, width = 1024, height = 1024
  )
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  pngs,
  "KleinFibonacciMoebius_cm4.gif",
  width = 512, height = 512,
  delay = 1/8
)

file.remove(pngs)

