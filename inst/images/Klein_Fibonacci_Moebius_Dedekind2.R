library(jacobi)
library(RcppColors)
library(PlaneGeometry)

# Dedekind tessellation
Circles <- list()
isInteger <- function(x) abs(x - floor(x)) < x * 1e-6
N <- 150L
for(n in 1L:N) {
  if(isInteger(n/2) && ((n/2L) %% 2L == 1L)) {
    next
  }
  for(p in 1:n) {
    q <- sqrt(n*n - p*p + 4L)
    cases <- (isInteger(q) && isInteger(q/2) && (n %% 2L == 1L)) ||
      (isInteger(q) && isInteger(q/4) && (n %% 4L == 0L))
    if(cases) {
      circ <- Circle$new(center = c(q, p)/n, radius = 2/n)
      Circles <- append(Circles, circ)
      circ <- Circle$new(center = c(-q, p)/n, radius = 2/n)
      Circles <- append(Circles, circ)
      circ <- Circle$new(center = c(q, -p)/n, radius = 2/n)
      Circles <- append(Circles, circ)
      circ <- Circle$new(center = c(-q, -p)/n, radius = 2/n)
      Circles <- append(Circles, circ)
    }
  }
}




# modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
Phi <- Mobius$new(rbind(c(1i,1), c(1,1i)))
Psi <- Phi$inverse()
PhiInv <- function(z) {
  1i + (2i*z) / (-z + 1i)
}

bkgcol <- "white"

T <- Mobius$new(rbind(c(0,-1), c(1,0)))
U <- Mobius$new(rbind(c(1,1), c(0,1)))
R <- U$compose(T)
# R^t, generalized power
Rt <- function(t) R$gpower(t)

t_ <- seq(0, 2, length.out = 121L)[-1L]
for(i in seq_along(t_)) {
  Mt <- Rt(t_[i])
  a <- Mt$a
  b <- Mt$b
  c <- Mt$c
  d <- Mt$d
  f <- function(x, y){
    z <- complex(real = x, imaginary = y)
    z <- (a*z + b) / (c*z + d)
    w <- PhiInv(z)
    ifelse(
      Mod(z) > 0.96, 
      NA_complex_,
      ifelse(
        y < 0, 
        -1/w, w
      )
    )
  }
  x <- seq(-1, 1, len = 2048)
  y <- seq(-1, 1, len = 2048)
  Z <- outer(x, y, f)
  K <- kleinj(Z) / 1728
  G <- K / (1 - K - K*K)
  image <- colorMap4(G, bkgcolor = bkgcol)
  #
  svg("x.svg", width = 15, height = 15)
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(
    c(-1, 1), c(-1, 1), type = "n", 
    xlab = NA, ylab = NA, axes = FALSE, asp = 1
  )
  rasterImage(image, -1, -1, 1, 1)
  # Dedekind tessellation
  Mob <- Mt
  CirclesImages <- vector("list", length(Circles))
  for(k in 1L:length(Circles)) {
    tcirc <- Mob$transformCircle(Circles[[k]])
    if(is(tcirc, "Circle")) {
      draw(tcirc, border = "white", lwd = 2)
    } else {
      draw(tcirc, col = "white", lwd = 2)
    }
  }
  par(opar)
  dev.off()
  pngf <- sprintf("zzpic%03d.png", i)
  rsvg::rsvg_png(
    "x.svg", pngf, width = 512, height = 512
  )
}

library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  pngs,
  "KleinMoebius_cm3.gif",
  width = 512, height = 512,
  delay = 1/11
)

file.remove(pngs)

