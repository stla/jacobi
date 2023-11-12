library(jacobi)
library(RcppColors)

# the modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i*z) / (1i - z)
}

# background color
bkgcol <- "#ffffff"


# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  tau <- PhiInv(z)
  ifelse(
    Mod(z) > 0.95, 
    NA_complex_,
    wzeta(z, omega=c(0.5, 0.5*tau))
  )
}
x <- y <- seq(-1, 1, length.out = 256L)
S <- outer(x, y, Vectorize(f))

img <- colorMap5(S)

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)








# animation ####
t_ <- seq(-2*pi, pi, length.out = 121)[-1L]

for(i in seq_along(t_)) {
  print(i)
  img <- colorMap5(Mob(R, s(t_[i])), bkgcolor = bkgcol)
  svg("x.svg", width = 16, height = 16)
  opar <- par(mar = c(0, 0, 0, 0))
  plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
       axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
  rasterImage(img, 0, 0, 1, 1)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "x.svg", file = sprintf("zzpic%03d.png", i), 
    width = 512, height = 512
  )
}

# mount animation
library(gifski)
pngs <- Sys.glob("zzpic*.png")
gifski(
  pngs,
  "RogersRamanujanMobius.gif",
  width = 512, height = 512,
  delay = 1/11
)

file.remove(pngs)
