library(jacobi)
library(RcppColors)

# the conformal map
Phi <- function(z) {
  -jtheta2(z, 2i) / jtheta3(z, 2i)
}

# background color
bkgcol <- "#ffffff"

# make the color mapping
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  ifelse(
    x %in% c(0, 1) || y %in% c(0,1), 
    NA_complex_,
    jacobi:::jtheta3(z,Phi(z*pi))
  )
}
x <- y <- seq(0, 1, length.out = 256L)
L <- outer(x, y, Vectorize(f))

img <- colorMap5(1/L)
img <- permuteRGB(img, "brg")

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)


# animation ####
# Möbius transformation of order 3
Mob <- function(z, t) {
  a <- pi*t/3
  ((sqrt(3)*cos(a) - sin(a)) * z - 2*sin(a))/
    (2*sin(a) * z + sqrt(3)*cos(a) + sin(a))
}

# smooth stair function
xmsinx <- function(x) x - sin(x)
s <- function(x) {
  (xmsinx(xmsinx(3*x)) + 6*pi) / (4*pi)  
}

# frames
t_ <- seq(-2*pi, 2*pi, length.out = 161)[-1L]
#t_ <- seq(0, 3, length.out = 61L)[-1L]
for(i in seq_along(t_)) {
  print(i)
  img <- colorMap5(Mob(L, s(t_[i])))#, bkgcolor = bkgcol)
  img <- permuteRGB(img, "brg")
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
  "lambdaOnSquare.gif",
  width = 512, height = 512,
  delay = 1/11
)

file.remove(pngs)
