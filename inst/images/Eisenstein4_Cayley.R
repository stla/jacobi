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
    Mod(z) >= 1, 
    NA_complex_,
    EisensteinE(4L, exp(1i*pi*tau))
  )
}
x <- y <- seq(-1, 1, length.out = 512L)
S <- outer(x, y, Vectorize(f))


img <- colorMap7(Im(S))
img <- colorMap5(S)

opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)

############################################################################
n <- 2048L
x <- y <- seq(0.0001, 0.9999, length.out = n)
Grid <- transform(
  expand.grid(X = x, Y = y), 
  Z = complex(real = X, imaginary = Y)
)
Tau <- square2H(Grid$Z)
Q <- exp(1i*pi*Tau)
K <- vapply(Q, function(q) EisensteinE(6L, q), complex(1L))
dim(K) <- c(n, n)
# plot
if(require("RcppColors")) {
  img <- colorMap5(K)
} else {
  img <- as.raster((Arg(K) + pi)/(2*pi))
}
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
  (xmsinx(xmsinx(2*x)) + 4*pi) / (2*pi)  
}


# frames
t_ <- seq(-2*pi, pi, length.out = 161)[-1L]
for(i in seq_along(t_)) {
  print(i)
  img <- colorMap5(Mob(K, s(t_[i])))#, bkgcolor = bkgcol)
  #img <- permuteRGB(img, "brg")
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
  "Eisenstein6OnSquare.gif",
  width = 512, height = 512,
  delay = 1/10
)

file.remove(pngs)
