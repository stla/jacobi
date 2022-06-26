library(jacobi)

# background color
bkgcol <- rgb(21, 25, 30, maxColorValue = 255)


modulo <- function(a, p) {
  a - p * ifelse(a > 0, floor(a/p), ceiling(a/p))
}

colormap1 <- function(z){
  if(is.na(z)) return(bkgcol)
  if(is.infinite(z) || is.nan(z)) return("#000000")
  x <- Re(z)
  y <- Im(z)
  a <- atan2(y, x)
  r <- modulo(Mod(z), 1)
  g <- abs(modulo(a, 0.5)) * 2
  b <- abs(modulo(x*y, 1))
  if(is.nan(b)){
    return("#000000")
  }
  rgb(
    (1.0 - cos(r-0.5))*8.0, 
    (1.0 - cos(g-0.5))*8.0, 
    (1.0 - cos(b-0.5))*8.0,
    maxColorValue = 1
  )
}

colormap2 <- function(z){
  if(is.na(z)) return(bkgcol)
  if(is.infinite(z) || is.nan(z)) return("#000000")
  arg <- Arg(z)
  if(arg < 0) arg <- 2*pi + arg
  h <- arg / 2 / pi
  x <- 2*pi*log(1+Mod(z))
  s <- sqrt((1 + sin(x))/2)
  v <- (1 + cos(x))/2
  hsv(h, s, v)
}

f <- Vectorize(function(x, y){
  w <- x + 1i*y
  z <- sigmaw(w, omega = c(1/2, 1i/2))
  colormap1(z)
})

x <- seq(-1.5, 1.5, len = 1000)
y <- seq(-1.5, 1.5, len = 1000)
z <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)

svg("sigma.svg")
opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("sigma.svg", "sigma.png", width = 512, height = 512)

