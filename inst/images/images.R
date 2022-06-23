library(jacobi)

modulo <- function(a, p) {
  a - p * ifelse(a > 0, floor(a/p), ceiling(a/p))
}

colormap1 <- function(z){
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
  if(is.infinite(z) || is.nan(z)) return("#000000")
  arg <- Arg(z)
  if(arg < 0) arg <- 2*pi + arg
  h <- arg / 2 / pi
  s <- sqrt((1 + sin(2*pi*log(1+Mod(z))))/2)
  v <- (1 + cos(2*pi*log(1+Mod(z))))/2
  hsv(h, s, v)
}

MyMatPow <- function(gamma, t){
  g2 <- Mod(gamma)^2
  h <- sqrt(1-g2)
  d2 <- h^t * (cos(t*pi/2) + 1i*sin(t*pi/2))
  d1 <- Conj(d2)
  H11 <- Re(d1) - 1i*Im(d1)/h
  H12 <- Im(d2) * gamma / h
  H21 <- Conj(H12)
  H22 <- Conj(H11)
  c(a = H11, b = H12, c = H21, d = H22)
}

# background color
bkgcol <- rgb(21, 25, 30, maxColorValue = 255)

M <- MyMatPow(0.3 + 0.3i, pi/5)
a <- M["a"]; b <- M["b"]; c <- M["c"]; d <- M["d"]; 
f <- Vectorize(function(x, y){
  q0 <- x + 1i*y
  q <- (a*q0 + b) / (c*q0 + d)
  if(Mod(q) >= 0.99 || (Im(q) == 0 && Re(q) <= 0)) return(bkgcol)
  z <- En(4, q)
  colormap2(z)
})

x <- y <- seq(-1, 1, len = 200)
image <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)




x <- y <- seq(-1, 1, len = 2000)
t_ <- head(seq(0, 2*pi, len = 120), -1L)
for(i in 1:length(t_)){
  M <- MyMatPow(0.3 + 0.3i, t_[i])
  a <- M["a"]; b <- M["b"]; c <- M["c"]; d <- M["d"]; 
  image <- outer(x, y, f)
  svg("zzz.svg")
  opar <- par(mar = c(0,0,0,0), bg = rgb(21,25,30,maxColorValue = 255))
  plot(c(-100, 100), c(-100, 100), type = "n", 
       xlab = "", ylab = "", axes = FALSE, asp = 1)
  rasterImage(image, -100, -100, 100, 100)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "zzz.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512
  )
}


a <- -0.5
p <- 3
modulo(a, p)
jacobi:::modulo(a, p)

fcolor <- function(z){
  if(is.infinite(z) || is.nan(z)) return("#000000")
  x <- Re(z)
  y <- Im(z)
  a = atan2(y, x)
  r <- jacobi:::modulo(Mod(z), 1)
  g <- abs(jacobi:::modulo(a, 0.5)) * 2
  b <- abs(jacobi:::modulo(x*y, 1))
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


# used for 1/Delta and E6
degrees <- function(z){
  arg <- Arg(z)
  if(arg < 0) arg <- 2*pi+arg
  arg / 2 / pi
}
fcolor <- function(z){
  if(is.infinite(z) || is.nan(z)) return("#000000")
  h = degrees(z)
  s = sqrt((1+sin(2*pi*log(1+Mod(z))))/2)
  v = (1+cos(2*pi*log(1+Mod(z))))/2
  hsv(h,s,v);
}

Delta <- function(tau){
  g <- jacobi:::g_from_omega(1, tau)
  g2 <- g[1L]
  g3 <- g[2L]
  g2^3 - 27*g3^2
}

palette <- c(
  "red", "orange", "yellow", "green", "turquoise",
  "cyan", "blue", "violet", "magenta", "red"
)
degrees <- function(z){
  arg <- Arg(z)
  if(arg < 0) arg <- 2*pi+arg
  arg / 2 / pi
}
tone <- function(z) colorRamp(palette)(degrees(z))
b <- function(r) atan(log(r))/pi + 1/2
sat <- function(z) 1 - b(Mod(z))^2
bt <- function(z) 1 - (1-b(Mod(z)))^2
fcolor <- function(z) hsv(degrees(z), sat(z), bt(z)^0.25)

bgcol <- rgb(21,25,30,maxColorValue = 255)
f <- Vectorize(function(x, y){
  q <- x + 1i*y
  if(Mod(q) >= 0.99 || (Im(q) == 0 && Re(q) <= 0)) return(bgcol)
  # tau <- -1i * log(q) / pi
  z <- En(2, q)
  colormap2(z)
})

x <- seq(-1, 1, len = 200)
y <- seq(-1, 1, len = 200)

z <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)

svg("Eisenstein6.svg")
opar <- par(mar = c(0,0,0,0), bg = rgb(21,25,30,maxColorValue = 255))
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("Eisenstein6.svg", "Eisenstein6.png", width = 512, height = 512)


image(z, axes = FALSE, col = hcl.colors(200))