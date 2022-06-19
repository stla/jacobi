library(jacobi)

modulo <- function(a, p) {
  a - (ifelse(a>0, floor(a/p), ceiling(a / p)) * p);
}

a <- -0.5
p <- 3
modulo(a, p)
jacobi:::modulo(a, p)

fcolor <- function(z){
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


# used for 1/Delta
degrees <- function(z){
  arg <- Arg(z)
  if(arg < 0) arg <- 2*pi+arg
  arg / 2 / pi
}
fcolor <- function(z){
  h = degrees(z);
  s = sqrt((1+sin(2*pi*log(1+Mod(z))))/2);
  v = (1+cos(2*pi*log(1+Mod(z))))/2;
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

f <- Vectorize(function(x, y){
  q <- x + 1i*y
  if(Mod(q) >= 0.99 || (Im(q) == 0 && Re(q) <= 0)) return("#000000")
  tau <- -1i * log(q) / pi
  z <- 1/kleinj(tau)
  if(is.infinite(z) || is.nan(z)) return("#000000")
  fcolor(z)
})

x <- seq(-1, 1, len = 2000)
y <- seq(-1, 1, len = 2000)

z <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)

svg("KleinInverse.svg")
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(c(-100, 100), c(-100, 100), type = "n", 
     xlab = "", ylab = "", axes = FALSE, asp = 1)
image <- z
rasterImage(image, -100, -100, 100, 100)
par(opar)
dev.off()

rsvg::rsvg_png("KleinInverse.svg", "KleinInverse.png", width = 512, height = 512)


image(z, axes = FALSE, col = hcl.colors(200))