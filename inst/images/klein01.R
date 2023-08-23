library(jacobi)


library(PlaneGeometry)
Phi <- Mobius$new(rbind(c(1i,1), c(1,1i)))

# modified Cayley transformation
Phi <- function(z) (1i*z + 1) / (z + 1i)
PhiInv <- function(z) {
  1i + (2i) / (-1 + 1i/z)
}



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


palette <- c(
  "red", "orange", "yellow", "green", "turquoise",
  "cyan", "blue", "violet", "magenta", "red"
)
degrees <- function(z){
  arg <- Arg(z)
  if(arg < 0) arg <- pi+arg
  arg / 2 / pi
}
tone <- function(z) colorRamp(palette)(degrees(z))
b <- function(r) atan(log(r))/pi + 1/2
sat <- function(z) 1 - b(Mod(z))^2
bt <- function(z) 1 - (1-b(Mod(z)))^2
fcolor <- function(z) hsv(degrees(z), sat(z)^0.05, bt(z))



# background color
bkgcol <- rgb(21, 25, 30, maxColorValue = 255)

f <- Vectorize(function(x, y){
  q <- x + 1i*y
  if(Mod(q) >= 0.99 || (Im(q) == 0 && Re(q) <= 0)) return(bkgcol)
  tau <- -1i * log(q) / pi
  z <- kleinj(tau)
  fcolor(1/z^2)
})

x <- y <- seq(-1, 1, len = 500)
image <- outer(x, y, f)

opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)




omega1 <- 1/2
omega2 <- tau/2
omega <- c(omega1, omega2)
e1 <- wp(omega1, omega = omega)
e2 <- wp(omega2, omega = omega)
e3 <- wp(omega1+omega2, omega = omega)
tau <- omega2 / omega1
lambda(tau)
(e3-e2)/(e1-e2)

lambda2 <- function(tau) {
  omega1 <- 1/2
  omega2 <- tau/2
  omega <- c(omega1, omega2)
  e1 <- wp(omega1, omega = omega)
  e2 <- wp(omega2, omega = omega)
  e3 <- wp(omega1+omega2, omega = omega)
  (e3-e2)/(e1-e2)
}
kleinj2 <- function(tau) 
{
  lbd <- lambda2(tau)
  x <- lbd * (1 - lbd)
  256 * (1 - x)^3/x^2
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
  if(Mod(q) >= 0.95){
    return(bkgcol) # || (Im(q) == 0 && Re(q) <= 0)) return(bkgcol)
#  tau <- -1i * log(q) / pi
  } else {
    z <- if(y < 0.8) kleinj3(-1/PhiInv(q)) else kleinj3(PhiInv(q))
    fcolor(z)
  }
})

x <- seq(-0.95, 0.95, len = 1024)
y <- seq(-0.95, 0.95, len = 1024)
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
