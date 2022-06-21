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

n <- 200L
x <- seq(-1, 1, len = n)
Z <- jacobi:::MOB(x, 0.3+0.3i, pi/3)
image <- apply(Z, c(1,2), colormap2)


n <- 2000L
x <- seq(-1, 1, len = n)
image <- jacobi:::Image(x, 0.3+0.3i, pi/3)

opar <- par(mar = c(0,0,0,0), bg = bkgcol)
plot(
  c(-100, 100), c(-100, 100), type = "n", 
  xlab = "", ylab = "", axes = FALSE, asp = 1
)
rasterImage(image, -100, -100, 100, 100)
par(opar)







x <- seq(-1, 1, len = 2000)
t_ <- head(seq(0, 2, len = 180), -1L)
for(i in 1:length(t_)){
  image <- jacobi:::Image_E6(x, 0.7-0.3i, t_[i])
  svg("zzz.svg")
  opar <- par(mar = c(0,0,0,0), bg = bkgcol)
  plot(c(-100, 100), c(-100, 100), type = "n", 
       xlab = "", ylab = "", axes = FALSE, asp = 1)
  rasterImage(image, -100, -100, 100, 100)
  par(opar)
  dev.off()
  rsvg::rsvg_png(
    "zzz.svg", sprintf("zzpic%03d.png", i), width = 512, height = 512
  )
}

for(i in 1:length(t_)){
  command <-  sprintf(
    "magick convert -density 300 zzpic%03d.png -bordercolor #15191e -border 10x10 -colorspace RGB -trim -quality 100 -resize 512x512! aapic%03d.png",
    i, i
  )
  system(command)
}

