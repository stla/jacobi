library(jacobi)
library(RcppColors)

npixls <- 512L
x <- seq(-4*2, 4*2, len = npixls)
y <- seq(-4*1, 4*1, len = npixls)

halfPeriods(c(0, 1))
1.529954 / 1.324979 # sqrt(4/3)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(1, 0.5 + 2i), derivative = 1L)
}

Z <- outer(x, y, f)
img <- colorMap6(Im(Z))


s <- 1/2
opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100)*s, c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100*s, -100, 100*s, 100)
par(opar)


library(imager)

im1 <- as.cimg(
  aperm(
    array(
      col2rgb(c(img)) / 255, 
      dim = c(3, npixls, npixls)
    ), c(3L, 2L, 1L)
  )
)

# s = sqrt(R²/r²-1) ; r=1
# s²+1 = R²
# R = sqrt(4/3 + 1) = sqrt(7/3) or sqrt(3/4 + 1) = sqrt(7/4) 

im2 <- im1
r <- c(squeeze(R(im2)))
g <- c(squeeze(G(im2)))
b <- c(squeeze(B(im2)))
clrs <- rgb(r, g, b)

##################################################################
library(rgl)
library(cgalMeshes)
mesh1 <- torusMesh(
  R = sqrt(5), r = 1, nu = npixls, nv = npixls, conformal = FALSE
)
mesh2 <- torusMesh(
  R = sqrt(5), r = 1, nu = npixls, nv = npixls, conformal = TRUE
)
#clrs <- c(img)
mesh1[["material"]] <- mesh2[["material"]] <- list("color" = clrs)
# plot
open3d(windowRect = 50 + c(0, 0, 640, 320))
clear3d(type = "lights") # destroy current lights
light3d(x = -50, y = 100, z = 100)
bg3d("#363940")
mfrow3d(1, 2)
view3d(-25, -25, zoom = 0.75)
shade3d(mesh1, specular = "gold")
next3d()
view3d(-25, -25, zoom = 0.75)
shade3d(mesh2, specular = "gold")





s <- sqrt(3/4)
svg("x.svg", width = 16*s, height = 16)
opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100)*s, c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100*s, -100, 100*s, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "wpprime_equianharmonic.png", 
               width = 512*s, height = 512)

