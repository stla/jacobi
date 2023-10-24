library(jacobi)
library(RcppColors)

x <- seq(-4*1.529954, 4*1.529954, len = 1024L)
y <- seq(-4*1.324979, 4*1.324979, len = 1024L)

halfPeriods(c(0, 1))
1.529954 / 1.324979 # sqrt(4/3)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, g = c(0, 1), derivative = 1L)
}

Z <- outer(x, y, f)
img <- colorMap6(Im(Z))

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)


library(imager)

im1 <- as.cimg(
  aperm(
    array(
      col2rgb(c(img)) / 255, 
      dim = c(3, 1024, 1024)
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
  R = sqrt(7/4), r = 1, nu = 1024L, nv = 1024L, conformal = FALSE
)
mesh2 <- torusMesh(
  R = sqrt(7/4), r = 1, nu = 1024L, nv = 1024L, conformal = TRUE
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






svg("x.svg", width = 16, height = 8)
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(
  c(-200, 200), c(-100, 100), type = "n", xaxs = "i", yaxs = "i",
  xlab = NA, ylab = NA, axes = FALSE, asp = 1
)
rasterImage(img, -200, -100, 200, 100)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "wp2x1.png", width = 1024, height = 512)

