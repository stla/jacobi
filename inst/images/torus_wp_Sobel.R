library(jacobi)
library(RcppColors)

x <- y <- seq(-4, 0, length.out = 1024L)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(0.5, 0.5 + 0.5i))
}

Z <- outer(x, y, f)
img <- colorMap5(Z)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-100, 100), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)


svg("x.svg", width = 16, height = 16)
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(c(-100, 100), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -100, -100, 100, 100)
par(opar)
dev.off()
rsvg::rsvg_png(
  "x.svg", "wp4x4_cm5.png", width = 1024, height = 1024
)

# Sobelize
library(imager)
# transform the raster to an 'imager' image
im1 <- as.cimg(
  aperm(
    array(
      col2rgb(img) / 255, 
      dim = c(3, 1024, 1024)
    ), c(2L, 3L, 1L)
  )
)
# Sobel transformation
Sobel <- function(im) {
  M <- rbind(
    c(-1, -2, -1),
    c( 0,  0,  0),
    c( 1,  2,  1)
  )
  imX <- convolve(im, as.cimg(M))
  imY <- convolve(im, as.cimg(t(M)))
  imXY <- enorm(list(imX, imY))
  pmax(pmin(imXY, 1), 0)
}
# apply Sobel transformation and get the colors as hex codes
im2 <- Sobel(im1)
r <- c(squeeze(R(im2)))
g <- c(squeeze(G(im2)))
b <- c(squeeze(B(im2)))
clrs <- rgb(r, g, b)

source("C:/SL/R/imagesProcessing/savePNG.R")
savePNG(im2, "wp4x4_cm5_Sobel.png", 1024L, 1024L)


##################################################################
library(rgl)
library(cgalMeshes)
mesh1 <- torusMesh(
  R = sqrt(2), r = 1, nu = 1024L, nv = 1024L, conformal = FALSE
)
mesh2 <- torusMesh(
  R = sqrt(2), r = 1, nu = 1024L, nv = 1024L, conformal = TRUE
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

snapshot3d("tori_square_wp_Sobel.png", webshot = FALSE)
