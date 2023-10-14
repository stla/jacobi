library(jacobi)
library(RcppColors)

x <- y <- seq(-4, 0, length.out = 512L)

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(0.5, 0.5 + 0.5i), derivative = 0L)
}

Z <- outer(x, y, f)
img <- colorMap5(Z)
img <- cbind(img, img)

opar <- par(mar = c(0,0,0,0), bg="black")
plot(c(-200, 200), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -200, -100, 200, 100)
par(opar)

library(imager)

im1 <- as.cimg(
  aperm(
    array(
      col2rgb(c(img)) / 255, 
      dim = c(3, 512, 1024)
    ), c(3L, 2L, 1L)
  )
)


svg("x.svg", width = 16, height = 8)
opar <- par(mar = c(0,0,0,0), bg = "black")
plot(c(-200, 200), c(-100, 100), type = "n", xaxs ="i", yaxs = "i",
     xlab = NA, ylab = NA, axes = FALSE, asp = 1)
rasterImage(img, -200, -100, 200, 100)
par(opar)
dev.off()
rsvg::rsvg_png(
  "x.svg", "wp8x4_cm5_Sobel.png", width = 512, height = 256
)

# Sobelize
sobel <- function(im) {
  M <- rbind(
    c(-1, -2, -1),
    c(0, 0, 0),
    c(1, 2, 1)
  )
  imX <- convolve(im, as.cimg(M))
  imY <- convolve(im, as.cimg(t(M)))
  imXY <- enorm(list(imX, imY))
  pmax(pmin(imXY, 1), 0)
}

im2 <- sobel(im1)
r <- c(squeeze(R(im2)))
g <- c(squeeze(G(im2)))
b <- c(squeeze(B(im2)))
clrs <- rgb(r, g, b)
img <- as.raster(im2)

source("C:/SL/R/imagesProcessing/savePNG.R")
savePNG(im2, "wp8x4_cm5_Sobel.png", 512L, 256L)


##################################################################
library(rgl)
library(cgalMeshes)
mesh1 <- torusMesh(
  R = sqrt(5), r = 1, nu = 1024L, nv = 512L, conformal = FALSE
)
mesh2 <- torusMesh(
  R = sqrt(5), r = 1, nu = 1024L, nv = 512L, conformal = TRUE
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

snapshot3d("tori_nonsquare_wp_Sobel.png", webshot = FALSE)
