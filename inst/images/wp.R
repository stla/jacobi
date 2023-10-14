library(jacobi)
library(RcppColors)

x <- seq(4, 8, len = 1024)-8
y <- seq(4, 8, len = 1024)-8

f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  wp(z, omega = c(0.5, 0.5+0.5i), derivative = 0L)
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

##################################################################
library(rgl)
library(cgalMeshes)
mesh <- torusMesh(
  R = sqrt(2), r = 1, nu = 1024, nv = 1024, conformal = FALSE, normals = FALSE
)

#clrs <- c(img)
mesh[["material"]] <- list("color" = clrs)

open3d(windowRect = 50 + c(0, 0, 512, 512))
bg3d("#363940")
view3d(-25, -25, zoom = 0.75)
clear3d(type = "lights") # destroy current lights
light3d(x = -50, y = 100, z = 100)
shade3d(mesh, specular = "forestgreen")
