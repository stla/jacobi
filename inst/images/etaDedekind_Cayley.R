library(jacobi)
library(RcppColors)

# background color
bkgcol <- "#ffffffff"

# map the disk to H and calculate kleinj
f <- function(x, y) {
  z <- complex(real = x, imaginary = y)
  K <- rep(NA_complex_, length(x))
  inDisk <- Mod(z) < 1
  K[inDisk] <- eta(disk2H(z[inDisk]))
  K
}
n <- 1024L
x <- y <- seq(-1, 1, length.out = n)
Grid <- expand.grid(X = x, Y = y)
K <- f(Grid$X, Grid$Y)^24
dim(K) <- c(n, n)
# plot
if(require("RcppColors")) {
  img <- colorMap5(K, bkgcolor = bkgcol)
} else {
  img <- as.raster(1 - abs(Im(K))/Mod(K))
}
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)

# isocuboids fork ####
library(isocuboids) 
library(ggplot2)
mat <- 1 - abs(Im(K))/Mod(K)
colors <- as.matrix(as.raster(mat))
colors <- colors[, ncol(colors):1L]
mat[is.na(mat)] <- 0
cuboid_matrix2(
  70*mat, a1 = 45, a2 = 45, show_axes = FALSE, verbose = TRUE, 
  cuboid_fill = c(colors)
) + theme_void() + 
  theme(
    panel.background = element_rect(fill = "white", color = "white")
  )
ggsave(
  "isocuboids_fork_etaDedekindOnCircle.png_gray",
  width = 5, height = 5, dpi = 600
)

# with colors ####
img <- colorMap5(K, bkgcolor = "white")
colors <- mat <- as.matrix(img)
colors <- colors[, ncol(colors):1L]
mat <- c(mat)
mat <- mat |> col2rgb() |> rgb2hsv()
mat <- matrix(mat[3L, ], nrow = n, ncol = n) 
mat[mat == 1] <- 0
gg <- cuboid_matrix2(
  70*mat, a1 = 45, a2 = 45, show_axes = FALSE, verbose = TRUE, 
  cuboid_fill = c(colors)
) + theme_void() + 
  theme(
    panel.background = 
      element_rect(fill = "white", color = "white")
  )
ggsave(
  "isocuboids_fork_etaDedekindOnCircle_colors.png",
  gg, width = 5, height = 5, dpi = 600
)

# isocuboids ####
library(isocuboids)
mat <- 1 - abs(Im(K))/Mod(K)
mat[is.na(mat)] <- 0
cuboid_matrix(
  70*mat, a1 = 45, a2 = 45, show_axes = FALSE, verbose = TRUE, 
  cuboid_fill = trekcolors::trek_pal("klingon")
)
ggplot2::ggsave(
  "isocuboids_etaDedekindOnCircle_klingon.png",
  width = 5, height = 5, dpi = 600
)

# elevated Delaunay ####
mat <- Arg(K) + pi
mat[is.na(mat)] <- 0
library(tessellation)
pts <- expand.grid(X = 1:n, Y = 1:n)
pts$Z <- c(mat)
pts <- as.matrix(pts)
del <- delaunay(pts, elevation = TRUE)

# plotting
library(rgl)
mesh <- del[["mesh"]]
open3d(windowRect = c(100, 100, 612, 356), zoom = 0.6)
aspect3d(1, 1, 20)
shade3d(mesh, color = "limegreen", polygon_offset = 1)



# svg to png ####
svg("x.svg", width = 15, height = 15)
opar <- par(mar = c(0, 0, 0, 0))
plot(NULL, xlim = c(0, 1), ylim = c(0, 1), asp = 1, 
     axes = FALSE, xaxs = "i", yaxs = "i", xlab = NA, ylab = NA)
rasterImage(img, 0, 0, 1, 1)
par(opar)
dev.off()

rsvg::rsvg_png("x.svg", "zzzetaDedekindOnCircle.png", 512, 512)
