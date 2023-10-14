# https://math.stackexchange.com/questions/107966/how-do-i-map-the-torus-to-a-plane
# https://static1.bridgesmathart.org/2011/cdrom/proceedings/134/paper_134.pdf

library(jacobi)
library(rgl)
library(Rvcg)
library(RcppColors)

s <- 1
fw <- function(v){
  sqrt(s*s + 1) - cos(2*pi*v)
}
fx <- function(u, w){
  s * cos(2*pi*u/s)/w
}
fy <- function(u, w){
  s * sin(2*pi*u/s)/w
}
fz <- function(v, w){
  sin(2*pi*v)/w
}

torusMesh <- function(nu = 50, nv = 30, rgl = TRUE){
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  tris1   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(-s/2, s/2, length.out = nu + 1L)[-1L]
  v_ <- seq(-1/2, 1/2, length.out = nv + 1L)[-1L]
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  w_ <- fw(v_)
  z_ <- fz(v_, w_)
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    rg <- (i_nv - nv + 1L):i_nv
    vs[, rg] <- rbind(
      fx(u_[i], w_), fy(u_[i], w_), z_
    )
    k1 <- i_nv - nv
    k_ <- k1 + j_
    l_ <- k1 + jp1_
    m_ <- i_nv + j_
    tris1[, k_] <- rbind(m_, l_, k_)
    tris2[, k_] <- rbind(m_, i_nv + jp1_, l_)
  }
  i_nv <- nunv
  rg <- (i_nv - nv + 1L):i_nv
  vs[, rg] <- rbind(
    fx(s/2, w_), fy(s/2, w_), z_
  )
  k1 <- i_nv - nv
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  tmesh <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE
  )
  vcgUpdateNormals(tmesh)
}


mesh <- Rvcg::vcgUpdateNormals(torusMesh(nu = 256, nv = 256))

coords <- function(xyz){
  x <- xyz[, 1L]
  y <- xyz[, 2L]
  z <- xyz[, 3L]
  u <- s/(2*pi) * atan(y/x)
  sgn <- ifelse(x*x + y*y >= s*s+1, 1, -1)
  v <- sign(z)/(2*pi)*acos((z*z*sqrt(s*s+1) + sgn*sqrt(1-z*z*s*s))/(z*z+1))
  cbind(u, v)
}
uv <- coords(t(mesh$vb[-4L, ]))
theta <- 2*uv[, 1L]
phi   <- uv[, 2L]
Z <- wsigma(theta + 1i * phi, omega=c(1/4, 1/4*(1/2+1i/2)))
color <- colorMap1(Z, reverse = c(FALSE, FALSE, TRUE))
# color <- ifelse(theta <= 0, ifelse(phi <= 0, "navy", "yellow"), ifelse(phi > 0, "navy", "yellow"))
mesh$material <- list(color = color)

open3d(windowRect = c(50, 50, 562, 562))
view3d(0, 0, zoom = 0.95)
bg3d("gainsboro")
shade3d(mesh)

# # -- if you want an animation
M <- par3d("userMatrix")
movie3d(
  par3dinterp(
    time = seq(0, 1, len = 9),
    userMatrix = list(
      M,
      rotate3d(M, pi, 1, 0, 0),
      rotate3d(M, pi, 1, 1, 0),
      rotate3d(M, pi, 1, 1, 1),
      rotate3d(M, pi, 0, 1, 1),
      rotate3d(M, pi, 0, 1, 0),
      rotate3d(M, pi, 1, 0, 1),
      rotate3d(M, pi, 0, 0, 1),
      M
    )
  ),
  fps = 120,
  duration = 1,
  dir = ".",
  movie = "zzpic",
  convert = FALSE,
  webshot = FALSE
)

command <- "gifski --fps=9 --frames=zzpic*.png -o SigmaTorus.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)
