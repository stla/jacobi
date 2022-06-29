# https://math.stackexchange.com/questions/107966/how-do-i-map-the-torus-to-a-plane

library(jacobi)
library(rgl)
library(RcppColors)

s <- 1
fx <- function(u,v){
  w <- sqrt(s*s + 1) - cos(2*pi*v)
  s * cos(2*pi*u/s)/w
}
fy <- function(u,v){
  w <- sqrt(s*s + 1) - cos(2*pi*v)
  s * sin(2*pi*u/s)/w
}
fz <- function(v){
  w <- sqrt(s*s + 1) - cos(2*pi*v)
  s * sin(2*pi*v)/w
}


torusMesh <- function(nu = 50, nv = 30, rgl = TRUE){
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  tris1   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(-1/2, 1/2, length.out = nu + 1L)[-1L]
  v_ <- seq(-1/2, 1/2, length.out = nv + 1L)[-1L]
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  color <- NULL
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    rg <- (i_nv - nv + 1L):i_nv
    vs[, rg] <- rbind(
      fx(u_[i], v_), fy(u_[i], v_), fz(v_)
    )
    # color <- c(
    #   color,
    #   wsigma((u_[i]) + 1i*(v_), tau = 2+2i)
    # )
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
    fx(1/2, v_), fy(1, v_), fz(v_)
  )
  # color <- c(
  #   color,
  #   wsigma(1 + 1i*(v_), tau = 2+2i)
  # )
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
  # tmesh$material <- list(color=colorMap1(color, reverse = c(T,T,T)))
  if(rgl){
    tmesh
  }else{
    out <- Mesh(mesh = tmesh)
    out[["normals"]] <- normals
    out
  }
}


mesh <- Rvcg::vcgUpdateNormals(torusMesh(nu = 600, nv = 600))

coords <- function(xyz){
  x <- xyz[, 1L]
  y <- xyz[, 2L]
  z <- xyz[, 3L]
  u <- s/(2*pi) * atan(y/x)
  sgn <- ifelse(x*x + y*y >= 2, 1, -1)
  v <- sign(z)/(2*pi)*acos((z*z*sqrt(s*s+1) + sgn*sqrt(1-z*z*s*s))/(z*z+1))
  cbind(u,v)
}
sphcoords <- coords(t(mesh$vb[-4L, ]))
theta <- sphcoords[, 1L]*2
phi   <- sphcoords[, 2L]
Z <- wsigma(theta + 1i * phi, omega=c(1/4, 1/4*(1/2+1i/2)))
color <- colorMap1(Z, reverse = c(FALSE, FALSE, TRUE))
mesh$material <- list(color = color)

open3d(windowRect = c(50, 50, 562, 562), zoom = 0.75)
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
