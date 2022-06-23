library(rgl)
library(Rvcg)
library(jacobi)


CostaMesh <- function(umin = 0.1, umax = 0.9, vmin = 0.1, vmax = 0.9, nu, nv){ 
  e1 <- Re(wp(1/2, omega = c(1/2, 1i/2)))
  c <- 4*e1^2
  zz1 <- function(u, v){
    w <- u + 1i*v
    zetaw(w, c(c,0))
  }
  zz2 <- function(u, v){
    w <- u + 1i*v
    zetaw(w-1/2, c(c,0))
  }
  zz3 <- function(u, v){
    w <- u + 1i*v
    zetaw(w-1i/2, c(c,0))
  }
  fx <- function(u, z1, z2, z3){
    Re(pi*(u+pi/4/e1) - z1 + pi/2/e1*(z2 - z3))/2
  }
  fy <- function(v, z1, z2, z3){
    Re(pi*(v+pi/4/e1) - 1i*z1 - 1i*pi/2/e1*(z2 - z3))/2
  }
  fz <- function(u,v){
    w <- u + 1i*v
    p <- wp(w, c(c,0))
    sqrt(pi/2)*log(Mod((p-e1)/(p+e1)))/2
  }
  #
  nu <- as.integer(nu + 1)
  nv <- as.integer(nv + 1)
  u_ <- seq(umin, umax, length.out = nu)
  u_ <- c(umin/2, u_, (umax+1)/2)
  v_ <- seq(vmin, vmax, length.out = nv)
  v_ <- c(vmin/2, v_, (vmax+1)/2)
  vsarray <- array(NA_real_, dim=c(3,nv+2,nu+2))
  for(j in 1L:(nv+2)){
    for(i in 1L:(nu+2)){
      z1 <- zz1(u_[i], v_[j])
      z2 <- zz2(u_[i], v_[j])
      z3 <- zz3(u_[i], v_[j])
      vsarray[, j, i] <- 
        c(fx(u_[i], z1, z2, z3), fy(v_[j], z1, z2, z3), fz(u_[i], v_[j]))
    }
  }
  vs <- matrix(vsarray[, -c(1L, nv+2L), -c(1L, nu+2L)], nrow = 3L, ncol = nu*nv)
  # triangles
  tris1 <- tris2 <- matrix(NA_integer_, nrow = 3L, ncol = (nu-1L)*(nv-1L))
  for(i in 1L:(nu-1L)){
    for(j in 1L:(nv-1L)){
      tris1[, (i-1L)*(nv-1)+j] <- c((i-1L)*nv+j, (i-1L)*nv+j+1L, i*nv+j)
      tris2[, (i-1L)*(nv-1)+j] <- c(i*nv+j, (i-1L)*nv+j+1L, i*nv+j+1L)
    }
  }
  # output
  mesh <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE
  )
  vcgUpdateNormals(mesh)
}


mesh <- CostaMesh(nu=100, nv=100)

open3d(windowRect = c(50, 50, 562, 562), zoom=0.95)
bg3d(rgb(21, 25, 30, maxColorValue = 255))
shade3d(mesh, color = "darkred", back="cull")
shade3d(mesh, color = "darkorange", front="cull")

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

command <- "gifski --fps=10 --frames=zzpic*.png -o Costa.gif"
system(command)

pngfiles <- list.files(pattern = "^zzpic?.*png$")
file.remove(pngfiles)

