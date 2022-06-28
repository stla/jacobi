#' @title Costa surface
#' @description Computes a mesh of the Costa surface.
#'
#' @param nu,nv numbers of subdivisions 
#'
#' @return A triangle \strong{rgl} mesh (object of class \code{mesh3d}).
#' @export
#' 
#' @importFrom rgl tmesh3d clipMesh3d
#' @importFrom Rvcg vcgClean
#'
#' @examples
#' library(jacobi)
#' library(rgl)
#' \donttest{
#' mesh <- CostaMesh(nu = 250, nv = 250)
#' open3d(windowRect = c(50, 50, 562, 562), zoom = 0.9)
#' bg3d("#15191E")
#' shade3d(mesh, color = "darkred", back = "cull")
#' shade3d(mesh, color = "orange", front = "cull")
#' }
CostaMesh <- function(nu = 50L, nv = 50L){ 
  e1 <- Re(wp(1/2, omega = c(1/2, 1i/2)))
  c <- 4*e1^2
  zf1 <- function(u, v){
    wzeta(u + 1i*v, omega = c(1/2, 1i/2))
  }
  zf2 <- function(u, v){
    z <- -v + 1i*(u - 1/2)
    1i * wzeta(z, omega = c(1/2, 1i/2))
  }
  zf3 <- function(u, v){
    z <- u + 1i*(v - 1/2)
    wzeta(z, omega = c(1/2, 1i/2))
  }
  wf <- function(u, v){
    pi/2/e1 * (zf2(u, v) - zf3(u, v))
  }
  fx <- function(u, z1, w){
    Re(pi*(u+pi/4/e1) - z1 + w)/2
  }
  fy <- function(v, z1, w){
    Re(pi*(v+pi/4/e1) - 1i*(z1 + w))/2
  }
  fz <- function(u, v){
    z <- u + 1i*v
    p <- wp(z, omega = c(1/2, 1i/2))
    sqrt(pi/2) * log(Mod((p-e1)/(p+e1)))/2
  }
  #
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  u_ <- seq(0, 1, length.out = nu)
  v_ <- seq(0, 1, length.out = nv)
  vsarray <- array(NA_real_, dim = c(3L, nv, nu))
  for(j in 1L:nv){
    v <- v_[j]
    z1 <- zf1(u_, v)
    w  <- wf(u_, v)
    z <- rbind(fx(u_, z1, w), fy(v, z1, w), fz(u_, v))
    z[is.nan(z)] <- Inf
    vsarray[, j, ] <- z
  }
  vs <- matrix(vsarray, nrow = 3L, ncol = nu*nv)
  # triangles
  tris1 <- tris2 <- matrix(NA_integer_, nrow = 3L, ncol = (nu-1L)*(nv-1L))
  for(i in 1L:(nu-1L)){
    for(j in 1L:(nv-1L)){
      a <- (i-1L)*nv + j
      b <- a + nv
      tris1[, (i-1L)*(nv-1)+j] <- c(a, a + 1L, b)
      tris2[, (i-1L)*(nv-1)+j] <- c(b, a + 1L, b + 1L)
    }
  }
  # output
  mesh0 <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    homogeneous = FALSE
  )
  mesh <- clipMesh3d(
    mesh0, 
    fn = function(vals) apply(vals, 1L, function(row) sum(row^2)), 
    bound = 25, greater = FALSE
  )
  vcgClean(mesh, sel = 6, tol = 0.001)
}
