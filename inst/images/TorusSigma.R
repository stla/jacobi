library(jacobi)
library(rgl)
library(RcppColors)

torusMesh <- function(R, r, nu = 50, nv = 30, rgl = TRUE){
  nu <- as.integer(nu)
  nv <- as.integer(nv)
  nunv <- nu*nv
  vs      <- matrix(NA_real_, nrow = 3L, ncol = nunv)
  normals <- matrix(NA_real_, nrow = nu*nv, ncol = 3L)
  tris1   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  tris2   <- matrix(NA_integer_, nrow = 3L, ncol = nunv)
  u_ <- seq(-pi, pi, length.out = nu + 1L)[-1L]
  cosu_ <- cos(u_)
  sinu_ <- sin(u_)
  v_ <- seq(-pi, pi, length.out = nv + 1L)[-1L]
  cosv_ <- cos(v_)
  sinv_ <- sin(v_)
  Rrcosv_ <- R + r*cosv_
  rsinv_ <- r*sinv_
  jp1_ <- c(2L:nv, 1L)
  j_ <- 1L:nv
  color <- NULL
  for(i in 1L:(nu-1L)){
    i_nv <- i*nv
    rg <- (i_nv - nv + 1L):i_nv
    cosu_i <- cosu_[i]
    sinu_i <- sinu_[i]
    vs[, rg] <- rbind(
      cosu_i * Rrcosv_,
      sinu_i * Rrcosv_,
      rsinv_
    )
    color <- c(
      color,
      wsigma((u_[i]) + 1i*(v_), tau = pi/2+1i*pi/2)
    )
    normals[rg, ] <- cbind(
      cosu_i * cosv_,
      sinu_i * cosv_,
      sinv_
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
    -Rrcosv_,
    0,
    rsinv_
  )
  color <- c(
    color,
    wsigma((pi) + 1i*(v_), tau = pi/2+1i*pi/2)
  )
  normals[rg, ] <- cbind(
    -cosv_,
    0,
    sinv_
  )
  k1 <- i_nv - nv
  l_ <- k1 + jp1_
  k_ <- k1 + j_
  tris1[, k_] <- rbind(j_, l_, k_)
  tris2[, k_] <- rbind(j_, jp1_, l_)
  tmesh <- tmesh3d(
    vertices = vs,
    indices = cbind(tris1, tris2),
    normals = normals,
    homogeneous = FALSE
  )
  tmesh$material <- list(color=colorMap1(color, reverse = c(T,T,T)))
  if(rgl){
    tmesh
  }else{
    out <- Mesh(mesh = tmesh)
    out[["normals"]] <- normals
    out
  }
}


mesh <- torusMesh(1, 0.4, nu = 700, nv = 700)
