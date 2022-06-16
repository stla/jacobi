isComplexPair <- function(x){
  (is.complex(x) || is.numeric(x)) && length(x) == 2L && !anyNA(x)
}

g_from_omega <- function(w1, w2){
  if(Im(w2)*Re(w1) < Im(w1)*Re(w2)){
    stop(
      "Invalid `omega` values. Do you want to exchange `omega1` and `omega2`?"
    )
  }
  ratio <- w2 / w1
  if(Im(ratio) <= 0){
    stop("The ratio `omega2/omega1` must have a positive imaginary part.")
  }
  q <- exp(1i * pi * ratio)
  j2 <- jtheta2(0, q)
  j3 <- jtheta3(0, q)
  g2 <- 4/3 * (pi/2/w1)**4 * (j2**8 - (j2*j3)**4 + j3**8) 
  g3 <- 8/27 * (pi/2/w1)**6 * (j2**12 - (
    (3/2 * j2**8 * j3**4) + (3/2 * j2**4 * j3**8) 
  ) + j3**12)
  c(g2, g3)
}

pweierstrass <- function(z, g = NULL, omega = NULL, derivative = 0L){
  if(!is.element(derivative, 0L:3L)){
    stop("`derivative` must be an integer between 0 and 3.") 
  }
  if(is.null(g) && is.null(omega)){
    stop("You must supply either `g` or `omega`.")
  }
  if(!is.null(tau) && !is.null(q)){
    stop("You must supply either `g` or `omega`, not both.")
  }
  if(!is.null(g)){
    stopifnot(isComplexPair(g))
  }
  if(!is.null(omega)){
    stopifnot(isComplexPair(omega))
    g <- g_from_omega(omega[1L], omega[2L])
  }
  g2 <- g[1L]
  g3 <- g[2L]
}
