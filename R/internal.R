isBoolean <- function(x){
  is.logical(x) && length(x) == 1L && !is.na(x)
}

isComplexNumber <- function(z){
  (is.complex(z) || is.numeric(z)) && length(z) == 1L && !is.na(z)
}

isComplexPair <- function(x){
  (is.complex(x) || is.numeric(x)) && length(x) == 2L && !anyNA(x)
}

isComplex <- function(x){
  (is.complex(x) || is.numeric(x)) && !anyNA(x)
}

isComplexObject <- function(x){
  is.complex(x) || is.numeric(x)
}

check_and_get_tau <- function(tau, q){
  if(is.null(tau) && is.null(q)){
    stop("You must supply either `tau` or `q`.")
  }
  if(!is.null(tau) && !is.null(q)){
    stop("You must supply either `tau` or `q`, not both.")
  }
  if(!is.null(tau)){
    stopifnot(isComplexNumber(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be strictly positive.")
    }
  }
  if(!is.null(q)){
    stopifnot(isComplexNumber(q))
    if(Mod(q) >= 1){
      stop("The modulus of `q` must be strictly less than one.")
    }
    if(Im(q) == 0 && Re(q) <= 0){
      stop("If `q` is real, it must be strictly positive.")
    }
    tau <- -1i * log(q) / pi
    # if(Im(tau) <= 0){
    #   stop("Invalid value of `q`.")
    # }
  }
  tau
}

isReal <- function(g){
  all(Im(g) == 0)
}

#' @importFrom Carlson elliptic_F
#' @noRd
tau_from_m <- function(m){
  1i * elliptic_F(pi/2, 1-m) / elliptic_F(pi/2, m) 
}

check_and_get_tau_from_m <- function(tau, m){
  if(is.null(tau) && is.null(m)){
    stop("You must supply either `tau` or `m`.")
  }
  if(!is.null(tau) && !is.null(m)){
    stop("You must supply either `tau` or `m`, not both.")
  }
  if(!is.null(tau)){
    stopifnot(isComplexNumber(tau))
    if(Im(tau) <= 0){
      stop("The imaginary part of `tau` must be strictly positive.")
    }
  }
  if(!is.null(m)){
    stopifnot(isComplexNumber(m))
    tau <- tau_from_m(m)
    if(Im(tau) <= 0){
      stop("Invalid value of `m`.")
    }
  }
  tau
}

# e3e2e1 <- function(g){
#   g2 <- g[1L]
#   g3 <- g[2L]
#   a <- 27*g3 + 3*sqrt(as.complex(-3*g2^3 + 81*g3^2))
#   b <- a^(1/3)
#   c <- g2/b
#   bp3c <- b + 3*c
#   bm3c <- sqrt(3)*(b - 3*c)
#   e1 <- bp3c/6
#   e2 <- -(bp3c + 1i*bm3c)/12
#   e3 <- -(bp3c - 1i*bm3c)/12
#   c(e3, e2, e1)
# }
