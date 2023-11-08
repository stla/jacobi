z <- 1i

circle2square <- function(z) {
  C <- 2*gamma(5/4)^2 / sqrt(pi)
  sqrt2half <- sqrt(2) / 2
  w <- complex(real = sqrt2half, imaginary = sqrt2half)
  -w * Carlson::elliptic_F(1i * asinh(w * z), -1) / C  
}

n <- 70L
r <- seq(0, 1, length.out = n)
theta <- seq(0, 2*pi, length.out = n+1L)[-1L]
Grid <- transform(
  expand.grid(R = r, Theta = theta),
  Z = R*exp(1i*Theta)
)
s <- vapply(Grid$Z, circle2square, complex(1L))
plot(Re(s), Im(s), pch = ".", asp = 1, cex = 2)
#
# a more insightful plot ####
r_ <- seq(0, 1, length.out = 10L)
theta_ <- seq(0, 2*pi, length.out = 33)[-1L]
plot(
  NULL, xlim = c(-1, 1), ylim = c(-1, 1), asp = 1, xlab = "x", ylab = "y"
)
for(r in r_) {
  theta <- sort(
    c(seq(0, 2, length.out = 200L), c(1/4, 3/4, 5/4, 7/4))
  )
  z <- r*(cospi(theta) + 1i*sinpi(theta))
  s <- vapply(z, circle2square, complex(1L))
  lines(Re(s), Im(s), col = "blue", lwd = 2)
}
for(theta in theta_) {
  r <- seq(0, 1, length.out = 30L)
  z <- r*exp(1i*theta)
  s <- vapply(z, circle2square, complex(1L))
  lines(Re(s), Im(s), col = "green", lwd = 2)
}
