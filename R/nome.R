#' @importFrom Carlson elliptic_F
#' @noRd
tau_from_m <- function(m){
  1i * elliptic_F(pi/2, 1-m) / elliptic_F(pi/2, m) 
}