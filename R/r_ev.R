#' Draw from the standard extreme-value distribution by inversion
#' @noRd
r_ev = function(n){
  u = runif(n)
  q_ev(u)
}
