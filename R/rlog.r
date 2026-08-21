#' Draw from the standard logistic distribution by inversion
#' @noRd
rlog = function(n){
  u = runif(n)
  qlog(u)
}
