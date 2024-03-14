#' rlog
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
rlog = function(n){
  u = runif(n)
  qlog(u)
}