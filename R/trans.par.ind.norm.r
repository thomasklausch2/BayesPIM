#' trans.par.ind.norm
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
trans.par.ind.norm = function(p, v, Z1){
  mu = Z1 %*% as.matrix(p)
  sd = exp(v)
  cbind(mu,sd)
}