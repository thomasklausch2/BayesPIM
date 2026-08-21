#' Map AFT coefficients to scale and shape parameters on the natural scale
#' @noRd
trans_par = function(x1, par){
  p = length(par)
  p1 = exp(x1 %*% as.matrix(par[1:(p-1)]))
  p2 = 1/exp(par[(p)])
  cbind(p1,p2)
}
 

