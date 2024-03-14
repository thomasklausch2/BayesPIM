#' pst.theta
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
pst.theta = function(C, a){
  x = table(C)
  alpha = a + x
  rdirichlet(1, alpha)
}
