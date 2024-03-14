#' logrob
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
logrob = function(x,tol){
  i = x==0
  x[i] = x[i]+tol
  log(x)
}