#' rloglog
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
rloglog = function(x,lambda, gamma){ rllogis(x, shape = gamma, scale= lambda ) }