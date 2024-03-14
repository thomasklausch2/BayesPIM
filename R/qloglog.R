#' qloglog
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
qloglog = function(x,lambda, gamma){ qllogis(x, shape = gamma, scale= lambda ) }