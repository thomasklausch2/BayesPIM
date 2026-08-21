#' Random generation from the log-logistic distribution
#' @noRd
rloglog = function(x,lambda, gamma){ rllogis(x, shape = gamma, scale= lambda ) }
