#' Quantile function of the log-logistic distribution
#' @noRd
qloglog = function(x,lambda, gamma){ qllogis(x, shape = gamma, scale= lambda ) }
