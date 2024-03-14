#' cor2cov
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
cor2cov <- function(R,S){
  diag(S) %*% R %*% diag(S)
}