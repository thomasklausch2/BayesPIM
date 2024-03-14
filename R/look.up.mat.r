#' look.up.mat
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
look.up.mat = function(L, a){
  b <- matrix(nrow=length(a), ncol=2)
  for(i in seq_along(a)) b[i,] <- c(L[[i]][a[i]], L[[i]][a[i]+1])
  b
}
# look.up = function(L, a){
#   b <- numeric(length(a))
#   for(i in seq_along(a)) b[i] <- L[[i]][a[i]]
#   b
# }
