#' Map generalized gamma AFT coefficients to distribution parameters
#' @noRd
trans_par_gengamma = function(x1, par){
  p = length(par)
  mu    <- drop(x1 %*% par[seq_len(p - 2L)])
  sigma <- exp(par[p - 1L])
  Q     <- par[p]
  cbind(mu = mu, sigma = sigma, Q = Q)
}
