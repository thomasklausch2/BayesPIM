#' Log posterior of the AFT parameters up to a normalizing constant
#' @noRd
log_pst <- function(
    eta,
    log_ll,
    log_prior_fun,
    dist,
    beta_prior,
    tau_t,
    sig_prior,
    q_prior_sd,
    ...
) {
  ll <- log_ll(
    eta = eta,
    dist = dist,
    ...
  )
  
  log_prior <- log_prior_fun(
    eta = eta,
    dist = dist,
    beta_prior = beta_prior,
    tau_t = tau_t,
    sig_prior = sig_prior,
    q_prior_sd = q_prior_sd
  )
  
  log_posterior <- ll + log_prior
  
  if (!is.finite(log_posterior)) {
    return(-Inf)
  }
  
  log_posterior
}
