#' Update AFT parameters with a random-walk Metropolis-Hastings step
#' @noRd
step_mh <- function(eta, data, x, dist,
                    log_ll, log_prior_fun, tau_t, sig_prior,
                    beta_prior = "norm", q_prior_sd, fix_sigma = FALSE, control) {
  rows <- data$rows
  times <- data$times[rows]
  x <- x[rows, , drop = FALSE]
  prop_var <- control$prop_var
  if (is.null(prop_var)) {
    stop("`control$prop_var` must be supplied for the MH sampler.", call. = FALSE)
  }

  eta_proposed <- MASS::mvrnorm(1, eta, prop_var)
  if (fix_sigma) eta_proposed[length(eta_proposed)] <- log(sig_prior)

  log_post_proposed <- log_pst(
    eta = eta_proposed, log_ll = log_ll, log_prior_fun = log_prior_fun,
    tau_t = tau_t, sig_prior = sig_prior, q_prior_sd = q_prior_sd, beta_prior = beta_prior,
    y = times, dist = dist, x = x
  )
  log_post_current <- log_pst(
    eta = eta, log_ll = log_ll, log_prior_fun = log_prior_fun,
    tau_t = tau_t, sig_prior = sig_prior, q_prior_sd = q_prior_sd, beta_prior = beta_prior,
    y = times, dist = dist, x = x
  )
  log_ratio <- log_post_proposed - log_post_current
  u <- stats::runif(1, 0, 1)
  if (is.nan(log_ratio)) {
    log_ratio <- -Inf
  }

  accepted <- log(u) <= log_ratio
  if (accepted) {
    return(list(
      eta = eta_proposed,
      log_post = log_post_proposed,
      accepted = TRUE
    ))
  } else {
    return(list(
      eta = eta,
      log_post = log_post_current,
      accepted = FALSE
    ))
  }
}
