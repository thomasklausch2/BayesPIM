#' Update AFT parameters with univariate slice sampling
#' @noRd
step_slice <- function(eta, data, x, dist,
                       log_ll, log_prior_fun, tau_t, sig_prior,
                       beta_prior = "norm", q_prior_sd, fix_sigma = FALSE, control) {

  width <- control$width
  if (!(is.numeric(width) && length(width) == 1L &&
        is.finite(width) && width > 0)) {
    stop("`control$width` must be a positive finite number.", call. = FALSE)
  }

  eta_updated <- slice_passthrough(
    eta = as.numeric(eta),
    w = width,
    log_ll = log_ll,
    fix_sigma = fix_sigma,
    fix_q = FALSE,
    has_q = FALSE,
    log_prior_fun = log_prior_fun,
    tau_t = tau_t,
    sig_prior = sig_prior,
    beta_prior = beta_prior,
    q_prior_sd = q_prior_sd,
    y = data$times[data$rows],
    x = x[data$rows, , drop = FALSE],
    dist = dist
  )

  list(
    eta = eta_updated,
    log_post = log_pst(
      eta = eta_updated,
      log_ll = log_ll,
      log_prior_fun = log_prior_fun,
      tau_t = tau_t,
      sig_prior = sig_prior,
      beta_prior = beta_prior,
      q_prior_sd = q_prior_sd,
      y = data$times[data$rows],
      x = x[data$rows, , drop = FALSE],
      dist = dist
    ),
    accepted = NA
  )
}
