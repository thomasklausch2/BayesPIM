#' Update AFT parameters after collapsing over exact incidence times
#' @noRd
step_slice_collapsed <- function(eta, data, x, dist,
                                 log_ll, log_prior_fun, tau_t, sig_prior,
                                 beta_prior = "norm", q_prior_sd, fix_sigma = FALSE, control) {

  width <- control$width
  if (!(is.numeric(width) && length(width) == 1L &&
        is.finite(width) && width > 0)) {
    stop("`control$width` must be a positive finite number.", call. = FALSE)
  }

  if (!is.list(data) || !all(c("L", "R", "rows") %in% names(data))) {
    stop("Collapsed slice sampling requires interval state with `L`, `R`, and `rows`.", call. = FALSE)
  }
  
  x_interval <- x[data$rows, , drop = FALSE]
  eta_updated <- slice_passthrough(
    eta = as.numeric(eta),
    w = width,
    log_ll = log_ll,
    fix_sigma = fix_sigma,
    fix_q = ifelse(dist == 'gengamma', control$fix_q, FALSE),
    has_q = dist == 'gengamma',
    log_prior_fun = log_prior_fun,
    tau_t = tau_t,
    sig_prior = sig_prior,
    beta_prior = beta_prior,
    q_prior_sd = q_prior_sd,
    L = data$L,
    R = data$R,
    x = x_interval,
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
      dist = dist,
      L = data$L,
      R = data$R,
      x = x_interval
    ),
    accepted = NA
  )
}
