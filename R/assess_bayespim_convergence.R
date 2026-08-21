#' Calculate BayesPIM convergence diagnostics
#' @noRd
assess_bayespim_convergence <- function(ret, max_rhat, min_ess, chains = NULL) {
  mats <- bayespim_convergence_chain_matrices(ret, chains = chains)
  parameter_metadata <- bayespim_parameter_metadata(ret, mats)
  sampled_metadata <- parameter_metadata[!parameter_metadata$fixed, , drop = FALSE]

  rhat <- bayespim_parameter_diagnostic(
    mats = mats,
    parameter_metadata = sampled_metadata,
    diagnostic = posterior::rhat,
    label = "R-hat"
  )
  ess <- bayespim_parameter_diagnostic(
    mats = mats,
    parameter_metadata = sampled_metadata,
    diagnostic = posterior::ess_mean,
    label = "ESS"
  )

  rhat_ok <- length(rhat) > 0L && all(is.finite(rhat)) && all(rhat <= max_rhat)
  ess_ok <- length(ess) > 0L && all(is.finite(ess)) && all(ess >= min_ess)

  list(
    rhat = rhat,
    ess = ess,
    table = format_bayespim_convergence_table(
      parameter_metadata = sampled_metadata,
      rhat = rhat,
      ess = ess
    ),
    max_rhat = max_rhat,
    min_ess = min_ess,
    criteria = bayespim_convergence_criteria(
      max_rhat = max_rhat,
      min_ess = min_ess,
      fixed_parameters = paste(
        parameter_metadata$block[parameter_metadata$fixed],
        parameter_metadata$parameter[parameter_metadata$fixed],
        sep = ":"
      )
    ),
    n_iter = ret$total_iterations,
    n_iter_diagnostic = nrow(mats[[1L]]),
    n_chains = length(mats),
    fixed_parameters = paste(
      parameter_metadata$block[parameter_metadata$fixed],
      parameter_metadata$parameter[parameter_metadata$fixed],
      sep = ":"
    ),
    converged = rhat_ok && ess_ok
  )
}
