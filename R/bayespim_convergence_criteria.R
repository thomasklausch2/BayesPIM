#' Describe the BayesPIM convergence criteria
#' @noRd
bayespim_convergence_criteria <- function(max_rhat, min_ess, fixed_parameters) {
  criteria <- sprintf(
    "Convergence criteria: R-hat <= %.3f and ESS >= %.1f",
    max_rhat,
    min_ess
  )

  if (length(fixed_parameters) > 0L) {
    criteria <- sprintf(
      "%s for sampled parameters; fixed parameters excluded: %s",
      criteria,
      paste(fixed_parameters, collapse = ", ")
    )
  }

  criteria
}

