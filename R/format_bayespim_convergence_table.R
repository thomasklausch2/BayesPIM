#' Format BayesPIM convergence diagnostics
#' @noRd
format_bayespim_convergence_table <- function(parameter_metadata, rhat, ess) {
  data.frame(
    block = parameter_metadata$block,
    parameter = parameter_metadata$parameter,
    R_hat = unname(round(rhat, 3L)),
    ESS = unname(round(ess, 1L)),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

