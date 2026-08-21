#' Safely calculate one posterior diagnostic for each parameter
#' @noRd
bayespim_parameter_diagnostic <- function(mats, parameter_metadata, diagnostic, label) {
  diagnostic_names <- paste(
    parameter_metadata$block,
    parameter_metadata$parameter,
    sep = ":"
  )
  values <- stats::setNames(
    rep(NA_real_, nrow(parameter_metadata)),
    make.unique(diagnostic_names)
  )

  for (parameter_row in seq_len(nrow(parameter_metadata))) {
    draws <- bayespim_parameter_draw_matrix(
      mats = mats,
      parameter_index = parameter_metadata$index[parameter_row]
    )
    value <- try(diagnostic(draws), silent = TRUE)

    if (inherits(value, "try-error") || length(value) != 1L) {
      warning(
        sprintf(
          "%s could not be calculated for parameter %s:%s.",
          label,
          parameter_metadata$block[parameter_row],
          parameter_metadata$parameter[parameter_row]
        ),
        call. = FALSE
      )
    } else {
      values[[parameter_row]] <- as.numeric(value)
    }
  }

  values
}

