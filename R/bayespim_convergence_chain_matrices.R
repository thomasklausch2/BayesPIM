#' Convert the BayesPIM chains to equally sized matrices
#' @noRd
bayespim_convergence_chain_matrices <- function(ret, chains = NULL) {
  if (is.null(chains)) {
    chains <- bayespim_analysis_chains(ret)
  }
  if (is.null(chains) || length(chains) < 1L) {
    stop(
      "Convergence assessment requires at least one parameter chain.",
      call. = FALSE
    )
  }

  mats <- vector("list", length(chains))
  for (chain in seq_along(chains)) {
    mats[[chain]] <- as.matrix(chains[[chain]])
  }

  n_iter <- vapply(mats, nrow, integer(1L))
  if (length(unique(n_iter)) != 1L) {
    stop("All chains must contain the same number of draws.", call. = FALSE)
  }

  parameter_names <- colnames(mats[[1L]])
  if (is.null(parameter_names) || length(parameter_names) == 0L) {
    stop("The MCMC chains must have parameter names.", call. = FALSE)
  }

  for (chain in seq_along(mats)) {
    if (!identical(colnames(mats[[chain]]), parameter_names)) {
      stop(
        "All chains must contain the same parameters in the same order.",
        call. = FALSE
      )
    }
  }

  mats
}

