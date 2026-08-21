#' Describe the parameter blocks in a BayesPIM fit
#' @noRd
bayespim_parameter_metadata <- function(ret, mats = NULL) {
  if (is.null(mats)) {
    mats <- bayespim_convergence_chain_matrices(ret)
  }

  parameter_names <- colnames(mats[[1L]])
  n_parameters <- length(parameter_names)
  n_t <- bayespim_design_ncol(ret$x_t) + 2L
  if (identical(ret$dist, "gengamma")) {
    n_t <- n_t + 1L
  }

  n_g <- 0L
  if (isTRUE(ret$prev)) {
    n_g <- bayespim_design_ncol(ret$x_g) + 1L
  }

  n_kappa <- as.integer(isTRUE(ret$update_kappa))
  if (n_t + n_g + n_kappa != n_parameters) {
    stop(
      paste0(
        "Could not match the stored parameters to the t, g, and kappa blocks. ",
        "The model metadata imply ", n_t + n_g + n_kappa,
        " parameters, but the chains contain ", n_parameters, "."
      ),
      call. = FALSE
    )
  }

  block <- c(
    rep("t", n_t),
    rep("g", n_g),
    rep("kappa", n_kappa)
  )
  fixed <- rep(FALSE, n_parameters)

  if (isTRUE(ret$fix_sigma)) {
    fixed[bayespim_design_ncol(ret$x_t) + 2L] <- TRUE
  }
  if (isTRUE(ret$fix_q) && identical(ret$dist, "gengamma")) {
    fixed[bayespim_design_ncol(ret$x_t) + 3L] <- TRUE
  }

  metadata <- data.frame(
    index = seq_len(n_parameters),
    block = block,
    parameter = parameter_names,
    fixed = fixed,
    stringsAsFactors = FALSE
  )

  if (!any(!metadata$fixed)) {
    stop("No sampled parameters are available for convergence assessment.", call. = FALSE)
  }

  metadata
}

