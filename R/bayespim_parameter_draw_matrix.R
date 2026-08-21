#' Arrange one parameter as iterations by chains
#' @noRd
bayespim_parameter_draw_matrix <- function(mats, parameter_index) {
  n_iter <- nrow(mats[[1L]])
  draws <- matrix(
    NA_real_,
    nrow = n_iter,
    ncol = length(mats),
    dimnames = list(NULL, paste0("chain", seq_along(mats)))
  )

  for (chain in seq_along(mats)) {
    draws[, chain] <- as.numeric(mats[[chain]][, parameter_index])
  }

  draws
}

