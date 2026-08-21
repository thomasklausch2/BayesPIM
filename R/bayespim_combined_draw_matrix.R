#' Combine all BayesPIM chains by row
#' @noRd
bayespim_combined_draw_matrix <- function(chains) {
  matrices <- vector("list", length(chains))
  for (chain in seq_along(chains)) {
    matrices[[chain]] <- as.matrix(chains[[chain]])
  }
  do.call(rbind, matrices)
}

