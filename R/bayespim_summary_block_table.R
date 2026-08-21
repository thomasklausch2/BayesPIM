#' Create one parameter-block table for a BayesPIM summary
#' @noRd
bayespim_summary_block_table <- function(
    chains,
    parameter_metadata,
    convergence_table,
    block,
    probs
) {
  block_metadata <- parameter_metadata[
    parameter_metadata$block == block,
    ,
    drop = FALSE
  ]
  if (nrow(block_metadata) == 0L) {
    return(NULL)
  }

  draws <- bayespim_combined_draw_matrix(chains)
  quantile_names <- names(stats::quantile(c(0, 1), probs = probs))
  quantiles <- matrix(
    NA_real_,
    nrow = nrow(block_metadata),
    ncol = length(probs),
    dimnames = list(NULL, quantile_names)
  )

  for (parameter_row in seq_len(nrow(block_metadata))) {
    quantiles[parameter_row, ] <- stats::quantile(
      draws[, block_metadata$index[parameter_row]],
      probs = probs,
      names = FALSE
    )
  }

  block_diagnostics <- convergence_table[
    convergence_table$block == block,
    ,
    drop = FALSE
  ]
  diagnostic_rows <- match(
    block_metadata$parameter,
    block_diagnostics$parameter
  )

  out <- as.data.frame(quantiles, check.names = FALSE)
  out$R_hat <- block_diagnostics$R_hat[diagnostic_rows]
  out$ESS <- block_diagnostics$ESS[diagnostic_rows]
  rownames(out) <- make.unique(block_metadata$parameter)
  out
}

