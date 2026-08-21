#' Subset stored BayesPIM chains by the raw warm-up cutoff
#' @noRd
bayespim_subset_chain_iterations <- function(
    chains,
    warmup = 0L,
    total_iterations
) {
  if (is.null(chains) || length(chains) < 1L) {
    stop("The BayesPIM object contains no posterior chains.", call. = FALSE)
  }

  n_stored <- nrow(as.matrix(chains[[1L]]))
  warmup <- bayespim_validate_warmup(warmup, total_iterations)
  subset_chains <- vector("list", length(chains))
  for (chain in seq_along(chains)) {
    current_chain <- chains[[chain]]
    chain_matrix <- as.matrix(current_chain)
    if (nrow(chain_matrix) != n_stored) {
      stop("All posterior chains must have the same length.", call. = FALSE)
    }
    mcpar <- attr(current_chain, "mcpar")
    iterations <- seq(mcpar[1L], mcpar[2L], by = mcpar[3L])
    if (length(iterations) != nrow(chain_matrix)) {
      stop("The stored MCMC iteration metadata are inconsistent.", call. = FALSE)
    }
    rows <- which(iterations > warmup)
    if (length(rows) == 0L) {
      stop("`warmup` must leave at least one stored draw in each chain.", call. = FALSE)
    }
    subset_chains[[chain]] <- coda::mcmc(
      chain_matrix[rows, , drop = FALSE],
      start = iterations[rows[1L]],
      end = iterations[rows[length(rows)]],
      thin = mcpar[3L]
    )
  }

  coda::mcmc.list(subset_chains)
}
