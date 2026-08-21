#' Calculate draw counts for a BayesPIM summary
#' @noRd
bayespim_summary_draw_counts <- function(
    object,
    retained_chains,
    warmup
) {
  all_chains <- object$par
  n_chains <- length(all_chains)
  n_stored <- nrow(as.matrix(all_chains[[1L]]))
  total_iterations <- object$total_iterations
  save_every <- object$save_every
  warmup <- bayespim_validate_warmup(warmup, total_iterations)

  total_generated <- total_iterations * n_chains
  total_stored <- n_stored * n_chains
  retained_stored_per_chain <- nrow(as.matrix(retained_chains[[1L]]))
  retained_stored <- retained_stored_per_chain * n_chains
  stored_iterations <- seq(attr(all_chains[[1L]], "mcpar")[1L],
                           attr(all_chains[[1L]], "mcpar")[2L],
                           by = attr(all_chains[[1L]], "mcpar")[3L])
  stored_warmup_per_chain <- sum(stored_iterations <= warmup)
  stored_warmup <- stored_warmup_per_chain * n_chains

  list(
    total_generated = total_generated,
    total_stored = total_stored,
    generated_per_chain = total_iterations,
    stored_per_chain = n_stored,
    save_every = save_every,
    warmup_per_chain = warmup,
    warmup_generated = warmup * n_chains,
    stored_warmup = stored_warmup,
    stored_warmup_per_chain = stored_warmup_per_chain,
    retained_stored = retained_stored,
    retained_stored_per_chain = retained_stored_per_chain
  )
}
