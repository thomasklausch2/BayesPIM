#' Format draw counts for a BayesPIM summary
#' @noRd
format_bayespim_summary_draw_counts <- function(draw_counts) {
  paste(
    sprintf(
      "MCMC iterations generated: %d (%d per chain)",
      draw_counts$total_generated,
      draw_counts$generated_per_chain
    ),
    sprintf(
      "Parameter draws stored: %d (%d per chain; save_every = %d)",
      draw_counts$total_stored,
      draw_counts$stored_per_chain,
      draw_counts$save_every
    ),
    sprintf(
      "Warm-up cutoff: %d generated iterations per chain",
      draw_counts$warmup_per_chain
    ),
    sprintf(
      "Stored warm-up draws omitted: %d (%d per chain)",
      draw_counts$stored_warmup,
      draw_counts$stored_warmup_per_chain
    ),
    sprintf(
      "Posterior draws used: %d (%d per chain)",
      draw_counts$retained_stored,
      draw_counts$retained_stored_per_chain
    ),
    sep = "\n"
  )
}
