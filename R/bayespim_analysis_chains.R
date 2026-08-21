#' Return every stored BayesPIM parameter draw after warm-up
#' @noRd
bayespim_analysis_chains <- function(object, warmup = NULL) {
  if (is.null(warmup)) warmup <- object$warmup

  bayespim_subset_chain_iterations(
    chains = object$par,
    warmup = warmup,
    total_iterations = object$total_iterations
  )
}
