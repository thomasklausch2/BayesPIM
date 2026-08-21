#' Summary method for bayespim objects
#'
#' @param object An object of class `"bayespim"`.
#' @param warmup Number of initial generated MCMC iterations to discard from
#'   each chain before calculating the summary from every stored draw. It is
#'   interpreted on the generated-iteration scale, not as a number of stored
#'   draws when `save_every > 1`. Defaults to the final warmup value stored by
#'   `bayespim()`, including any update or explicit override. Supplying this
#'   argument overrides the stored value for this summary only.
#' @param probs Numeric vector of posterior quantiles.
#' @param ... Additional arguments, currently unused.
#'
#' @return Invisibly returns a list of class `"summary.bayespim"` containing
#'   the fitted latent-time distribution and incidence sampler, posterior and
#'   convergence tables for the latent-time model and prevalence model, and
#'   kappa when estimated, together with the convergence criteria, formatted
#'   draw information, and the underlying draw counts. The same information is
#'   printed to the console.
#'
#' @examples
#' data(mod)
#' summary(mod)
#'
#' @method summary bayespim
#' @export
summary.bayespim <- function(
    object,
    warmup = object$warmup,
    probs = c(0.025, 0.5, 0.975),
    ...
) {
  if (!inherits(object, "bayespim")) {
    stop("`object` must be of class 'bayespim'.", call. = FALSE)
  }
  probs <- bayespim_validate_probs(probs)
  if (is.null(warmup)) warmup <- 0L
  chains <- bayespim_analysis_chains(
    object = object,
    warmup = warmup
  )

  max_rhat <- object$convergence$max_rhat
  min_ess <- object$convergence$min_ess
  convergence <- assess_bayespim_convergence(
    ret = object,
    max_rhat = max_rhat,
    min_ess = min_ess,
    chains = chains
  )
  parameter_metadata <- bayespim_parameter_metadata(
    object,
    lapply(chains, as.matrix)
  )

  table_t <- bayespim_summary_block_table(
    chains = chains,
    parameter_metadata = parameter_metadata,
    convergence_table = convergence$table,
    block = "t",
    probs = probs
  )
  table_g <- bayespim_summary_block_table(
    chains = chains,
    parameter_metadata = parameter_metadata,
    convergence_table = convergence$table,
    block = "g",
    probs = probs
  )
  table_kappa <- bayespim_summary_block_table(
    chains = chains,
    parameter_metadata = parameter_metadata,
    convergence_table = convergence$table,
    block = "kappa",
    probs = probs
  )

  draw_counts <- bayespim_summary_draw_counts(
    object = object,
    retained_chains = chains,
    warmup = warmup
  )
  draws_info <- format_bayespim_summary_draw_counts(
    draw_counts = draw_counts
  )

  printed <- c(
    paste0("Latent-time distribution: ", object$dist),
    paste0("Incidence sampler: ", object$sampler),
    "",
    "Parameters of the latent-time model",
    capture.output(print(round(table_t, 3)))
  )
  if (!is.null(table_g)) {
    printed <- c(
      printed,
      "",
      "Parameters of the prevalence model",
      capture.output(print(round(table_g, 3)))
    )
  }
  if (!is.null(table_kappa)) {
    printed <- c(
      printed,
      "",
      "Test sensitivity",
      capture.output(print(round(table_kappa, 3)))
    )
  }
  printed <- c(
    printed,
    "",
    convergence$criteria,
    "",
    draws_info
  )
  cat(paste(printed, collapse = "\n"), "\n", sep = "")

  out <- structure(
    list(
      distribution = object$dist,
      sampler = object$sampler,
      table_t = table_t,
      table_g = table_g,
      table_kappa = table_kappa,
      convergence_criteria = convergence$criteria,
      draws_info = draws_info,
      draws = draw_counts
    ),
    class = "summary.bayespim"
  )
  invisible(out)
}
