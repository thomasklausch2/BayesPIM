#' Plot method for bayespim objects
#'
#' @param x An object of class `"bayespim"`.
#' @param warmup Number of initial generated iterations to discard from each
#'   chain before optional plot-only thinning. Defaults to the final warmup value stored by
#'   `bayespim()`, including any update or explicit override. Supplying this
#'   argument overrides the stored value for this plot
#'   only.
#' @param thinning Positive integer thinning interval used only for this plot.
#'   The default is no additional thinning. When the argument is omitted and
#'   more than 20,000 post-warm-up draws are stored per chain, it is increased
#'   automatically to limit plotting cost.
#' @param ... Additional arguments passed to the MCMC plotting method.
#'
#' @return Invisibly returns `x`. The method is called for its side effect of
#'   producing trace and density plots for the requested parameter blocks.
#'
#' @examples
#' data(mod)
#' plot(mod, thinning = 20)
#'
#' @method plot bayespim
#' @export
plot.bayespim <- function(
    x,
    warmup = x$warmup,
    thinning = 1L,
    ...
) {
  if (!inherits(x, "bayespim")) {
    stop("`x` must be of class 'bayespim'.", call. = FALSE)
  }
  chains <- bayespim_analysis_chains(x, warmup = warmup)
  n_stored <- nrow(as.matrix(chains[[1L]]))
  if (missing(thinning) && n_stored > 20000L) {
    thinning <- ceiling(n_stored / 20000L)
  }
  chains <- trim_mcmc(chains, burnin = 0L, thinning = thinning)

  graphics::plot(chains, ...)
  invisible(x)
}
