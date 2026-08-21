#' Trim and thin an `mcmc.list`
#'
#' Convenience function for trimming burn-in iterations and applying thinning to
#' an object of class `mcmc.list`.
#'
#' @param obj An object of class `mcmc.list`.
#' @param burnin Non-negative integer giving the number of stored rows to
#'   discard. Defaults to `0`.
#' @param end Integer; final stored row to retain. Defaults to the number of rows
#'   in the first chain.
#' @param thinning Integer; thinning interval. Defaults to `1`.
#'
#' @return An object of class `mcmc.list` containing the trimmed and thinned
#'   chains.
#'
#' @details
#' The function selects stored rows
#' `seq(burnin + 1, end, by = thinning)` from each chain and reconstructs the
#' result as an `mcmc.list`, preserving its original iteration numbering and
#' multiplying its existing storage interval by `thinning`.
#'
#' @examples
#' data(mod)
#' trimmed <- trim_mcmc(
#'   mod$par,
#'   burnin = mod$warmup,
#'   thinning = 20
#' )
#' nrow(as.matrix(trimmed[[1]]))
#'
#' @export
trim_mcmc <- function(obj, burnin = 0L, end = NULL, thinning = 1L) {
  if (!inherits(obj, "mcmc.list")) {
    stop("obj must be an mcmc.list.", call. = FALSE)
  }

  burnin <- as.integer(round(burnin))
  thinning <- as.integer(round(thinning))

  if (is.na(burnin) || burnin < 0L) {
    stop("burnin must be a non-negative integer.", call. = FALSE)
  }

  if (is.na(thinning) || thinning < 1L) {
    stop("thinning must be a positive integer.", call. = FALSE)
  }

  n_iter <- nrow(as.matrix(obj[[1L]]))

  if (is.null(end)) {
    end <- n_iter
  } else {
    end <- as.integer(round(end))
  }

  if (is.na(end) || end < 1L || end > n_iter) {
    stop("end must be between 1 and the number of stored draws.", call. = FALSE)
  }

  rows <- seq.int(from = burnin + 1L, to = end, by = thinning)

  if (length(rows) == 0L) {
    stop("No draws remain after applying burnin/end/thinning.", call. = FALSE)
  }

  coda::mcmc.list(lapply(obj, function(x) {
    x_mat <- as.matrix(x)
    mcpar <- attr(x, "mcpar")
    iterations <- seq(mcpar[1L], mcpar[2L], by = mcpar[3L])

    coda::mcmc(
      x_mat[rows, , drop = FALSE],
      start = iterations[rows[1L]],
      end = iterations[rows[length(rows)]],
      thin = mcpar[3L] * thinning
    )
  }))
}
