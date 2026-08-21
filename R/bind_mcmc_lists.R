#' Concatenate the draws of two mcmc.lists chain by chain
#' @noRd
bind_mcmc_lists <- function(pr, r) {
  if (length(pr) != length(r)) {
    stop("MCMC lists must contain the same number of chains.", call. = FALSE)
  }

  bound <- lapply(seq_along(pr), function(i) {
    previous_chain <- as.mcmc(pr[[i]])
    current_chain <- as.mcmc(r[[i]])
    previous_mcpar <- attr(previous_chain, "mcpar")
    current_mcpar <- attr(current_chain, "mcpar")

    if (previous_mcpar[3L] != current_mcpar[3L]) {
      stop("MCMC chains must use the same storage interval.", call. = FALSE)
    }
    if (current_mcpar[1L] != previous_mcpar[2L] + previous_mcpar[3L]) {
      stop("MCMC chain iteration indices must be consecutive.", call. = FALSE)
    }

    mcmc(
      rbind(as.matrix(previous_chain), as.matrix(current_chain)),
      start = previous_mcpar[1L],
      end = current_mcpar[2L],
      thin = previous_mcpar[3L]
    )
  })

  mcmc.list(bound)
}
