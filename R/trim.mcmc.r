#' Subset MCMC draws (burn-in and thinning)
#'
#' Takes an \code{mcmc.list} object (or a list of MCMC chains) and returns
#' a new \code{mcmc.list} containing only the specified subset of iterations (from
#' \code{burnin} to \code{end}) with the specified thinning interval.
#'
#' @param obj An object of class \code{mcmc.list} (or a list of matrices) 
#'            containing MCMC draws.
#' @param burnin A numeric scalar giving the starting iteration of the MCMC 
#'               sample to keep. Defaults to \code{1}.
#' @param end A numeric scalar giving the last iteration of the MCMC sample 
#'            to keep. Defaults to the number of rows in the first chain of 
#'            \code{obj}.
#' @param thining A numeric scalar for the thinning interval. Defaults to \code{1}.
#'
#' @details
#' This function subsets each chain of the input \code{obj} to the 
#' specified iteration indices and creates a new \code{mcmc.list}.
#' If you have multiple MCMC chains, each chain is trimmed in the same way.
#'
#' @return An object of class \code{mcmc.list}, representing the trimmed subset 
#'         of the original MCMC draws.
#'
#' @examples
#' # Example with a toy mcmc.list
#' set.seed(123)
#' x1 <- matrix(rnorm(2000), ncol = 2)
#' x2 <- matrix(rnorm(2000), ncol = 2)
#' mcmc_list <- mcmc.list(mcmc(x1), mcmc(x2))
#'
#' # Trim and thin the chains
#' trimmed_mcmc <- trim.mcmc(mcmc_list, burnin = 100, end = 800, thining = 5)
#' summary(trimmed_mcmc)
#'
#' @export
trim.mcmc <- function(obj, burnin = 1, end = nrow(as.matrix(obj[[1]])), thining = 1) {
  mcmc.list(
    lapply(obj, function(x) {
      mcmc(
        x[seq(burnin, end, by = thining), ],
        start = burnin, end = end, thin = thining
      )
    })
  )
}
