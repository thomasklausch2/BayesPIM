#' Handle one BayesPIM convergence check
#' @noRd
handle_bayespim_convergence <- function(
    ret,
    min_effss,
    maxit,
    ndraws_update,
    max_rhat = 1.01,
    update_till_converge = TRUE,
    silent = FALSE
) {
  convergence <- assess_bayespim_convergence(
    ret = ret,
    max_rhat = max_rhat,
    min_ess = min_effss
  )
  ret$convergence <- convergence

  if (!isTRUE(update_till_converge)) {
    print_bayespim_convergence(
      silent = silent,
      convergence = convergence,
      status = sprintf(
        paste0(
          "Convergence diagnostics after %d iterations per chain ",
          "(%d stored post-warm-up draws used)."
        ),
        convergence$n_iter,
        convergence$n_iter_diagnostic
      )
    )
    return(list(ret = ret, needs_update = FALSE, ndraws_next = NULL))
  }

  if (isTRUE(convergence$converged)) {
    print_bayespim_convergence(
      silent = silent,
      convergence = convergence,
      status = sprintf(
        paste0(
          "Converged after %d iterations per chain ",
          "(%d stored post-warm-up draws used)."
        ),
        convergence$n_iter,
        convergence$n_iter_diagnostic
      )
    )
    return(list(ret = ret, needs_update = FALSE, ndraws_next = NULL))
  }

  if (convergence$n_iter >= maxit) {
    print_bayespim_convergence(
      silent = silent,
      convergence = convergence,
      status = sprintf(
        "Stopped at maxit = %g without convergence.",
        maxit
      )
    )
    return(list(ret = ret, needs_update = FALSE, ndraws_next = NULL))
  }

  ndraws_next <- ndraws_update
  if (is.finite(maxit)) {
    ndraws_next <- min(ndraws_next, maxit - convergence$n_iter)
  }

  if (ndraws_next < 2L) {
    print_bayespim_convergence(
      silent = silent,
      convergence = convergence,
      status = sprintf(
        paste0(
          "Stopped after %d iterations per chain without convergence: ",
          "fewer than two draws remain before maxit."
        ),
        convergence$n_iter
      )
    )
    return(list(ret = ret, needs_update = FALSE, ndraws_next = NULL))
  }

  print_bayespim_convergence(
    silent = silent,
    convergence = convergence,
    status = sprintf(
      "Not converged after %d iterations per chain; updating with %d iterations.",
      convergence$n_iter,
      as.integer(ndraws_next)
    )
  )

  list(
    ret = ret,
    needs_update = TRUE,
    ndraws_next = as.integer(ndraws_next)
  )
}
