#' Print the current BayesPIM convergence assessment
#' @noRd
print_bayespim_convergence <- function(convergence, status, silent = FALSE) {
  if (isTRUE(silent)) {
    return(invisible(NULL))
  }

  table_text <- paste(
    capture.output(print(convergence$table, row.names = FALSE)),
    collapse = "\n"
  )

  message(
    sprintf(
      "%s\n%s\n%s",
      status,
      convergence$criteria,
      table_text
    )
  )
}

