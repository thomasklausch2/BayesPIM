#' Solve for beta shape parameters matching a target mean and standard deviation
#' @noRd
find_ab = function(m, s){
  if (!(is.numeric(m) && length(m) == 1L && is.finite(m) &&
        m > 0 && m < 1)) {
    stop("The Beta prior mean `m` must be a single finite number strictly between 0 and 1.", call. = FALSE)
  }
  if (!(is.numeric(s) && length(s) == 1L && is.finite(s) && s > 0)) {
    stop("The Beta prior standard deviation `s` must be a single positive finite number.", call. = FALSE)
  }
  
  max_sd = sqrt(m * (1 - m))
  if (s >= max_sd) {
    stop(
      sprintf(
        "The requested Beta prior is not feasible: `s` must be smaller than sqrt(m * (1 - m)) = %.4f.",
        max_sd
      ),
      call. = FALSE
    )
  }
  
  concentration = m * (1 - m) / s^2 - 1
  a = m * concentration
  b = (1 - m) * concentration
  
  if (!all(is.finite(c(a, b))) || any(c(a, b) <= 0)) {
    stop(
      paste0(
        "The requested mean and standard deviation imply Beta shape parameters ",
        "that are not finite and positive. Choose a less extreme `kappa_prior`."
      ),
      call. = FALSE
    )
  }
  
  list(
    a = a,
    b = b,
    check_mean = a / (a + b),
    check_sd = sqrt(a * b / ((a + b)^2 * (a + b + 1)))
  )
}
