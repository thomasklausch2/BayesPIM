#' @keywords internal
"_PACKAGE"

# `idx` and `j` are foreach() iteration variables, bound by non-standard
# evaluation inside %dopar% in get_ic() and bayespim(). Declaring them avoids
# a spurious "no visible binding for global variable" note.
utils::globalVariables(c("idx", "j"))

# The following block is used by usethis to automatically manage
# roxygen namespace tags. Modify with care!
## usethis namespace: start
#' @useDynLib BayesPIM, .registration = TRUE
#' @importFrom stats dexp dgamma dlnorm dnorm dt dweibull ecdf median optim pexp pgamma plnorm pnorm pweibull qexp qgamma qlnorm qnorm quantile qweibull rbeta rbinom rchisq rexp rgamma rlnorm rmultinom rnorm runif rweibull var
#' @importFrom coda gelman.diag effectiveSize mcmc
#' @importFrom actuar dllogis pllogis qllogis rllogis
#' @importFrom utils capture.output
#' @import Rcpp coda MASS doParallel foreach parallel
## usethis namespace: end
NULL
