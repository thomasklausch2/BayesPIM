#' Distribution function of the log-logistic distribution
#' @noRd
ploglog = function(x,lambda, gamma, lower_tail = TRUE, log_p = FALSE){ 
  pllogis(x, shape = gamma, scale= lambda, 
          lower.tail = lower_tail, log.p = log_p ) 
}
