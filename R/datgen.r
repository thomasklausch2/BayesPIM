#' gen.dat
#'
#' @param x Numeric; Description of x
#' @return Numeric; Description of the return value
#' @export
gen.dat <- function(kappa = 0.7,
                    n = 1000,
                    p = 2,
                    p.discrete = 0,
                    r = 0,
                    s = 1,
                    sigma.X = 1/2,
                    mu.X = 4,
                    beta.X = NULL,
                    beta.W = NULL,
                    theta  = 0.15,
                    Tmax = 20,
                    v.min = 1,
                    v.max = 6,
                    mean.rc = 40,
                    dist.X = 'weibull',
                    k = 1,
                    sel.mod = 'probit',
                    prob.r  = 0) {

  # Sim Z
  R <- matrix(r, p, p)
  diag(R) <- 1
  S <- rep(s, p)
  Sigma <- cor2cov(R, S)
  if (p > 0) {
    Z <- mvrnorm(n, mu = rep(0, p), Sigma)
  } else {
    Z <- NULL
  }

  # Sim discrete Z
  if (p.discrete == 1) {
    Z.discrete <- rbinom(n, 1, 0.5)
    Z <- cbind(Z, Z.discrete)
    colnames(Z) <- paste(1:ncol(Z))
  }

  Z1 <- cbind(as.matrix(rep(1, n)), Z)

  # Sim.X
  if (dist.X == 'weibull' | dist.X == 'loglog' | dist.X == 'lognormal') {
    if (dist.X == 'weibull') {
      e.X <- r.ev(n)
    }
    if (dist.X == 'loglog') {
      e.X <- rlog(n)
    }
    if (dist.X == 'lognormal') {
      e.X <- rnorm(n)
    }
    X <- Z1 %*% c(mu.X, beta.X) + sigma.X * e.X  # log surv times
    X <- exp(X)  # surv times
    X <- as.numeric(X)
  }

  # Simulate X using the generalized gamma distribution
  if (dist.X == 'gengamma') {
    # Parameters for generalized gamma distribution
    lambda <- exp(-Z1 %*% c(mu.X, beta.X))
    gamma <- 1 / sigma.X

    # Generate samples from the generalized gamma distribution
    X <- rggamma(n, a = lambda^(-1), b = gamma, k = k)
  }

  # Screening times
  # Generate screening sequences
  V <- list()
  t.rc <- numeric()
  for (i in 1:n) {
    v_i <- runif(1, v.min, v.max)
    t.rc[i] <- v_i + rexp(1, 1 / mean.rc)
    j <- 1
    while (v_i[j] < t.rc[i]) {
      v_i[j+1] <- runif(1, v_i[j] + v.min, v_i[j] + v.max)
      j <- j + 1
    }
    v_i <- v_i[-length(v_i)]
    v_i <- c(0, v_i, Inf)
    V[[i]] <- v_i
  }

  if(sel.mod == 'probit'){
  mu.W   = Z1 %*% as.matrix(c(qnorm(theta), beta.W)) + rnorm(n)
  g <- as.numeric(mu.W > 0)}
  if(sel.mod == 'logit'){
  g <-rbinom(n,1, 1/(1+exp(-(Z1 %*% as.matrix(c( log(theta/(1-theta)), beta.W))) )))
  }

  r = rbinom(n, 1, prob.r)
  Vobs <- v_to_vobs(V, X, kappa, C=g, baseline.test = r )
  X.true <- X

  ret <- list(Vobs = Vobs, X.true = X.true, Z = Z, C = g, r = r, p.W = pnorm(Z1 %*% as.matrix(c(qnorm(theta), beta.W))))
  ret
}
