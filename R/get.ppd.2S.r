#' Posterior Predictive Cumulative Incidence Function
#'
#' Computes the posterior predictive cumulative incidence function (CIF) from a
#' \link{Bayes.2S} prevalence-incidence mixture model. The function can return
#' \emph{quantiles} corresponding to user-specified \emph{percentiles} (i.e., time
#' points at which the cumulative probability reaches certain thresholds), or
#' vice versa (\emph{percentiles} at user-specified \emph{quantiles}). Additionally,
#' it allows for marginal or conditional CIFs of either the mixture population
#' (including prevalence as a point-mass at time zero) or the non-prevalent (healthy)
#' subpopulation.
#'
#' @details
#' For a prevalence-incidence mixture model, some fraction of the population may
#' already have experienced the event (prevalent cases) at baseline, while the
#' remaining (healthy) fraction has not. This function estimates the CIF in two ways:
#' \itemize{
#'   \item \code{type = "xstar"} (mixture CIF): Includes a point-mass at time zero
#'     representing baseline prevalence, with incidence beginning thereafter.
#'   \item \code{type = "x"} (non-prevalent CIF): Excludes prevalent cases,
#'     so it only shows incidence among the initially healthy subpopulation.
#' }
#'
#' You may request a \emph{marginal} CIF by setting both \code{fix_Z.X = NULL} and 
#' \code{fix_Z.W = NULL}, thus integrating over all covariates. Alternatively, a 
#' \emph{conditional} CIF can be obtained by partially or fully specifying fixed 
#' covariate values in \code{fix_Z.X} (and optionally \code{fix_Z.W}) while integrating
#' out the unspecified covariates (\code{NA} entries). 
#'
#' The function operates in two main modes:
#' \itemize{
#'   \item \code{ppd.type = "quantiles"}: Given a set of \code{perc} (cumulative
#'     probabilities), returns corresponding \emph{quantiles} (time points).
#'   \item \code{ppd.type = "percentiles"}: Given a set of \code{quant} (time points),
#'     returns corresponding \emph{percentiles} (cumulative probabilities).
#' }
#'
#' @param mod A fitted prevalence-incidence mixture model of class \code{bayes.2S}.
#' @param fix_Z.X Either \code{NULL} for a marginal CIF or a numeric vector of length
#'   \code{ncol(Z.X)} to request a conditional CIF. Numeric entries fix those covariates
#'   at the given value, whereas \code{NA} entries are integrated out. See \emph{Details}.
#' @param fix_Z.W Same as \code{fix_Z.X} but for the prevalence model covariates;
#'   must be \code{NULL} to obtain a marginal CIF.
#' @param pst.samples Integer; number of posterior samples to draw when computing the
#'   posterior predictive CIF. Must not exceed the total available posterior samples
#'   in \code{mod}. Larger values can improve precision but increase computation time.
#' @param perc A numeric vector of cumulative probabilities (i.e., percentiles in (0,1))
#'   for which time points are returned when \code{ppd.type = "quantiles"}.
#' @param type Character; \code{"xstar"} for the mixture CIF (prevalence + incidence),
#'   or \code{"x"} for the non-prevalent (healthy) population CIF.
#' @param ppd.type Character; \code{"percentiles"} to return cumulative probabilities
#'   at times \code{quant}, or \code{"quantiles"} to return time points at cumulative
#'   probabilities \code{perc}.
#' @param quant A numeric vector of time points for which the function returns
#'   cumulative probabilities when \code{ppd.type = "percentiles"}.
#'
#' @return A \code{list} with some or all of the following elements:
#' \describe{
#'   \item{\code{med.cdf}}{
#'     \itemize{
#'       \item If \code{ppd.type = "quantiles"}, \code{med.cdf} contains the median
#'         quantiles (time points) across posterior samples for each percentile in
#'         \code{perc}.
#'       \item If \code{ppd.type = "percentiles"}, \code{med.cdf} contains the
#'         median cumulative probabilities across posterior samples for each time
#'         point in \code{quant}.
#'     }
#'   }
#'   \item{\code{med.cdf.ci}}{
#'     A 2-row matrix with the 2.5% and 97.5% posterior quantiles for the
#'     estimated \code{med.cdf}, reflecting uncertainty.
#'   }
#'   \item{\code{quant}}{
#'     If \code{ppd.type = "percentiles"}, this is a copy of the input \code{quant}
#'     (time points).
#'   }
#'   \item{\code{perc}}{
#'     If \code{ppd.type = "quantiles"}, this is a copy of the input \code{perc}
#'     (cumulative probabilities).
#'   }
#' }
#'
#' @examples
#' \dontrun{
#' # Generate data according to the Klausch et al. (2024) PIM
#' set.seed(2025)
#' dat <- gen.dat(kappa = 0.7, n = 1e3, theta = 0.2,
#'                p = 1, p.discrete = 1,
#'                beta.X = c(0.2, 0.2), beta.W = c(0.2, 0.2),
#'                v.min = 20, v.max = 30, mean.rc = 80,
#'                sigma.X = 0.2, mu.X = 5, dist.X = "weibull",
#'                prob.r  = 1)
#'
#' # Fit a Bayes.2S model (example with moderate ndraws = 2e4)
#' mod <- bayes.2S(Vobs = dat$Vobs, Z.X = dat$Z, Z.W = dat$Z, r = dat$r,
#'                 kappa = 0.7, update.kappa = FALSE, ndraws = 1e4,
#'                 chains = 4, prop.sd.X = 0.008, parallel = TRUE,
#'                 dist.X = "weibull")
#'
#' ###################
#' # (1) Provide percentiles, get back quantiles (times)
#' ###################
#' cif_nonprev <- get.ppd.2S(mod, pst.samples = 1e3, type = "x",
#'                           ppd.type = "quantiles", perc = seq(0, 1, 0.01))
#' cif_mix     <- get.ppd.2S(mod, pst.samples = 1e3, type = "xstar",
#'                           ppd.type = "quantiles", perc = seq(0, 1, 0.01))
#'
#' # Plot: Non-prevalent stratum CIF vs. mixture CIF (marginal)
#' par(mfrow = c(1,2))
#' plot(cif_nonprev$med.cdf, cif_nonprev$perc, type = "l", xlim = c(0,300), ylim = c(0,1),
#'      xlab = "Time", ylab = "Cumulative Incidence")
#' lines(cif_nonprev$med.cdf.ci[1,], cif_nonprev$perc, lty = 2)
#' lines(cif_nonprev$med.cdf.ci[2,], cif_nonprev$perc, lty = 2)
#'
#' plot(cif_mix$med.cdf, cif_mix$perc, type = "l", xlim = c(0,300), ylim = c(0,1),
#'      xlab = "Time", ylab = "Cumulative Incidence")
#' lines(cif_mix$med.cdf.ci[1,], cif_mix$perc, lty = 2)
#' lines(cif_mix$med.cdf.ci[2,], cif_mix$perc, lty = 2)
#'
#' ###################
#' # (2) Provide quantiles (times), get back percentiles (cumulative probabilities)
#' ###################
#' cif2_nonprev <- get.ppd.2S(mod, pst.samples = 1e3, type = "x",
#'                            ppd.type = "percentiles", quant = 1:300)
#' cif2_mix     <- get.ppd.2S(mod, pst.samples = 1e3, type = "xstar",
#'                            ppd.type = "percentiles", quant = 1:300)
#'
#' # Plot: Non-prevalent vs. mixture CIF using times on the x-axis
#' plot(cif2_nonprev$quant, cif2_nonprev$med.cdf, type = "l", xlim = c(0,300), ylim = c(0,1),
#'      xlab = "Time", ylab = "Cumulative Incidence")
#' lines(cif2_nonprev$quant, cif2_nonprev$med.cdf.ci[1,], lty = 2)
#' lines(cif2_nonprev$quant, cif2_nonprev$med.cdf.ci[2,], lty = 2)
#'
#' plot(cif2_mix$quant, cif2_mix$med.cdf, type = "l", xlim = c(0,300), ylim = c(0,1),
#'      xlab = "Time", ylab = "Cumulative Incidence")
#' lines(cif2_mix$quant, cif2_mix$med.cdf.ci[1,], lty = 2)
#' lines(cif2_mix$quant, cif2_mix$med.cdf.ci[2,], lty = 2)
#'
#' ###################
#' # (3) Conditional CIFs by fixing some covariates
#' ###################
#' cif_mix_m1 <- get.ppd.2S(mod, fix_Z.X = c(-1, NA), pst.samples = 1e3,
#'                          type = "xstar", ppd.type = "quantiles", perc = seq(0,1,0.01))
#' cif_mix_0  <- get.ppd.2S(mod, fix_Z.X = c(0, NA),  pst.samples = 1e3,
#'                          type = "xstar", ppd.type = "quantiles", perc = seq(0,1,0.01))
#' cif_mix_p1 <- get.ppd.2S(mod, fix_Z.X = c(1, NA),  pst.samples = 1e3,
#'                          type = "xstar", ppd.type = "quantiles", perc = seq(0,1,0.01))
#'
#' # Plot: mixture CIF for three different values of the first covariate
#' par(mfrow = c(1,1))
#' plot(cif_mix_m1$med.cdf, cif_mix_m1$perc, type = "l", xlim = c(0,300), ylim = c(0,1),
#'      xlab = "Time", ylab = "Cumulative Incidence", col=1)
#' lines(cif_mix_0$med.cdf,  cif_mix_m1$perc, col=2)
#' lines(cif_mix_p1$med.cdf, cif_mix_m1$perc, col=3)
#' }
#'
#' @export
get.ppd.2S = function (mod, fix_Z.X = NULL, fix_Z.W = NULL, pst.samples = 1e3, perc = seq(0, 1, 0.01), type = "x", 
                       ppd.type = "percentiles", quant = NULL) 
{ 
  
  Vobs = mod$Vobs
  vanilla = mod$vanilla
  max.fu = max(sapply(Vobs, function(x) if (is.infinite(x[length(x)])) return(x[length(x) - 1]) else return(x[length(x)])))
  if (is.null(quant)) 
    quant = seq(0, max.fu, max.fu/1000)
  par.list.X = (mod$par.X.bi)
  Z.X = as.matrix(mod$Z.X)
  Z.W = as.matrix(mod$Z.W)
  dist.X = mod$dist.X
  g.fixed = numeric(length = length(Vobs))
  for (i in 1:length(mod$Vobs)) if (length(mod$Vobs[[i]]) == 1) g.fixed[i] = 1
  
  if( !is.null(fix_Z.X) ){
    if(length(fix_Z.X) != ncol(Z.X)) stop('Length of fix_Z.X has to be equal to ncol(Z.X)')
    Z.X[,!is.na(fix_Z.X)] = fix_Z.X[!is.na(fix_Z.X)]
    if(is.null(fix_Z.W)) fix_Z.W = fix_Z.X
  }
  if( !is.null(fix_Z.X) ){
    if(length(fix_Z.W) != ncol(Z.W)) stop('Length of fix_Z.W has to be equal to ncol(Z.W)')
    Z.W[,!is.na(fix_Z.W)] = fix_Z.W[!is.na(fix_Z.W)]
  }
  
  s = sample(1:(length(par.list.X) * nrow(as.matrix(par.list.X[1]))), 
             pst.samples, replace = F)
  if (vanilla)  ppd = sample.ppd.vanilla(par.list = par.list.X, Z.X = Z.X, 
                             dist.X = dist.X, s = s)
  if (!vanilla) ppd = sample.ppd.xstar(par.list = par.list.X, Z.X = Z.X, 
                           Z.W = Z.W, dist.X = dist.X, s = s, 
                           g.fixed = g.fixed, type = type)
  degenerate <- apply(ppd, 2, function(x) sum(is.na(x)) == length(x))
  if (sum(degenerate) > 0) {
    ppd <- ppd[, !degenerate]
    warning("The prevalence model produced posterior predicted values equal to 1 for all n. This probably signifies a convergence issue.")
  }
  ret = list()
  if (ppd.type == "percentiles") {
    perc = apply(ppd, 2, function(x) {
      e = ecdf(x)
      e(quant)
    })
    ret$med.cdf = apply(perc, 1, median)
    ret$med.cdf.ci = apply(perc, 1, quantile, c(0.025, 0.975))
    ret$quant = quant
  }
  if (ppd.type == "quantiles") {
    q = apply(ppd, 2, quantile, perc, na.rm = T)
    ret$med.cdf = apply(q, 1, median)
    ret$med.cdf.ci = apply(q, 1, quantile, c(0.025, 0.975))
    ret$perc = perc
  }
  ret
}