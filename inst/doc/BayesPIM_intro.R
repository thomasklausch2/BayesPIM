## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## -----------------------------------------------------------------------------
library(BayesPIM)

## -----------------------------------------------------------------------------
# Generate data according to the Klausch et al. (2024) PIM
set.seed(2025)
dat = gen.dat(
  kappa = 0.7, n = 1e3, theta = 0.2,
  p = 1, p.discrete = 1,
  beta.X = c(0.2, 0.2), beta.W = c(0.2, 0.2),
  v.min = 20, v.max = 30, mean.rc = 80,
  sigma.X = 0.2, mu.X = 5, dist.X = "weibull",
  prob.r = 1
)

## -----------------------------------------------------------------------------
head(dat$Vobs)

## -----------------------------------------------------------------------------
# An initial model fit with fixed test sensitivity kappa (approx. 1-3 minutes, depending on machine)
mod = bayes.2S( Vobs = dat$Vobs,
                Z.X = dat$Z,
                Z.W = dat$Z,
                r= dat$r,
                kappa = 0.7,
                update.kappa = F,
                ndraws= 1e4,
                chains = 4,
                prop.sd.X = 0.008,
                parallel = T,
                dist.X = 'weibull'
)

## ----eval=FALSE---------------------------------------------------------------
# # Inspect results
# mod$runtime # runtime of Gibbs sampler
# plot( trim.mcmc( mod$par.X.all, thining = 10) ) # MCMC chains including burn-in
# plot( trim.mcmc( mod$par.X.bi, thining = 10) ) # MCMC chains excluding burn-in
# summary( mod$par.X.bi)
# apply(mod$ac.X, 2, mean) # Acceptance rates per chain
# gelman.diag(mod$par.X.bi) # Gelman convergence diagnostics

## -----------------------------------------------------------------------------
# An initial model fit with a moderate number of ndraws (here 1e3)
mod.ini = bayes.2S( Vobs = dat$Vobs,
                Z.X = dat$Z,
                Z.W = dat$Z,
                r= dat$r,
                kappa = 0.7,
                update.kappa = F,
                ndraws= 1e3,
                chains = 4,
                prop.sd.X = 0.005,
                parallel = T,
                dist.X = 'weibull'
)

# Running the automated search
search.sd <- search.prop.sd(m = mod.ini)
print(search.sd$prop.sd.X)

## ----eval=FALSE---------------------------------------------------------------
# # Model updating
# mod_update = bayes.2S( prev.run = mod ) # ndraws additional MCMC draws
# mod_update = bayes.2S( prev.run = mod, ndraws.update = 1e3 ) # ndraws.update additional MCMC draws

## -----------------------------------------------------------------------------
# Get posterior predictive marginal CIF, we provide percentiles and obtain quantiles
cif_nonprev <- get.ppd.2S(mod, pst.samples = 1e3, type = 'x', 
                          ppd.type = "quantiles", perc = seq(0, 1, 0.01))
cif_mix     <- get.ppd.2S(mod, pst.samples = 1e3, type = 'xstar', 
                          ppd.type = "quantiles", perc = seq(0, 1, 0.01))

## ----comparison_plots1, fig.width=8, fig.height=4-----------------------------
# Comparison plot non-prevalent stratum CIF vs. mixture CIF (marginal)
par(mfrow = c(1,2))
plot(cif_nonprev$med.cdf, cif_nonprev$perc, ty = 'l', ylim=c(0,1), 
     xlim=c(0,300), xlab = 'Time', ylab ='Cum. prob.', main = 'Non-prevalence CIF')
lines(cif_nonprev$med.cdf.ci[1,], cif_nonprev$perc, lty=2, col=2)
lines(cif_nonprev$med.cdf.ci[2,], cif_nonprev$perc, lty=2, col=2)
plot(cif_mix$med.cdf, cif_nonprev$perc, ty = 'l', ylim=c(0,1), 
     xlim=c(0,300), xlab = 'Time', ylab ='Cum. prob.', main = 'Mixture CIF')
lines(cif_mix$med.cdf.ci[1,], cif_nonprev$perc, lty=2, col=2)
lines(cif_mix$med.cdf.ci[2,], cif_nonprev$perc, lty=2, col=2)

## ----eval=FALSE---------------------------------------------------------------
# # Alternatively, we can provide quantiles and obtain percentiles
# cif2_nonprev <- get.ppd.2S(mod, pst.samples = 1e3, type = 'x',
#                            ppd.type = "percentiles", quant = 1:300)
# cif2_mix     <- get.ppd.2S(mod, pst.samples = 1e3, type = 'xstar',
#                            ppd.type = "percentiles", quant = 1:300)

## -----------------------------------------------------------------------------
get.ppd.2S(mod, pst.samples = 1e3, type = 'xstar', ppd.type = "percentiles", quant = c(0,100))

## -----------------------------------------------------------------------------
# Conditional CIFs, example conditional mixture CIF integrating over one covariate
cif_mix_m1 <- get.ppd.2S(mod, fix_Z.X = c(-1,NA), pst.samples = 1e3, type = 'xstar', 
                         ppd.type = "quantiles", perc = seq(0, 1, 0.01))
cif_mix_0  <- get.ppd.2S(mod, fix_Z.X = c(0,NA), pst.samples = 1e3, type = 'xstar', 
                         ppd.type = "quantiles", perc = seq(0, 1, 0.01))
cif_mix_p1 <- get.ppd.2S(mod, fix_Z.X = c(1,NA), pst.samples = 1e3, type = 'xstar', 
                         ppd.type = "quantiles", perc = seq(0, 1, 0.01))

## ----comparison_plots2, fig.width=4, fig.height=4-----------------------------
# Plot of CIFs for three levels of the first covariate
par(mfrow = c(1,1))
plot(cif_mix_m1$med.cdf, cif_mix_m1$perc, ty = 'l', ylim=c(0,1), 
     xlim=c(0,300), xlab = 'Time', ylab ='Cum. prob.', main = 'Conditional mixture CIF')
lines(cif_mix_0$med.cdf, cif_mix_m1$perc, col=2)
lines(cif_mix_p1$med.cdf, cif_mix_m1$perc, col=3)

## -----------------------------------------------------------------------------
# An exponential model
mod_exp = bayes.2S( Vobs = dat$Vobs,
                Z.X = dat$Z,
                Z.W = dat$Z,
                r= dat$r,
                kappa = 0.7,
                update.kappa = F,
                ndraws= 1e4,
                chains = 4,
                prop.sd.X = 0.008,
                parallel = T,
                fix.sigma.X = T,
                sig.prior.X = 1,
                dist.X = 'weibull'
)

# Get information criteria
get.IC_2S(mod, samples = 1e3)
get.IC_2S(mod_exp, samples = 1e3)

