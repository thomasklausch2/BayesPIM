# BayesPIM 2.0

## Breaking changes

* Version 2.0 is a substantial rewrite: the user-facing API, the sampler, and
  the post-estimation tooling have all changed, and code written for 1.0.1 will
  not run unchanged.
* All exported functions and arguments were renamed to snake_case.
  `bayes.2S()` is now `bayespim()`, `gen.dat()` is `gen_data()`, `get.IC_2S()`
  is `get_ic()`, `trim.mcmc()` is `trim_mcmc()`, and `search.prop.sd()` is
  `search_prop_sd()`. Arguments follow the same convention (`Vobs` -> `v_obs`,
  `Z.X` -> `x_t`, `Z.W` -> `x_g`, `dist.X` -> `dist`, `tau.w` -> `tau_g`).
* `bayes.2S_seq()` and `search.prop.sd_seq()` are removed; the parallel and
  sequential code paths are unified.
* `get.ppd.2S()` is replaced by `ppCIF()`, which computes the mixture and the
  non-prevalent cumulative incidence function in a single call and has a
  `plot()` method.
* The `thining`, `conv.crit`, `parallel`, `vanilla`, and `ndraws.naive`
  arguments are gone. Draw storage is controlled by `save_every` and the
  warm-up cutoff by `warmup`.
* In `gen_data()`, the covariate-correlation argument is renamed from `r` to
  `rho`, so that `r` unambiguously denotes the baseline-test indicator, as it
  does in `bayespim()` and in the returned `$r`.
* Two defaults changed: `kappa` no longer defaults to `0.5` and must be given
  explicitly when `update_kappa = FALSE`, and the effective-sample-size target
  `min_effss` rose from `chains * 10` to `chains * 100`.
* The generalized gamma is reparameterised to the Prentice form (location,
  scale, signed shape `Q`) using **flexsurv**, replacing the **ggamma**
  parameterisation of 1.0.1. Fitted shape values are not comparable across
  versions, and the model is now available only with the collapsed sampler.

## New features

* **Two slice samplers.** `sampler = "slice_collapsed"`, the new default,
  augments only the latent screening interval and updates the incidence
  parameters from the interval-censored likelihood; `sampler = "slice"`
  augments exact event times. Both show lower autocorrelation and faster
  convergence than the Metropolis-Hastings sampler of 1.0.1, which remains
  available as `sampler = "mh"`.
* **Gamma incidence times** (`dist = "gamma"`), parameterised through the
  conditional mean and coefficient of variation.
* **`summary()` and `plot()` methods** for fitted models, reporting posterior
  quantiles and convergence diagnostics for each parameter block.
* **Revised convergence assessment** using rank-normalised split R-hat and
  effective sample size from the **posterior** package, replacing
  `coda::gelman.diag()`. `update_till_converge = TRUE` extends sampling
  automatically until `max_rhat` and `min_effss` are met.
* **Reproducibility.** `seed_chains` sets one seed per chain and the
  end-of-chain RNG state is stored, so a run continued through `prev_run` is
  identical to an uninterrupted run of the same length.
* **Covariate standardisation** (`standardize_covariates`) and internal time
  rescaling (`rescale_times`), both enabled by default, with returned
  coefficients on the original scale.
* **User-supplied priors** through `log_prior_fun`; the default is exported as
  `log_aft_prior()`.
* **Input validation.** Malformed inputs produce a single, specific error
  before sampling begins.
* Further additions: `save_every` to limit memory, `silent` to suppress
  progress output, and `fix_q` to hold the generalized-gamma shape fixed.
* **Example fit.** `data(mod)` provides a converged model so post-estimation
  examples run without refitting.
* **Rewritten vignette** (`vignette("BayesPIM_intro")`): a user guide covering
  estimation, convergence, model comparison, and posterior CIFs.
* **Test suite.** The package now ships contract tests covering the samplers,
  distributions, draw storage, covariate scaling, and post-estimation
  functions.

## Other

* Internal code was consolidated and reduced: duplicate parameterisations,
  unreachable branches, and redundant state were removed, and naming is
  consistently snake_case throughout.
* Dependencies: **posterior**, **survival**, and **flexsurv** added;
  **ggamma** and **mvtnorm** removed; **coda** moved to `Depends`.
* Smaller bug fixes and behind-the-scenes efficiency updates and leaner code.
