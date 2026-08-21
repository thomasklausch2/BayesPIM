## Submission

This is a major update of BayesPIM, from 1.0.1 to 2.0.

Version 2.0 is a substantial rewrite. All exported functions and arguments were
renamed to snake_case, two slice samplers were added (one of which is the new
default), the generalized-gamma model was reparameterised to the Prentice form,
and post-estimation was reorganised around `summary()`, `plot()`, and `ppCIF()`
methods. Because this breaks code written for 1.0.1, the changes are listed in
NEWS.md.

Examples, tests, and vignette code use at most two parallel workers, in line
with CRAN policy. The vignette does not fit models when built; its results are
precomputed and included as text and figures.

## Test environments

* local: macOS 26.5, R 4.6.0 (aarch64-apple-darwin23)
* win-builder: R-devel

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies.
