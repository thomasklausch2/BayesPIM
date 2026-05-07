## Resubmission

This is a resubmission of BayesPIM 1.0.1 after archival of 1.0.0 from CRAN.

In this version I have:

* Fixed an issue in the starting-value search for `bayes.2S()` and `bayes.2S_seq()`, where the
  preliminary naive run could spawn more parallel workers than intended.
* Revised examples so that CRAN-executed do not spawn more than 2 parallel workers.
* Moved long-running MCMC examples to `\dontrun{}`.
* Removed non-portable output redirection to `NUL`.
* Added an input validation routine

## R CMD check results

0 errors | 0 warnings | 0 notes

## Downstream dependencies

There are no downstream dependencies.