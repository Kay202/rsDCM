# rsDCM 0.2.0

## Robust and sparse group inference

- Replaced mode-only inclusion summaries with integrated posterior inclusion
  probabilities estimated by a registered C++ parallel-tempering sampler.
- Added batch-means Monte Carlo standard errors, chain-specific PIPs,
  acceptance rates, swap rates, and maximum between-chain PIP ranges.
- Reworked mode fitting as a multistart generalized EM algorithm with
  conditional maximization, monotonicity safeguards, stationarity checks, and
  an `nlminb()` fallback for covariance-component optimization.
- Added conditional Laplace covariance matrices and standard deviations for
  group coefficients.
- Added selectable precision-Cholesky, covariance-Cholesky, and symmetric
  whitening implementations. Returned weights are explicitly identified as
  whitening-coordinate diagnostics.
- Fixed the residual scale by default for scale identification and added an
  opt-in free-scale mode with a warning.
- Added diagnostics for multistart convergence, objectives, Q-function
  monotonicity, and fitted covariance components.
- Replaced the bundled preliminary NARPS summaries with corrected within-
  subject PEB summaries for the same 48 participants.

## Package organization

- Consolidated the novel method in `R/rsDCM.R`.
- Consolidated MATLAB-ported DCM and PEB functions in `R/DCM.R`.
- Moved the C++ sampler to the standard `src/` directory with registered
  native routines.
- Updated manuals, examples, tests, citation information, and the vignette.

# rsDCM 0.1.0

Initial release.
