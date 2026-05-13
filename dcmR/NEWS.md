# dcmR 0.1.0

Initial release.

## Highlights

* `spm_dcm_estimate()`, `spm_nlsi_GN()`, `spm_int()` ported from SPM12
  with optimizations: propagator caching in the bilinear integrator,
  identity-projection fast path in `spm_diff()`, diagonal fast paths in
  `spm_inv()` / `spm_logdet()`. Roughly 20x faster end-to-end than a
  literal port on representative two-region DCMs.
* Synthetic `toy_dcm` dataset for documentation, smoke-testing, and
  examples.
* Vignette `introduction` with a worked end-to-end example.

## CRAN-friendly behavioural changes vs. the MATLAB original

* Progress messages from `spm_nlsi_GN()` use `message()` and can be
  silenced with `suppressMessages()`.
* `spm_dcm_estimate()` no longer writes to disk by default. Pass
  `save = TRUE` together with a file path to recover the SPM12
  behaviour.
* The previous `.spm_env` global is replaced by a package-internal
  environment whose mutable settings are exposed through
  `dcm_options()`.
