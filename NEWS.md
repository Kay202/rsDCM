# rsDCM 0.1.0

First release.

## Robust and sparse group-level DCM

* `rsdcm()` fits a robust, sparse group-level model to subject-level DCM
  parameters: Student-t weighting of subjects for robustness to outliers,
  a nonlocal product-moment (pMOM) spike-and-slab prior for sparse
  selection of group effects, and ReML-estimated between-subject variance
  components. It
  returns the group effects, their posterior inclusion probabilities, the
  per-subject robustness weights, and the learned spike/slab scales. See
  <arXiv:2609.06379>.
* `rsdcm_fit()` is a convenience wrapper that assembles the inputs from a
  list of `dcm_estimate()` fits (the posterior `Ep`/`Cp` restricted to a
  chosen parameter field) and a between-subject design.
* `narps_dcm`: subject-level DCM summaries for 48 subjects (with
  group/gender/age covariates), derived from the openly shared NARPS
  dataset. It is the runnable real-data example for `rsdcm()`.
* Progress is reported via `message()` and is silenceable, matching the
  rest of the package.

## Dynamic Causal Modelling for fMRI (R port of SPM25)

* `dcm_estimate()`, `dcm_nlsi_GN()`, `dcm_int()` ported from SPM25 with
  optimizations: propagator caching in the bilinear integrator,
  identity-projection fast path in `dcm_diff()`, diagonal fast paths in
  `dcm_inv()` / `dcm_logdet()`. Roughly 20x faster end-to-end than a
  literal port on representative DCMs.
* `toy_dcm` dataset (three regions, 482 scans) for documentation,
  smoke-testing, and examples.
* Vignette `introduction` with a worked end-to-end example.
* Help pages for all exported functions.

## Parametric Empirical Bayes (PEB)

* Hierarchical (group-level) modelling of DCM parameters across subjects,
  ported from SPM25's `spm_dcm_peb`:
  * `dcm_peb_prepare()` / `dcm_peb_run()`: prepare and fit a second- or
    third-level PEB by variational Laplace.
  * `dcm_peb_of_pebs()`: third-level PEB-of-PEBs over a directory of
    subject PEBs, with a between-subject design.
  * `dcm_peb_design()`: build the group-level design matrix.
  * `dcm_peb_files()` / `dcm_peb_load()`: locate and load subject PEBs
    from `.rds` or MATLAB `.mat` files (`.mat` needs the suggested
    `R.matlab` package).
* `dcm_peb_of_pebs()` does not write to disk unless a `save_path` is supplied.

## Naming

Ported SPM routines keep the `dcm_` prefix (from `spm_`), to avoid
colliding with an installed SPM-derived package and to make it obvious
which functions come from here: `spm_dcm_estimate` becomes `dcm_estimate`,
`spm_nlsi_GN` becomes `dcm_nlsi_GN`, `spm_int` becomes `dcm_int`, and so on.

## Supported model variants

* This release implements only the **deterministic, single-state** fMRI DCM
  (non-linear DCM, via a non-empty `d` array, is also supported).
* The **two-state**, **stochastic**, and **spectral / cross-spectral-density
  (`induced`)** variants are not yet implemented. `dcm_estimate()`
  **stops with an informative error** if any of these options is switched
  on, instead of silently running the deterministic model. They are planned
  for a future update.

## Behavioural changes vs. the MATLAB original

* Progress messages from `dcm_nlsi_GN()` use `message()` and can be
  silenced with `suppressMessages()`.
* `dcm_estimate()` does not write to disk by default. Pass `save = TRUE`
  together with a file path to recover the SPM25 behaviour.
* The previous `.spm_env` global is replaced by a package-internal
  environment whose mutable settings are exposed through `rsdcm_options()`.
