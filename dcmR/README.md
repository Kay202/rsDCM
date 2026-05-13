# dcmR

<!-- badges: start -->
[![R-CMD-check](https://github.com/Kay202/dcmR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Kay202/dcmR/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

`dcmR` is an R package for Dynamic Causal Modelling (DCM) of functional
MRI data. It is an R port of the corresponding routines in the MATLAB
SPM12 toolbox, with optimizations for the integration and Gauss-Newton
inner loops.

## Installation

From GitHub (development version):

```r
# install.packages("remotes")
remotes::install_github("Kay202/dcmR")
```

From CRAN, once published:

```r
install.packages("dcmR")
```

## Quick start

```r
library(dcmR)

# Synthetic 2-region toy DCM shipped with the package
data(toy_dcm)

# Invert
fit <- dcm_estimate(toy_dcm)

# Posterior connectivity
round(fit$Ep$A, 3)

# Variational free energy (log-evidence proxy)
fit$F
```

The introductory vignette has more detail:

```r
vignette("introduction", package = "dcmR")
```

## What's exported

User-facing entry points:

- `dcm_estimate()` — full DCM inversion.
- `dcm_fmri_priors()` — build priors from adjacency arrays.
- `dcm_nlsi_GN()` — variational Laplace optimisation (lower-level).
- `dcm_int()` — bilinear-system integrator (lower-level).
- `dcm_fx_fmri()`, `dcm_gx_fmri()` — neural and BOLD equations.
- `dcm_evidence()` — AIC / BIC summary of a fit.
- `dcm_options()` — adjust runtime options (e.g. FD step).

Lower-level utilities (`dcm_vec`, `dcm_inv`, `dcm_logdet`, ...) are
documented but kept as internal helpers; access them via
`dcmR:::dcm_vec` if needed.

## Acknowledgment

`dcmR` is a derivative work of the MATLAB SPM12 toolbox, distributed
under GPL-2 by the Wellcome Centre for Human Neuroimaging. See
[`LICENSE.note`](LICENSE.note) for the list of ported routines and
the references.

## License

GPL-2.

## References

- Friston, K.J., Harrison, L., Penny, W. (2003). Dynamic causal
  modelling. *NeuroImage*, **19**(4), 1273-1302.
- Friston, K.J., Mattout, J., Trujillo-Barreto, N., Ashburner, J.,
  Penny, W. (2007). Variational free energy and the Laplace
  approximation. *NeuroImage*, **34**(1), 220-234.
