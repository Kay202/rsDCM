#' rsDCM: Robust and Sparse Dynamic Causal Modelling for Functional MRI
#'
#' rsDCM provides a robust and sparse method for group-level Dynamic Causal
#' Modelling (DCM) of fMRI data, together with an R implementation of
#' single-subject DCM ported from the MATLAB SPM25 toolbox.
#'
#' The package's main method is \code{\link{rsdcm}} (with the convenience
#' wrapper \code{\link{rsdcm_fit}}): a group-level model that weights subjects
#' by a Student-t likelihood (robustness to outliers) and selects group-level
#' effects with a nonlocal product-moment (pMOM) spike-and-slab prior
#' (sparsity), with ReML-estimated between-subject variance components. The
#' \dQuote{rs} is
#' \emph{robust and sparse}; it does not denote resting-state fMRI.
#'
#' For single-subject inversion, the main entry point is
#' \code{\link{dcm_estimate}}.
#'
#' @section Group-level modelling (robust and sparse):
#' \describe{
#'   \item{\code{\link{rsdcm}}}{Robust, sparse group DCM (Student-t + pMOM).}
#'   \item{\code{\link{rsdcm_fit}}}{Assemble \code{rsdcm} inputs from a list of fitted DCMs.}
#'   \item{\code{\link{rsdcm_options}}}{Get or set package numerical options.}
#' }
#'
#' @section Single-subject DCM estimation (SPM25 port):
#' \describe{
#'   \item{\code{\link{dcm_estimate}}}{Full DCM inversion (main entry point).}
#'   \item{\code{\link{dcm_nlsi_GN}}}{Variational Laplace / Gauss-Newton inversion.}
#'   \item{\code{\link{dcm_int}}}{Bilinear-system integrator (forward model).}
#'   \item{\code{\link{dcm_fmri_priors}}}{Construct fMRI DCM priors.}
#'   \item{\code{\link{dcm_fx_fmri}}, \code{\link{dcm_gx_fmri}}}{Neural-state and BOLD observation equations.}
#'   \item{\code{\link{dcm_bireduce}}, \code{\link{dcm_kernels}}}{Bilinear reduction and Volterra kernels.}
#'   \item{\code{\link{dcm_evidence}}, \code{\link{dcm_log_evidence}}, \code{\link{dcm_log_evidence_reduce}}}{Model evidence and Bayesian model reduction.}
#' }
#'
#' @section Parametric Empirical Bayes (Gaussian group modelling):
#' \describe{
#'   \item{\code{\link{dcm_peb_prepare}}, \code{\link{dcm_peb_run}}}{Prepare and fit a second- or third-level PEB.}
#'   \item{\code{\link{dcm_peb_of_pebs}}}{Third-level PEB-of-PEBs over a directory of subject PEBs.}
#'   \item{\code{\link{dcm_peb_design}}, \code{\link{dcm_peb_files}}, \code{\link{dcm_peb_load}}}{Design-matrix and subject-PEB loading helpers.}
#' }
#'
#' @section Datasets:
#' \describe{
#'   \item{\code{\link{toy_dcm}}}{Three-region single-subject DCM specification.}
#'   \item{\code{\link{narps_dcm}}}{48-subject NARPS DCM summaries used by \code{rsdcm}.}
#' }
#'
#' Lower-level numerical helpers (matrix and vector utilities, numerical
#' differentiation, matrix exponentials, and error-covariance bases) are also
#' exported and documented; see, for example, \code{\link{dcm_vec}},
#' \code{\link{dcm_inv}}, \code{\link{dcm_diff}}, and \code{\link{dcm_Ce}}.
#'
#' @section Acknowledgment:
#' The single-subject DCM routines are a derivative work of the MATLAB SPM25
#' (version 25.01.02) toolbox, distributed under GPL-2 by the Wellcome Centre
#' for Human
#' Neuroimaging. See the \code{LICENSE.note} file in the package source for
#' the list of ported routines.
#'
#' @references
#' Arhin, G., Sanyal, N. (2026). Robust and Sparse Group Dynamic Causal
#' Modeling via Student-t Parametric Empirical Bayes and Nonlocal Priors.
#' arXiv:2609.06379. \doi{10.48550/arXiv.2609.06379}
#'
#' Botvinik-Nezer, R., Holzmeister, F., Camerer, C.F., et al. (2020).
#' Variability in the analysis of a single neuroimaging dataset by many teams.
#' Nature, 582(7810), 84-88. \doi{10.1038/s41586-020-2314-9}
#'
#' Friston, K.J., Harrison, L., Penny, W. (2003). Dynamic causal modelling.
#' NeuroImage, 19(4), 1273-1302. \doi{10.1016/S1053-8119(03)00202-7}
#'
#' Friston, K.J., Mattout, J., Trujillo-Barreto, N., Ashburner, J., Penny,
#' W. (2007). Variational free energy and the Laplace approximation.
#' NeuroImage, 34(1), 220-234. \doi{10.1016/j.neuroimage.2006.08.035}
#'
#' Valerio, D., Peres, A., Bergstrom, F., Seidel, P., Almeida, J. (2025).
#' Neural and behavioral similarity-driven tuning curves for manipulable
#' objects. Imaging Neuroscience, 3. \doi{10.1162/imag_a_00482}
#'
#' @aliases rsDCM-package
#'
#' @importFrom stats lm.fit optimize pnorm toeplitz var
#' @importFrom methods as
#'
# The Matrix S4 class/method imports below are load-bearing, not cosmetic: they
# let t(), %*%, solve(), diag() and friends dispatch on Matrix objects from
# inside rsDCM. Do not prune them because they look unused; nothing references
# them by name, and removing them breaks dispatch at runtime, not build time.
#' @importClassesFrom Matrix Matrix sparseMatrix denseMatrix diagonalMatrix
#' @importClassesFrom Matrix dgCMatrix dgTMatrix dgeMatrix ddiMatrix
#' @importMethodsFrom Matrix t solve diag kronecker crossprod tcrossprod
#' @importMethodsFrom Matrix norm determinant
#' @rawNamespace importMethodsFrom(Matrix, "%*%")
"_PACKAGE"

# Package-internal environment used to hold mutable runtime options
# (e.g. the finite-difference step size used by dcm_diff).
.rsDCM_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Force Matrix's namespace to load so its S4 methods become discoverable
  # to dispatch from within rsDCM. The importMethodsFrom directives in
  # NAMESPACE handle this in principle, but loadNamespace() makes it
  # explicit and tolerant of devtools::load_all() scenarios where the
  # NAMESPACE machinery is partially active.
  loadNamespace("Matrix")
  loadNamespace("expm")

  # Default finite-difference step for dcm_diff
  assign("GLOBAL_DX", exp(-8), envir = .rsDCM_env)
  invisible(NULL)
}

#' Get or set rsDCM runtime options
#'
#' Read or modify package-internal numerical options. Currently the only
#' option is \code{GLOBAL_DX}, the finite-difference step used by
#' \code{\link{dcm_diff}} when computing numerical Jacobians.
#'
#' @param ... Named arguments to set. Call with no arguments to retrieve
#'   the current option list.
#' @return The current option list (invisibly when setting).
#' @examples
#' rsdcm_options()                 # show current options
#' rsdcm_options(GLOBAL_DX = 1e-4) # widen the FD step
#' rsdcm_options(GLOBAL_DX = exp(-8)) # restore default
#' @export
rsdcm_options <- function(...) {
  args <- list(...)
  if (length(args) == 0L) {
    return(as.list(.rsDCM_env))
  }
  if (is.null(names(args)) || any(names(args) == "")) {
    stop("All arguments to rsdcm_options() must be named.")
  }
  for (nm in names(args)) assign(nm, args[[nm]], envir = .rsDCM_env)
  invisible(as.list(.rsDCM_env))
}
