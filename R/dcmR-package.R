#' dcmR: Dynamic Causal Modelling for Functional MRI
#'
#' Tools for estimating Dynamic Causal Models (DCM) for fMRI data using
#' variational Laplace inversion. The package is an R port of the
#' corresponding routines in the MATLAB SPM12 toolbox, with optimizations
#' for the integration and Gauss-Newton inner loops.
#'
#' The main user-facing entry point is \code{\link{dcm_estimate}},
#' which inverts a DCM specification given BOLD time series and an input
#' design matrix. Lower-level routines (\code{\link{dcm_nlsi_GN}},
#' \code{\link{dcm_int}}) are also exported for advanced use.
#'
#' @section Acknowledgment:
#' This package is a derivative work of the MATLAB SPM12 toolbox,
#' distributed under GPL-2 by the Wellcome Centre for Human Neuroimaging.
#' See the \code{LICENSE.note} file in the package source for details.
#'
#' @references
#' Friston, K.J., Harrison, L., Penny, W. (2003). Dynamic causal modelling.
#' NeuroImage, 19(4), 1273-1302.
#'
#' Friston, K.J., Mattout, J., Trujillo-Barreto, N., Ashburner, J., Penny,
#' W. (2007). Variational free energy and the Laplace approximation.
#' NeuroImage, 34(1), 220-234.
#'
#' @keywords internal
#' @aliases dcmR-package
"_PACKAGE"

# Package-internal environment used to hold mutable runtime options
# (e.g. the finite-difference step size used by dcm_diff).
.dcmR_env <- new.env(parent = emptyenv())

.onLoad <- function(libname, pkgname) {
  # Force Matrix's namespace to load so its S4 methods become discoverable
  # to dispatch from within dcmR. The importMethodsFrom directives in
  # NAMESPACE handle this in principle, but loadNamespace() makes it
  # explicit and tolerant of devtools::load_all() scenarios where the
  # NAMESPACE machinery is partially active.
  loadNamespace("Matrix")
  loadNamespace("expm")

  # Default finite-difference step for dcm_diff
  assign("GLOBAL_DX", exp(-8), envir = .dcmR_env)
  invisible(NULL)
}

#' Get or set dcmR runtime options
#'
#' Read or modify package-internal numerical options. Currently the only
#' option is \code{GLOBAL_DX}, the finite-difference step used by
#' \code{\link{dcm_diff}} when computing numerical Jacobians.
#'
#' @param ... Named arguments to set. Call with no arguments to retrieve
#'   the current option list.
#' @return The current option list (invisibly when setting).
#' @examples
#' dcm_options()                 # show current options
#' dcm_options(GLOBAL_DX = 1e-4) # widen the FD step
#' dcm_options(GLOBAL_DX = exp(-8)) # restore default
#' @export
dcm_options <- function(...) {
  args <- list(...)
  if (length(args) == 0L) {
    return(as.list(.dcmR_env))
  }
  if (is.null(names(args)) || any(names(args) == "")) {
    stop("All arguments to dcm_options() must be named.")
  }
  for (nm in names(args)) assign(nm, args[[nm]], envir = .dcmR_env)
  invisible(as.list(.dcmR_env))
}
