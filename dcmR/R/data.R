#' Example DCM for fMRI inversion
#'
#' A two-region Dynamic Causal Model used in the documentation, vignette,
#' and tests. The object is an input specification (no posterior fields)
#' so that running \code{\link{dcm_estimate}(toy_dcm)} performs the full
#' inversion from scratch.
#'
#' @format A list with the standard SPM-style DCM fields:
#' \describe{
#'   \item{a, b, c, d}{Connectivity, modulatory, driving, and non-linear
#'     adjacency arrays.}
#'   \item{Y}{Observed BOLD data. \code{Y$y} is the time-series matrix
#'     (rows = scans, columns = regions); \code{Y$dt} is the sampling
#'     interval in seconds; \code{Y$X0} is the confound design.}
#'   \item{U}{Input design. \code{U$u} is the input matrix at the
#'     microtime resolution; \code{U$dt} is the microtime step.}
#'   \item{n, v}{Number of regions and number of scans, respectively.}
#'   \item{TE}{Echo time in seconds.}
#'   \item{options}{Estimation options: \code{nonlinear}, \code{two_state},
#'     \code{stochastic}, \code{centre}, \code{induced}, \code{maxit},
#'     \code{maxnodes}, \code{hE}, \code{hC}.}
#' }
#'
#' @source The dataset was derived from a real fMRI Dynamic Causal Modelling
#'   exercise carried out during package development. The build script that
#'   prepared it from the raw source is in \code{data-raw/make-toy-dcm.R}.
#'
#' @examples
#' data(toy_dcm)
#' str(toy_dcm, max.level = 1)
#' \dontrun{
#' fit <- dcm_estimate(toy_dcm)
#' round(fit$Ep$A, 3)
#' }
"toy_dcm"
