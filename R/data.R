#' Example DCM for fMRI inversion
#'
#' A three-region Dynamic Causal Model used in the documentation, vignette,
#' and tests. It is an input specification (no posterior fields), ready to
#' invert with \code{\link{dcm_estimate}}.
#'
#' The model has 3 regions and 482 scans, with a single driving input
#' entering region 3 and a single modulatory input. All entries of the
#' \code{a} matrix are enabled, so every directed connection between the
#' three regions is estimated.
#'
#' @format A list with the standard SPM-style DCM fields:
#' \describe{
#'   \item{a}{3 x 3 binary matrix of endogenous connections to estimate.}
#'   \item{b}{3 x 3 x 1 binary array of modulatory connections.}
#'   \item{c}{3 x 1 binary matrix of driving inputs; here the input enters
#'     region 3 only.}
#'   \item{Y}{Observed BOLD data. \code{Y$y} is the 482 x 3 time-series
#'     matrix (rows = scans, columns = regions); \code{Y$dt} is the
#'     sampling interval in seconds; \code{Y$X0} is the confound design;
#'     \code{Y$Q} is the error-covariance component list; \code{Y$name}
#'     holds the region names.}
#'   \item{U}{Input design. \code{U$u} is the 482 x 1 input matrix at the
#'     microtime resolution; \code{U$dt} is the microtime step;
#'     \code{U$name} holds the input names.}
#'   \item{n, v}{Number of regions (3) and number of scans (482).}
#'   \item{TE}{Echo time in seconds.}
#'   \item{name}{Name of the model.}
#'   \item{xY}{Region (VOI) structure carried over from the source DCM.}
#'   \item{options}{Estimation options: \code{nonlinear}, \code{two_state},
#'     \code{stochastic}, \code{centre}, \code{mode}, \code{maxit},
#'     \code{induced}, \code{maxnodes}, \code{hE}, \code{hC}.}
#' }
#'
#' @source Derived from the fMRI data of Valerio et al. (2025), "Neural and
#'   behavioral similarity-driven tuning curves for manipulable objects",
#'   Imaging Neuroscience, 3, \doi{10.1162/imag_a_00482}, via a Dynamic Causal
#'   Modelling exercise carried out with that data during package development.
#'   The build script that prepared the shipped object is in
#'   \code{data-raw/make-toy-dcm.R}. See the package \code{CITATION} for the
#'   full author list.
#'
#' @examples
#' data(toy_dcm)
#' str(toy_dcm, max.level = 1)
#' # Invert it with dcm_estimate(toy_dcm); see ?dcm_estimate for a runnable
#' # (donttest) inversion example.
"toy_dcm"

#' Subject-level DCM summaries for a group analysis (NARPS)
#'
#' Subject-level DCM posterior summaries for 48 subjects, used to demonstrate
#' the robust and sparse group-level model \code{\link{rsdcm}}. The summaries
#' are the inputs a group analysis needs: a posterior mean and covariance per
#' subject.
#'
#' Each subject has 22 parameters: the 16 intrinsic connections of a
#' four-region DCM (A) plus 6 task-modulatory connections (B), across the
#' regions vmPFC, vStr (ventral striatum), amygdala, and anterior insula.
#'
#' @format A list with:
#' \describe{
#'   \item{eta}{48 x 22 matrix of subject-level posterior means (subjects in
#'     rows, parameters in columns).}
#'   \item{Cp}{Length-48 list of 22 x 22 posterior covariance matrices.}
#'   \item{parameter_names}{Length-22 character vector, e.g. \code{"A(1,1)"},
#'     \code{"B(2,1,1)"}.}
#'   \item{subject}{Length-48 character vector of de-identified subject ids.}
#'   \item{covariates}{Data frame of between-subject covariates: \code{group},
#'     \code{gender}, \code{age}.}
#'   \item{regions}{The four region labels, in order.}
#' }
#'
#' @source Subject-level Dynamic Causal Modelling summaries computed by the
#'   package author from the openly shared NARPS dataset: Botvinik-Nezer, R.,
#'   Holzmeister, F., Camerer, C. F., et al. (2020), "Variability in the
#'   analysis of a single neuroimaging dataset by many teams", Nature,
#'   582(7810), 84-88, \doi{10.1038/s41586-020-2314-9}. The build script that
#'   prepared the shipped object is in \code{data-raw/make-narps-dcm.R}.
#'
#' @seealso \code{\link{rsdcm}}, which this dataset is the example
#'   input for.
#' @examples
#' data(narps_dcm)
#' dim(narps_dcm$eta)               # 48 subjects x 22 parameters
#' length(narps_dcm$Cp)             # one posterior covariance per subject
#' table(narps_dcm$covariates$group)
"narps_dcm"
