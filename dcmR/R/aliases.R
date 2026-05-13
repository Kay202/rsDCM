# =============================================================================
# Internal aliases (not exported).
#
# These mirror the helpers used in the standalone working version of the code
# (helpers_final_optimized.R / exact_core_final_2_optimized.R). Keeping them
# in the package — even though the .has_* flags are always TRUE here because
# Matrix, expm, and MASS are hard Imports — preserves the line-for-line
# structural parity with the upstream working files. That makes it easier
# to port bug fixes or new features (e.g. stochastic DCM, two-state DCM)
# between the two codebases without having to re-translate.
#
# These are evaluated at package-namespace-load time, so they're available
# to every package function from the moment the package is attached.
# =============================================================================

# Runtime presence flags. In the package these are always TRUE because the
# corresponding packages are listed in Imports:, but the conditionals are
# preserved throughout the package source for parity with the working file.
.has_Matrix <- requireNamespace("Matrix", quietly = TRUE)
.has_expm   <- requireNamespace("expm",   quietly = TRUE)
.has_MASS   <- requireNamespace("MASS",   quietly = TRUE)

# Constructor shortcuts. These are direct references — `.Diag(n)` is exactly
# `Matrix::Diagonal(n)` — and exist purely to keep the package source
# visually aligned with the upstream helpers_final_optimized.R file.
if (.has_Matrix) {
  .Diag    <- Matrix::Diagonal
  .spmat   <- Matrix::sparseMatrix
  .spzero  <- function(m, n) Matrix::Matrix(0, nrow = m, ncol = n,
                                            sparse = TRUE)
}
