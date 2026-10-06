# Internal aliases (not exported).
#
# These mirror the helpers in the standalone working files
# (helpers_final_optimized.R / exact_core_final_2_optimized.R). The .has_*
# flags are always TRUE here (Matrix, expm and MASS are hard Imports) but are
# kept for parity with those files. They load at namespace-load time, so every
# package function can use them.

.has_Matrix <- requireNamespace("Matrix", quietly = TRUE)
.has_expm   <- requireNamespace("expm",   quietly = TRUE)
.has_MASS   <- requireNamespace("MASS",   quietly = TRUE)

# Constructor shortcuts: .Diag is Matrix::Diagonal, etc.
if (.has_Matrix) {
  .Diag    <- Matrix::Diagonal
  .spmat   <- Matrix::sparseMatrix
  .spzero  <- function(m, n) Matrix::Matrix(0, nrow = m, ncol = n,
                                            sparse = TRUE)
}
