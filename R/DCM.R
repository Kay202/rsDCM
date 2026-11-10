# SPDX-License-Identifier: GPL-2.0-only
# MATLAB-ported DCM and PEB functions. See LICENSE.note.

# ---- aliases ----

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

# ---- utils-vec ----

# Vector / list flattening utilities (internal). Performance-critical.

#' Flatten a nested numeric structure into a vector
#'
#' Recursively walks a numeric, logical, or nested list structure and returns
#' a single flat numeric vector. Internal helper; mirrors SPM25's \code{spm_vec}.
#'
#' @param X A numeric, logical, or list (possibly nested).
#' @param ... Additional structures to flatten and concatenate.
#' @return A numeric vector.
#' @examples
#' # A DCM parameter structure flattens in field order
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' dcm_vec(P)
#'
#' # Nested lists are walked recursively
#' dcm_vec(list(a = 1, b = list(c = 2:3, d = 4)))
#' @seealso \code{\link{dcm_unvec}} for the inverse operation.
#' @keywords internal
#' @export
dcm_vec <- function(X, ...) {
  args <- list(...)
  if (length(args) > 0) X <- c(list(X), args)

  if (is.numeric(X))  return(as.vector(as.numeric(X)))
  if (is.logical(X))  return(as.numeric(X))

  .fill <- function(X, buf, off) {
    if (is.numeric(X) || is.logical(X)) {
      v <- as.numeric(X)
      n <- length(v)
      if (n > 0) buf[off + seq_len(n)] <- v
      return(list(buf = buf, off = off + n))
    }
    if (is.list(X)) {
      items <- if (!is.null(names(X)) && all(names(X) != "")) X[names(X)] else X
      for (elem in items) {
        r   <- .fill(elem, buf, off)
        buf <- r$buf
        off <- r$off
      }
      return(list(buf = buf, off = off))
    }
    list(buf = buf, off = off)
  }

  total <- dcm_length(X)
  if (total == 0L) return(numeric(0))
  buf <- numeric(total)
  r   <- .fill(X, buf, 0L)
  r$buf
}

#' Total number of numeric entries in a nested structure
#'
#' Returns the length \code{dcm_vec(X)} would produce, without allocating it.
#'
#' @param X A numeric, logical, or list (possibly nested).
#' @return Integer length.
#' @examples
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' dcm_length(P)            # 6
#' length(dcm_vec(P))       # same, but allocates the vector
#' @keywords internal
#' @export
dcm_length <- function(X) {
  if (is.numeric(X) || is.logical(X)) return(length(X))
  if (!is.list(X)) return(0L)

  stack <- X
  total <- 0L
  while (length(stack) > 0) {
    item  <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (is.numeric(item) || is.logical(item)) {
      total <- total + length(item)
    } else if (is.list(item)) {
      stack <- c(stack, item)
    }
  }
  total
}

#' Reshape a flat vector back into a template structure
#'
#' Inverse of \code{\link{dcm_vec}}: splits a flat numeric vector across the
#' shape implied by one or more template arguments.
#'
#' @param vX Flat numeric vector.
#' @param ... One or more templates whose structure is copied; \code{vX} is
#'   distributed across them in order.
#' @return A structure (or list of structures) matching the templates.
#' @examples
#' # Round-trip: flatten a structure and rebuild it
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' identical(dcm_unvec(dcm_vec(P), P), P)
#'
#' # The template supplies the shape; the vector supplies the values
#' tmpl <- list(A = matrix(0, 2, 2), C = numeric(2))
#' dcm_unvec(1:6, tmpl)
#' @seealso \code{\link{dcm_vec}} for the forward operation.
#' @keywords internal
#' @export
dcm_unvec <- function(vX, ...) {
  vX <- as.numeric(vX)
  args <- list(...)

  if (length(args) > 1) {
    out <- vector("list", length(args))
    off <- 0L
    for (k in seq_along(args)) {
      nk <- dcm_length(args[[k]])
      out[[k]] <- dcm_unvec(if (nk > 0) vX[off + seq_len(nk)] else numeric(0),
                            args[[k]])
      off <- off + nk
    }
    return(out)
  }

  X <- if (length(args) == 1L) args[[1L]] else args
  if (!(is.numeric(vX) || is.logical(vX))) vX <- dcm_vec(vX)

  if (is.numeric(X) || is.logical(X)) {
    X[] <- vX
    return(X)
  }

  env <- new.env(parent = emptyenv())
  env$off <- 0L

  .fill_unvec <- function(X) {
    if (is.numeric(X) || is.logical(X)) {
      n <- length(X)
      if (n > 0) {
        X[] <- vX[env$off + seq_len(n)]
        env$off <- env$off + n
      }
      return(X)
    }
    if (is.list(X)) {
      for (i in seq_along(X)) X[[i]] <- .fill_unvec(X[[i]])
      return(X)
    }
    X
  }

  .fill_unvec(X)
}

#' Zero-fill a structure
#'
#' Returns a copy of \code{X} with the same shape but every numeric entry
#' replaced by 0.
#'
#' @param X A nested structure.
#' @return A structure with zeros, same shape as \code{X}.
#' @examples
#' dcm_zeros(list(A = matrix(1:4, 2), C = c(5, 6)))
#' @keywords internal
#' @export
dcm_zeros <- function(X) dcm_unvec(rep(0, dcm_length(X)), X)

# ---- utils-matrix ----

# Matrix utilities (internal): dense/sparse-agnostic wrappers over Matrix/MASS.

#' Inverse of an ill-conditioned matrix
#'
#' Computes \code{solve(A + tol*I)} with an automatically chosen tolerance,
#' matching SPM25's \code{spm_inv}. Includes a fast path for diagonal
#' \code{Matrix} objects.
#'
#' @param A Square numeric matrix (dense or \code{Matrix}).
#' @param TOL Optional tolerance. If \code{NULL}, set automatically.
#' @return Inverse matrix.
#' @examples
#' A <- matrix(c(2, 1, 1, 2), 2, 2)
#' dcm_inv(A)
#'
#' # Unlike solve(), a singular matrix is regularised rather than an error
#' singular <- matrix(1, 2, 2)
#' dcm_inv(singular)
#' @keywords internal
#' @export
dcm_inv <- function(A, TOL = NULL) {
  m <- nrow(A); n <- ncol(A)
  if (length(A) == 0) {
    return(if (.has_Matrix) .spzero(n, m) else matrix(0, nrow = n, ncol = m))
  }
  # Fast path: diagonal Matrix stays in Matrix world.
  if (.has_Matrix && inherits(A, "diagonalMatrix") && m == n) {
    d <- Matrix::diag(A)
    if (is.null(TOL)) TOL <- max(.Machine$double.eps * max(abs(d)) * m, exp(-32))
    return(.Diag(n, x = 1 / (d + TOL)))
  }
  # Defensive coercion: from here on, work in base R so we don't depend on
  # Matrix S4 dispatch being registered for the session.
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  m <- nrow(A); n <- ncol(A)
  if (is.null(TOL)) {
    TOL <- max(.Machine$double.eps * norm(A, "I") * max(m, n), exp(-32))
  }
  solve(A + diag(TOL, m, n))
}

#' Pseudo-inverse via SVD with truncation
#'
#' @param A Numeric matrix.
#' @param TOL Singular-value tolerance.
#' @return Pseudo-inverse of \code{A}.
#' @examples
#' # Left inverse of a tall matrix
#' A <- matrix(c(1, 2, 3, 4, 5, 7), nrow = 3)
#' round(dcm_pinv(A) %*% A, 8)   # ~ 2 x 2 identity
#' @keywords internal
#' @export
dcm_pinv <- function(A, TOL = NULL) {
  if (inherits(A, "Matrix")) A <- as.matrix(A)
  m <- nrow(A); n <- ncol(A)
  if (!length(A)) return(if (.has_Matrix) .spzero(n, m) else matrix(0, n, m))
  if (is.null(TOL)) {
    X <- suppressWarnings(dcm_inv(t(A) %*% A))
    if (all(is.finite(X))) return(X %*% t(A))
  }
  svd_result <- dcm_svd(A, 0)
  U <- as.matrix(svd_result$U); S <- svd_result$S; V <- as.matrix(svd_result$V)
  S <- as.vector(Matrix::diag(S))
  if (is.null(TOL)) TOL <- max(m, n) * .Machine$double.eps * max(S)
  r <- sum(abs(S) > TOL)
  if (r == 0) return(if (.has_Matrix) .spzero(n, m) else matrix(0, n, m))
  i <- seq_len(r)
  S_inv <- diag(1 / S[i], r, r)
  V[, i, drop = FALSE] %*% S_inv %*% t(U[, i, drop = FALSE])
}

#' Log-determinant of a (semi-)definite matrix
#'
#' Robust log-determinant suitable for sparse, dense, or rank-deficient
#' covariance matrices. Mirrors SPM25's \code{spm_logdet}, with fast paths
#' for 1x1 and diagonal inputs.
#'
#' @param C Square matrix.
#' @return Numeric log-determinant (or \code{NaN} if non-positive).
#' @examples
#' C <- diag(c(1, 2, 4))
#' dcm_logdet(C)          # log(1 * 2 * 4)
#' log(prod(c(1, 2, 4)))  # same
#'
#' # Zero rows/columns are dropped rather than sending the result to -Inf
#' dcm_logdet(diag(c(1, 2, 0)))
#' @keywords internal
#' @export
dcm_logdet <- function(C) {
  if (is.null(dim(C))) stop("dcm_logdet: C must be a matrix, not a vector")

  if (prod(dim(C)) == 1L) {
    v <- as.numeric(C)
    if (v <= 0) return(NaN)
    return(log(v))
  }

  if (inherits(C, "diagonalMatrix")) {
    d <- Matrix::diag(C)
    d <- d[d != 0]
    if (!length(d)) return(0)
    if (any(d <= 0 | is.nan(d))) return(NaN)
    return(sum(log(d)))
  }

  diagv <- if (inherits(C, "Matrix")) Matrix::diag(C) else diag(C)
  idx   <- which(diagv != 0)
  if (!length(idx)) return(0)
  C <- C[idx, idx, drop = FALSE]

  if (inherits(C, "Matrix")) {
    if (any(is.nan(Matrix::summary(C)$x))) return(NaN)
  } else if (any(is.nan(C))) return(NaN)

  TOL <- 1e-16
  svd_logdet <- function(M) {
    d <- svd(as.matrix(M), nu = 0, nv = 0)$d
    sum(log(d[d > TOL & d < 1/TOL]))
  }
  is_sym <- function(M) {
    Nm <- if (inherits(M, "Matrix")) Matrix::norm(M - Matrix::t(M), "I")
          else norm(M - t(M), "I")
    !is.na(Nm) && Nm <= TOL
  }

  if (!is_sym(C)) return(svd_logdet(C))

  if (inherits(C, "Matrix")) {
    Csym <- Matrix::forceSymmetric(C, uplo = "U")
    chol_ok <- TRUE
    tryCatch(Matrix::Cholesky(Csym, LDL = FALSE, perm = TRUE),
             error = function(e) { chol_ok <<- FALSE })
    if (chol_ok) {
      detv <- Matrix::determinant(Csym, logarithm = TRUE)
      if (detv$sign <= 0) return(NaN)
      return(as.numeric(detv$modulus))
    }
    return(svd_logdet(Csym))
  } else {
    Csym <- (C + t(C)) / 2
    R <- tryCatch(chol(Csym), error = function(e) NULL)
    if (!is.null(R)) return(2 * sum(log(diag(R))))
    return(svd_logdet(Csym))
  }
}

#' Truncated SVD with sparsity awareness
#'
#' Computes a thin SVD and drops singular components below a relative
#' tolerance \code{U}. Mirrors SPM25's \code{spm_svd}.
#'
#' @param X Numeric matrix.
#' @param U Relative tolerance for retaining singular values.
#' @return List with \code{U}, \code{S}, \code{V}.
#' @keywords internal
#' @export
dcm_svd <- function(X, U = NULL) {
  if (is.null(U)) U <- 1e-6
  if (U >= 1) U <- U - 1e-6
  if (U <= 0) U <- 64 * .Machine$double.eps

  M <- nrow(X); N <- ncol(X)
  p <- which(apply(X, 1, function(row) any(row != 0)))
  q <- which(apply(X, 2, function(col) any(col != 0)))
  X <- X[p, q, drop = FALSE]

  if (inherits(X, "sparseMatrix")) {
    temp <- Matrix::summary(X); i <- temp$i; j <- temp$j; s <- temp$x
  } else {
    indices <- which(X != 0, arr.ind = TRUE)
    if (!length(indices)) {
      i <- integer(0); j <- integer(0); s <- numeric(0)
    } else {
      i <- indices[, 1]; j <- indices[, 2]; s <- X[indices]
    }
  }

  m <- nrow(X); n <- ncol(X)
  if (any(i - j != 0)) {
    Xd <- as.matrix(X)
    if (m > n) {
      sv <- svd(t(Xd) %*% Xd)
      v <- sv$u; S <- diag(sv$d); sv_d <- diag(S)
      jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      v <- v[, jj, drop = FALSE]; u <- dcm_en(Xd %*% v)
      S <- sqrt(S[jj, jj, drop = FALSE])
    } else if (m < n) {
      sv <- svd(Xd %*% t(Xd))
      u <- sv$u; S <- diag(sv$d); sv_d <- diag(S)
      jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      u <- u[, jj, drop = FALSE]; v <- dcm_en(t(Xd) %*% u)
      S <- sqrt(S[jj, jj, drop = FALSE])
    } else {
      sv <- svd(Xd); u <- sv$u; S <- diag(sv$d); v <- sv$v
      sv_d <- diag(S)^2
      jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
      v <- v[, jj, drop = FALSE]; u <- u[, jj, drop = FALSE]
      S <- S[jj, jj, drop = FALSE]
    }
  } else {
    S <- .spmat(i = seq_len(n), j = seq_len(n), x = s,
                              dims = c(m, n))
    u <- .Diag(m)[, seq_len(n), drop = FALSE]
    v <- .Diag(m)[, seq_len(n), drop = FALSE]
    ord <- sort(-s, index.return = TRUE)$ix
    S <- S[ord, ord, drop = FALSE]
    v <- v[, ord, drop = FALSE]
    u <- u[, ord, drop = FALSE]
    sv_d <- Matrix::diag(S)^2
    jj <- which(sv_d * length(sv_d) / sum(sv_d) > U)
    v <- v[, jj, drop = FALSE]; u <- u[, jj, drop = FALSE]
    S <- S[jj, jj, drop = FALSE]
  }

  jlen <- ncol(S)
  U_out <- .spzero(M, jlen)
  V_out <- .spzero(N, jlen)
  if (jlen > 0) { U_out[p, ] <- u; V_out[q, ] <- v }
  list(U = U_out, S = S, V = V_out)
}

#' Euclidean column normalization
#'
#' @param X Matrix.
#' @param p Optional polynomial detrend order applied first.
#' @return \code{X} with each non-zero column scaled to unit Euclidean norm.
#' @examples
#' X <- matrix(c(3, 4, 0, 0, 5, 12), nrow = 2)
#' dcm_en(X)
#' sqrt(colSums(dcm_en(X)^2))  # 1, 0, 1 (the all-zero column is left alone)
#' @keywords internal
#' @export
dcm_en <- function(X, p = NULL) {
  if (!is.null(p)) X <- dcm_detrend(X, p)
  for (i in seq_len(ncol(X)))
    if (any(X[, i])) X[, i] <- X[, i] / sqrt(sum(X[, i]^2))
  X
}

#' Polynomial detrend (column-wise)
#'
#' @param x Matrix or vector.
#' @param p Polynomial order. \code{p = 0} just centres the columns.
#' @return Detrended matrix of the same shape.
#' @examples
#' # p = 0 centres each column
#' x <- cbind(1:10, (1:10) * 2 + 5)
#' round(colMeans(dcm_detrend(x, 0)), 10)
#'
#' # p = 1 removes a linear trend, so perfectly linear columns go to ~0
#' round(dcm_detrend(x, 1), 8)
#' @keywords internal
#' @export
dcm_detrend <- function(x, p = 0) {
  if (is.list(x)) {
    y <- x; for (i in seq_along(x)) y[[i]] <- dcm_detrend(x[[i]], p)
    return(y)
  }
  if (is.data.frame(x)) x <- as.matrix(x)
  dim_x <- dim(x)
  if (is.null(dim_x)) {
    m <- length(x); n <- 1; x <- matrix(x, ncol = 1)
  } else {
    m <- dim_x[1]; n <- dim_x[2]
  }
  if (m == 0 || n == 0) return(matrix(0, m, n))
  if (p == 0) return(x - matrix(colMeans(x), nrow = m, ncol = n, byrow = TRUE))

  G <- matrix(0, m, p + 1)
  for (i in 0:p) G[, i + 1] <- (seq_len(m))^i

  if (p == 1) {
    coef <- stats::lm.fit(G, x)$coefficients
    return(x - G %*% coef)
  }
  x - G %*% (MASS::ginv(G) %*% x)
}

#' Sparse identity-like matrix
#'
#' Build an identity-like sparse matrix of any shape, with optional shifted
#' diagonals. Mirrors SPM25's \code{spm_speye}.
#'
#' @param m,n Output dimensions.
#' @param k Diagonal offset (\code{0} for main).
#' @param c Cyclic-fill flag (0/1/2) as in SPM25.
#' @return A sparse \code{Matrix}.
#' @examples
#' dcm_speye(3)           # 3 x 3 sparse identity
#' dcm_speye(3, 3, 1)     # ones on the first super-diagonal
#' dcm_speye(2, 4)        # non-square is fine
#' @keywords internal
#' @export
dcm_speye <- function(m, n = m, k = 0, c = 0) {
  D <- Matrix::bandSparse(m, n, k, diagonals = list(rep(1, m)))
  if (c == 1) {
    if (k < 0) D <- D + dcm_speye(m, n, min(n, m) + k)
    else if (k > 0) D <- D + dcm_speye(m, n, k - min(n, m))
  } else if (c == 2) {
    col_sums <- Matrix::colSums(D); i <- which(col_sums == 0)
    if (length(i) > 0) {
      D <- D + .spmat(i = i, j = i, x = 1, dims = c(n, m))
    }
  }
  D
}

#' Concatenate a list of matrix blocks
#'
#' @param x A list (or matrix of lists) of numeric blocks.
#' @param d Optional dimension along which to concatenate.
#' @return A single matrix.
#' @keywords internal
#' @export
dcm_cat <- function(x, d = NULL) {
  if (!is.list(x)) return(x)

  if (!is.null(dim(x)) && !is.null(d)) {
    dims <- dim(x)
    if (d == 1) {
      y <- vector("list", dims[2])
      for (i in seq_len(dims[2])) y[[i]] <- dcm_cat(x[, i, drop = FALSE])
    } else if (d == 2) {
      y <- vector("list", dims[1])
      for (i in seq_len(dims[1])) y[[i]] <- dcm_cat(x[i, , drop = FALSE])
    } else stop("unknown option")
    return(y)
  }

  if (is.null(dim(x)) && all(vapply(x, is.numeric, TRUE))) {
    mats <- lapply(x, function(e) if (is.vector(e)) matrix(e, ncol = 1) else e)
    return(do.call(cbind, mats))
  }

  if (is.null(dim(x))) x <- matrix(x, nrow = 1)
  dims <- dim(x); n <- dims[1]; m <- dims[2]

  I <- matrix(0, n, m); J <- matrix(0, n, m)
  for (i in seq_len(n)) for (j in seq_len(m)) {
    elem <- x[[i, j]]
    if (is.list(elem)) { elem <- dcm_cat(elem); x[[i, j]] <- elem }
    if (!is.null(elem) && length(elem) > 0) {
      de <- dim(elem)
      if (is.null(de)) {
        I[i, j] <- if (length(elem) == 1) 1 else length(elem); J[i, j] <- 1
      } else { I[i, j] <- de[1]; J[i, j] <- de[2] }
    }
  }
  I_max <- apply(I, 1, max); J_max <- apply(J, 2, max)

  for (i in seq_len(n)) for (j in seq_len(m)) {
    elem <- x[[i, j]]
    if (is.null(elem) || (is.numeric(elem) && !length(elem))) {
      x[[i, j]] <- matrix(0, I_max[i], J_max[j])
    } else if (is.numeric(elem) && length(elem) == 1) {
      x[[i, j]] <- matrix(elem, I_max[i], J_max[j])
    } else if (is.vector(elem) && length(elem) > 1) {
      x[[i, j]] <- if (I[i, j] == I_max[i]) matrix(elem, ncol = 1)
                   else matrix(elem, nrow = 1)
    } else if (is.matrix(elem)) {
      cd <- dim(elem)
      if (cd[1] < I_max[i] || cd[2] < J_max[j]) {
        newm <- matrix(0, I_max[i], J_max[j])
        newm[seq_len(cd[1]), seq_len(cd[2])] <- elem
        x[[i, j]] <- newm
      }
    }
  }

  row_mats <- vector("list", n)
  for (i in seq_len(n)) row_mats[[i]] <- do.call(cbind, as.list(x[i, ]))
  do.call(rbind, row_mats)
}

#' Block-diagonal layout from a list
#'
#' @param ... A list and an optional offset \code{K}.
#' @return A list (or matrix) with the inputs on the diagonal.
#' @keywords internal
#' @export
dcm_diag <- function(...) {
  args <- list(...)
  if (!length(args)) stop("Not enough input arguments.")
  X <- args[[1]]
  if (!is.list(X)) return(diag(...))

  K <- if (length(args) < 2) 0 else args[[2]]
  if (is.null(dim(X))) {
    size <- length(X) + abs(K)
    D <- vector("list", size * size); dim(D) <- c(size, size)
    diag_idx <- which(diag(1, length(X)) == 1, arr.ind = TRUE)
    if (K != 0) {
      if (K > 0) diag_idx[, 2] <- diag_idx[, 2] + K
      else       diag_idx[, 1] <- diag_idx[, 1] - K
    }
    valid <- diag_idx[, 1] <= size & diag_idx[, 2] <= size
    diag_idx <- diag_idx[valid, , drop = FALSE]
    for (i in seq_len(nrow(diag_idx)))
      D[[diag_idx[i, 1], diag_idx[i, 2]]] <- X[[i]]
    return(D)
  } else {
    m <- nrow(X); n <- ncol(X)
    if (K != 0) {
      tmp <- matrix(0, max(m, n), max(m, n))
      if (K > 0) for (i in seq_len(n - K)) tmp[i, i + K] <- 1
      else       for (i in seq_len(m + K)) tmp[i - K, i] <- 1
      diag_idx <- which(tmp == 1, arr.ind = TRUE)
    } else {
      diag_idx <- which(diag(1, max(m, n)) == 1, arr.ind = TRUE)
    }
    D <- list()
    for (i in seq_len(nrow(diag_idx))) {
      r <- diag_idx[i, 1]; cc <- diag_idx[i, 2]
      if (r <= m && cc <= n) D[[length(D) + 1]] <- X[[r, cc]]
    }
    return(D)
  }
}

#' Stack a list of equally-shaped block matrices into a column-vectorized matrix
#'
#' @param blocks List of numeric matrices, all the same shape.
#' @param row_major If \code{TRUE}, vectorize by rows.
#' @return A matrix whose columns are the vectorized blocks.
#' @keywords internal
#' @export
dcm_blocks_to_mat <- function(blocks, row_major = FALSE) {
  stopifnot(is.list(blocks), length(blocks) > 0)
  B1 <- blocks[[1]]
  if (!is.matrix(B1)) stop("All blocks must be matrices.")
  rb <- nrow(B1); cb <- ncol(B1)
  for (k in seq_along(blocks)) {
    bk <- blocks[[k]]
    if (!is.matrix(bk) || nrow(bk) != rb || ncol(bk) != cb) {
      stop(sprintf("Block %d has shape %s; expected %dx%d.", k,
                   paste0(dim(bk), collapse = "x"), rb, cb))
    }
  }
  n_cols <- length(blocks); n_rows <- rb * cb
  out <- matrix(0, nrow = n_rows, ncol = n_cols)
  if (!row_major) for (k in seq_len(n_cols)) out[, k] <- as.vector(blocks[[k]])
  else            for (k in seq_len(n_cols)) out[, k] <- as.vector(t(blocks[[k]]))
  out
}

#' Trace of a matrix product
#'
#' Computes \code{trace(A \%*\% B)} efficiently as \code{sum(t(A) * B)}.
#'
#' @param A,B Numeric matrices.
#' @return Numeric scalar.
#' @keywords internal
#' @export
dcm_trace <- function(A, B) sum(t(A) * B)

# ---- utils-misc ----

# Miscellaneous utilities (internal).

#' Coerce a function-like object to a function
#'
#' Accepts a function, a function name (as character), or a one-line body and
#' returns a real \code{function}. Mirrors SPM25's \code{spm_funcheck}.
#'
#' @param f A function or character string.
#' @return A \code{function}.
#' @keywords internal
#' @export
dcm_funcheck <- function(f) {
  if (is.function(f) || is.null(f)) return(f)
  if (is.character(f)) {
    if (exists(f, mode = "function")) return(get(f, mode = "function"))
    tryCatch(eval(parse(text = paste0("function(...) { ", f, " }"))),
             error = function(e) stop(paste("Cannot convert", f, "to a function")))
  } else stop("Input must be a function or character string")
}

#' Compute a stable identifier for a numeric structure
#'
#' Used by SPM25 to fingerprint a dataset; returned in the DCM result.
#'
#' @param ... Numeric data (any shape, possibly nested).
#' @return Numeric scalar.
#' @keywords internal
#' @export
dcm_data_id <- function(...) {
  X <- list(...); if (length(X) == 1) X <- X[[1]]
  ID <- 0
  if (is.character(X)) { ID <- 0
  } else if (is.numeric(X)) {
    Y <- as.numeric(X)
    ID <- sum(abs(Y[!is.nan(Y) & !is.infinite(Y)]))
  } else if (is.list(X)) {
    for (elem in X) ID <- ID + dcm_data_id(elem)
  }
  if (ID > 0) ID <- 10^(-(floor(log10(ID)) - 2)) * ID
  ID
}

#' Univariate Normal CDF (SPM-style)
#'
#' Vectorized Normal CDF with \code{NaN} for non-positive variances.
#'
#' @param x Quantile.
#' @param u Mean.
#' @param v Variance.
#' @return Numeric vector of CDF values.
#' @examples
#' dcm_Ncdf(0)                     # 0.5
#' round(dcm_Ncdf(c(-1.96, 0, 1.96)), 4)
#'
#' # Note the third argument is a VARIANCE, not a standard deviation
#' dcm_Ncdf(2, u = 0, v = 4)       # == pnorm(2, mean = 0, sd = 2)
#' @keywords internal
#' @export
dcm_Ncdf <- function(x, u = 0, v = 1) {
  if (missing(x)) return(numeric(0))

  lx <- length(x); lu <- length(u); lv <- length(v)
  target_len <- max(lx, lu, lv)

  if (lx < target_len) x <- rep_len(x, target_len)
  if (lu < target_len) u <- rep_len(u, target_len)
  if (lv < target_len) v <- rep_len(v, target_len)

  F  <- numeric(target_len)
  md <- v > 0
  if (any(!md)) {
    F[!md] <- NaN
    warning("Returning NaN for variance <= 0.")
  }
  if (any(md)) F[md] <- stats::pnorm(x[md], mean = u[md], sd = sqrt(v[md]))
  F
}

#' Locate or describe entries in a structured parameter
#'
#' Given a structured parameter object \code{X} (with named fields), returns
#' the linear indices of one or more named fields, or the structural location
#' of a given linear index.
#'
#' @param X A structured parameter (list with numeric fields).
#' @param ... One or more field names or numeric indices.
#' @return Integer indices or a character description.
#' @keywords internal
#' @export
dcm_fieldindices <- function(X, ...) {
  varargin <- list(...)
  if (length(varargin) > 0 && is.numeric(varargin[[1]]) &&
      length(varargin[[1]]) > 1) {
    return(lapply(varargin[[1]], function(j) dcm_fieldindices(X, j)))
  }
  X0 <- dcm_zeros(X); ix <- dcm_vec(X0)
  for (i in seq_along(varargin)) {
    if (is.character(varargin[[i]])) {
      field <- varargin[[i]]; x <- X0
      if (field %in% names(x)) {
        f <- x[[field]]; f <- dcm_unvec(dcm_vec(f) + 1, f); x[[field]] <- f
        ix <- ix + dcm_vec(x)
      } else warning(paste("Field", field, "not found"))
    } else if (is.numeric(varargin[[i]])) {
      ind <- varargin[[i]]; x <- ix; x[ind] <- 1; x <- dcm_unvec(x, X)
      for (nm in names(x)) {
        if (any(dcm_vec(x[[nm]]) != 0)) {
          v <- x[[nm]]
          if (is.numeric(v) || is.matrix(v)) {
            idx <- which(v != 0, arr.ind = TRUE)
            if (is.null(dim(v)) || min(dim(v)) == 1) {
              return(sprintf("%s[%i]", nm, which(v != 0)[1]))
            } else if (length(dim(v)) < 3) {
              return(sprintf("%s[%i,%i]", nm, idx[1, 1], idx[1, 2]))
            } else {
              return(sprintf("%s[%i,%i,%i]", nm, idx[1, 1], idx[1, 2], idx[1, 3]))
            }
          }
        }
      }
    }
  }
  which(ix != 0)
}

#' Find indices of free parameters in a DCM
#'
#' @param ... Either \code{(DCM, fields)} or \code{(pC, pE, fields)}.
#' @return List with \code{i} (free indices), plus \code{pC}, \code{pE}, \code{Np}.
#' @keywords internal
#' @export
dcm_find_pC <- function(...) {
  varargin <- list(...)
  pC <- NULL; pE <- NULL; fields <- NULL
  if (length(varargin) > 2) {
    pC <- varargin[[1]]; pE <- varargin[[2]]; fields <- varargin[[3]]
  } else if (length(varargin) > 1) {
    DCM <- varargin[[1]]; fields <- varargin[[2]]
  } else { DCM <- varargin[[1]]; fields <- NULL }

  if (is.null(pC) || is.null(pE)) {
    if (is.character(DCM)) DCM <- readRDS(DCM)
    if (any(c("options", "M") %in% names(DCM))) {
      tryCatch({
        result <- dcm_find_rC(DCM); pC <- result$pC; pE <- result$pE
      }, error = function(e) {})
    }
  }

  q <- if (is.list(pC) && !is.matrix(pC)) dcm_vec(pC) else diag(pC)
  threshold <- mean(q[q < 1024]) / 1024
  i <- which(q > threshold); Np <- length(q)

  if (!is.null(fields)) {
    if (is.character(fields) && length(fields) == 1) fields <- list(fields)
    if (is.list(pE)) {
      j <- do.call(dcm_fieldindices, c(list(pE), fields))
      if (!length(j) && !(length(fields) > 0 && fields[[1]] == "none")) {
        warning(sprintf("%s not found. Returning all fields",
                        paste(unlist(fields), collapse = ",")))
      } else i <- j[j %in% i]
    }
  }
  list(i = i, pC = pC, pE = pE, Np = Np)
}

#' Recover priors from a DCM container
#'
#' @param DCM A DCM list.
#' @return List with \code{pC}, \code{pE}.
#' @keywords internal
#' @export
dcm_find_rC <- function(DCM) {
  if (!is.null(DCM$M$pC) && !is.null(DCM$M$pE)) {
    return(list(pC = DCM$M$pC, pE = DCM$M$pE))
  }
  if (is.null(DCM$d)) DCM$d <- array(0, dim = c(nrow(DCM$a), nrow(DCM$a), 0))
  result <- dcm_fmri_priors(DCM$a, DCM$b, DCM$c, DCM$d, DCM$options)
  list(pC = result$pC, pE = result$pE)
}

# ---- cov-basis ----

# Error covariance basis functions (internal): the Q matrices for dcm_nlsi_GN.

# Internal: sparse identity of size n
.speye_Ce <- function(n) .Diag(n)

# Internal: place block M into an n_total x n_total sparse matrix at offset
# k_offset on both axes.
.sparse_offset <- function(M, n_total, k_offset) {
  if (inherits(M, "diagonalMatrix")) {
    d <- Matrix::diag(M); nz <- which(d != 0)
    if (length(nz) == 0L) {
      return(.spmat(i = integer(0), j = integer(0),
                                  x = numeric(0),
                                  dims = c(n_total, n_total)))
    }
    return(.spmat(i = nz + k_offset, j = nz + k_offset,
                                x = d[nz], dims = c(n_total, n_total)))
  }
  if (inherits(M, "sparseMatrix")) {
    Mt <- methods::as(M, "TsparseMatrix")
    i <- Mt@i + 1L; j <- Mt@j + 1L; x <- Mt@x
    if (!length(x)) {
      return(.spmat(i = integer(0), j = integer(0),
                                  x = numeric(0),
                                  dims = c(n_total, n_total)))
    }
    return(.spmat(i = i + k_offset, j = j + k_offset,
                                x = x, dims = c(n_total, n_total)))
  }
  if (is.matrix(M)) {
    Ms <- Matrix::Matrix(M, sparse = TRUE)
    Mt <- methods::as(Ms, "TsparseMatrix")
    i <- Mt@i + 1L; j <- Mt@j + 1L; x <- Mt@x
    if (!length(x)) {
      return(.spmat(i = integer(0), j = integer(0),
                                  x = numeric(0),
                                  dims = c(n_total, n_total)))
    }
    return(.spmat(i = i + k_offset, j = j + k_offset,
                                x = x, dims = c(n_total, n_total)))
  }
  stop("Unsupported matrix type: cannot offset block.")
}

#' Error-covariance basis (AR or FAST)
#'
#' Construct a list of covariance basis matrices for the noise model in
#' variational inversion. Mirrors SPM25's \code{spm_Ce}.
#'
#' @param t Either a numeric vector of session lengths (defaults to AR basis)
#'   or a string \code{"ar"} / \code{"fast"} selecting the basis type.
#' @param v Session lengths (when \code{t} is a string).
#' @param a AR coefficient (or sampling interval for FAST).
#' @return A list of sparse \code{Matrix} objects.
#' @keywords internal
#' @export
dcm_Ce <- function(t, v = NULL, a = NULL) {
  if (!is.character(t)) {
    if (!missing(v)) a <- v else a <- NULL
    v <- t; t <- "ar"
  } else {
    if (missing(a)) a <- NULL
  }
  t <- tolower(t)

  if (t == "ar") {
    C <- list(); l <- length(v); n <- sum(v); k <- 0L
    if (l > 1L) {
      for (i in seq_len(l)) {
        dCda <- dcm_Ce(v[i], a)
        for (j in seq_along(dCda))
          C[[length(C) + 1L]] <- .sparse_offset(dCda[[j]], n_total = n,
                                                k_offset = k)
        k <- k + v[i]
      }
      return(C)
    }
    if (!is.null(a) && length(a) > 0L) {
      Q    <- dcm_Q(a, v)
      dQda <- dcm_diff("dcm_Q", a, v, 1)
      dQ   <- if (is.list(dQda$J)) dQda$J[[1L]] else dQda$J
      C[[1L]] <- Matrix::Matrix(Q - dQ * a, sparse = TRUE)
      C[[2L]] <- Matrix::Matrix(Q + dQ * a, sparse = TRUE)
    } else {
      C[[1L]] <- .speye_Ce(v)
    }
    return(C)
  }

  if (t == "fast") {
    dt <- a; C <- list(); n <- sum(v); k <- 0L
    for (m in seq_along(v)) {
      T <- (0:(v[m] - 1)) * dt
      if (!is.numeric(dt) || length(dt) != 1L || !is.finite(dt) || dt <= 0) {
        d <- numeric(0)
      } else {
        start_exp <- floor(log2(dt / 4)); end_exp <- log2(64)
        d <- if (is.na(start_exp) || start_exp > end_exp) numeric(0)
             else 2^seq(from = start_exp, to = end_exp, by = 1)
      }
      if (length(d) > 0L) {
        for (i in seq_len(min(6L, length(d)))) {
          for (j in 0:2) {
            vtoe <- (T^j) * exp(-T / d[i])
            C[[length(C) + 1L]] <- .sparse_offset(stats::toeplitz(vtoe),
                                                  n_total = n, k_offset = k)
          }
        }
      }
      k <- k + v[m]
    }
    return(C)
  }

  stop("Unknown error covariance constraints.")
}

#' AR(p) autocorrelation / precision matrix
#'
#' Returns either a banded precision (for \code{q != 0}) or a Toeplitz
#' autocorrelation derived from the AR coefficients \code{a}.
#'
#' @param a AR coefficient vector.
#' @param n Matrix dimension.
#' @param q Output flag: \code{0} for autocorrelation, otherwise precision.
#' @return Sparse \code{Matrix} or numeric matrix.
#' @keywords internal
#' @export
dcm_Q <- function(a, n, q = 0) {
  if (q != 0) {
    A <- c(-a[1], (1 + a[1]^2), -a[1])
    return(Matrix::bandSparse(n, n, (-1):1,
                              list(rep(A[1], n), rep(A[2], n), rep(A[3], n))))
  }
  p <- length(a); A <- c(1, -a)
  diag_list <- lapply(0:p, function(k) rep(A[k + 1], n))
  P <- Matrix::bandSparse(n, n, -(0:p), diag_list)
  K <- solve(P)
  K <- K * (abs(K) > 1e-4)
  Q <- K %*% t(K)
  stats::toeplitz(Q[, 1])
}

# ---- diff ----

# Numerical differentiation (internal). dcm_diff is a hot-path function:
# called once per parameter direction in every Gauss-Newton iteration.

#' High-order numerical Jacobian
#'
#' Forward-difference numerical Jacobian of a (possibly nested-output)
#' function with respect to one or more of its arguments. Mirrors SPM25's
#' \code{dcm_diff}.
#'
#' @param ... Function, its arguments, the index (or vector of indices) of
#'   the argument(s) to differentiate with respect to, and optionally a list
#'   of projection matrices.
#' @return A list containing the Jacobian \code{J} (or higher-order
#'   derivatives) and the function value \code{f0}.
#' @keywords internal
#' @export
dcm_diff <- function(...) {
  args <- list(...)
  if (length(args) < 2) stop("Improper call: need at least f and one x argument")

  dx <- tryCatch({
    if (exists("GLOBAL_DX", envir = .rsDCM_env, inherits = FALSE) &&
        !is.null(get("GLOBAL_DX", envir = .rsDCM_env)))
      get("GLOBAL_DX", envir = .rsDCM_env)
    else exp(-8)
  }, error = function(e) exp(-8))

  f <- args[[1]]
  if (is.character(f)) {
    if (!exists(f, mode = "function")) stop("Function name not found: ", f)
    f <- get(f, mode = "function")
  } else if (!is.function(f)) {
    stop("First argument must be a function or name of a function")
  }
  if (!is.function(f)) f <- dcm_funcheck(f)

  q_flag <- TRUE; V_list <- NULL; n <- NULL
  last <- args[[length(args)]]

  if (is.list(last)) {
    V_list <- last; n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]; q_flag <- TRUE
  } else if (is.character(last)) {
    n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]
    V_list <- vector("list", length(xs)); q_flag <- FALSE
  } else if (is.numeric(last)) {
    n <- as.integer(last); xs <- args[2:(length(args) - 1)]
    V_list <- vector("list", length(xs)); q_flag <- TRUE
  } else stop("Improper call.")

  if (any(n < 1L) || any(n > length(xs))) stop("n contains invalid argument indices")
  if (is.null(V_list)) V_list <- vector("list", length(xs))
  if (length(V_list) < length(xs))
    V_list <- c(V_list, vector("list", length(xs) - length(V_list)))

  for (i in seq_along(xs)) {
    if (is.null(V_list[[i]]) && any(n == i))
      V_list[[i]] <- diag(1, dcm_length(xs[[i]]))
  }

  m  <- n[length(n)]
  xm <- dcm_vec(xs[[m]])
  Vm_obj <- V_list[[m]]
  ncols_m <- if (!is.null(Vm_obj)) ncol(as.matrix(Vm_obj)) else 0L
  Jcols <- vector("list", ncols_m)

  V_is_identity <- FALSE
  if (ncols_m > 0L) {
    if (is.matrix(Vm_obj) && nrow(Vm_obj) == ncols_m && ncols_m == length(xm)) {
      V_is_identity <- isTRUE(all(diag(Vm_obj) == 1)) &&
                       isTRUE(sum(abs(Vm_obj)) == ncols_m)
    }
  }

  if (length(n) == 1L) {
    f0 <- do.call(f, xs)
    if (ncols_m > 0L) {
      if (V_is_identity) {
        for (i in seq_len(ncols_m)) {
          xm_p <- xm; xm_p[i] <- xm_p[i] + dx
          xi <- xs; xi[[m]] <- dcm_unvec(xm_p, xs[[m]])
          fi <- do.call(f, xi)
          Jcols[[i]] <- dcm_dfdx(fi, f0, dx)
        }
      } else {
        Vm <- as.matrix(Vm_obj)
        for (i in seq_len(ncols_m)) {
          xi <- xs
          xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
          fi <- do.call(f, xi)
          Jcols[[i]] <- dcm_dfdx(fi, f0, dx)
        }
      }
    }
    fv <- dcm_vec(f0)
    if (length(xm) == 0L) {
      J <- matrix(numeric(0), nrow = length(fv), ncol = 0)
    } else if (length(fv) == 0L) {
      J <- matrix(numeric(0), nrow = 0, ncol = length(xm))
    } else if (is.numeric(f0) && ncols_m > 0L && q_flag) {
      J <- dcm_dfdx_cat(Jcols)
    } else {
      J <- Jcols
    }
    return(list(J = J, f0 = f0))
  }

  if (length(n) == 2L) {
    d1 <- do.call(dcm_diff, c(list(f), xs, list(n[1], V_list)))
    if (ncols_m > 0L) {
      Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
      for (i in seq_len(ncols_m)) {
        xi <- xs
        xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
        fi <- do.call(dcm_diff, c(list(f), xi, list(n[1], V_list)))
        Jcols[[i]] <- dcm_dfdx(fi$J, d1$J, dx)
        p_numeric <- p_numeric && is.numeric(Jcols[[i]])
      }
      d2 <- if (p_numeric && q_flag) dcm_dfdx_cat(Jcols) else Jcols
    } else {
      d2 <- if (q_flag) matrix(0, nrow = length(dcm_vec(d1$f0)), ncol = 0)
            else list()
    }
    return(list(d2 = d2, d1 = d1$J, f0 = d1$f0))
  }

  lower <- do.call(dcm_diff, c(list(f), xs,
                               list(n[1:(length(n) - 1)], V_list)))
  if (ncols_m > 0L) {
    Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
    for (i in seq_len(ncols_m)) {
      xi <- xs
      xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
      fi <- do.call(dcm_diff, c(list(f), xi,
                                list(n[1:(length(n) - 1)], V_list)))
      Jcols[[i]] <- dcm_dfdx(fi$J, lower$J, dx)
      p_numeric <- p_numeric && is.numeric(Jcols[[i]])
    }
    J <- if (p_numeric && q_flag) dcm_dfdx_cat(Jcols) else Jcols
  } else {
    J <- if (q_flag) matrix(0, nrow = length(dcm_vec(lower$f0)), ncol = 0)
         else list()
  }
  list(J = J, f0 = lower$f0)
}

# Internal: forward-difference quotient that handles list-valued outputs.
dcm_dfdx <- function(f, f0, dx) {
  if (is.list(f)) {
    dfdx <- f
    for (i in seq_along(f)) dfdx[[i]] <- dcm_dfdx(f[[i]], f0[[i]], dx)
    return(dfdx)
  } else if (is.list(f0) && !is.data.frame(f0)) {
    return((dcm_vec(f) - dcm_vec(f0)) / dx)
  } else {
    return((f - f0) / dx)
  }
}

# Internal: concatenate per-column Jacobian pieces into a single matrix.
dcm_dfdx_cat <- function(J) {
  if (is.vector(J[[1]])) {
    if (is.null(dim(J[[1]])) || ncol(as.matrix(J[[1]])) == 1)
      return(dcm_cat(J))
    else
      return(t(dcm_cat(lapply(J, function(x) t(as.matrix(x))))))
  }
  J
}

# ---- expm ----

# Matrix exponential utilities and the local-linearization integration step
# dcm_dx (internal).

#' Local linearisation integration step
#'
#' Computes a one-step Bayesian update \code{dx = (expm(t*J) - I) * J^{-1} * f}
#' via an augmented matrix exponential. Mirrors SPM25's \code{spm_dx}.
#'
#' @param dfdx Jacobian \code{df/dx} of the model.
#' @param f Current residual / gradient.
#' @param t Step length (or list with a single regulariser).
#' @param Q Optional skew matrix.
#' @return Update \code{dx} matching the structure of \code{f}.
#' @keywords internal
#' @export
dcm_dx <- function(dfdx, f, t = Inf, Q = NULL) {
  nmax <- 512; xf <- f
  f <- dcm_vec(f); n <- length(f)

  if (is.list(t)) {
    tval <- t[[1]]
    if (length(tval) == 1L) {
      t <- suppressWarnings(exp(tval - dcm_logdet(dfdx) / n))
    } else {
      t <- suppressWarnings(exp(tval - log(diag(-dfdx))))
    }
  }

  if (!missing(Q) && !is.null(Q)) {
    L <- as.matrix(dfdx); L[upper.tri(L)] <- 0; Q <- L - t(L)
    Qn <- norm(Q, type = "2")
    Q <- if (is.finite(Qn) && Qn > 0) Q / Qn / 8 else Q * 0
    f <- f - Q %*% f; dfdx <- dfdx - Q %*% dfdx
  }

  if (all(t > exp(16))) {
    pinv_dfdx <- dcm_pinv(dfdx)
    pinv_dfdx <- if (inherits(pinv_dfdx, "Matrix")) as.matrix(pinv_dfdx)
                 else pinv_dfdx
    dx <- -pinv_dfdx %*% matrix(as.numeric(f), ncol = 1)
    return(dcm_unvec(as.numeric(dx), xf))
  }

  Tmat <- if (is.vector(t) && length(t) > 1L) diag(as.numeric(t))
          else as.numeric(t)[1L]

  if (length(Tmat) == 1L) {
    tf <- Tmat * as.matrix(f); tJ <- Tmat * as.matrix(dfdx)
  } else {
    tf <- Tmat %*% f; tJ <- Tmat %*% dfdx
  }

  J <- rbind(c(0, rep(0, n)), cbind(tf, tJ))

  if (n <= nmax) {
    E <- expm::expm(J)
    dx <- E[2:(n + 1), 1, drop = FALSE]
  } else {
    dx <- matrix(expv(1, J, c(1, rep(0, n)))[2:(n + 1)], ncol = 1)
  }
  dcm_unvec(as.numeric(dx), xf)
}

#' Krylov-based action of matrix exponential
#'
#' Computes \code{w = expm(t*A) \%*\% v} without forming \code{expm(t*A)}.
#' Used as a fallback in \code{\link{dcm_dx}} for very large systems.
#'
#' @param t Step length.
#' @param A Square matrix.
#' @param v Vector.
#' @param tol Tolerance.
#' @param m Krylov subspace dimension.
#' @return Vector \code{expm(t*A) \%*\% v}.
#' @keywords internal
#' @export
expv <- function(t, A, v, tol = 1.0e-7, m = NULL) {
  A <- as.matrix(A); n <- nrow(A); if (is.null(m)) m <- min(n, 30L)
  anorm <- norm(A, type = "I"); mxrej <- 10; btol <- 1.0e-7
  gamma <- 0.9; delta <- 1.2; mb <- m; t_out <- abs(t)
  nstep <- 0; t_new <- 0; t_now <- 0; s_error <- 0
  rndoff <- anorm * .Machine$double.eps
  k1 <- 2; xm <- 1 / m; normv <- norm(as.matrix(v), type = "2"); beta <- normv
  fact <- (((m + 1) / exp(1))^(m + 1)) * sqrt(2 * pi * (m + 1))
  t_new <- (1 / anorm) * ((fact * tol) / (4 * beta * anorm))^xm
  s <- 10^(floor(log10(t_new)) - 1); t_new <- ceiling(t_new / s) * s
  sgn <- sign(t); w <- v; hump <- normv

  while (t_now < t_out) {
    nstep <- nstep + 1; t_step <- min(t_out - t_now, t_new)
    V <- matrix(0, n, m + 1); H <- matrix(0, m + 2, m + 2)
    V[, 1] <- (1 / beta) * w; k1 <- 2
    for (j in seq_len(m)) {
      p <- A %*% V[, j]
      for (i in seq_len(j)) {
        H[i, j] <- crossprod(V[, i], p)[1, 1]
        p <- p - H[i, j] * V[, i]
      }
      s <- norm(as.matrix(p), type = "2")
      if (s < btol) {
        k1 <- 0; mb <- j; t_step <- t_out - t_now; break
      }
      H[j + 1, j] <- s; V[, j + 1] <- (1 / s) * p
    }
    if (k1 != 0) {
      H[m + 2, m + 1] <- 1; avnorm <- norm(A %*% V[, m + 1], type = "2")
    }
    ireject <- 0
    repeat {
      mx <- mb + k1; Hsub <- H[1:mx, 1:mx, drop = FALSE]
      F <- expm::expm(sgn * t_step * Hsub)
      if (k1 == 0) { err_loc <- btol; break }
      phi1 <- abs(beta * F[m + 1, 1])
      phi2 <- abs(beta * F[m + 2, 1] * avnorm)
      if      (phi1 > 10 * phi2) { err_loc <- phi2; xm <- 1 / m }
      else if (phi1 > phi2)      { err_loc <- (phi1 * phi2) / (phi1 - phi2); xm <- 1 / m }
      else                       { err_loc <- phi1; xm <- 1 / (m - 1) }
      if (err_loc <= delta * t_step * tol) break
      t_step <- gamma * t_step * (t_step * tol / err_loc)^xm
      s <- 10^(floor(log10(t_step)) - 1); t_step <- ceiling(t_step / s) * s
      if (ireject == mxrej) stop("The requested tolerance is too high.")
      ireject <- ireject + 1
    }
    mx <- mb + max(0, k1 - 1)
    w <- V[, 1:mx, drop = FALSE] %*% (beta * F[1:mx, 1, drop = FALSE])
    beta <- norm(as.matrix(w), type = "2"); hump <- max(hump, beta)
    t_now <- t_now + t_step
    t_new <- gamma * t_step * (t_step * tol / err_loc)^xm
    s <- 10^(floor(log10(t_new)) - 1); t_new <- ceiling(t_new / s) * s
    err_loc <- max(err_loc, rndoff); s_error <- s_error + err_loc
  }
  as.vector(w)
}

#' Pade approximation of the matrix exponential
#'
#' Pure-R Pade approximation, used as a fallback. See also
#' \code{\link[expm]{expm}}.
#'
#' @param A Square matrix.
#' @param p Order.
#' @return Matrix exponential of \code{A}.
#' @keywords internal
#' @export
padm <- function(A, p = 6) {
  A <- as.matrix(A); n <- nrow(A)
  c_vec <- numeric(p + 1); c_vec[1] <- 1
  if (p >= 1) {
    for (k in seq_len(p))
      c_vec[k + 1] <- c_vec[k] * ((p + 1 - k) / (k * (2 * p + 1 - k)))
  }
  s <- norm(A, type = "I")
  if (s > 0.5) {
    s <- max(0, floor(log(s) / log(2)) + 2); A <- A / (2^s)
  } else s <- 0
  I <- diag(1, n); A2 <- A %*% A
  Q <- c_vec[p + 1] * I; P <- c_vec[p] * I; odd <- TRUE
  for (k in (p - 1):1) {
    if (odd) Q <- Q %*% A2 + c_vec[k] * I
    else     P <- P %*% A2 + c_vec[k] * I
    odd <- !odd
  }
  if (odd) {
    Q <- Q %*% A; Q <- Q - P; E <- -(I + 2 * solve(Q, diag(n)))
  } else {
    P <- P %*% A; Q <- Q - P; E <- I + 2 * solve(Q, diag(n))
  }
  if (s > 0) for (k in seq_len(s)) E <- E %*% E
  E
}

#' Scaled-and-squared matrix exponential
#'
#' SPM25-style matrix exponential used as a fallback to \code{\link[expm]{expm}}.
#'
#' @param J Square matrix.
#' @param x Optional vector to multiply by \code{expm(J)} on the right.
#' @return Matrix or vector.
#' @examples
#' # expm of a diagonal matrix is just exp() of the diagonal
#' J <- diag(c(-1, -2))
#' round(dcm_expm(J), 8)
#' round(diag(exp(c(-1, -2))), 8)
#'
#' # Supplying x returns expm(J) %*% x without forming the product yourself
#' round(dcm_expm(J, c(1, 1)), 8)
#' @keywords internal
#' @export
dcm_expm <- function(J, x = NULL) {
  I <- .Diag(nrow(J))
  e <- floor(log2(norm(J, "I"))); s <- max(0, e + 1)
  J <- J / 2^s; X <- J; c_v <- 1 / 2; E <- I + c_v * J; D <- I - c_v * J
  q <- 6; p <- 1
  for (k in 2:q) {
    c_v <- c_v * (q - k + 1) / (k * (2 * q - k + 1))
    X <- J %*% X; cX <- c_v * X
    E <- E + cX
    if (p == 1) D <- D + cX else D <- D - cX
    p <- 1 - p
  }
  E <- solve(as.matrix(D), as.matrix(E))
  for (k in seq_len(s)) E <- E %*% E
  if (!is.null(x)) E %*% x else E
}

# ---- fmri-model ----

# fMRI generative model: state equation, observation equation, priors, modes.

#' Logistic sigmoid
#'
#' Standard SPM logistic function \code{1/(1 + exp(-x))}, vectorized.
#'
#' @param x Numeric input.
#' @return Numeric output of the same shape.
#' @keywords internal
#' @export
dcm_phi <- function(x) 1 / (1 + exp(-x))

#' DCM mode generator
#'
#' Maps a vector of mode parameters to a connectivity matrix. Mirrors
#' SPM25's \code{spm_dcm_fmri_mode_gen}.
#'
#' @param Ev Numeric vector of mode parameters.
#' @param modes Numeric matrix of mode columns.
#' @param Cv Optional covariance of \code{Ev}.
#' @return If \code{Cv} is \code{NULL}, the connectivity matrix \code{Ep}.
#'   Otherwise a list with \code{Ep} and propagated covariance \code{Cp}.
#' @keywords internal
#' @export
dcm_fmri_mode_gen <- function(Ev, modes, Cv = NULL) {
  Ev <- as.numeric(Ev); modes <- as.matrix(modes)
  Dv <- diag(-exp(-Ev), nrow = length(Ev))
  Ep <- modes %*% Dv %*% t(modes)
  if (is.null(Cv)) return(Ep)
  dAdv <- dcm_diff(dcm_fmri_mode_gen, Ev, modes, 1)
  Jpart <- if (!is.null(dAdv$J)) dAdv$J
           else if (!is.null(dAdv$d1)) dAdv$d1 else dAdv
  n <- nrow(Ep); p <- length(Ev); G <- matrix(0, n * n, p)
  if (is.list(Jpart)) {
    if (length(Jpart) != p) stop("Derivative list length mismatch.")
    for (i in seq_len(p)) {
      Gi <- as.matrix(Jpart[[i]]); G[, i] <- as.vector(Gi)
    }
  } else if (is.matrix(Jpart)) {
    G <- Jpart
  } else stop("Unrecognized derivative container from dcm_diff.")
  Cp <- G %*% as.matrix(Cv) %*% t(G)
  list(Ep = Ep, Cp = Cp)
}

#' fMRI neural and hemodynamic state equation
#'
#' Computes \code{dx/dt} for the bilinear neural state equation coupled to
#' the Buxton-Friston hemodynamic model. Mirrors SPM25's \code{spm_fx_fmri}.
#'
#' @param x State matrix (rows = regions, cols = state variables).
#' @param u Driving inputs at the current time.
#' @param P List of model parameters (\code{A}, \code{B}, \code{C}, \code{D},
#'   \code{transit}, \code{decay}, \code{epsilon}, ...).
#' @param M Model structure list (optional).
#' @return Matrix \code{dx/dt} of the same shape as \code{x}.
#' @export
dcm_fx_fmri <- function(x, u, P, M) {
  if (missing(M) || is.null(M)) M <- list()
  symmetry <- if (!is.null(M$symmetry)) M$symmetry else 0

  A <- P$A; B <- P$B; C <- P$C / 16; D <- P$D
  n <- nrow(x); mcols <- ncol(x); u <- as.numeric(u)
  f <- x

  if (mcols == 5) {
    if (is.vector(A) && length(A) > 1) {
      EE <- dcm_fmri_mode_gen(A, M$modes)
      if (length(dim(B)) == 3 && dim(B)[3] > 0) {
        for (i in seq_len(dim(B)[3])) EE <- EE + u[i] * B[, , i]
      }
      if (length(dim(D)) == 3 && dim(D)[3] > 0) {
        for (i in seq_len(dim(D)[3])) EE <- EE + x[i, 1] * D[, , i]
      }
    } else {
      Aeff <- A
      if (length(dim(B)) == 3 && dim(B)[3] > 0) {
        for (i in seq_len(dim(B)[3])) Aeff <- Aeff + u[i] * B[, , i]
      }
      if (length(dim(D)) == 3 && dim(D)[3] > 0) {
        for (i in seq_len(dim(D)[3])) Aeff <- Aeff + x[i, 1] * D[, , i]
      }
      if (length(dim(Aeff)) == 3 && dim(Aeff)[3] > 1) {
        Aeff <- exp(Aeff[, , 1]) - exp(Aeff[, , 2])
      }
      SE <- diag(Aeff)
      EE <- Aeff - diag(exp(SE) / 2 + SE)
      if (isTRUE(symmetry)) EE <- (EE + t(EE)) / 2
    }
    f[, 1] <- as.vector(EE %*% x[, 1] + C %*% u)
  } else {
    Aeff <- A
    if (length(dim(B)) == 3 && dim(B)[3] > 0) {
      for (i in seq_len(dim(B)[3])) Aeff[, , 1] <- Aeff[, , 1] + u[i] * B[, , i]
    }
    if (length(dim(D)) == 3 && dim(D)[3] > 0) {
      for (i in seq_len(dim(D)[3])) Aeff[, , 1] <- Aeff[, , 1] + x[i, 1] * D[, , i]
    }
    nreg <- n; EE <- exp(Aeff[, , 1]) / 8; IE <- diag(diag(EE)); EE <- EE - IE
    EI <- diag(nreg); SE <- diag(nreg) / 2; SI <- diag(nreg)
    if (length(dim(Aeff)) > 2 && dim(Aeff)[3] > 1) {
      phi <- dcm_phi(Aeff[, , 2] * 2)
      EI <- EI + EE * (1 - phi); EE <- EE * phi - SE
    } else EE <- EE - SE
    f[, 1] <- as.vector(EE %*% x[, 1] - IE %*% x[, 6] + C %*% u)
    f[, 6] <- as.vector(EI %*% x[, 1] - SI %*% x[, 6])
  }

  H <- c(0.64, 0.32, 2.00, 0.32, 0.4)

  z3 <- exp(x[, 3]); z4 <- exp(x[, 4]); z5 <- exp(x[, 5])

  sd <- H[1] * exp(P$decay); tt <- H[3] * exp(P$transit)
  fv <- z4^(1 / H[4])
  ff <- (1 - (1 - H[5])^(1 / z3)) / H[5]

  f[, 2] <- x[, 1] - as.numeric(sd) * x[, 2] - H[2] * (z3 - 1)
  f[, 3] <- x[, 2] / z3
  f[, 4] <- (z3 - fv) / (as.numeric(tt) * z4)
  f[, 5] <- (ff * z3 - fv * z5 / z4) / (as.numeric(tt) * z5)
  f
}

#' fMRI state equation with analytic Jacobians
#'
#' Same as \code{\link{dcm_fx_fmri}} but additionally returns analytic
#' Jacobians \code{dfdx} and \code{dfdu}. Mirrors SPM25's \code{spm_fx_fmri}
#' when called with \code{nargout > 1}.
#'
#' @inheritParams dcm_fx_fmri
#' @return List with \code{f}, \code{dfdx}, \code{D = 1}, \code{dfdu}.
#' @keywords internal
#' @export
dcm_fx_fmri2 <- function(x, u, P, M) {
  if (missing(M) || is.null(M)) M <- list()
  symmetry <- if (!is.null(M$symmetry)) M$symmetry else 0

  P$A <- as.matrix(P$A)
  if (!is.array(P$B)) P$B <- as.array(P$B)
  P$C <- as.matrix(P$C) / 16
  if (!is.array(P$D)) P$D <- as.array(P$D)

  n <- nrow(x); mcols <- ncol(x); f <- x

  if (mcols == 5) {
    if (is.vector(P$A) && length(P$A) > 1) {
      EE <- dcm_fmri_mode_gen(P$A, M$modes)
      if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0) {
        for (i in seq_len(dim(P$B)[3]))
          EE <- EE + as.numeric(u[i]) * P$B[, , i]
      }
      if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0) {
        for (i in seq_len(dim(P$D)[3]))
          EE <- EE + x[i, 1] * P$D[, , i]
      }
    } else {
      Aeff <- P$A
      if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0) {
        for (i in seq_len(dim(P$B)[3]))
          Aeff <- Aeff + as.numeric(u[i]) * P$B[, , i]
      }
      if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0) {
        for (i in seq_len(dim(P$D)[3]))
          Aeff <- Aeff + x[i, 1] * P$D[, , i]
      }
      if (length(dim(P$A)) == 3 && dim(P$A)[3] > 1)
        Aeff <- exp(P$A[, , 1]) - exp(P$A[, , 2])
      SE <- diag(Aeff); EE <- Aeff - diag(exp(SE) / 2 + SE)
      if (isTRUE(symmetry)) EE <- (EE + t(EE)) / 2
    }
    f[, 1] <- as.vector(EE %*% x[, 1] + P$C %*% as.numeric(u))
  } else {
    Aeff <- P$A
    if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0) {
      for (i in seq_len(dim(P$B)[3]))
        Aeff[, , 1] <- Aeff[, , 1] + as.numeric(u[i]) * P$B[, , i]
    }
    if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0) {
      for (i in seq_len(dim(P$D)[3]))
        Aeff[, , 1] <- Aeff[, , 1] + x[i, 1] * P$D[, , i]
    }
    nreg <- n; EE <- exp(Aeff[, , 1]) / 8; IE <- diag(diag(EE)); EE <- EE - IE
    EI <- diag(nreg); SE <- diag(nreg) / 2; SI <- diag(nreg)
    if (length(dim(P$A)) > 2 && dim(P$A)[3] > 1) {
      phi <- dcm_phi(P$A[, , 2] * 2)
      EI <- EI + EE * (1 - phi); EE <- EE * phi - SE
    } else EE <- EE - SE
    f[, 1] <- as.vector(EE %*% x[, 1] - IE %*% x[, 6] + P$C %*% as.numeric(u))
    f[, 6] <- as.vector(EI %*% x[, 1] - SI %*% x[, 6])
  }

  H <- c(0.64, 0.32, 2.00, 0.32, 0.4)

  ex3 <- exp(x[, 3]); ex4 <- exp(x[, 4]); ex5 <- exp(x[, 5])
  x[, 3] <- ex3; x[, 4] <- ex4; x[, 5] <- ex5

  sd <- H[1] * exp(P$decay); tt <- H[3] * exp(P$transit)
  fv <- x[, 4]^(1 / H[4])
  ff <- (1 - (1 - H[5])^(1 / x[, 3])) / H[5]

  f[, 2] <- x[, 1] - as.numeric(sd) * x[, 2] - H[2] * (x[, 3] - 1)
  f[, 3] <- x[, 2] / x[, 3]
  f[, 4] <- (x[, 3] - fv) / (as.numeric(tt) * x[, 4])
  f[, 5] <- (ff * x[, 3] - fv * x[, 5] / x[, 4]) / (as.numeric(tt) * x[, 5])

  f <- as.vector(t(f))
  if (nargs() < 2) return(f)

  dfdx <- vector("list", 0); dfdu <- vector("list", 0)
  set_block <- function(L, i, j, M) {
    if (length(L) < i) length(L) <- i
    if (is.null(L[[i]])) L[[i]] <- list()
    L[[i]][[j]] <- M; L
  }

  I_n <- .Diag(n)

  if (mcols == 5) {
    d11 <- EE
    if (length(dim(P$D)) == 3 && dim(P$D)[3] > 0) {
      for (i in seq_len(dim(P$D)[3])) {
        D_i <- P$D[, , i] + diag((diag(EE) - 1) * diag(P$D[, , i]))
        d11[, i] <- d11[, i] + D_i %*% x[, 1]
      }
    }
    dfdx <- set_block(dfdx, 1, 1, d11)
  } else {
    dfdx <- set_block(dfdx, 1, 1, EE)
    dfdx <- set_block(dfdx, 1, 6, -IE)
    dfdx <- set_block(dfdx, 6, 1, EI)
    dfdx <- set_block(dfdx, 6, 6, -SI)
  }

  Bu <- P$C
  if (length(dim(P$B)) == 3 && dim(P$B)[3] > 0) {
    for (i in seq_len(dim(P$B)[3])) {
      Bm <- P$B[, , i] + diag((diag(EE) - 1) * diag(P$B[, , i]))
      Bu[, i] <- Bu[, i] + Bm %*% x[, 1]
    }
  }
  dfdu <- set_block(dfdu, 1, 1, Bu)

  dfdx <- set_block(dfdx, 2, 1, I_n)
  dfdx <- set_block(dfdx, 2, 2, (-as.numeric(sd)) * I_n)
  dfdx <- set_block(dfdx, 2, 3, .Diag(n, -H[2] * x[, 3]))
  dfdx <- set_block(dfdx, 3, 2, .Diag(n, 1 / x[, 3]))
  dfdx <- set_block(dfdx, 3, 3, .Diag(n, -x[, 2] / x[, 3]))

  d44 <- .Diag(n, -x[, 4]^(1 / H[4] - 1) /
                                (as.numeric(tt) * H[4]) -
                              (x[, 3] - x[, 4]^(1 / H[4])) /
                                as.numeric(tt) / x[, 4])
  dfdx <- set_block(dfdx, 4, 3,
                    .Diag(n, x[, 3] / (as.numeric(tt) * x[, 4])))
  dfdx <- set_block(dfdx, 4, 4, d44)

  d53 <- .Diag(n,
    (x[, 3] + log(1 - H[5]) * (1 - H[5])^(1 / x[, 3]) -
     x[, 3] * (1 - H[5])^(1 / x[, 3])) /
      (as.numeric(tt) * x[, 5] * H[5]))
  d54 <- .Diag(n,
    (x[, 4]^(1 / H[4] - 1) * (H[4] - 1)) / (as.numeric(tt) * H[4]))
  d55 <- .Diag(n,
    (x[, 3] / x[, 5]) * ((1 - H[5])^(1 / x[, 3]) - 1) /
      (as.numeric(tt) * H[5]))

  dfdx <- set_block(dfdx, 5, 3, d53)
  dfdx <- set_block(dfdx, 5, 4, d54)
  dfdx <- set_block(dfdx, 5, 5, d55)

  dfdx <- dcm_cat(dfdx); dfdu <- dcm_cat(dfdu)
  list(f = f, dfdx = dfdx, D = 1, dfdu = dfdu)
}

#' fMRI BOLD observation equation
#'
#' Computes the BOLD signal from the hemodynamic state variables. Mirrors
#' SPM25's \code{spm_gx_fmri}.
#'
#' @inheritParams dcm_fx_fmri
#' @return List with the BOLD prediction \code{g} and Jacobian \code{dgdx}.
#' @export
dcm_gx_fmri <- function(x, u, P, M) {
  x <- as.matrix(x); n <- nrow(x); m <- ncol(x)
  TE <- 0.04; V0 <- 4; ep <- as.numeric(exp(P$epsilon))
  r0 <- 25; nu0 <- 40.3; E0 <- 0.4
  k1 <- as.numeric(4.3 * nu0 * E0 * TE)
  k2 <- as.numeric(ep * r0 * E0 * TE)
  k3 <- as.numeric(1 - ep)
  v <- exp(x[, 4]); q <- exp(x[, 5])
  g <- V0 * (k1 - k1 * q + k2 - k2 * q / v + k3 - k3 * v)

  zero_block <- function() matrix(0, n, n)
  dgdx_blocks <- vector("list", m)
  for (j in seq_len(m)) dgdx_blocks[[j]] <- zero_block()
  dgdx_blocks[[4]] <- diag(-V0 * (k3 * v - k2 * q / v), n, n)
  dgdx_blocks[[5]] <- diag(-V0 * (k1 * q + k2 * q / v), n, n)
  dGdX <- dcm_cat(list(dgdx_blocks))
  list(g = g, dgdx = dGdX)
}

#' Construct fMRI DCM priors
#'
#' Builds the prior expectations \code{pE} and prior covariance \code{pC}
#' for a deterministic fMRI DCM. Mirrors SPM25's \code{dcm_fmri_priors}.
#'
#' @param A Connectivity adjacency matrix.
#' @param B Modulatory adjacency array.
#' @param C Driving-input adjacency matrix.
#' @param D Non-linear adjacency array.
#' @param options List of model options (\code{stochastic}, \code{induced},
#'   \code{two_state}, \code{backwards}, \code{precision}, \code{decay}).
#' @return List with \code{pE} (prior expectation), \code{x} (initial state
#'   template), and \code{pC} (prior covariance).
#' @section Model variants:
#' Only the deterministic, single-state model is supported end-to-end in this
#' release. Branches for the two-state (\code{options$two_state}) and spectral
#' (\code{options$induced}) variants exist but are experimental and are not
#' wired through \code{\link{dcm_estimate}}, which rejects them. They are
#' planned for a future update.
#' @examples
#' # Priors for the bundled three-region (deterministic) model
#' data(toy_dcm)
#' pr <- dcm_fmri_priors(toy_dcm$a, toy_dcm$b, toy_dcm$c,
#'                       D = NULL, options = toy_dcm$options)
#' names(pr)
#' pr$pE$A          # prior expectation of the endogenous connections
#' dim(pr$x)        # 3 regions x 5 hemodynamic states
#' @seealso \code{\link{dcm_estimate}}, which builds these priors for you.
#' @export
dcm_fmri_priors <- function(A, B, C, D, options = list()) {
  A_in <- A; B_in <- B; C_in <- C; D_in <- D
  n <- nrow(A_in)
  if (is.null(options$stochastic)) options$stochastic <- 0
  if (is.null(options$induced))    options$induced    <- 0
  if (is.null(options$two_state))  options$two_state  <- 0
  if (is.null(options$backwards))  options$backwards  <- 0
  if (is.null(D_in))               D_in <- array(0, dim = c(n, n, 0))

  pE <- list(); pC <- list()

  if (isTRUE(options$two_state)) {
    x <- matrix(0, n, 6)
    pA <- if (!is.null(options$precision)) exp(options$precision) else 16
    A_mask <- A_in != 0
    pE$A <- A_mask * 32 - 32; pE$B <- B_in * 0
    pE$C <- C_in * 0; pE$D <- D_in * 0

    Ad <- if (length(dim(A_in)) < 3) array(A_mask, dim = c(n, n, 1)) else A_mask
    pC$A <- Ad / pA
    pC$B <- B_in / 4; pC$C <- C_in * 4; pC$D <- D_in / 4

    if (isTRUE(options$backwards)) {
      if (length(dim(pE$A)) < 3) pE$A <- array(pE$A, dim = c(n, n, 2))
      if (length(dim(pC$A)) < 3) pC$A <- array(pC$A, dim = c(n, n, 2))
      pE$A[, , 2] <- 0; pC$A[, , 2] <- A_mask / pA
    }
  } else {
    x <- matrix(0, n, 5)
    pA <- if (!is.null(options$precision)) exp(options$precision) else 64
    dA <- if (!is.null(options$decay)) options$decay else 1
    A_mask <- A_in != 0

    if (is.null(dim(A_in))) {
      pE$A <- (as.numeric(A_mask) - 1) * dA; pC$A <- as.numeric(A_mask)
    } else {
      pE$A <- A_mask / 128
      Ad <- if (length(dim(A_in)) < 3) array(A_mask, dim = c(n, n, 1)) else A_mask
      pC$A <- Ad / pA
    }
    pE$B <- B_in * 0; pE$C <- C_in * 0; pE$D <- D_in * 0
    pC$B <- B_in; pC$C <- C_in; pC$D <- D_in
  }

  pE$transit <- matrix(0, n, 1); pC$transit <- matrix(1 / 256, n, 1)
  pE$decay   <- matrix(0, 1, 1); pC$decay   <- matrix(1 / 256, 1, 1)
  pE$epsilon <- matrix(0, 1, 1); pC$epsilon <- matrix(1 / 256, 1, 1)

  if (isTRUE(options$induced)) {
    pE$a <- matrix(0, 2, 1); pC$a <- matrix(1 / 64, 2, 1)
    pE$b <- matrix(0, 2, 1); pC$b <- matrix(1 / 64, 2, 1)
    pE$c <- matrix(0, 1, n); pC$c <- matrix(1 / 64, 1, n)
  }

  Cvec <- dcm_vec(pC)
  Cmat <- if (length(Cvec) > 0) diag(as.numeric(Cvec))
          else matrix(numeric(0), 0, 0)
  pC <- Cmat
  list(pE = pE, x = x, pC = pC)
}

# ---- integrate ----

# Bilinear-system integrator (exported): the performance-critical inner loop,
# called O(np) times per Gauss-Newton iteration via dcm_diff.
#
# Key optimizations vs SPM25:
#   1. M0, M1 densified once up front (no per-step sparse->dense).
#   2. expm::expm bound directly (avoids the slower pure-R fallback).
#   3. The propagator E_dt = expm(J*dt) is cached and reused while J and dt
#      are unchanged (the dominant speedup).

#' Integrate a bilinear DCM
#'
#' Forward-integrates the model defined by \code{M$f} (state equation) and
#' \code{M$g} (observation equation) under the input \code{U}, returning the
#' predicted observations at the requested sample points.
#'
#' @param P Parameter structure.
#' @param M Model list. Required fields: \code{f}, \code{x}; optional
#'   \code{g}, \code{ns}, \code{delays}, \code{l}.
#' @param U List with input matrix \code{u} and time step \code{dt}.
#' @return Numeric matrix of predicted observations (rows = samples).
#' @details
#' \code{M$f} and \code{M$g} may be given as functions or as the names of
#' functions. \code{M$n} must be the length of the flattened state
#' (\code{length(dcm_vec(M$x))}, i.e. regions x hidden states), while
#' \code{M$l} is the number of observed outputs (regions) and \code{M$m}
#' the number of inputs.
#' @examples
#' # Simulate the BOLD response of a two-region model to a boxcar input.
#' n <- 2L
#' pri <- dcm_fmri_priors(A = matrix(1, n, n),
#'                        B = array(0, c(n, n, 1)),
#'                        C = matrix(c(1, 0), n, 1),
#'                        D = array(0, c(n, n, 0)),
#'                        options = list())
#'
#' U <- list(u = matrix(c(rep(1, 16), rep(0, 16)), ncol = 1), dt = 1)
#' M <- list(f = "dcm_fx_fmri", g = "dcm_gx_fmri", x = pri$x,
#'           m = ncol(U$u), n = length(pri$x), l = nrow(pri$x), ns = 32)
#'
#' # The priors put C at zero, so start from them and switch on a driving
#' # input to region 1 and a connection from region 1 to region 2.
#' P <- pri$pE
#' P$C[1, 1] <- 1
#' P$A[2, 1] <- 0.4
#'
#' y <- dcm_int(P, M, U)
#' dim(y)                  # 32 samples x 2 regions
#' round(y[seq(1, 32, 4), ], 3)
#' @seealso \code{\link{dcm_estimate}}, which calls this as its forward model.
#' @export
dcm_int <- function(P, M, U) {
  if (!is.list(U)) U <- list(u = U)
  if (is.null(U$u)) stop("U must contain U$u")
  if (is.null(U$dt)) U$dt <- 1

  u_bins <- nrow(U$u)
  v_out  <- if (!is.null(M$ns)) M$ns else u_bins
  x0 <- dcm_vec(M$x); x <- c(1, x0)

  f_fun <- if (!is.null(M$f)) {
    if (is.function(M$f)) M$f
    else if (is.character(M$f) && exists(M$f, mode = "function"))
      get(M$f, mode = "function")
    else { M$x <- numeric(0); function(x, u, P, M) numeric(0) }
  } else { M$x <- numeric(0); function(x, u, P, M) numeric(0) }

  g_fun <- if (!is.null(M$g)) {
    if (is.function(M$g)) M$g
    else if (is.character(M$g) && exists(M$g, mode = "function"))
      get(M$g, mode = "function")
    else function(x, u, P, M) x
  } else function(x, u, P, M) x

  bi <- dcm_bireduce(M, P)
  M0 <- as.matrix(bi$M0)
  M1 <- lapply(bi$M1, as.matrix)
  m_inputs <- length(M1)
  l_out <- if (!is.null(M$l)) M$l
           else nrow(as.matrix(g_fun(M$x, rep(0, m_inputs), P, M)))

  D <- tryCatch({
    if (!is.null(M$delays)) pmax(round(M$delays / U$dt), 1)
    else rep(round(u_bins / v_out), l_out)
  }, error = function(e) rep(round(u_bins / v_out), l_out))

  if (ncol(U$u) == 0) {
    change_rows <- integer(0)
  } else {
    dU <- diff(U$u)
    change_rows <- which(apply(dU, 1, function(r) any(r != 0)))
  }
  i_idx <- c(1, 1 + change_rows)
  su_full <- logical(u_bins); su_full[i_idx] <- TRUE

  s_idx <- ceiling((0:(v_out - 1)) * (u_bins / v_out))
  s_idx[s_idx < 1] <- 1; s_idx[s_idx > u_bins] <- u_bins

  sy_full <- matrix(0L, l_out, u_bins)
  for (j in seq_len(l_out)) {
    ii <- s_idx + D[j]; ii[ii < 1] <- 1; ii[ii > u_bins] <- u_bins
    sy_full[j, ii] <- seq_len(v_out)
  }

  t_mask <- su_full | apply(sy_full != 0L, 2, any)
  t_idx  <- which(t_mask); su <- su_full[t_idx]
  sy     <- sy_full[, t_idx, drop = FALSE]
  dt_vec <- c(diff(t_idx), 0) * U$dt

  any_sy   <- .colSums(sy != 0L, l_out, ncol(sy)) > 0
  Uu_dense <- as.matrix(U$u)

  y <- matrix(0, v_out, l_out); J <- M0
  u_row <- numeric(m_inputs)

  E_dt <- NULL; E_dt_dt <- NA_real_; J_dirty <- TRUE

  for (ii in seq_along(t_idx)) {

    if (su[ii]) {
      if (m_inputs > 0) {
        u_row <- as.numeric(Uu_dense[t_idx[ii], seq_len(m_inputs), drop = TRUE])
        bad <- !is.finite(u_row)
        if (any(bad)) u_row[bad] <- 0
        J <- M0
        for (j in seq_len(m_inputs)) {
          uj <- u_row[j]
          if (uj != 0) J <- J + uj * M1[[j]]
        }
      } else {
        J <- M0
      }
      J_dirty <- TRUE
    }

    if (any_sy[ii]) {
      x_state <- dcm_unvec(x[-1], M$x)
      q_out   <- dcm_vec(g_fun(x_state, u_row, P, M))
      js <- which(sy[, ii] != 0L)
      for (j in js) y[sy[j, ii], j] <- q_out[j]
    }

    dt <- dt_vec[ii]
    if (dt != 0) {
      if (J_dirty || dt != E_dt_dt) {
        E_dt    <- expm::expm(J * dt)
        E_dt_dt <- dt
        J_dirty <- FALSE
      }
      x <- as.numeric(E_dt %*% x)
    }
    x_norm <- sum(abs(x))
    if (!is.finite(x_norm) || x_norm > 1e6) break
  }
  Re(y)
}

# ---- bireduce ----

# Bilinear reduction (M0/M1 form) and Volterra kernels (internal).

#' Bilinear reduction of a non-linear DCM
#'
#' Reduces a non-linear state equation to bilinear form \code{dx/dt = M0 x +
#' sum_i u_i M1[[i]] x}, also returning the lead-field expansion. Mirrors
#' SPM25's \code{spm_bireduce}.
#'
#' @param M Model list (with \code{f}, \code{g}, \code{x}, \code{u}, ...).
#' @param P Parameter structure.
#' @return List with \code{M0}, \code{M1}, \code{L1}, \code{L2}.
#' @keywords internal
#' @export
dcm_bireduce <- function(M, P) {
  funx <- NULL
  if (!is.null(M$f)) {
    if (is.function(M$f)) funx <- M$f
    else if (is.character(M$f) && exists(M$f, mode = "function"))
      funx <- get(M$f, mode = "function")
  }
  if (is.null(funx)) {
    funx <- function(x, u, P, M) numeric(0)
    M$n <- 0; M$x <- numeric(0)
  }

  if (!is.null(M$u)) u <- dcm_vec(M$u)
  else {
    m_dim <- if (!is.null(M$m)) as.integer(M$m) else 0L
    u <- if (m_dim > 0L) matrix(0, m_dim, 1) else numeric(0)
  }

  tmp12 <- dcm_diff(funx, M$x, u, P, M, c(1, 2))
  if (!is.null(tmp12$d2) && !is.null(tmp12$d1)) {
    dfdxu <- tmp12$d2; dfdx <- tmp12$d1
  } else stop("dcm_diff(funx, ..., c(1,2)) did not return the expected {d2,d1}.")

  tmp2 <- dcm_diff(funx, M$x, u, P, M, 2)
  dfdu <- tmp2$J; f0 <- dcm_vec(tmp2$f0); n <- length(f0)

  normalize_dfdu <- function(dfdu_raw, n, m, n_block_cols_hint = NULL) {
    if (is.matrix(dfdu_raw) && nrow(dfdu_raw) == n && ncol(dfdu_raw) == m)
      return(dfdu_raw)
    if (is.list(dfdu_raw) && length(dfdu_raw) == m &&
        all(vapply(dfdu_raw, function(z)
          is.numeric(z) && length(as.vector(z)) == n, TRUE))) {
      out <- matrix(0, n, m)
      for (j in seq_len(m)) out[, j] <- as.vector(dfdu_raw[[j]])
      return(out)
    }
    if (is.list(dfdu_raw) && length(dfdu_raw) > 0 &&
        is.matrix(dfdu_raw[[1]])) {
      big <- dcm_blocks_to_mat(dfdu_raw)
      if (m == 1 && ncol(big) != 1)
        big <- matrix(as.vector(big), nrow = n, ncol = 1)
      if (nrow(big) != n) stop("dfdu stitching produced wrong row count.")
      if (ncol(big) != m) big <- big[, seq_len(min(ncol(big), m)), drop = FALSE]
      return(big)
    }
    if (is.numeric(dfdu_raw) && length(dfdu_raw) == n && m == 1)
      return(matrix(dfdu_raw, nrow = n, ncol = 1))
    stop("Unexpected dfdu shape.")
  }

  m <- if (is.list(dfdxu)) length(dfdxu) else {
    if (!is.null(M$m)) as.integer(M$m) else max(1L, ncol(as.matrix(u)))
  }
  dfdu <- normalize_dfdu(dfdu, n = n, m = m)

  if (is.list(dfdx) && length(dfdx) > 0 && is.matrix(dfdx[[1]]))
    dfdx <- dcm_blocks_to_mat(dfdx)
  else dfdx <- as.matrix(dfdx)

  if (is.list(dfdxu)) {
    dfdxu <- lapply(dfdxu, function(slc) {
      if (is.list(slc) && length(slc) > 0 && is.matrix(slc[[1]]))
        dcm_blocks_to_mat(slc)
      else as.matrix(slc)
    })
  } else if (length(dim(dfdxu)) == 3L && all(dim(dfdxu)[1:2] == c(n, n))) {
    dfdxu <- lapply(seq_len(dim(dfdxu)[3]), function(i) dfdxu[, , i])
  } else stop("Unexpected shape for dfdxu.")

  n <- length(f0)
  x <- dcm_vec(M$x)
  if (!is.null(M$D)) {
    D <- M$D; f0 <- D %*% f0; dfdx <- D %*% dfdx; dfdu <- D %*% dfdu
    if (is.list(dfdxu)) dfdxu <- lapply(dfdxu, function(A) D %*% A)
  }

  M0 <- dcm_cat(as.matrix(list(list(0, NULL),
                                list(f0 - dfdx %*% x, dfdx)), 2, 2,
                           byrow = FALSE))
  M1 <- vector("list", m)
  if (m > 0) {
    for (i in seq_len(m)) {
      dfdxu_i <- if (is.list(dfdxu)) dfdxu[[i]] else dfdxu[, , i]
      rhs_i   <- dfdu[, i, drop = FALSE] - dfdxu_i %*% x
      M1[[i]] <- dcm_cat(as.matrix(list(list(0, NULL),
                                         list(rhs_i, dfdxu_i)), 2, 2))
    }
  }

  fung <- NULL
  if (!is.null(M$g)) {
    if (is.function(M$g)) fung <- M$g
    else if (is.character(M$g) && exists(M$g, mode = "function"))
      fung <- get(M$g, mode = "function")
  }
  if (is.null(fung)) {
    fung <- function(x, u, P, M) dcm_vec(x)
    if (is.null(M$l)) M$l <- n
  }

  tmpg1 <- dcm_diff(fung, M$x, u, P, M, 1)
  if (is.list(tmpg1$f0) && !is.null(tmpg1$f0$g) && !is.null(tmpg1$f0$dgdx)) {
    g0 <- as.numeric(tmpg1$f0$g); dgdx <- as.matrix(tmpg1$f0$dgdx)
  } else stop("dcm_diff(fung, ..., 1) didn't return f0$g and f0$dgdx as expected.")

  l <- length(g0); x <- dcm_vec(M$x); n <- length(x)
  if (!is.matrix(dgdx) || nrow(dgdx) != l || ncol(dgdx) != n) {
    stop(sprintf("dgdx must be l x n; got %dx%d, expected %dx%d",
                 nrow(dgdx), ncol(dgdx), l, n))
  }

  g0_col <- matrix(g0, ncol = 1); x_col <- matrix(as.numeric(x), ncol = 1)
  L1 <- dcm_cat(list(list(g0_col - dgdx %*% x_col, dgdx)))

  .ensure_mat <- function(M, ny, nx) {
    if (is.null(M)) return(matrix(0, ny, nx))
    if (is.matrix(M) || (is.array(M) && length(dim(M)) == 2L)) {
      M <- as.matrix(M)
      if (nrow(M) == ny && ncol(M) == nx) return(M)
      Z <- matrix(0, ny, nx); nr <- min(nrow(M), ny); nc <- min(ncol(M), nx)
      if (nr > 0 && nc > 0)
        Z[seq_len(nr), seq_len(nc)] <- M[seq_len(nr), seq_len(nc), drop = FALSE]
      return(Z)
    }
    if (is.atomic(M) && is.null(dim(M))) {
      v <- as.numeric(M)
      if (!length(v)) return(matrix(0, ny, nx))
      if (length(v) < ny * nx) v <- c(v, rep(0, ny * nx - length(v)))
      return(matrix(v[seq_len(ny * nx)], nrow = ny, ncol = nx))
    }
    if (is.list(M)) {
      v <- suppressWarnings(as.numeric(unlist(M, use.names = FALSE)))
      v[is.na(v)] <- 0
      if (length(v) < ny * nx) v <- c(v, rep(0, ny * nx - length(v)))
      return(matrix(v[seq_len(ny * nx)], nrow = ny, ncol = nx))
    }
    matrix(0, ny, nx)
  }

  tmpg2 <- dcm_diff(fung, M$x, u, P, M, c(1, 1), "nocat")
  dgdxx <- if (!is.null(tmpg2$d2)) tmpg2$d2
           else if (!is.null(tmpg2$J)) tmpg2$J else tmpg2
  if (!is.list(dgdxx) || length(dgdxx) < n)
    stop("dcm_diff did not return expected second-derivative structure.")

  ny <- nrow(dgdx); nx <- ncol(dgdx)
  dgdxx <- lapply(dgdxx, .ensure_mat, ny = ny, nx = nx)

  L2 <- vector("list", l)
  for (i in seq_len(l)) {
    Di <- matrix(0, n, n)
    for (j in seq_len(n))
      Di[j, ] <- .ensure_mat(dgdxx[[j]], ny, nx)[i, , drop = FALSE]
    L2[[i]] <- dcm_cat(dcm_diag(list(0, Di)))
  }
  list(M0 = M0, M1 = M1, L1 = L1, L2 = L2)
}

#' Volterra kernels of a bilinear system
#'
#' Computes the first- and second-order Volterra kernels of a bilinear
#' system specified in M0/M1 form. Mirrors SPM25's \code{spm_kernels}.
#'
#' @param ... Either \code{(M0, M1, N, dt)}, \code{(M0, M1, L1, N, dt)}, or
#'   \code{(M0, M1, L1, L2, N, dt)}.
#' @return List with kernels \code{K0}, \code{K1}, \code{K2}, and the
#'   intermediate response \code{H1}.
#' @keywords internal
#' @export
dcm_kernels <- function(...) {
  args <- list(...); na <- length(args)
  if (na == 4) {
    M0 <- args[[1]]; M1 <- args[[2]]; N <- as.integer(args[[3]])
    dt <- as.numeric(args[[4]]); L1 <- NULL; L2 <- NULL
  } else if (na == 5) {
    M0 <- args[[1]]; M1 <- args[[2]]; L1 <- args[[3]]
    N <- as.integer(args[[4]]); dt <- as.numeric(args[[5]]); L2 <- NULL
  } else if (na == 6) {
    M0 <- args[[1]]; M1 <- args[[2]]; L1 <- args[[3]]; L2 <- args[[4]]
    N <- as.integer(args[[5]]); dt <- as.numeric(args[[6]])
  } else stop("Improper call to dcm_kernels")

  if (is.list(M0) && !is.matrix(M0)) {
    bi <- dcm_bireduce(M0, M1); M0 <- bi$M0; M1 <- bi$M1
    L1 <- bi$L1; L2 <- bi$L2
  }
  if (is.null(L1)) {
    nn <- nrow(as.matrix(M0)); L1 <- diag(nn)[-1, , drop = FALSE]
  }
  if (is.null(L2)) L2 <- list()

  N <- max(1L, N); M0 <- as.matrix(M0)
  n <- nrow(M0); l <- nrow(as.matrix(L1))
  if (!is.list(M1)) stop("M1 must be a list of n x n matrices")
  m <- length(M1)

  H1 <- array(0, dim = c(N, n, m))
  K1 <- array(0, dim = c(N, l, m))
  K2 <- array(0, dim = c(N, N, l, m, m))

  e1 <- expm::expm(dt * M0); e2 <- expm::expm(-dt * M0)
  Mip <- vector("list", N)
  Mip[[1]] <- lapply(seq_len(m),
                     function(p) e1 %*% as.matrix(M1[[p]]) %*% e2)
  if (N >= 2) {
    for (i in 2:N) {
      Mip[[i]] <- lapply(seq_len(m),
                         function(p) e1 %*% Mip[[i - 1]][[p]] %*% e2)
    }
  }

  X0 <- matrix(0, n, 1); X0[1, 1] <- 1
  ei <- e1; if (N >= 2) for (i in 2:N) ei <- e1 %*% ei
  H0 <- ei %*% X0; K0 <- L1 %*% H0

  L1_mat <- as.matrix(L1)
  for (p in seq_len(m)) for (i in seq_len(N)) {
    v <- Mip[[i]][[p]] %*% H0
    H1[i, , p] <- as.numeric(v)
    K1[i, , p] <- as.numeric(t(v) %*% t(L1_mat))
  }

  if (m > 0) {
    for (p in seq_len(m)) for (q in seq_len(m)) for (j in seq_len(N)) {
      H1_slice_t <- t(matrix(H1[j:N, , p], nrow = N - j + 1, ncol = n))
      Hmat <- L1_mat %*% Mip[[j]][[q]] %*% H1_slice_t
      Ht   <- t(Hmat)
      K2[j, j:N, , q, p] <- array(Ht, dim = c(N - j + 1, l))
      K2[j:N, j, , p, q] <- array(Ht, dim = c(N - j + 1, l))
    }
  }

  if (length(L2) > 0) {
    for (i in seq_len(m)) for (j in seq_len(m)) {
      Hi <- matrix(H1[, , i], nrow = N, ncol = n)
      Hj <- matrix(H1[, , j], nrow = N, ncol = n)
      for (pout in seq_len(l)) {
        K2[, , pout, i, j] <- K2[, , pout, i, j] +
                              Hi %*% as.matrix(L2[[pout]]) %*% t(Hj)
      }
    }
  }

  list(K0 = K0, K1 = K1, K2 = K2, H1 = H1)
}

# ---- nlsi-gn ----

# Variational Laplace / Gauss-Newton inversion (exported).

# Internal: NULL/empty fallback operator
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b

#' Variational Laplace inversion of a non-linear system
#'
#' Performs Gauss-Newton optimisation of the variational free energy for a
#' non-linear forward model with Gaussian priors. Mirrors SPM25's
#' \code{dcm_nlsi_GN}.
#'
#' @param M Model specification (list with \code{IS}/\code{f}/\code{g},
#'   priors \code{pE}, \code{pC}, hyperpriors \code{hE}, \code{hC}, ...).
#' @param U Input structure passed through to \code{M$IS}.
#' @param Y Data (or list with \code{$y}, \code{$Q}, \code{$X0}, \code{$dt}).
#' @return List with posterior expectation \code{Ep}, covariance \code{Cp},
#'   log-precision estimate \code{Eh}, free energy \code{F}, and components.
#' @details
#' Most users should call \code{\link{dcm_estimate}} instead, which assembles
#' \code{M}, \code{U} and \code{Y} from a DCM specification and calls this
#' function. Use \code{dcm_nlsi_GN} directly only to invert a non-linear
#' model that is not an fMRI DCM.
#'
#' Progress is reported per Gauss-Newton iteration as
#' \code{EM:(+) k  F: ...}, where \code{(+)} marks an accepted step and
#' \code{(-)} a rejected one. Set \code{M$noprint <- 1} to silence it.
#' @examples
#' # dcm_nlsi_GN is the Gauss-Newton optimiser that dcm_estimate() calls after
#' # assembling M, U and Y from a DCM specification. For a runnable inversion
#' # see the example in ?dcm_estimate; dcm_nlsi_GN returns the same posterior
#' # fields (Ep, Cp, Eh, F).
#' @seealso \code{\link{dcm_estimate}} for the user-facing entry point.
#' @export
dcm_nlsi_GN <- function(M, U, Y) {
  if (is.null(M$nograph)) M$nograph <- 1L
  if (is.null(M$noprint)) M$noprint <- 0L
  if (is.null(M$Nmax))    M$Nmax    <- 128L

  if (is.null(M$IS)) M$IS <- if (!is.null(M$G)) M$G else "dcm_int"
  if (is.null(M$FS)) M$FS <- function(x, ...) x

  y <- tryCatch(Y$y, error = function(e) Y)
  if (is.function(M$IS) && is.function(M$FS)) {
    IS <- function(P, M, U) M$FS(M$IS(P, M, U))
    y  <- M$FS(y, M)
  } else {
    y <- tryCatch(M$FS(y, M),
                  error = function(e)
                    tryCatch(M$FS(y), error = function(e2) y))
    IS <- tryCatch(
      dcm_funcheck(function(P, M, U) M$FS(do.call(M$IS, list(P, M, U)), M)),
      error = function(e) tryCatch(
        dcm_funcheck(function(P, M, U) M$FS(do.call(M$IS, list(P, M, U)))),
        error = function(e2) dcm_funcheck(M$IS)))
  }
  IS <- dcm_funcheck(IS)
  if (!is.null(M$f)) M$f <- dcm_funcheck(M$f)
  if (!is.null(M$g)) M$g <- dcm_funcheck(M$g)
  if (!is.null(M$h)) M$h <- dcm_funcheck(M$h)

  ns <- if (is.list(y)) nrow(y[[1]]) else nrow(as.matrix(y))
  ny <- length(dcm_vec(y)); nr <- ny / ns; M$ns <- ns

  if (is.null(M$x)) {
    if (is.null(M$n)) M$n <- 0L
    M$x <- .spzero(M$n, 1)
  }
  if (missing(U) || is.null(U)) U <- list()
  if (is.null(M$P)) M$P <- M$pE
  dt <- tryCatch(Y$dt, error = function(e) 1)

  # Build error-covariance basis Q.
  .valid_Qblock <- function(M) {
    if (is.null(M)) return(FALSE)
    nr_ <- tryCatch(nrow(M), error = function(e) NULL)
    nc_ <- tryCatch(ncol(M), error = function(e) NULL)
    if (is.null(nr_) || is.null(nc_)) return(FALSE)
    if (length(nr_) != 1L || length(nc_) != 1L) return(FALSE)
    if (!is.finite(nr_) || !is.finite(nc_)) return(FALSE)
    nr_ > 0L && nc_ > 0L
  }
  Q <- tryCatch({
    Q0 <- Y$Q
    if (is.null(Q0))                              dcm_Ce(rep(ns, nr))
    else if (is.numeric(Q0) && !is.matrix(Q0))    list(Matrix::Matrix(Q0, sparse = TRUE))
    else if (inherits(Q0, "Matrix") || is.matrix(Q0)) list(Q0)
    else if (is.list(Q0) && length(Q0) > 0L)      Q0
    else                                          dcm_Ce(rep(ns, nr))
  }, error = function(e) dcm_Ce(rep(ns, nr)))

  if (length(Q) == 0L || !.valid_Qblock(Q[[1]])) {
    Q <- dcm_Ce(rep(ns, nr))
  }
  nh <- length(Q); nq <- ny / nrow(Q[[1]])

  pE <- M$pE
  pC <- tryCatch(M$pC,
                 error = function(e)
                   .Diag(dcm_length(M$pE), x = exp(16)))

  dfdu <- tryCatch({
    nb <- nrow(Y$X0); nx <- ny / nb
    Matrix::kronecker(.Diag(nx), Y$X0)
  }, error = function(e) .spzero(ny, 0))
  if (!length(dfdu)) dfdu <- .spzero(ny, 0)
  # Coerce to a plain numeric matrix so downstream t(), %*%, and solve()
  # don't depend on Matrix S4 dispatch being registered for the session.
  if (inherits(dfdu, "Matrix")) dfdu <- as.matrix(dfdu)

  hE <- tryCatch({
    he <- as.numeric(M$hE)
    if (length(he) != nh) he + rep(0, nh) else he
  }, error = function(e) rep(0, nh) - log(stats::var(dcm_vec(y))) + 4)
  h <- hE

  ihC <- tryCatch({
    H <- dcm_inv(M$hC)
    if (length(H) != nh) H * .Diag(nh) else H
  }, error = function(e) .Diag(nh, x = exp(4)))

  if (is.list(pC)) pC <- dcm_diag(dcm_vec(pC))

  # Coerce pC to a plain numeric matrix so downstream %*%, t(), and solve()
  # use base R dispatch rather than relying on S4 method registration.
  if (inherits(pC, "Matrix")) pC <- as.matrix(pC)
  if (!is.matrix(pC)) {
    stop(sprintf("dcm_nlsi_GN: pC must be a square numeric matrix, got %s.",
                 paste(class(pC), collapse = "/")))
  }
  if (nrow(pC) != ncol(pC)) {
    stop(sprintf("dcm_nlsi_GN: pC must be square; got %dx%d.",
                 nrow(pC), ncol(pC)))
  }

  # Extract V from the SVD and coerce to a plain matrix. This is what fixes
  # the "Error in t.default(V) : argument is not a matrix" failure that can
  # arise when Matrix's S4 dispatch isn't registered (e.g. when the package
  # is source()'d rather than installed).
  svd_pC <- dcm_svd(pC, 0)
  if (is.null(svd_pC$V)) {
    stop("dcm_nlsi_GN: dcm_svd(pC, 0) returned no V; pC may be degenerate.")
  }
  V    <- as.matrix(svd_pC$V)
  if (length(dim(V)) != 2L) {
    stop(sprintf(
      "dcm_nlsi_GN: V from dcm_svd is not a 2D matrix (class %s, length %d).",
      paste(class(V), collapse = "/"), length(V)
    ))
  }
  nu   <- ncol(dfdu); np <- ncol(V)
  ip   <- seq_len(np); iu <- np + seq_len(nu)
  pC_r <- t(V) %*% pC %*% V
  uC   <- .Diag(nu, x = 1e8)
  ipC  <- solve(as.matrix(dcm_cat(dcm_diag(list(pC_r, uC)))))

  Eu  <- if (nu > 0) as.matrix(dcm_pinv(dfdu)) %*% dcm_vec(y)
         else matrix(0, 0, 1)
  p   <- rbind(t(V) %*% (dcm_vec(M$P) - dcm_vec(M$pE)), Eu)
  Ep  <- dcm_unvec(as.numeric(dcm_vec(pE)) +
                     as.numeric(V %*% p[ip, , drop = FALSE]), pE)

  criterion <- c(FALSE, FALSE, FALSE, FALSE); C <- list(F = -Inf); v <- -4
  dFdh  <- matrix(0, nh, 1); dFdhh <- matrix(0, nh, nh)

  h_prev       <- h + Inf
  iS_big_cache <- NULL
  iS_cache     <- NULL
  S_cache      <- NULL

  for (k in seq_len(M$Nmax)) {
    revert <- FALSE; dfdp <- NULL; f_out <- NULL

    tryCatch({
      dfdp_f <- dcm_diff(IS, Ep, M, U, 1, list(V))
      dfdp   <- matrix(dcm_vec(dfdp_f[[1]]), ny, np)
      f_out  <- dfdp_f[[2]]
      normdfdp <- norm(as.matrix(dfdp), type = "I")
      revert <- !is.finite(normdfdp) || normdfdp > 1e32
    }, error = function(e) { revert <<- TRUE })

    if (revert && k > 1) {
      for (i in seq_len(4)) {
        v  <- min(v - 2, -4)
        p  <- C$p + dcm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(v))
        Ep <- dcm_unvec(as.numeric(dcm_vec(pE)) +
                         as.numeric(as.matrix(V) %*% p[ip, , drop = FALSE]), pE)
        tryCatch({
          dfdp_f   <- dcm_diff(IS, Ep, M, U, 1, list(V))
          dfdp     <- matrix(dcm_vec(dfdp_f[[1]]), ny, np)
          f_out    <- dfdp_f[[2]]
          normdfdp <- norm(as.matrix(dfdp), type = "I")
          revert   <- !is.finite(normdfdp) || normdfdp > exp(32)
        }, error = function(e) { revert <<- TRUE })
        if (!revert) break
      }
    }
    if (revert) stop("rsDCM:dcm_nlsi_GN Convergence failure.")

    e  <- dcm_vec(y) - dcm_vec(f_out) - dfdu %*% p[iu, , drop = FALSE]
    J  <- -cbind(dfdp, dfdu)
    tJ <- t(J)

    for (m_step in seq_len(8)) {

      h_diff   <- abs(h - h_prev)
      h_changed <- !all(is.finite(h_diff)) || any(h_diff > 1e-10)
      if (h_changed || is.null(iS_big_cache)) {
        iS <- .spzero(nrow(Q[[1]]), ncol(Q[[1]]))
        for (i in seq_len(nh)) iS <- iS + Q[[i]] * (exp(-32) + exp(h[i]))
        if (nh > 1) S_cache <- dcm_inv(iS)
        iS_big_cache <- Matrix::kronecker(.Diag(nq), iS)
        iS_cache     <- iS
        h_prev       <- h
      }
      iS     <- iS_cache
      iS_big <- iS_big_cache
      S      <- S_cache

      tJ_iSbig <- tJ %*% iS_big
      Pp <- Re(as.numeric(tJ_iSbig %*% J))
      Cp <- as.matrix(dcm_inv(Pp + ipC))

      if (nh > 1) {
        P_list <- vector("list", nh); PS <- vector("list", nh)
        JPJ <- vector("list", nh)
        for (i in seq_len(nh)) {
          Pi       <- Q[[i]] * exp(h[i])
          PS[[i]]  <- Pi %*% S
          P_big    <- Matrix::kronecker(.Diag(nq), Pi)
          P_list[[i]] <- P_big
          JPJ[[i]] <- Re(as.numeric(tJ %*% P_big %*% J))
        }
        for (i in seq_len(nh)) {
          dFdh[i, 1] <- sum(PS[[i]] * .Diag(nrow(PS[[i]]))) *
                          nq / 2 -
                       as.numeric(t(e) %*% P_list[[i]] %*% e) / 2 -
                       sum(t(Cp) * JPJ[[i]]) / 2
          for (j in i:nh) {
            val <- -sum(t(PS[[i]]) * PS[[j]]) * nq / 2
            dFdhh[i, j] <- val; dFdhh[j, i] <- val
          }
        }
      } else {
        dFdh[1, 1]  <- ny / 2 -
                       as.numeric(t(e) %*% iS_big %*% e) / 2 -
                       sum(t(Cp) * Pp) / 2
        dFdhh[1, 1] <- -ny / 2
      }

      dFdhh <- dFdhh + diag(as.numeric(dFdh))
      d <- h - hE; dFdh <- dFdh - ihC %*% d; dFdhh <- dFdhh - ihC

      H_mat <- Re(as.numeric(-dFdhh)); H_mat <- matrix(H_mat, nh, nh)
      H_mat <- 0.5 * (H_mat + t(H_mat))
      Ch <- if (length(H_mat) == 1L) matrix(1 / H_mat, 1, 1) else dcm_inv(H_mat)

      dh <- dcm_dx(as.matrix(dFdhh), as.matrix(dFdh), list(4))
      dh <- pmin(pmax(dh, -1), 1); h <- h + dh
      if (as.numeric(t(dFdh) %*% dh) < 1e-2) break
    }

    A <- 0; B <- 0
    L1 <- dcm_logdet(iS) * nq / 2 -
          as.numeric(t(e) %*% iS_big %*% e) / 2 -
          ny * log(8 * atan(1)) / 2
    L2 <- dcm_logdet(as.matrix(ipC %*% Cp)) / 2 -
          as.numeric(t(p) %*% ipC %*% p) / 2
    L3 <- dcm_logdet(as.matrix(ihC %*% Ch)) / 2 -
          as.numeric(t(d) %*% ihC %*% d) / 2
    F_val <- L1 + L2 + L3; L <- c(L1, L2, L3)

    if (!exists("F0", inherits = FALSE)) F0 <- F_val

    if (F_val > C$F || k < 3) {
      C$p <- p; C$h <- h; C$F <- F_val; C$L <- L; C$Cp <- Cp
      dFdp  <- -Re(as.numeric(tJ_iSbig %*% e)) - ipC %*% p + A / 2
      dFdpp <- -Re(as.numeric(tJ_iSbig %*% J)) - ipC + B / 2
      v <- min(v + 0.5, 4); str <- "EM:(+)"
    } else {
      p <- C$p; h <- C$h; Cp <- C$Cp; v <- min(v - 2, -4); str <- "EM:(-)"
    }

    dp <- dcm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(v), "Q")
    p  <- p + dp
    Ep <- dcm_unvec(as.numeric(dcm_vec(pE)) +
                     as.numeric(as.matrix(V) %*% p[ip, , drop = FALSE]), pE)

    dF_pred <- as.numeric(t(dFdp) %*% dp)
    if (!isTRUE(M$noprint == 1)) {
      message(sprintf("%-6s: %d  F: %-10.3e  dF predicted: %.3e",
                      str, k, C$F - F0, dF_pred))
    }
    criterion <- c((dF_pred < 1e-1), criterion[1:3])
    if (all(criterion)) break
  }

  Ep <- dcm_unvec(as.numeric(dcm_vec(pE)) +
                   as.numeric(as.matrix(V) %*% C$p[ip, , drop = FALSE]), pE)
  Cp <- V %*% C$Cp[ip, ip, drop = FALSE] %*% t(V)
  list(Ep = Ep, Cp = Cp, Eh = C$h, F = C$F, L = C$L,
       dFdp = dFdp, dFdpp = dFdpp)
}

# ---- estimate ----

# Top-level DCM estimation entry point and Bayesian model evidence helpers.

#' Estimate a Dynamic Causal Model for fMRI
#'
#' Inverts a fully-specified DCM using variational Laplace inversion. This
#' is the main user-facing entry point of the package. Mirrors SPM25's
#' \code{dcm_estimate}. Progress is reported via \code{message()} and can be
#' silenced with \code{suppressMessages()}.
#'
#' @param P A DCM list, or a path to an \code{.RData}/\code{.rds} file
#'   containing one.
#' @param save Logical. If \code{TRUE} and \code{P} is a file path, the
#'   estimated DCM is saved back to that path.
#' @return The estimated DCM (a list) with posterior fields populated. The
#'   main fields of interest are \code{Ep} (posterior expectations, with
#'   \code{Ep$A}, \code{Ep$B}, \code{Ep$C} the connectivity estimates),
#'   \code{Cp} (posterior covariance), \code{Pp} (posterior probabilities),
#'   and \code{F} (the negative free energy, used for model comparison).
#' @section Supported model variants:
#' This release implements only the \strong{deterministic, single-state} fMRI
#' DCM. The two-state (\code{options$two_state}), stochastic
#' (\code{options$stochastic}) and spectral / cross-spectral-density
#' (\code{options$induced}) variants are \emph{not} yet implemented: switching
#' any of them on causes \code{dcm_estimate} to stop with an informative error.
#' These variants are planned for a future update. Non-linear DCM (a non-empty
#' \code{d} array) is supported.
#' @examples
#' data(toy_dcm)
#' str(toy_dcm, max.level = 1)
#'
#' \donttest{
#' # Full inversion of the bundled three-region model (takes about a minute).
#' # Progress reporting goes through message(), so it can be silenced.
#' fit <- suppressMessages(dcm_estimate(toy_dcm))
#' round(fit$Ep$A, 3)   # posterior connectivity estimates
#' round(fit$Pp$A, 3)   # posterior probability each connection is non-zero
#' fit$F                # negative free energy, for model comparison
#' }
#' @seealso \code{\link{dcm_nlsi_GN}} for the underlying inversion,
#'   \code{\link{dcm_fmri_priors}} for the priors it constructs.
#' @export
dcm_estimate <- function(P, save = FALSE) {
  if (missing(P)) stop("DCM structure or filename required.")
  loaded_from_file <- FALSE
  file_path <- NULL
  env <- NULL
  if (is.list(P)) {
    DCM <- P
  } else if (is.character(P) && length(P) == 1) {
    file_path <- P
    env <- new.env(parent = emptyenv())
    load(P, envir = env)
    if (!exists("DCM", envir = env))
      stop("File does not contain object 'DCM'.")
    DCM <- get("DCM", envir = env)
    loaded_from_file <- TRUE
  } else stop("Unsupported input for P.")

  if (is.null(DCM$options$two_state))  DCM$options$two_state  <- 0
  if (is.null(DCM$options$stochastic)) DCM$options$stochastic <- 0
  if (is.null(DCM$options$nonlinear))  DCM$options$nonlinear  <- 0
  if (is.null(DCM$options$centre))     DCM$options$centre     <- 0
  if (is.null(DCM$options$hidden))     DCM$options$hidden     <- integer(0)
  if (is.null(DCM$options$hE))         DCM$options$hE         <- 6
  if (is.null(DCM$options$hC))         DCM$options$hC         <- 1 / 128

  # Reject DCM variants that are not yet implemented. Only the deterministic,
  # single-state model is wired through the estimator; two-state / stochastic /
  # induced are checked against the user-supplied options (before the internal
  # defaults below, e.g. the stochastic flag auto-set for null inputs) so we
  # stop loudly instead of silently running the deterministic model.
  .dcm_unsupported <- list(
    two_state  = "the two-state neural model",
    stochastic = "stochastic DCM",
    induced    = "spectral (cross-spectral density) DCM"
  )
  for (.opt in names(.dcm_unsupported)) {
    .val <- DCM$options[[.opt]]
    .on  <- isTRUE(.val) ||
      (is.numeric(.val) && length(.val) == 1L && !is.na(.val) && .val != 0)
    if (.on) {
      stop(sprintf(
        paste0("options$%s = %s requests %s, which is not yet implemented in ",
               "rsDCM. This release supports only the deterministic, ",
               "single-state fMRI DCM; %s is planned for a future update. ",
               "Set options$%s = 0 (or remove it) to run the deterministic ",
               "model."),
        .opt, format(.val), .dcm_unsupported[[.opt]],
        .dcm_unsupported[[.opt]], .opt),
        call. = FALSE)
    }
  }

  if (is.null(DCM$n)) DCM$n <- nrow(DCM$a)
  if (is.null(DCM$v)) DCM$v <- nrow(DCM$Y$y)

  M <- list()
  M$nograph <- if (!is.null(DCM$options$nograph)) DCM$options$nograph else 1L
  M$noprint <- 0L

  if (is.null(DCM$options$maxit)) {
    if (!is.null(DCM$options$nN)) {
      DCM$options$maxit <- DCM$options$nN
      warning("options$nN is deprecated; please use options$maxit")
    } else DCM$options$maxit <- if (DCM$options$stochastic) 32 else 128
  }
  M$Nmax <- if (!is.null(DCM$M$Nmax)) DCM$M$Nmax else DCM$options$maxit

  if (is.null(DCM$options$maxnodes)) {
    if (!is.null(DCM$options$nmax)) {
      DCM$options$maxnodes <- DCM$options$nmax
      warning("options$nmax is deprecated; please use options$maxnodes")
    } else DCM$options$maxnodes <- 8
  }
  DCM$options$induced <- 0

  U <- DCM$U; Y <- DCM$Y; n <- DCM$n; v <- DCM$v
  Y$y <- dcm_detrend(Y$y)
  if (DCM$options$centre == 1) U$u <- dcm_detrend(U$u)

  rng <- max(Y$y) - min(Y$y); scale <- 4 / max(rng, 4)
  Y$y <- Y$y * scale; Y$scale <- scale

  if (is.null(Y$X0)) Y$X0 <- matrix(1, v, 1)
  if (ncol(Y$X0) == 0) Y$X0 <- matrix(1, v, 1)

  if (!is.null(DCM$delays)) M$delays <- DCM$delays
  else M$delays <- matrix(1, n, 1)
  if (!is.null(DCM$TE)) M$TE <- DCM$TE

  if (is.null(DCM$d)) {
    DCM$d <- array(0, dim = c(n, n, 0))
    DCM$options$nonlinear <- 0
  } else DCM$options$nonlinear <- as.integer(dim(DCM$d)[3] > 0)

  if (DCM$options$nonlinear) {
    M$IS <- "dcm_int"
    M$nsteps <- round(max(Y$dt, 1)); M$states <- seq_len(n)
  } else M$IS <- "dcm_int"

  if (is.null(DCM$c) || is.null(U$u)) {
    DCM$c <- matrix(0, n, 1); DCM$b <- array(0, dim = c(n, n, 1))
    U$u <- matrix(0, v, 1); U$name <- list("null")
  }
  if (all(dcm_vec(U$u) == 0) || all(dcm_vec(DCM$c) == 0))
    DCM$options$stochastic <- 1

  pri <- dcm_fmri_priors(DCM$a, DCM$b, DCM$c, DCM$d, DCM$options)
  pE  <- pri$pE; pC <- pri$pC; x0 <- pri$x

  if (!is.null(DCM$options$P))  M$P  <- DCM$options$P
  if (!is.null(DCM$options$pE)) pE   <- DCM$options$pE
  if (!is.null(DCM$options$pC)) pC   <- DCM$options$pC
  if (!is.null(DCM$M$P))        M$P  <- DCM$M$P
  if (!is.null(DCM$M$pE))       pE   <- DCM$M$pE
  if (!is.null(DCM$M$pC))       pC   <- DCM$M$pC

  if (n > DCM$options$maxnodes) {
    y0 <- Y$y - Y$X0 %*% (solve(t(Y$X0) %*% Y$X0) %*% t(Y$X0) %*% Y$y)
    Vn <- as.matrix(dcm_svd(t(y0))$V)
    Vn <- Vn[, seq_len(DCM$options$maxnodes), drop = FALSE]
    j  <- seq_len(n * n); VV <- kronecker(Vn %*% t(Vn), Vn %*% t(Vn))
    pC[j, j] <- VV %*% pC[j, j] %*% t(VV)
  }

  hE <- matrix(DCM$options$hE, n, 1); hC <- diag(DCM$options$hC, n, n)
  i  <- DCM$options$hidden
  if (length(i) > 0) { hE[i] <- -4; diag(hC)[i] <- exp(-16) }

  M$f  <- "dcm_fx_fmri"; M$g  <- "dcm_gx_fmri"; M$x  <- x0
  M$pE <- pE; M$pC <- pC; M$hE <- hE; M$hC <- hC
  M$m  <- ncol(U$u); M$n  <- length(x0); M$l  <- nrow(x0)
  M$N  <- 64; M$dt <- 32 / M$N; M$ns <- v

  fit  <- dcm_nlsi_GN(M, U, Y)
  Ep <- fit$Ep; Cp <- fit$Cp; Eh <- fit$Eh; F <- fit$F

  yhat <- do.call(M$IS, list(Ep, M, U))
  R    <- Y$y - yhat
  X0tX0_inv <- dcm_inv(t(Y$X0) %*% Y$X0)
  R    <- R - Y$X0 %*% (X0tX0_inv %*% t(Y$X0) %*% R)
  Ce   <- exp(-Eh)

  bi <- dcm_bireduce(M, Ep)
  M0 <- bi$M0; M1 <- bi$M1; L1 <- bi$L1; L2 <- bi$L2
  H  <- dcm_kernels(M0, M1, L1, L2, M$N, M$dt)
  H1 <- H$H1; H0 <- H$H0

  L_mat <- .spmat(i = seq_len(n), j = seq_len(n) + 1,
                                x = 1, dims = c(n, ncol(M0)))
  K  <- dcm_kernels(M0, M1, L_mat, M$N, M$dt); K1 <- K$H1

  Tvec <- as.numeric(dcm_vec(pE))
  Pp_vec <- suppressWarnings(
    1 - dcm_Ncdf(Tvec, abs(as.numeric(dcm_vec(Ep))), diag(Cp)))
  Pp <- dcm_unvec(Pp_vec, Ep); Vp <- dcm_unvec(diag(Cp), Ep)

  DCM$M <- M; DCM$Y <- Y; DCM$U <- U; DCM$Ce <- Ce
  DCM$Ep <- Ep; DCM$Cp <- Cp; DCM$Pp <- Pp; DCM$Vp <- Vp
  DCM$H1 <- H1; DCM$K1 <- K1; DCM$R <- R; DCM$y <- yhat; DCM$T <- 0

  y_id <- if (!is.null(M$FS))
            tryCatch(do.call(M$FS, list(Y$y, M)),
                     error = function(e)
                       tryCatch(do.call(M$FS, list(Y$y)),
                                error = function(e2) Y$y))
          else Y$y
  DCM$ID <- dcm_data_id(y_id); DCM$F <- F

  if (loaded_from_file && isTRUE(save)) {
    assign("DCM", DCM, envir = env)
    save(list = "DCM", file = file_path)
  }
  DCM
}

#' Bayesian model reduction (full-rank)
#'
#' Computes the change in log-evidence and the reduced posterior when
#' replacing the original prior with a reduced prior. Mirrors SPM25's
#' \code{dcm_log_evidence}.
#'
#' @param qE Posterior expectation under the original priors.
#' @param qC Posterior covariance.
#' @param pE Original prior expectation.
#' @param pC Original prior covariance.
#' @param rE Reduced prior expectation.
#' @param rC Reduced prior covariance.
#' @param ... Passed to \code{rE} if it is a function.
#' @return List with \code{F}, \code{sE}, \code{sC}.
#' @keywords internal
#' @export
dcm_log_evidence <- function(qE, qC, pE, pC, rE = NULL, rC = NULL, ...) {
  if (!is.null(rE) && is.function(rE)) {
    priors <- tryCatch(rE(...), error = function(e) NULL)
    if (!is.null(priors)) {
      rE_val <- priors[[1]]; rC_val <- priors[[2]]
    } else { rE_val <- rE; rC_val <- rC }
  } else { rE_val <- rE; rC_val <- rC }

  if (is.null(rE_val) || is.null(rC_val)) {
    n <- if (is.matrix(qC)) nrow(qC) else dcm_length(qC)
    rE_val <- .spzero(n, 1)
    rC_val <- .spzero(n, n)
  }

  if (is.list(pC) && !is.matrix(pC)) pC <- diag(dcm_vec(pC))
  if (is.list(qC) && !is.matrix(qC)) qC <- diag(dcm_vec(qC))
  if (is.list(rC_val) && !is.matrix(rC_val)) rC_val <- diag(dcm_vec(rC_val))

  qE_vec <- dcm_vec(qE); pE_vec <- dcm_vec(pE); rE_vec <- dcm_vec(rE_val)
  TOL <- exp(-16); i <- which(diag(pC) > TOL)
  if (!length(i)) return(list(F = 0, sE = qE, sC = qC))

  qP <- dcm_inv(qC[i, i, drop = FALSE], TOL)
  pP <- dcm_inv(pC[i, i, drop = FALSE], TOL)
  rP <- dcm_inv(rC_val[i, i, drop = FALSE], TOL)
  sP <- qP + rP - pP
  sC_r <- dcm_inv(sP, TOL); pC_r <- dcm_inv(pP, TOL)
  sE_vec <- as.matrix(qP) %*% qE_vec[i] +
            as.matrix(rP) %*% rE_vec[i] -
            as.matrix(pP) %*% pE_vec[i]

  term1 <- dcm_logdet(rP %*% qP %*% sC_r %*% pC_r)
  term2 <- as.numeric(t(qE_vec[i]) %*% qP %*% qE_vec[i]) +
           as.numeric(t(rE_vec[i]) %*% rP %*% rE_vec[i]) -
           as.numeric(t(pE_vec[i]) %*% pP %*% pE_vec[i]) -
           as.numeric(t(sE_vec) %*% sC_r %*% sE_vec)
  F <- (term1 - term2) / 2

  rE_vec[i] <- as.numeric(sC_r %*% sE_vec)
  rC_full <- rC_val; rC_full[i, i] <- sC_r
  list(F = F, sE = dcm_unvec(rE_vec, qE), sC = rC_full)
}

#' Bayesian model reduction (subspace projection)
#'
#' Reduced-rank version of \code{\link{dcm_log_evidence}}. Mirrors SPM25's
#' \code{dcm_log_evidence_reduce}.
#'
#' @inheritParams dcm_log_evidence
#' @param TOL Tolerance.
#' @return List with \code{F}, \code{sE}, \code{sC}.
#' @keywords internal
#' @export
dcm_log_evidence_reduce <- function(qE, qC, pE, pC, rE, rC, TOL = 1e-8) {
  pC <- as.matrix(pC)
  if (is.list(pC) && !is.matrix(pC)) {
    v <- dcm_vec(pC); pC <- diag(v[v != 0])
  }
  if (is.list(rC) && !is.matrix(rC)) rC <- diag(dcm_vec(rC))

  RE_orig <- rE; SE_orig <- qE
  svd_result <- dcm_svd(pC, 1e-6); U <- as.matrix(svd_result$U)

  qE_vec <- as.numeric(t(U) %*% dcm_vec(qE))
  pE_vec <- as.numeric(t(U) %*% dcm_vec(pE))
  rE_vec <- as.numeric(t(U) %*% dcm_vec(rE))
  qC_red <- t(U) %*% as.matrix(qC) %*% U
  pC_red <- t(U) %*% pC %*% U
  rC_red <- t(U) %*% rC %*% U

  qP <- dcm_inv(qC_red, TOL)
  pP <- dcm_inv(pC_red, TOL)
  rP <- dcm_inv(rC_red, TOL)
  sP <- qP + rP - pP
  sC_red <- dcm_inv(sP, TOL); pC_inv <- dcm_inv(pP, TOL)
  sE_vec <- as.matrix(qP) %*% qE_vec +
            as.matrix(rP) %*% rE_vec -
            as.matrix(pP) %*% pE_vec

  term1 <- dcm_logdet(as.matrix(rP %*% qP %*% sC_red %*% pC_inv))
  term2 <- as.numeric(t(qE_vec) %*% qP %*% qE_vec) +
           as.numeric(t(rE_vec) %*% rP %*% rE_vec) -
           as.numeric(t(pE_vec) %*% pP %*% pE_vec) -
           as.numeric(t(sE_vec) %*% sC_red %*% sE_vec)
  F <- (term1 - term2) / 2

  pE_orig_vec <- dcm_vec(RE_orig)
  rE_restored <- as.numeric(sC_red %*% sE_vec)
  sE_final_vec <- as.numeric(U %*% rE_restored + pE_orig_vec -
                               U %*% (t(U) %*% pE_orig_vec))
  sC_final <- U %*% sC_red %*% t(U)
  list(F = F, sE = dcm_unvec(sE_final_vec, SE_orig), sC = sC_final)
}

#' Approximate model evidence (AIC, BIC)
#'
#' AIC and BIC penalties for an estimated DCM, plus per-region cost terms.
#' Mirrors SPM25's \code{spm_dcm_evidence}.
#'
#' @param DCM An estimated DCM.
#' @return List with per-region cost, AIC penalty, BIC penalty, and overall
#'   AIC and BIC.
#' @export
dcm_evidence <- function(DCM) {
  v <- DCM$v; n <- DCM$n
  Cp_diag <- if (inherits(DCM$Cp, c("Matrix", "sparseMatrix", "denseMatrix")))
                Matrix::diag(DCM$Cp)
             else if (is.matrix(DCM$Cp)) diag(DCM$Cp) else DCM$Cp
  wsel <- which(Cp_diag != 0)
  evidence <- list(region_cost = numeric(n))
  for (i in seq_len(n)) {
    lambda_i <- tryCatch({
      if (is.matrix(DCM$Ce) || inherits(DCM$Ce, "Matrix"))
        DCM$Ce[i * v, i * v]
      else if (is.vector(DCM$Ce) && length(DCM$Ce) > 1) DCM$Ce[i]
      else DCM$Ce
    }, error = function(e)
        if (length(DCM$Ce) >= i) DCM$Ce[i] else DCM$Ce[1])
    R_i <- DCM$R[, i]
    evidence$region_cost[i] <- -0.5 * v * log(lambda_i) -
                                0.5 * sum(R_i * (1 / lambda_i) * R_i)
  }
  evidence$aic_penalty <- length(wsel)
  evidence$bic_penalty <- 0.5 * length(wsel) * log(v)
  evidence$aic_overall <- sum(evidence$region_cost) - evidence$aic_penalty
  evidence$bic_overall <- sum(evidence$region_cost) - evidence$bic_penalty
  evidence
}

# ---- peb ----

# Parametric Empirical Bayes (PEB): hierarchical modelling of DCM parameters.
#
# R-native estimator. MATLAB .mat subject PEBs are normalised to native R-list
# shape at load time (.normalize_mat_peb) so a single prepare/run path serves
# both .rds and .mat inputs. %||% is not redefined here; the package-level
# operator from nlsi-gn.R (NULL or zero-length falls through) is used throughout.

# Normalise a MATLAB-struct PEB (from R.matlab::readMat) to native shape.
# readMat wraps structs as 3-D arrays. The estimator only needs the top-level
# struct and its $M sub-struct unwrapped; all other fields are plain after that.
.normalize_mat_peb <- function(x) {
  if (!is.null(dim(x)) && length(dim(x)) >= 3) x <- x[, , 1]
  if (!is.null(x$M) && !is.null(dim(x$M)) && length(dim(x$M)) >= 3) x$M <- x$M[, , 1]
  x
}

#' Locate subject PEB files in a directory
#'
#' Finds files named like \code{PEB_sub*.rds} or \code{PEB_sub*.mat} and
#' extracts the numeric subject identifier from each filename.
#'
#' @param peb_dir Directory to search.
#' @return List with \code{files} (full paths) and \code{subjects} (integer
#'   identifiers), in the order returned by \code{\link{list.files}}.
#' @examples
#' # Set up a directory with two dummy subject PEB files
#' d <- tempfile(); dir.create(d)
#' saveRDS(list(), file.path(d, "PEB_sub-01.rds"))
#' saveRDS(list(), file.path(d, "PEB_sub-02.rds"))
#'
#' found <- dcm_peb_files(d)
#' basename(found$files)
#' found$subjects
#'
#' unlink(d, recursive = TRUE)
#' @keywords internal
#' @export
dcm_peb_files <- function(peb_dir) {
  files <- list.files(peb_dir,
                      pattern = "PEB_sub.*\\.(rds|mat)$",
                      full.names = TRUE,
                      ignore.case = TRUE)

  if (!length(files)) stop("No PEB files found in: ", peb_dir)

  # Extract subject IDs (flexible: sub01, sub-01, sub_01)
  subjects <- sub(".*sub[-_]?([0-9]+).*", "\\1", basename(files))
  subjects <- suppressWarnings(as.integer(subjects))

  valid <- !is.na(subjects)
  if (!any(valid))
    stop("Found PEB files but no subject number could be parsed from any ",
         "filename. Expected something like 'PEB_sub-01.rds'.")

  list(
    files    = files[valid],
    subjects = subjects[valid]
  )
}

#' Load subject PEBs from .rds or .mat files
#'
#' Loads every subject PEB found by \code{\link{dcm_peb_files}} and normalises
#' MATLAB structs to native R-list shape.
#'
#' @param peb_dir Directory containing the subject PEB files.
#' @param subjects Optional integer vector; keep only these subject numbers.
#' @return List with \code{pebs} (named list of PEB structures),
#'   \code{subjects}, and \code{files}.
#' @details
#' Reading \code{.mat} files requires the \pkg{R.matlab} package. It is listed
#' in \code{Suggests}, so install it only if you have MATLAB PEBs to read;
#' \code{.rds} input needs nothing extra.
#' @seealso \code{\link{dcm_peb_of_pebs}}, which calls this.
#' @keywords internal
#' @export
dcm_peb_load <- function(peb_dir, subjects = NULL) {

  D <- dcm_peb_files(peb_dir)

  if (!is.null(subjects)) {
    keep <- D$subjects %in% subjects
    D$files    <- D$files[keep]
    D$subjects <- D$subjects[keep]
  }

  if (!length(D$files)) stop("No matching subjects found.")

  needs_matlab <- any(grepl("\\.mat$", D$files, ignore.case = TRUE))
  if (needs_matlab && !requireNamespace("R.matlab", quietly = TRUE)) {
    stop("Reading .mat PEB files requires the 'R.matlab' package. ",
         "Install it with install.packages(\"R.matlab\"), or supply .rds files.")
  }

  pebs <- lapply(D$files, function(f) {
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      readRDS(f)                                    # native R PEB
    } else {
      m <- R.matlab::readMat(f, sparseMatrix = "matrix")
      m$PEB                                         # MATLAB PEB struct
    }
  })

  # unify both formats to native R-list shape for the single prepare/run path
  pebs <- lapply(pebs, .normalize_mat_peb)

  names(pebs) <- paste0("sub", sprintf("%03d", D$subjects))
  list(pebs = pebs, subjects = D$subjects, files = D$files)
}

#' Prepare a PEB (second- or third-level) model
#'
#' Gathers the first-level posterior densities, selects the parameter indices
#' implied by \code{field}, projects them onto a rank-reduced subspace, and
#' builds the second-level priors, hyperpriors and precision components.
#' Mirrors the preparation half of SPM25's \code{spm_dcm_peb}.
#'
#' @param P List of estimated DCMs (or, for a PEB-of-PEBs, a list of PEBs).
#'   Each element must carry \code{M$pE}, \code{M$pC}, \code{Ep}, \code{Cp}
#'   and \code{F}.
#' @param M Second-level model specification. Recognised fields: \code{X}
#'   (design matrix, subjects x covariates), \code{W}, \code{Q} (covariance
#'   component option, see Details), \code{bE}, \code{bC}, \code{pC},
#'   \code{hE}, \code{hC}, \code{alpha}, \code{beta}, \code{Xnames},
#'   \code{maxit}. Defaults to an intercept-only design.
#' @param field Character vector of parameter blocks to model (e.g.
#'   \code{c("A", "B")}), the string \code{"all"} for every named block, or
#'   numeric indices into \code{dcm_vec(M$pE)}.
#' @return A list ("prep") consumed by \code{\link{dcm_peb_run}}, holding the
#'   projected densities, design, priors, precision components and sizes.
#' @details
#' \code{M$Q} selects the form of the second-level precision components:
#' \code{"single"} (one component, the default), \code{"fields"} (one per
#' requested field), \code{"all"} (one per parameter), \code{"none"} (no
#' components), or a list of numeric matrices for manual specification.
#' @seealso \code{\link{dcm_peb_run}} to fit the prepared model.
#' @export
dcm_peb_prepare <- function(P, M = list(), field = c("A", "B")) {

  if (is.null(M)) {
    M <- list(X = matrix(1, nrow = length(P), ncol = 1))
  }
  if (is.numeric(M)) {
    M <- list(X = M)
  }
  if (!("X" %in% names(M)) || is.null(M$X)) {
    M$X <- matrix(1, nrow = length(P), ncol = 1)
  }

  # Normalize field
  if (missing(field) || is.null(field)) field <- c("A", "B")
  if (is.list(field)) field <- unlist(field, use.names = FALSE)

  if (!(is.character(field) || is.numeric(field))) {
    stop("'field' must be a character vector (e.g. c('A','B')), the string ",
         "'all', or numeric indices")
  }

  Ns  <- length(P)
  DCM <- P[[1]]
  if (is.null(DCM$M) || is.null(DCM$M$pE))
    stop("P[[1]]$M$pE is required")

  pE <- DCM$M$pE

  # Resolve 'field' names if needed
  # "all" means every parameter. For a first-level DCM, M$pE has named blocks
  # (A/B/C/...) and we expand to those names. For a PEB input (a PEB-of-PEBs),
  # M$pE is an unnamed kronecker product, so "all" resolves positionally to
  # every element instead; this is the default path of dcm_peb_of_pebs.
  if (is.character(field) && length(field) == 1L && tolower(field) == "all") {
    field <- if (is.null(names(pE))) {
      seq_len(length(dcm_vec(pE)))
    } else {
      names(pE)
    }
  }

  maxit <- if ("maxit" %in% names(M)) M$maxit else 64

  if ("bC" %in% names(M) && Ns > 1) {
    q <- dcm_find_pC(M$bC, M$bE, field)$i
  } else if (is.numeric(field)) {
    q <- field
  } else {
    q <- dcm_find_pC(DCM, field)$i
  }

  is_peb_of_pebs <- !is.null(DCM$Pnames)
  if (is_peb_of_pebs) {
    Pstr <- character(0)
    for (i in seq_along(DCM$Xnames)) {
      Pstr <- c(Pstr, paste0(DCM$Xnames[[i]], ": ", DCM$Pnames))
    }
  } else {
    Pstr <- tryCatch(dcm_fieldindices(DCM$M$pE, field),
                     error = function(e) paste0("P", q))
    if (length(q) == 1) Pstr <- list(Pstr)
  }

  pE <- vector("list", Ns)
  pC <- vector("list", Ns)
  qE <- vector("list", Ns)
  qC <- vector("list", Ns)
  iF <- numeric(Ns)
  Ne <- NULL

  for (i in seq_len(Ns)) {
    D <- P[[i]]
    pC[[i]] <- if (is.list(D$M$pC)) diag(dcm_vec(D$M$pC)) else as.matrix(D$M$pC)
    pE[[i]] <- dcm_vec(D$M$pE)
    qE[[i]] <- dcm_vec(D$Ep)
    qC[[i]] <- as.matrix(D$Cp)
    iF[i]   <- as.numeric(D$F)
    if (i == 1) Ne <- length(pE[[i]])
    if (length(pE[[i]]) != Ne) stop("All DCMs must have same parameterisation.")
  }

  # Average priors over subjects (on selected indices q)
  PE <- numeric(length(q))
  PC <- matrix(0, length(q), length(q))
  for (i in seq_len(Ns)) {
    PE <- PE + pE[[i]][q]
    PC <- PC + pC[[i]][q, q, drop = FALSE]
  }
  PE <- PE / Ns
  PC <- PC / Ns

  # Rank-reduction (projector U)
  U <- if (Ns > 1) as.matrix(dcm_svd(PC)$U) else diag(length(q))

  # Project densities and apply shrinkage
  for (i in seq_len(Ns)) {
    pE[[i]] <- as.vector(t(U) %*% pE[[i]][q])
    pC[[i]] <- t(U) %*% pC[[i]][q, q, drop = FALSE] %*% U
    qE[[i]] <- as.vector(t(U) %*% qE[[i]][q])
    qC[[i]] <- t(U) %*% qC[[i]][q, q, drop = FALSE] %*% U
    if (Ns > 1) {
      qC[[i]] <- dcm_inv(dcm_inv(qC[[i]]) + dcm_inv(pC[[i]]) / 16)
    }
  }

  if (Ns > 1) {
    W <- t(U) %*% (M$W %||% diag(length(q))) %*% U
    X <- M$X
  } else {
    W <- M$W %||% diag(length(q))
    X <- matrix(1, 1, 1)
  }

  # Covariance component option
  OPTION <- "single"; Qcells <- list()
  if (!is.null(M$Q)) {
    if (is.list(M$Q) && is.numeric(M$Q[[1]])) { OPTION <- "manual"; Qcells <- M$Q }
    else { OPTION <- M$Q }
  }

  alpha <- M$alpha %||% 1
  beta  <- M$beta  %||% 16
  Nx <- ncol(X)
  # Prefer explicit M$Xnames, then the design's own column names (which
  # dcm_peb_design sets from the covariate data frame), and only fall back to
  # positional labels if neither is available.
  Xnames <- M$Xnames %||% as.list(colnames(X)) %||%
    lapply(seq_len(Nx), function(i) sprintf("Covariate %d", i))

  Np <- length(q)
  if (!is.null(M$bE)) {
    bE0 <- dcm_vec(M$bE); if (Ns > 1 && length(bE0) > Np) bE0 <- bE0[q]
  } else bE0 <- PE
  if (!is.null(M$bC)) {
    bC0 <- if (is.list(M$bC)) diag(dcm_vec(M$bC)) else as.matrix(M$bC)
    if (Ns > 1 && nrow(bC0) > Np) bC0 <- bC0[q, q, drop = FALSE]
  } else bC0 <- PC / alpha
  if (!is.null(M$pC)) {
    pC2 <- if (is.list(M$pC)) diag(dcm_vec(M$pC)) else as.matrix(M$pC)
    if (Ns > 1 && nrow(pC2) > Np) pC2 <- pC2[q, q, drop = FALSE]
  } else if (isTRUE(beta == 0)) {
    qEmat <- do.call(cbind, qE)
    pC2 <- diag(apply(qEmat, 1, stats::var))
  } else {
    pC2 <- PC / beta
  }

  # Base precision in projected space
  pQ <- dcm_inv(t(U) %*% pC2 %*% U)

  # Build Q components in projected space
  Q <- list()
  if (OPTION == "single") {
    Q <- list(pQ)
  } else if (OPTION == "fields") {
    pq <- dcm_inv(pC2)
    k <- 1; off <- 0L
    for (nm in names(DCM$M$pE)) {
      len <- length(dcm_vec(DCM$M$pE[[nm]]))
      if (nm %in% field) {
        ids <- off + seq_len(len)
        ids <- which(ids %in% q)
        if (length(ids)) {
          Qi <- matrix(0, Np, Np)
          Qi[ids, ids] <- pq[ids, ids, drop = FALSE]
          Q[[k]] <- t(U) %*% Qi %*% U; k <- k + 1
        }
      }
      off <- off + len
    }
  } else if (OPTION == "all") {
    pq <- dcm_inv(pC2)
    k <- 1
    for (i in seq_len(Np)) {
      Qi <- matrix(0, Np, Np); Qi[i, i] <- pq[i, i]
      qk <- t(U) %*% Qi %*% U
      if (any(qk != 0)) { Q[[k]] <- qk; k <- k + 1 }
    }
  } else if (OPTION == "manual") {
    pq <- dcm_inv(pC2)
    k <- 1
    for (i in seq_along(Qcells)) {
      diag_idx <- which(diag(Qcells[[i]]) != 0)
      jj <- which(seq_len(length(dcm_vec(DCM$M$pE))) %in% q[diag_idx])
      Qi <- matrix(0, Np, Np)
      if (length(jj)) Qi[jj, jj] <- pq[jj, jj, drop = FALSE]
      qk <- t(U) %*% Qi %*% U
      if (any(qk != 0)) { Q[[k]] <- qk; k <- k + 1 }
    }
  } # 'none' -> Q empty

  Ng <- length(Q); Nw <- ncol(W); Nb <- Nw * Nx

  gE <- rep(M$hE %||% 0, Ng)
  gC <- diag(M$hC %||% (1 / 16), nrow = Ng, ncol = Ng)

  # X scaling (Xc)
  Xc <- tryCatch(diag(nrow(X) / colSums(X^2)), error = function(e) diag(1, Nx))

  # Prior expectations and covariances for second-level parameters
  Xe <- matrix(1, nrow = Nx, ncol = 1)   # column of ones (dcm_speye(Nx, 1))
  bE <- as.vector(kronecker(Xe, t(U) %*% bE0))
  bC <- kronecker(Xc, t(U) %*% bC0 %*% U)
  bP <- dcm_inv(bC)
  gP <- dcm_inv(gC)

  b <- bE
  g <- gE

  list(
    # Inputs
    P = P, M = M, field = field,

    # Core sets
    Ns = Ns, X = X, W = W, U = U, q = q, Pstr = Pstr,
    pE = pE, pC = pC, qE = qE, qC = qC, iF = iF,

    # Sizes
    Nx = ncol(X), Nw = ncol(W), Ng = Ng, Np = length(q), Nb = Nb,

    # Priors/components
    bE0 = bE0, bC0 = bC0, pC2 = pC2, pQ = pQ, Q = Q,
    gE = gE, gC = gC, Xnames = Xnames, Xc = Xc,

    # Second-level priors (projected)
    bE = bE, bC = bC, bP = bP, gP = gP,

    # Initial parameter vectors
    b = b, g = g,

    # Control
    maxit = maxit,
    is_peb_of_pebs = is_peb_of_pebs
  )
}

#' Fit a prepared PEB model by variational Laplace
#'
#' Runs the variational Laplace (Fisher scoring) loop on the structure built by
#' \code{\link{dcm_peb_prepare}}, and returns the group-level PEB together with
#' the first-level DCMs updated under the empirical priors. Mirrors the
#' estimation half of SPM25's \code{spm_dcm_peb}.
#'
#' @param prep The list returned by \code{\link{dcm_peb_prepare}}.
#' @param verbose Logical. Report free energy per iteration via
#'   \code{message()}. Silence with \code{suppressMessages()} or
#'   \code{verbose = FALSE}.
#' @return List with two elements: \code{PEB}, the group-level result
#'   (\code{Ep} group parameter expectations, \code{Cp} their covariance,
#'   \code{Eh}/\code{Ch} log-precision estimates, \code{F} free energy, and the
#'   \code{Pnames}/\code{Xnames}/\code{Snames} labels), and \code{P}, the input
#'   DCMs with \code{Ep}, \code{Cp}, \code{M$pE}, \code{M$pC} and \code{F}
#'   updated under the empirical priors.
#' @details
#' Progress is reported as \code{VL Iter k: F=... dF=... [t=...]}, where
#' \code{t} is the log step size of the Fisher scoring update. The loop stops
#' when the step size collapses or the free-energy increase falls below 1e-4.
#' @seealso \code{\link{dcm_peb_prepare}}, \code{\link{dcm_peb_of_pebs}}.
#' @export
dcm_peb_run <- function(prep, verbose = TRUE) {
  P   <- prep$P; M <- prep$M
  Ns  <- prep$Ns; X <- prep$X; W <- prep$W; U <- prep$U
  q   <- prep$q; Pstr <- prep$Pstr
  pE  <- prep$pE; pC <- prep$pC; qE <- prep$qE; qC <- prep$qC; iF <- prep$iF
  Nx  <- prep$Nx; Nw <- prep$Nw; Ng <- prep$Ng; Nb <- prep$Nb
  bE0 <- prep$bE0; bC0 <- prep$bC0; pC2 <- prep$pC2; pQ <- prep$pQ; Q <- prep$Q
  gE  <- prep$gE; gC <- prep$gC; Xc <- prep$Xc; Xnames <- prep$Xnames
  bE  <- prep$bE; bC <- prep$bC; bP <- prep$bP; gP <- prep$gP
  b   <- prep$b; g <- prep$g
  maxit <- prep$maxit
  is_peb_of_pebs <- prep$is_peb_of_pebs

  # Joint prior precision over (b, g). With Q = "none" there are no covariance
  # components, hence no hyperparameters, and the g block is absent entirely --
  # building it would ask dcm_cat to reconcile a 1x4 scalar block against a
  # 0x0 gC.
  ipC <- if (Ng > 0) {
    dcm_cat(matrix(list(bP, NULL, 0, dcm_inv(gC)), 2, 2, byrow = TRUE))
  } else {
    as.matrix(bP)
  }

  # VL loop
  t_fs <- -4
  F0 <- -Inf
  tmp <- NULL
  dF <- NA_real_
  dFdgg <- NULL; dFdbb <- NULL

  for (niter in seq_len(maxit)) {
    # rP = floor + sum exp(g_i) Q_i
    rP <- pQ * exp(-8)
    if (Ng > 0) {
      for (i in seq_len(Ng)) rP <- rP + exp(g[i]) * Q[[i]]
    }
    rC <- dcm_inv(rP)

    # Gradients and curvature. Keep the g blocks as base matrices with zero
    # extent when Ng == 0 so the rbind/cbind assembly below stays conformable.
    F <- 0
    dFdb  <- -bP %*% (b - bE)
    dFdbb <- -bP
    if (Ng > 0) {
      dFdg  <- -dcm_inv(gC) %*% (g - gE)
      dFdgg <- as.matrix(-dcm_inv(gC))
    } else {
      dFdg  <- matrix(0, nrow = 0, ncol = 1)
      dFdgg <- matrix(0, nrow = 0, ncol = 0)
    }
    dFdbg <- matrix(0, nrow = Nb, ncol = Ng)

    for (i in seq_len(Ns)) {
      Xi <- kronecker(matrix(X[i, ], nrow = 1), W)
      rE_i <- as.vector(Xi %*% b)

      red <- dcm_log_evidence_reduce(qE[[i]], qC[[i]], pE[[i]], pC[[i]], rE_i, rC)
      Fi <- red$F; sE <- red$sE; sC <- red$sC
      F <- F + Fi + iF[i]
      dE <- sE - rE_i

      dFdb  <- dFdb  + t(Xi) %*% (rP %*% dE)
      dFdbb <- dFdbb + t(Xi) %*% ((rP %*% sC %*% rP - rP)) %*% Xi

      if (Ng > 0) {
        for (j in seq_len(Ng)) {
          dFdgj <- exp(g[j]) * ((sum(diag((rC - sC) %*% Q[[j]])) -
                                   t(dE) %*% Q[[j]] %*% dE) / 2)
          dFdg[j] <- dFdg[j] + as.numeric(dFdgj)
          dFdgg[j, j] <- dFdgg[j, j] + as.numeric(dFdgj)

          dFdbgj <- exp(g[j]) * t((Xi - sC %*% rP %*% Xi)) %*% Q[[j]] %*% dE
          dFdbg[, j] <- dFdbg[, j] + as.vector(dFdbgj)

          for (k in seq_len(Ng)) {
            dFdggj <- exp(g[j] + g[k]) * (
              (sum(diag((rC %*% Q[[k]] %*% rC - sC %*% Q[[k]] %*% sC) %*% Q[[j]])) / 2)
              - t(dE) %*% Q[[k]] %*% sC %*% Q[[j]] %*% dE
            )
            dFdgg[j, k] <- dFdgg[j, k] - as.numeric(dFdggj)
          }
        }
      }
    }

    # Block matrices and covariances
    dFdp  <- rbind(dFdb, dFdg)
    dFdpp <- rbind(cbind(as.matrix(dFdbb), dFdbg),
                   cbind(t(dFdbg), dFdgg))
    Cp    <- dcm_inv(-dFdpp)

    # Complexity of second level
    Fb <- as.numeric(t(b) %*% bP %*% b)
    Fg <- if (Ng > 0) as.numeric(t(g) %*% dcm_inv(gC) %*% g) else 0
    Fc <- Fb / 2 + Fg / 2 - dcm_logdet(as.matrix(ipC) %*% as.matrix(Cp)) / 2
    F  <- F - Fc

    if (niter == 1) F0 <- F
    if (F >= F0) {
      dF <- F - F0; F0 <- F
      tmp <- list(b = b, g = g, F0 = F0, dFdb = dFdb, dFdbb = dFdbb,
                  dFdg = dFdg, dFdgg = dFdgg)
      t_fs <- min(t_fs + 0.25, 2)
    } else {
      t_fs <- max(t_fs - 1, -4)
      b    <- tmp$b; g <- tmp$g; F0 <- tmp$F0
      dFdb <- tmp$dFdb; dFdbb <- tmp$dFdbb; dFdg <- tmp$dFdg; dFdgg <- tmp$dFdgg
    }

    # Fisher scoring step (Gauss-Newton update via dcm_dx)
    dp <- dcm_dx(as.matrix(dFdpp), as.matrix(dFdp), list(t_fs))
    if (sum(abs(dp)) >= 8) {
      dFdpp_blk <- rbind(cbind(dFdbb, matrix(0, nrow(dFdbb), ncol(dFdgg))),
                         cbind(matrix(0, nrow(dFdgg), ncol(dFdbb)), dFdgg))
      dp <- dcm_dx(dFdpp_blk, dFdp, list(t_fs))
    }
    db <- dp[seq_len(length(b))]
    b  <- b + db
    if (Ng > 0) {
      # NB: guard the slice; with Ng == 0 the expression (length(b)+1):length(dp)
      # would count *backwards* and index past the end of dp.
      dg <- dp[(length(b) + 1):length(dp)]
      g  <- g + tanh(dg)
    }

    if (isTRUE(verbose)) {
      message(sprintf("VL Iter %3d: F=%8.3f  dF=%8.4f  [t=%+.2f]",
                      niter, F, if (is.na(dF)) 0 else dF, t_fs))
    }
    if (niter > 4 && (t_fs <= -4 || (!is.na(dF) && dF < 1e-4))) break
  }

  # Assemble PEB
  Ub <- kronecker(diag(Nx), U)
  PEB <- list()
  PEB$Snames <- character(Ns)
  PEB$Pnames <- unlist(Pstr)
  PEB$Pind   <- q
  PEB$Xnames <- unlist(Xnames)
  if (is_peb_of_pebs) {
    PEB$Pind0 <- do.call(rbind, replicate(ncol(M$X), P[[1]]$Pind0, simplify = FALSE))
  } else {
    PEB$Pind0 <- PEB$Pind
  }

  PEB$M <- list()
  PEB$M$X  <- X
  PEB$M$W  <- W
  PEB$M$U  <- U
  PEB$M$pE <- kronecker(matrix(1, nrow = Nx, ncol = 1), M$bE %||% bE0)
  PEB$M$pC <- kronecker(Xc, M$bC %||% bC0)
  PEB$M$hE <- gE
  PEB$M$hC <- gC

  rP_final <- pQ * exp(-8)
  if (Ng > 0) for (i in seq_len(Ng)) rP_final <- rP_final + exp(g[i]) * Q[[i]]

  PEB$Ep <- U %*% matrix(b, nrow = Nw, ncol = Nx)
  PEB$Eh <- g
  PEB$Ch <- if (Ng > 0) dcm_inv(-dFdgg) else matrix(0, nrow = 0, ncol = 0)
  PEB$Cp <- Ub %*% dcm_inv(-dFdbb) %*% t(Ub)
  PEB$Ce <- U %*% dcm_inv(rP_final) %*% t(U)
  PEB$F  <- F
  PEB$M$Q <- lapply(Q, function(Qi) U %*% Qi %*% t(U))

  # Update first-level DCMs with the empirical priors
  P_out <- vector("list", Ns)
  for (i in seq_len(Ns)) {
    D <- P[[i]]
    D_pC <- if (is.list(D$M$pC)) diag(dcm_vec(D$M$pC)) else as.matrix(D$M$pC)
    RP <- dcm_inv(D_pC)
    RP[q, q] <- U %*% rP_final %*% t(U)
    RC <- dcm_inv(RP)

    RE <- dcm_vec(D$M$pE)
    # Use first-row covariate for reconstruction (matches the MATLAB simple path)
    RE[q] <- U %*% (matrix(b, nrow = Nw, ncol = Nx) %*% matrix(X[1, ], ncol = 1))
    RE <- dcm_unvec(RE, D$M$pE)

    red <- dcm_log_evidence_reduce(dcm_vec(D$Ep), D$Cp, dcm_vec(D$M$pE),
                                   D_pC, dcm_vec(RE), RC)

    D$M$pE <- RE
    D$M$pC <- RC
    D$Ep   <- dcm_unvec(red$sE, D$Ep)
    D$Cp   <- red$sC
    # share of complexity (simple split, as in the MATLAB code)
    Fb <- as.numeric(t(b) %*% bP %*% b)
    Fg <- as.numeric(t(g) %*% dcm_inv(gC) %*% g)
    Fc <- Fb / 2 + Fg / 2
    D$F    <- as.numeric(red$F + iF[i] - Fc / Ns)
    P_out[[i]] <- D

    nm <- NA_character_
    if (!is.null(D$name) && length(D$name) >= 1 && nzchar(as.character(D$name[1]))) {
      nm <- as.character(D$name[1])
    } else if (!is.null(D$Snames) && length(D$Snames) >= 1 &&
               nzchar(as.character(D$Snames[1]))) {
      nm <- paste0("PEB_", i, ":", as.character(D$Snames[1]))
    } else {
      nm <- sprintf("Subject %d", i)
    }
    PEB$Snames[i] <- nm
  }

  list(PEB = PEB, P = P_out)
}

#' Build a group-level (between-subject) design matrix
#'
#' @param n_subj Number of subjects (rows).
#' @param covariates Optional data frame of between-subject covariates with
#'   \code{n_subj} rows. \code{NULL} gives an intercept-only design.
#' @param Xnames Optional character vector of column names for the design.
#' @return A numeric design matrix with \code{n_subj} rows.
#' @examples
#' # Intercept only: tests the group mean
#' dcm_peb_design(4)
#'
#' # With a between-subject covariate
#' dcm_peb_design(4, covariates = data.frame(age = c(21, 34, 46, 58)))
#' @keywords internal
#' @export
dcm_peb_design <- function(n_subj, covariates = NULL, Xnames = NULL) {
  if (is.null(covariates)) {
    X <- matrix(1, nrow = n_subj, ncol = 1)
    colnames(X) <- "Intercept"
  } else {
    stopifnot(is.data.frame(covariates), nrow(covariates) == n_subj)
    # model.matrix includes an intercept by default
    X <- stats::model.matrix(~ ., data = covariates)
  }
  if (!is.null(Xnames)) colnames(X) <- Xnames
  X
}

#' Run a third-level PEB-of-PEBs over a directory of subject PEBs
#'
#' Loads subject PEBs (\code{.rds} or \code{.mat}), builds the third-level
#' design, runs the PEB-of-PEBs, and returns a tidy summary. Optionally saves
#' the group PEB.
#'
#' @param peb_dir Directory containing the subject PEB files.
#' @param subjects Optional integer vector; restrict to these subject numbers.
#' @param covariates Optional data frame of between-subject covariates, one row
#'   per subject in the order returned by \code{\link{dcm_peb_files}}.
#'   \code{NULL} gives an intercept-only design (the group mean).
#' @param Xnames Optional custom design column names.
#' @param save_path Optional path to save the group PEB to, as \code{.rds}.
#'   \code{NULL} (the default) does not write anything to disk.
#' @param M Third-level model options passed to \code{\link{dcm_peb_prepare}}.
#' @param field Parameter blocks to model; \code{"all"} is typical for PEB
#'   inputs.
#' @param verbose Logical. Report progress via \code{message()}.
#' @return List with \code{group_PEB}, \code{save_path} (\code{NULL} if not
#'   saved), \code{subjects}, \code{X} (the design), and \code{summary} (a tidy
#'   data frame of parameter estimates).
#' @details
#' Nothing is written to disk unless \code{save_path} is supplied.
#' @seealso \code{\link{dcm_peb_prepare}}, \code{\link{dcm_peb_run}}.
#' @export
dcm_peb_of_pebs <- function(peb_dir,
                            subjects   = NULL,
                            covariates = NULL,
                            Xnames     = NULL,
                            save_path  = NULL,
                            M     = list(Q = "single", alpha = 1, beta = 16,
                                         maxit = 64),
                            field = "all",
                            verbose = TRUE) {

  # 1) Load subject PEBs (.rds or .mat, normalised to native shape)
  L      <- dcm_peb_load(peb_dir, subjects = subjects)
  # Carry the file-derived subject label onto each PEB so that dcm_peb_run
  # reports "sub001..." in PEB$Snames rather than falling through to the
  # inner PEB's own first-subject name (which makes every input look alike).
  pebs   <- Map(function(peb, nm) { if (is.null(peb$name)) peb$name <- nm; peb },
                L$pebs, names(L$pebs))
  pebs   <- unname(pebs)
  subj   <- L$subjects
  n_subj <- length(pebs)
  if (verbose) message(sprintf("Loaded %d subject PEBs.", n_subj))

  # 2) Third-level design
  X <- dcm_peb_design(n_subj, covariates = covariates, Xnames = Xnames)
  if (verbose) {
    message("Third-level design X (rows = subjects, cols = covariates):")
    message(paste(utils::capture.output(print(X)), collapse = "\n"))
  }

  # 3) Third-level model options
  M2 <- utils::modifyList(list(X = X), M)

  # 4) Prepare and run the third-level PEB-of-PEBs
  prep3 <- dcm_peb_prepare(P = pebs, M = M2, field = field)
  run3  <- dcm_peb_run(prep3, verbose = verbose)

  # 5) Save group PEB only if asked
  if (!is.null(save_path)) {
    dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
    saveRDS(run3$PEB, save_path)
    if (verbose) message(sprintf("Saved group PEB: %s", save_path))
  }

  # 6) Tidy summary
  Ep <- as.numeric(run3$PEB$Ep)
  Eh <- if (length(run3$PEB$Eh)) run3$PEB$Eh else NA_real_
  F  <- run3$PEB$F
  Pnames <- run3$PEB$Pnames
  if (!is.null(dim(Pnames))) Pnames <- as.vector(Pnames)

  summary_df <- data.frame(
    param = if (length(Pnames)) Pnames else paste0("P", seq_along(Ep)),
    Ep    = Ep,
    Eh    = if (length(Eh)) Eh[1] else NA_real_,
    F     = F,
    stringsAsFactors = FALSE
  )

  list(
    group_PEB = run3$PEB,
    save_path = save_path,
    subjects  = subj,
    X         = X,
    summary   = summary_df
  )
}

# ---- bundled DCM example ----

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
