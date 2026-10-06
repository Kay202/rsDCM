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
