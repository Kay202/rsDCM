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
