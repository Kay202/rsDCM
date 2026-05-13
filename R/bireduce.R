# =============================================================================
# Bilinear reduction (M0/M1 form) and Volterra kernels (internal).
# =============================================================================

#' Bilinear reduction of a non-linear DCM
#'
#' Reduces a non-linear state equation to bilinear form \code{dx/dt = M0 x +
#' sum_i u_i M1[[i]] x}, also returning the lead-field expansion. Mirrors
#' SPM12's \code{spm_bireduce}.
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
#' system specified in M0/M1 form. Mirrors SPM12's \code{spm_kernels}.
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
