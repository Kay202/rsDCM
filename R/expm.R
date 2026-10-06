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
