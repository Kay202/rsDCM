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
