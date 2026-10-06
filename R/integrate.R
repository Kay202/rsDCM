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
