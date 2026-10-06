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
