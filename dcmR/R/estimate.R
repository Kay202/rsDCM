# =============================================================================
# Top-level DCM estimation entry point and Bayesian model evidence helpers.
# =============================================================================

#' Estimate a Dynamic Causal Model for fMRI
#'
#' Inverts a fully-specified DCM using variational Laplace inversion. This
#' is the main user-facing entry point of the package. Mirrors SPM12's
#' \code{dcm_estimate}, with the following CRAN-friendly behavioural
#' changes: progress messages are emitted via \code{message()} (silenceable
#' with \code{suppressMessages()}), and writing the result back to disk is
#' opt-in via the \code{save} argument.
#'
#' @param P A DCM list, or a path to an \code{.RData}/\code{.rds} file
#'   containing one.
#' @param save Logical. If \code{TRUE} and \code{P} is a file path, the
#'   estimated DCM is saved back to that path. Defaults to \code{FALSE} to
#'   comply with CRAN policy on writing to user files without consent.
#' @return The estimated DCM (a list) with posterior fields populated.
#' @examples
#' \dontrun{
#'   data(toy_dcm)
#'   fit <- dcm_estimate(toy_dcm)
#'   round(fit$Ep$A, 3)
#' }
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
  sw   <- options(warn = -1); on.exit(options(sw), add = TRUE)
  Pp_vec <- 1 - dcm_Ncdf(Tvec, abs(as.numeric(dcm_vec(Ep))), diag(Cp))
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
#' replacing the original prior with a reduced prior. Mirrors SPM12's
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
#' Reduced-rank version of \code{\link{dcm_log_evidence}}. Mirrors SPM12's
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
#' Mirrors SPM12's \code{spm_dcm_evidence}.
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
