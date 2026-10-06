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
