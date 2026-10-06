# Robust + sparse group DCM: Student-t weighting, pMOM spike-and-slab prior, ReML.
#
# A group-level model for DCM parameters across subjects that is robust to
# outlier subjects (Student-t weighting) and performs Bayesian variable
# selection on the group effects (a nonlocal product-moment (pMOM)
# spike-and-slab prior),
# with ReML-estimated between-subject variance components. This is an
# alternative to the Gaussian PEB layer (dcm_peb_*); the numerics are the
# author's own, not an SPM port.

# pMOM spike-and-slab prior: normalized spike N(0, tau0^2) and first-order pMOM
# slab, their mixture m, and the first/second derivatives of log m in beta.
# base::pi is written explicitly because the `pi` argument (inclusion
# probability) shadows the constant 3.14159... inside this function.
.mom_prior_components <- function(beta, tau0, tau1, pi) {
  f0 <- (2 * base::pi * tau0^2)^(-0.5) * exp(-0.5 * (beta^2 / tau0^2))
  f1 <- (2 * base::pi)^(-0.5) * tau1^(-3) * (beta^2) * exp(-0.5 * (beta^2 / tau1^2))

  m  <- (1 - pi) * f0 + pi * f1

  f0_prime <- -(beta / tau0^2) * f0

  # Numerically robust treatment near zero: sign(0) == 0 would zero the
  # denominator, so map exact zeros (and sub-eps magnitudes) to +/- eps.
  eps <- 1e-8
  beta_safe <- ifelse(abs(beta) < eps, ifelse(beta >= 0, eps, -eps), beta)

  f1_prime <- f1 * (2 / beta_safe - beta_safe / tau1^2)
  m_prime  <- (1 - pi) * f0_prime + pi * f1_prime

  f0_second    <- ((beta^2 / tau0^4) - (1 / tau0^2)) * f0
  logf1_second <- -2 / (beta_safe^2) - 1 / tau1^2
  logf1_prime  <-  2 / beta_safe - beta_safe / tau1^2
  f1_second    <- f1 * (logf1_prime^2 + logf1_second)
  m_second     <- (1 - pi) * f0_second + pi * f1_second

  logm_prime  <- m_prime / (m + 1e-16)
  logm_second <- (m_second * m - m_prime^2) / (m^2 + 1e-16)

  list(f0 = f0, f1 = f1, m = m,
       logm_prime = logm_prime, logm_second = logm_second)
}

# Update the spike (tau0) and slab (tau1) scales given the current beta and
# Inv-Gamma hyperpriors. Works on the log(tau^2) scale with a 1D optimize().
.mom_update_tau <- function(beta, pi, tau0, tau1, a0, b0, a1, b1,
                            log_tau2_lower = -20, log_tau2_upper = 5) {
  tiny <- 1e-16

  f_tau0 <- function(log_tau2) {
    tau2 <- exp(log_tau2)
    f0 <- (2 * base::pi * tau2)^(-0.5) * exp(-0.5 * (beta^2 / tau2))
    f1 <- (2 * base::pi)^(-0.5) * (tau1^2)^(-1.5) * (beta^2) *
      exp(-0.5 * (beta^2 / (tau1^2)))
    m  <- (1 - pi) * f0 + pi * f1
    loglik    <- sum(log(m + tiny))
    log_prior <- a0 * log(b0) - lgamma(a0) - (a0 + 1) * log(tau2) - b0 / tau2
    -(loglik + log_prior)
  }
  opt0 <- stats::optimize(f_tau0, lower = log_tau2_lower, upper = log_tau2_upper)
  tau0_new <- sqrt(exp(opt0$minimum))

  f_tau1 <- function(log_tau2) {
    tau2 <- exp(log_tau2)
    f0 <- (2 * base::pi * (tau0^2))^(-0.5) * exp(-0.5 * (beta^2 / (tau0^2)))
    f1 <- (2 * base::pi)^(-0.5) * (tau2)^(-1.5) * (beta^2) *
      exp(-0.5 * (beta^2 / tau2))
    m  <- (1 - pi) * f0 + pi * f1
    loglik    <- sum(log(m + tiny))
    log_prior <- a1 * log(b1) - lgamma(a1) - (a1 + 1) * log(tau2) - b1 / tau2
    -(loglik + log_prior)
  }
  opt1 <- stats::optimize(f_tau1, lower = log_tau2_lower, upper = log_tau2_upper)
  tau1_new <- sqrt(exp(opt1$minimum))

  list(tau0 = tau0_new, tau1 = tau1_new)
}

# Build Sigma_b(alpha) = sum_k alpha_k V_k  (each V_k is p x p).
.mom_build_sigma_b <- function(alpha, V_list, p) {
  K <- length(V_list)
  if (is.null(alpha) || K == 0) return(matrix(0, nrow = p, ncol = p))
  if (length(alpha) != K) stop("Length of alpha must match length of V_list.")
  Sigma_b <- matrix(0, nrow = p, ncol = p)
  for (k in seq_len(K)) {
    if (alpha[k] != 0) Sigma_b <- Sigma_b + alpha[k] * V_list[[k]]
  }
  Sigma_b
}

# Pre-whiten. For each subject n: M_n = C_n + Sigma_b, S_n = solve(M_n),
# L_n = chol(S_n) so L_n^T L_n = S_n. Then y_star = B(alpha) vec(eta) and
# X_star = B(alpha) (X_G kron I_p), giving an identity-error whitened model.
.mom_whiten <- function(eta_theta_y, C_theta_y_list, X_G, Sigma_b) {
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G <- as.matrix(X_G)
  N <- nrow(eta_theta_y); p <- ncol(eta_theta_y); r <- ncol(X_G)

  if (length(C_theta_y_list) != N)
    stop("C_theta_y_list must be a list of length N with p x p matrices.")

  Np <- N * p; q <- p * r
  y_star <- numeric(Np)
  X_star <- matrix(0, nrow = Np, ncol = q)
  S_list <- vector("list", N)

  row_start <- 1L
  for (n in seq_len(N)) {
    Cn <- C_theta_y_list[[n]]
    if (!all(dim(Cn) == c(p, p)))
      stop("Each element of C_theta_y_list must be p x p.")

    M_n <- Cn + Sigma_b + 1e-8 * diag(p)   # jitter for stability
    S_n <- solve(M_n)
    L_n <- chol(S_n)                        # upper-triangular, L_n^T L_n = S_n
    S_list[[n]] <- S_n

    idx <- row_start:(row_start + p - 1L)
    Xn_tilde <- kronecker(matrix(X_G[n, ], nrow = 1), diag(p))  # p x (p*r)
    y_star[idx]   <- as.numeric(L_n %*% eta_theta_y[n, ])
    X_star[idx, ] <- L_n %*% Xn_tilde
    row_start <- row_start + p
  }

  list(y_star = y_star, X_star = X_star, S_list = S_list)
}

# ReML gradient/Hessian for the variance components alpha, plus the separated
# trace and quadratic terms Tr_k = sum_n tr(S_n V_k) and
# Qd_k = sum_n e_n^T S_n V_k S_n e_n that drive the multiplicative update.
.mom_alpha_grad_hess <- function(S_list, V_list, X_G, beta_vec, eta_theta_y) {
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G <- as.matrix(X_G)
  N <- nrow(eta_theta_y); p <- ncol(eta_theta_y); r <- ncol(X_G)
  q <- length(beta_vec)
  if (q != p * r) stop("Length of beta_vec must be p * r.")

  K <- length(V_list)
  eta_vec <- as.vector(t(eta_theta_y))
  X_big   <- kronecker(X_G, diag(p))
  e_vec   <- eta_vec - as.numeric(X_big %*% beta_vec)

  g  <- numeric(K); Tr <- numeric(K); Qd <- numeric(K)
  H  <- matrix(0, nrow = K, ncol = K)

  offset <- 0L
  for (n in seq_len(N)) {
    S_n <- S_list[[n]]
    idx <- (offset + 1L):(offset + p)
    e_n <- e_vec[idx]

    for (k in seq_len(K)) {
      tr_k  <- sum(diag(S_n %*% V_list[[k]]))
      Tr[k] <- Tr[k] + tr_k
      g[k]  <- g[k] + 0.5 * tr_k
    }
    for (k in seq_len(K)) {
      SK <- S_n %*% V_list[[k]] %*% S_n
      q_k   <- as.numeric(t(e_n) %*% SK %*% e_n)
      Qd[k] <- Qd[k] + q_k
      g[k]  <- g[k] - 0.5 * q_k
      for (l in k:K) {
        H_add <- 0.5 * sum(diag(SK %*% V_list[[l]] %*% S_n))
        H[k, l] <- H[k, l] + H_add
        if (l != k) H[l, k] <- H[l, k] + H_add
      }
    }
    offset <- offset + p
  }

  list(g = g, H = H, Tr = Tr, Qd = Qd)
}

#' Robust and sparse group-level DCM (Student-t + pMOM)
#'
#' Fits a group-level model to subject-level DCM parameter estimates that is
#' robust to outlier subjects and performs Bayesian variable selection on the
#' group effects. The model combines three ingredients: Student-t weighting of
#' subjects (robustness), a nonlocal product-moment (pMOM) spike-and-slab prior
#' on the group effects (sparsity / inclusion probabilities), and
#' ReML-estimated between-subject variance components. Estimation is by EM.
#'
#' This is an alternative to the Gaussian Parametric Empirical Bayes layer
#' (\code{\link{dcm_peb_run}}); the numerics are original to the package, not a
#' port of SPM. Most users will call the convenience wrapper
#' \code{\link{rsdcm_fit}}, which assembles the arguments below from a
#' list of fitted DCMs.
#'
#' @param eta_theta_y N x p matrix of subject-level posterior means (subjects
#'   in rows, parameters in columns).
#' @param C_theta_y_list Length-N list of p x p subject-level posterior
#'   covariance matrices.
#' @param X_G N x r group (between-subject) design matrix.
#' @param V_list List of p x p basis matrices for the between-subject
#'   covariance \code{Sigma_b = sum_k alpha_k V_k}. A per-parameter diagonal
#'   basis is a common default.
#' @param nu Student-t degrees of freedom (smaller = heavier-tailed, more robust).
#' @param tau0,tau1 Initial spike and slab standard deviations.
#' @param pi Prior inclusion probability (slab weight).
#' @param a0,b0,a1,b1 Inverse-Gamma hyperprior parameters on \code{tau0^2} and
#'   \code{tau1^2}.
#' @param min_slab_spike_ratio Identifiability guard: enforce
#'   \code{tau1 >= min_slab_spike_ratio * tau0}.
#' @param max_iter Maximum EM iterations.
#' @param tol Convergence tolerance on the change in \code{beta} and
#'   \code{alpha}.
#' @param inner_sweeps Coordinate-Newton sweeps per M-step for \code{beta}.
#' @param verbose Logical. Report progress via \code{message()}; silence with
#'   \code{suppressMessages()} or \code{verbose = FALSE}.
#' @param max_beta_step Per-coordinate cap on the \code{beta} Newton step.
#' @param min_alpha,max_alpha Bounds on the variance components.
#' @return A list with the group-level effects \code{beta_mat} (p x r) and
#'   \code{beta_vec}, posterior inclusion probabilities \code{inclusion}
#'   (p x r), Student-t \code{weights}, variance components \code{alpha},
#'   learned scales \code{tau0}/\code{tau1}, residual variance \code{sigma2},
#'   fitted means \code{mu_hat} (N x p), and the sizes \code{N}, \code{p},
#'   \code{r}.
#' @examples
#' # Real 48-subject NARPS group analysis. The inputs are precomputed
#' # subject-level DCM summaries (a posterior mean and covariance per subject).
#' data(narps_dcm)
#' N <- nrow(narps_dcm$eta)
#' p <- ncol(narps_dcm$eta)
#' # Per-parameter variance-component basis, and an intercept-only design.
#' V_list <- lapply(seq_len(p), function(k) { V <- matrix(0, p, p); V[k, k] <- 1; V })
#' X_G <- matrix(1, N, 1, dimnames = list(NULL, "intercept"))
#'
#' \donttest{
#' fit <- rsdcm(narps_dcm$eta, narps_dcm$Cp, X_G, V_list, verbose = FALSE)
#' rownames(fit$beta_mat) <- rownames(fit$inclusion) <- narps_dcm$parameter_names
#'
#' # Group-level connections selected by the spike-and-slab prior (PIP > 0.5)
#' sel <- fit$inclusion[, 1] > 0.5
#' round(cbind(estimate = fit$beta_mat[sel, 1],
#'             PIP      = fit$inclusion[sel, 1]), 3)
#'
#' # Subjects the Student-t weighting down-weights most
#' w <- rowMeans(matrix(fit$weights, N, p, byrow = TRUE))
#' narps_dcm$subject[order(w)][1:5]
#' }
#' @seealso \code{\link{rsdcm_fit}} for the wrapper over fitted DCMs,
#'   \code{\link{dcm_peb_run}} for the Gaussian PEB alternative,
#'   \code{\link{narps_dcm}} for the example dataset.
#' @references
#' Arhin, G., Sanyal, N. (2026). Robust and sparse group dynamic causal
#' modeling via Student-t parametric empirical Bayes and nonlocal priors.
#' arXiv:2609.06379. \doi{10.48550/arXiv.2609.06379}
#' @export
rsdcm <- function(eta_theta_y, C_theta_y_list, X_G, V_list,
                            nu = 3, tau0 = 0.05, tau1 = 1.0, pi = 0.5,
                            a0 = 2, b0 = 1e-2, a1 = 2, b1 = 1.0,
                            min_slab_spike_ratio = 1,
                            max_iter = 500, tol = 1e-6, inner_sweeps = 3,
                            verbose = TRUE, max_beta_step = 1.0,
                            min_alpha = 1e-8, max_alpha = 100) {
  eta_theta_y <- as.matrix(eta_theta_y)
  X_G <- as.matrix(X_G)
  if (is.null(colnames(X_G)))
    colnames(X_G) <- if (ncol(X_G) == 1) "intercept"
                     else paste0("col", seq_len(ncol(X_G)))
  N <- nrow(eta_theta_y); p <- ncol(eta_theta_y); r <- ncol(X_G)
  q <- p * r; Np <- N * p

  if (nrow(X_G) != N)
    stop("Number of rows of X_G must match the number of subjects in eta_theta_y.")
  if (length(C_theta_y_list) != N)
    stop("C_theta_y_list must be a list of length N.")

  K <- length(V_list)
  alpha <- if (K > 0) rep(0.1, K) else numeric(0)

  # Initial beta via a simple ridge-regularised GLS step.
  X_big <- kronecker(X_G, diag(p))
  y_big <- as.vector(t(eta_theta_y))
  beta  <- solve(crossprod(X_big) + 1e-6 * diag(q), crossprod(X_big, y_big))
  sigma2 <- stats::var(y_big - as.numeric(X_big %*% beta))

  if (verbose)
    message("Starting EM for robust + sparse group DCM (Student-t + pMOM + ReML alpha)")

  for (iter in seq_len(max_iter)) {
    beta_old <- beta; alpha_old <- alpha

    # Pre-whiten for the current alpha.
    Sigma_b <- .mom_build_sigma_b(alpha, V_list, p)
    wd <- .mom_whiten(eta_theta_y, C_theta_y_list, X_G, Sigma_b)
    y_star <- wd$y_star; X_star <- wd$X_star; S_list <- wd$S_list

    # E-step: Student-t weights on whitened residuals.
    r_star <- as.numeric(y_star - X_star %*% beta)
    W <- (nu + 1) / (nu + (r_star^2) / (sigma2 + 1e-12))

    # M-step: beta by coordinate-Newton with the pMOM prior.
    for (sweep in seq_len(inner_sweeps)) {
      r_star <- as.numeric(y_star - X_star %*% beta)
      g_lik  <- as.numeric(crossprod(X_star, W * r_star)) / sigma2
      h_lik_diag <- numeric(q)
      for (j in seq_len(q)) h_lik_diag[j] <- -sum(W * X_star[, j]^2) / sigma2

      for (j in seq_len(q)) {
        bj <- beta[j]
        pc <- .mom_prior_components(beta = bj, tau0 = tau0, tau1 = tau1, pi = pi)
        g_j <- g_lik[j] + pc$logm_prime
        h_j <- h_lik_diag[j] + pc$logm_second
        # Ensure a valid ascent direction; near the pMOM singularity at 0 the
        # combined curvature can be non-negative, so fall back to the
        # always-negative likelihood curvature rather than freeze the coord.
        if (!is.finite(h_j) || h_j > -1e-8) h_j <- -abs(h_lik_diag[j]) - 1e-8
        step_j <- -g_j / h_j
        if (is.finite(step_j)) {
          if (abs(step_j) > max_beta_step) step_j <- sign(step_j) * max_beta_step
          beta[j] <- bj + step_j
        }
      }
    }

    # M-step: sigma^2 (MAP under a Jeffreys prior, hence Np + 2).
    r_star <- as.numeric(y_star - X_star %*% beta)
    sigma2 <- sum(W * r_star^2) / (Np + 2)

    # M-step: multiplicative ReML/EM update for alpha. alpha_k <- alpha_k * Q_k / T_k
    # is monotone, keeps alpha_k >= 0, and avoids boundary oscillation; the
    # per-iteration factor is damped to a bounded range.
    if (K > 0) {
      ah <- .mom_alpha_grad_hess(S_list, V_list, X_G, beta_vec = beta,
                                 eta_theta_y = eta_theta_y)
      ratio <- pmin(pmax(ah$Qd / (ah$Tr + 1e-12), 0.2), 5)
      alpha <- alpha * ratio
      alpha[alpha < min_alpha] <- min_alpha
      alpha[alpha > max_alpha] <- max_alpha
    }

    # M-step: spike/slab scales, with the slab kept at least as wide as the spike.
    ts <- .mom_update_tau(beta, pi, tau0, tau1, a0, b0, a1, b1)
    tau0 <- ts$tau0
    tau1 <- max(ts$tau1, min_slab_spike_ratio * ts$tau0)

    d_beta  <- max(abs(beta - beta_old))
    d_alpha <- if (K > 0) max(abs(alpha - alpha_old)) else 0
    if (verbose && (iter %% 20 == 0L || iter == 1L)) {
      message(sprintf(
        "Iter %3d: max|d.beta| = %.3e, max|d.alpha| = %.3e, sigma2 = %.3f, tau0 = %.4f, tau1 = %.4f",
        iter, d_beta, d_alpha, sigma2, tau0, tau1))
    }
    if (max(d_beta, d_alpha) < tol) {
      if (verbose) message("Converged at iter ", iter)
      break
    }
  }

  # Posterior inclusion probabilities from the final beta.
  pc <- .mom_prior_components(beta = beta, tau0 = tau0, tau1 = tau1, pi = pi)
  inclusion_vec <- (pi * pc$f1) / (pc$m + 1e-16)

  # beta matches kronecker(X_G, diag(p)) (regressor-major, parameter-minor),
  # so the p x r reshape is column-major (byrow = FALSE).
  beta_mat <- matrix(beta, nrow = p, ncol = r, byrow = FALSE)
  colnames(beta_mat) <- colnames(X_G); rownames(beta_mat) <- colnames(eta_theta_y)
  inclusion_mat <- matrix(inclusion_vec, nrow = p, ncol = r, byrow = FALSE)
  colnames(inclusion_mat) <- colnames(X_G); rownames(inclusion_mat) <- colnames(eta_theta_y)

  mu_hat_vec <- as.numeric(kronecker(X_G, diag(p)) %*% beta)
  mu_hat_mat <- t(matrix(mu_hat_vec, nrow = p, ncol = N))
  colnames(mu_hat_mat) <- colnames(eta_theta_y)

  list(beta_mat = beta_mat, beta_vec = beta, inclusion = inclusion_mat,
       sigma2 = sigma2, weights = W, alpha = alpha, tau0 = tau0, tau1 = tau1,
       mu_hat = mu_hat_mat, N = N, p = p, r = r,
       a0 = a0, b0 = b0, a1 = a1, b1 = b1)
}

#' Robust sparse group DCM from a list of fitted DCMs
#'
#' Convenience wrapper around \code{\link{rsdcm}} that assembles its
#' inputs from a list of estimated DCMs (as returned by
#' \code{\link{dcm_estimate}}). For each subject it takes the posterior mean
#' \code{Ep} and covariance \code{Cp} restricted to the requested parameter
#' field, builds a between-subject design, and fits the robust + sparse group
#' model.
#'
#' @param P List of estimated DCMs; each must carry \code{Ep}, \code{Cp}, and
#'   \code{M$pE}/\code{M$pC} (or \code{options}) so the field indices can be
#'   resolved.
#' @param field Parameter block(s) to model at the group level, named as in the
#'   DCM parameter structure \code{Ep}. For an fMRI DCM these are \code{"A"},
#'   \code{"B"}, \code{"C"}, \code{"D"}, \code{"transit"}, \code{"decay"} and
#'   \code{"epsilon"}; pass one (e.g. \code{"A"}) or several (e.g.
#'   \code{c("A", "B")}). Case-sensitive: the connectivity blocks are uppercase
#'   (\code{"A"}, not \code{"a"}). Passed to \code{\link{dcm_find_pC}}.
#' @param covariates Optional data frame of between-subject covariates, one row
#'   per subject; \code{NULL} gives an intercept-only design (the group mean).
#'   Ignored if \code{X_G} is supplied.
#' @param X_G Optional group design matrix (N x r); overrides \code{covariates}.
#' @param V_list Optional list of p x p variance-component bases; defaults to a
#'   per-parameter diagonal basis.
#' @param ... Further arguments passed to \code{\link{rsdcm}} (e.g.
#'   \code{nu}, \code{pi}, \code{verbose}).
#' @return The list returned by \code{\link{rsdcm}}, with the group
#'   effects and inclusion probabilities labelled by the selected parameters and
#'   design columns, plus the selected \code{param_index}.
#' @details
#' All subjects must share the same parameterisation, so the selected field
#' indices are required to match across \code{P}.
#' @examples
#' # rsdcm_fit() consumes a list of fitted DCMs (each from dcm_estimate), then
#' # fits the group model on a chosen field:
#' #   fits <- lapply(dcm_files, function(f) dcm_estimate(readRDS(f)))
#' #   grp  <- rsdcm_fit(fits, field = "A", covariates = my_covariates)
#' # For a runnable group-level example on precomputed subject summaries, see
#' # ?rsdcm and the narps_dcm dataset.
#' @seealso \code{\link{rsdcm}}, \code{\link{dcm_estimate}},
#'   \code{\link{dcm_peb_design}}.
#' @export
rsdcm_fit <- function(P, field = "A", covariates = NULL,
                                X_G = NULL, V_list = NULL, ...) {
  if (!is.list(P) || length(P) < 2L)
    stop("P must be a list of at least two fitted DCMs.")
  N <- length(P)

  if (!is.character(field) || !length(field))
    stop("`field` must be one or more parameter-block names (e.g. \"A\" or ",
         "c(\"A\", \"B\")). Numeric indices are not supported.")

  # Resolve the requested field(s) to free-parameter indices, validated against
  # the DCM parameter structure. An unknown or miscased name (e.g. "a" for "A")
  # is rejected here, and a field with no free parameters (e.g. "D" in a
  # deterministic DCM) is caught below, rather than silently expanding to every
  # parameter as dcm_find_pC() does on its own.
  rc1   <- dcm_find_pC(P[[1]])
  valid <- names(rc1$pE)
  if (is.null(valid))
    stop("Could not recover the parameter structure of P[[1]]; each DCM needs ",
         "M$pE (or a/b/c and options to rebuild the priors).")
  bad <- setdiff(field, valid)
  if (length(bad))
    stop("Unknown field(s): ", paste(paste0("\"", bad, "\""), collapse = ", "),
         ". Valid fields are: ", paste(paste0("\"", valid, "\""), collapse = ", "),
         " (case-sensitive).")

  field_index <- function(rc)
    intersect(suppressWarnings(
      unlist(do.call(dcm_fieldindices, c(list(rc$pE), as.list(field))))), rc$i)

  q <- field_index(rc1)
  if (!length(q))
    stop("Field ", paste(paste0("\"", field, "\""), collapse = ", "),
         " has no free parameters to model in these DCMs.")
  p <- length(q)

  eta_theta_y    <- matrix(0, nrow = N, ncol = p)
  C_theta_y_list <- vector("list", N)
  for (n in seq_len(N)) {
    qn <- if (n == 1L) q else field_index(dcm_find_pC(P[[n]]))
    if (!identical(as.integer(qn), as.integer(q)))
      stop("Subject ", n, " has a different parameterisation than subject 1; ",
           "all DCMs must match.")
    ev <- dcm_vec(P[[n]]$Ep)
    eta_theta_y[n, ] <- ev[q]
    Cp <- as.matrix(P[[n]]$Cp)
    C_theta_y_list[[n]] <- Cp[q, q, drop = FALSE]
  }

  # Column labels for the selected parameters, if available.
  labs <- tryCatch(unlist(lapply(q, function(j) dcm_fieldindices(P[[1]]$Ep, j))),
                   error = function(e) NULL)
  if (length(labs) == p) colnames(eta_theta_y) <- labs
  else colnames(eta_theta_y) <- paste0("P", q)

  # Design matrix: explicit X_G wins, else build from covariates (intercept default).
  if (is.null(X_G)) X_G <- dcm_peb_design(N, covariates = covariates)
  if (nrow(X_G) != N)
    stop("X_G / covariates must have one row per subject (", N, ").")

  # Default variance-component basis: one diagonal entry per parameter.
  if (is.null(V_list))
    V_list <- lapply(seq_len(p), function(k) { V <- matrix(0, p, p); V[k, k] <- 1; V })

  out <- rsdcm(eta_theta_y = eta_theta_y,
                         C_theta_y_list = C_theta_y_list,
                         X_G = X_G, V_list = V_list, ...)
  out$param_index <- q
  out$field <- field
  out
}
