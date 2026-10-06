# Tests for the robust + sparse group DCM (rsdcm). These are pure
# linear algebra on simulated subject summaries (no inversion), so they are
# fast and CRAN-runnable.

# Simulate N subjects with a known p x r group effect, heavy-tailed noise, and
# two injected outliers. Returns the pieces plus the truth.
gm_fixture <- function(N = 24L, p = 4L, seed = 1L) {
  set.seed(seed)
  X_G <- cbind(intercept = 1, group = stats::rbinom(N, 1, 0.5))
  beta_true <- matrix(0, p, 2, dimnames = list(NULL, c("intercept", "group")))
  beta_true[1, "intercept"] <-  0.7
  beta_true[3, "group"]     <-  1.1
  eta <- X_G %*% t(beta_true) + matrix(stats::rnorm(N * p, sd = 0.3), N, p)
  o1 <- min(5L, N)
  o2 <- if (N >= 17L) 17L else max(1L, N - 1L)
  c1 <- min(2L, p)
  c2 <- min(4L, p)
  outliers <- c(o1, o2)
  eta[o1, c1] <- eta[o1, c1] + 6
  eta[o2, c2] <- eta[o2, c2] - 7
  C_list <- replicate(N, diag(0.2, p), simplify = FALSE)
  V_list <- lapply(seq_len(p), function(k) { V <- matrix(0, p, p); V[k, k] <- 1; V })
  list(eta = eta, C_list = C_list, X_G = X_G, V_list = V_list,
       beta_true = beta_true, outliers = outliers, N = N, p = p)
}

test_that("rsdcm recovers a planted group effect", {
  fx  <- gm_fixture()
  fit <- rsdcm(fx$eta, fx$C_list, fx$X_G, fx$V_list, verbose = FALSE)

  expect_equal(dim(fit$beta_mat), c(fx$p, 2L))
  expect_true(all(is.finite(fit$beta_mat)))
  expect_gt(stats::cor(as.vector(fit$beta_mat), as.vector(fx$beta_true)), 0.9)

  # The two true effects should dominate the estimate.
  expect_equal(which.max(fit$beta_mat[, "intercept"]), 1L)
  expect_equal(which.max(fit$beta_mat[, "group"]), 3L)
})

test_that("posterior inclusion is higher at true effects than at zeros", {
  fx  <- gm_fixture()
  fit <- rsdcm(fx$eta, fx$C_list, fx$X_G, fx$V_list, verbose = FALSE)

  nz <- fx$beta_true != 0
  expect_gt(mean(fit$inclusion[nz]), mean(fit$inclusion[!nz]))
  expect_true(all(fit$inclusion >= 0 & fit$inclusion <= 1))
})

test_that("Student-t weights down-weight the injected outliers", {
  fx  <- gm_fixture()
  fit <- rsdcm(fx$eta, fx$C_list, fx$X_G, fx$V_list, verbose = FALSE)

  # weights are in subject-major order (p per subject); average within subject.
  subj_w <- vapply(seq_len(fx$N), function(n) {
    mean(fit$weights[((n - 1L) * fx$p + 1L):(n * fx$p)])
  }, numeric(1))
  # Both outlier subjects sit below the median subject weight.
  expect_true(all(subj_w[fx$outliers] < stats::median(subj_w)))
})

test_that("verbose progress is a silenceable message()", {
  fx <- gm_fixture(N = 12L, p = 3L)
  expect_message(
    rsdcm(fx$eta, fx$C_list, fx$X_G, fx$V_list,
                    max_iter = 2, verbose = TRUE),
    "robust \\+ sparse group DCM")
  expect_silent(
    suppressMessages(rsdcm(fx$eta, fx$C_list, fx$X_G, fx$V_list,
                                     max_iter = 2, verbose = TRUE)))
})

test_that("rsdcm_fit assembles from DCM fits and matches the core", {
  set.seed(3)
  n <- 2L
  pri <- dcm_fmri_priors(matrix(1, n, n), array(0, c(n, n, 1)),
                         matrix(c(1, 0), n, 1), array(0, c(n, n, 0)),
                         options = list())
  q  <- dcm_find_pC(list(M = list(pE = pri$pE, pC = pri$pC)), "A")$i
  p  <- length(q)
  N  <- 15L
  Cp <- diag(0.1, length(dcm_vec(pri$pE)))

  # Build N synthetic "fitted DCMs" with a group effect on A.
  A_group <- matrix(c(0.4, 0, -0.5, 0), n, n)
  P <- lapply(seq_len(N), function(i) {
    Ep <- pri$pE
    Ep$A <- A_group + matrix(rnorm(n * n, sd = 0.15), n, n)
    list(M = list(pE = pri$pE, pC = pri$pC), Ep = Ep, Cp = Cp, F = -50 + rnorm(1))
  })

  fit <- rsdcm_fit(P, field = "A", verbose = FALSE)
  expect_equal(dim(fit$beta_mat), c(p, 1L))          # intercept-only design
  expect_true(all(is.finite(fit$beta_mat)))
  expect_equal(fit$param_index, q)
  expect_true(all(is.finite(fit$inclusion)))

  # Same result when the identical inputs are assembled by hand and passed to
  # the core fitter directly.
  eta <- t(vapply(P, function(d) dcm_vec(d$Ep)[q], numeric(p)))
  C_list <- replicate(N, Cp[q, q, drop = FALSE], simplify = FALSE)
  X_G <- dcm_peb_design(N)
  V_list <- lapply(seq_len(p), function(k) { V <- matrix(0, p, p); V[k, k] <- 1; V })
  core <- rsdcm(eta, C_list, X_G, V_list, verbose = FALSE)
  expect_equal(unname(fit$beta_mat), unname(core$beta_mat), tolerance = 1e-8)
})
