group_fixture <- function(N = 12L, p = 2L, seed = 1L) {
  set.seed(seed)
  X <- cbind(intercept = 1, group = rep(c(0, 1), length.out = N))
  beta <- matrix(c(0.6, 0, 0, 0.8), p, 2)
  eta <- X %*% t(beta) + matrix(rnorm(N * p, sd = 0.2), N, p)
  colnames(eta) <- paste0("P", seq_len(p))
  Cp <- replicate(N, diag(0.1, p), simplify = FALSE)
  V <- lapply(seq_len(p), function(j) {
    z <- matrix(0, p, p)
    z[j, j] <- 1
    z
  })
  list(eta = eta, Cp = Cp, X = X, V = V, beta = beta)
}

quick_group_fit <- function(x, verbose = FALSE) {
  rsdcm(
    x$eta, x$Cp, x$X, x$V,
    n_starts = 1L,
    max_iter = 60L,
    tol = 1e-4,
    gradient_tol = 1e-3,
    pip_mcmc_iter = 100L,
    pip_mcmc_burn = 40L,
    pip_mcmc_thin = 3L,
    pip_mcmc_chains = 2L,
    pip_mcmc_temperatures = 3L,
    pip_mcmc_min_inverse_temp = 0.2,
    start_seed = 1L,
    pip_mcmc_seed = 2L,
    verbose = verbose
  )
}

test_that("rsdcm returns integrated PIPs and revised diagnostics", {
  x <- group_fixture()
  fit <- quick_group_fit(x)

  expect_equal(dim(fit$beta_mat), c(2L, 2L))
  expect_equal(dim(fit$beta_sd_approx), c(2L, 2L))
  expect_equal(dim(fit$inclusion), c(2L, 2L))
  expect_equal(dim(fit$inclusion_mcse), c(2L, 2L))
  expect_true(all(is.finite(fit$beta_mat)))
  expect_true(all(is.finite(fit$beta_sd_approx)))
  expect_true(all(fit$beta_sd_approx >= 0))
  expect_true(all(fit$inclusion >= 0 & fit$inclusion <= 1))
  expect_true(all(fit$inclusion_mcse >= 0))
  expect_true(is.finite(fit$pip_max_chain_range))
  expect_true(is.logical(fit$objective_monotone))
  expect_true(is.logical(fit$q_monotone))
  expect_match(fit$algorithm, "generalized EM")
  expect_identical(fit$scale_identification, "sigma2 fixed")
})

test_that("the compiled parallel-tempering sampler is registered", {
  expect_true(is.loaded("_rsDCM_mom_parallel_tempering_cpp"))
})

test_that("verbose output reports generalized EM and PIP diagnostics", {
  x <- group_fixture(N = 8L)
  expect_output(
    quick_group_fit(x, verbose = TRUE),
    "Starting generalized EM"
  )
})

test_that("linearly dependent covariance bases are rejected", {
  x <- group_fixture()
  expect_error(
    rsdcm(x$eta, x$Cp, x$X, c(x$V, x$V[1])),
    "linearly dependent"
  )
})
