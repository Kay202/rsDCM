# =============================================================================
# Smoke tests — fast (< 5 s), runnable on CRAN. They check that the building
# blocks return objects of the right shape, without running a full estimation.
# =============================================================================

test_that("dcm_inv handles diagonal and dense matrices", {
  D  <- diag(c(2, 4, 8))
  Di <- as.matrix(dcm_inv(D))
  expect_equal(diag(Di), 1 / c(2, 4, 8), tolerance = 1e-8)

  M  <- matrix(c(2, 1, 1, 3), 2, 2)
  Mi <- as.matrix(dcm_inv(M))
  expect_equal(M %*% Mi, diag(2), tolerance = 1e-6)
})

test_that("dcm_logdet matches log(det(.)) for a positive-definite matrix", {
  set.seed(1)
  A <- matrix(rnorm(16), 4, 4); S <- A %*% t(A) + diag(4)
  expect_equal(dcm_logdet(S), as.numeric(determinant(S)$modulus),
               tolerance = 1e-8)
})

test_that("dcm_logdet returns 0 for an empty all-zero matrix", {
  Z <- matrix(0, 3, 3)
  expect_equal(dcm_logdet(Z), 0)
})

test_that("dcm_phi is the standard logistic", {
  expect_equal(dcm_phi(0), 0.5)
  expect_lt(abs(dcm_phi(10) - 1), 1e-3)
  expect_lt(abs(dcm_phi(-10)),     1e-3)
})

test_that("dcm_fmri_priors returns expected shapes", {
  n <- 2L
  A <- matrix(1, n, n); B <- array(0, dim = c(n, n, 1))
  C <- matrix(c(1, 0), n, 1); D <- array(0, dim = c(n, n, 0))
  pri <- dcm_fmri_priors(A, B, C, D, options = list())

  expect_named(pri, c("pE", "x", "pC"), ignore.order = TRUE)
  expect_equal(dim(pri$pE$A), c(n, n))
  expect_equal(dim(pri$x), c(n, 5))
  expect_equal(nrow(pri$pC), dcm_length(pri$pE))
  expect_equal(ncol(pri$pC), dcm_length(pri$pE))
})

test_that("dcm_fx_fmri returns same-shape state derivative", {
  n <- 2L
  pri <- dcm_fmri_priors(matrix(1, n, n), array(0, c(n, n, 1)),
                             matrix(c(1, 0), n, 1), array(0, c(n, n, 0)),
                             options = list())
  M <- list(); P <- pri$pE; x <- pri$x; u <- 0
  f <- dcm_fx_fmri(x, u, P, M)
  expect_equal(dim(f), dim(x))
  expect_true(all(is.finite(f)))
})

test_that("dcm_gx_fmri returns BOLD prediction of the right length", {
  n <- 2L
  pri <- dcm_fmri_priors(matrix(1, n, n), array(0, c(n, n, 1)),
                             matrix(c(1, 0), n, 1), array(0, c(n, n, 0)),
                             options = list())
  out <- dcm_gx_fmri(pri$x, 0, pri$pE, list())
  expect_named(out, c("g", "dgdx"), ignore.order = TRUE)
  expect_equal(length(out$g), n)
  expect_true(all(is.finite(out$g)))
})

test_that("rsdcm_options round-trips a value", {
  old <- rsdcm_options()$GLOBAL_DX
  rsdcm_options(GLOBAL_DX = 1e-6)
  expect_equal(rsdcm_options()$GLOBAL_DX, 1e-6)
  rsdcm_options(GLOBAL_DX = old)
  expect_equal(rsdcm_options()$GLOBAL_DX, old)
})
