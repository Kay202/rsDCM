# =============================================================================
# Full estimation smoke test on the toy dataset. Skipped on CRAN because it
# runs the Gauss-Newton inner loop (a few seconds).
# =============================================================================

test_that("dcm_estimate runs end-to-end on toy_dcm", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("rsDCM")

  data(toy_dcm, package = "rsDCM")

  fit <- suppressMessages(dcm_estimate(toy_dcm))

  # Required output fields
  for (nm in c("Ep", "Cp", "F", "y", "R", "ID")) {
    expect_true(nm %in% names(fit), info = paste("missing field:", nm))
  }

  # Posterior expectation has the SPM-style structured field names
  expect_true(all(c("A", "B", "C", "transit", "decay", "epsilon") %in%
                    names(fit$Ep)))
  expect_equal(dim(fit$Ep$A), c(toy_dcm$n, toy_dcm$n))
  expect_true(all(is.finite(dcm_vec(fit$Ep))))

  # Posterior covariance is square and matches parameter length
  expect_equal(nrow(fit$Cp), dcm_length(fit$Ep))
  expect_equal(ncol(fit$Cp), dcm_length(fit$Ep))

  # Free energy is finite
  expect_true(is.finite(fit$F))

  # Predicted BOLD has the data shape
  expect_equal(dim(fit$y), dim(toy_dcm$Y$y))
})

# Helper: tolerate missing %||% in older R
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
