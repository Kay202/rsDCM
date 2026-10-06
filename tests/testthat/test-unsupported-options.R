# =============================================================================
# dcm_estimate must fail loudly on model variants that are not yet implemented
# (two_state, stochastic, induced), rather than silently running the
# deterministic model. These checks are fast: the error is raised during option
# validation, before any inversion runs, so no Gauss-Newton loop executes.
# =============================================================================

test_that("dcm_estimate rejects unsupported variants set as numeric 1", {
  data(toy_dcm, package = "rsDCM")

  for (opt in c("two_state", "stochastic", "induced")) {
    d <- toy_dcm
    d$options[[opt]] <- 1                      # numeric flag, as DCMs store it
    expect_error(dcm_estimate(d),
                 regexp = "not yet implemented",
                 info   = paste("option:", opt))
    # The message should name the offending option and point to the future.
    expect_error(dcm_estimate(d), regexp = opt, info = opt)
    expect_error(dcm_estimate(d), regexp = "future update", info = opt)
  }
})

test_that("the guard also catches logical TRUE, not just numeric 1", {
  data(toy_dcm, package = "rsDCM")
  d <- toy_dcm
  d$options$two_state <- TRUE
  expect_error(dcm_estimate(d), "not yet implemented")
})

test_that("explicit zero / FALSE does not trip the guard", {
  # We can't run a full (75 s) inversion here, so just confirm that option
  # validation passes for the deterministic settings and execution proceeds
  # past the guard. A deliberately malformed DCM makes it fail *later* (not on
  # the unsupported-option check), which is what we assert.
  d <- list(options = list(two_state = 0, stochastic = FALSE, induced = 0))
  err <- tryCatch(dcm_estimate(d), error = function(e) conditionMessage(e))
  expect_false(grepl("not yet implemented", err))
  expect_false(grepl("future update", err))
})

test_that("a resting-state DCM (null inputs) is not blocked by the guard", {
  # dcm_estimate internally sets stochastic <- 1 for null-input data; that
  # internal default must NOT trigger the user-facing guard, which only checks
  # the options the caller supplied. Again we only need to get *past* the
  # guard, so we assert the error (from the incomplete DCM) is not the guard's.
  d <- list(options = list())                  # nothing switched on by the user
  err <- tryCatch(dcm_estimate(d), error = function(e) conditionMessage(e))
  expect_false(grepl("not yet implemented", err))
})
