test_that("dcm_length is consistent with dcm_vec on numerics", {
  expect_equal(dcm_length(1:10), 10L)
  expect_equal(dcm_length(matrix(1:6, 2, 3)), 6L)
  expect_equal(dcm_length(numeric(0)), 0L)
})

test_that("dcm_vec / dcm_unvec round-trip a flat numeric", {
  x  <- 1:5
  vx <- dcm_vec(x)
  expect_equal(length(vx), 5L)
  expect_equal(dcm_unvec(vx, x), x)
})

test_that("dcm_vec / dcm_unvec round-trip a matrix", {
  x  <- matrix(seq_len(12), nrow = 3)
  vx <- dcm_vec(x)
  expect_equal(length(vx), 12L)
  expect_equal(dcm_unvec(vx, x), x)
})

test_that("dcm_vec / dcm_unvec round-trip a nested list", {
  x <- list(
    A = matrix(rnorm(6), 2, 3),
    B = array(rnorm(8), dim = c(2, 2, 2)),
    C = matrix(rnorm(4), 2, 2),
    transit = matrix(rnorm(2), 2, 1),
    epsilon = matrix(rnorm(1), 1, 1)
  )
  vx <- dcm_vec(x)
  expect_equal(length(vx), dcm_length(x))
  out <- dcm_unvec(vx, x)
  for (nm in names(x)) {
    expect_equal(out[[nm]], x[[nm]], info = nm)
    expect_equal(dim(out[[nm]]), dim(x[[nm]]), info = nm)
  }
})

test_that("dcm_zeros preserves shape", {
  x  <- list(A = matrix(1:4, 2, 2), b = c(5, 6, 7))
  z  <- dcm_zeros(x)
  expect_equal(dim(z$A), dim(x$A))
  expect_equal(length(z$b), length(x$b))
  expect_true(all(dcm_vec(z) == 0))
})

test_that("dcm_vec on logicals returns numeric vector", {
  x  <- c(TRUE, FALSE, TRUE)
  vx <- dcm_vec(x)
  expect_type(vx, "double")
  expect_equal(vx, c(1, 0, 1))
})
