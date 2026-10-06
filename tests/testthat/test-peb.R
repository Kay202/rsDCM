# =============================================================================
# PEB tests. These build subject-level densities directly (no inversion), so
# they stay fast enough to run on CRAN.
#
# Subject posteriors are constructed by Bayes from a known likelihood, i.e.
#   qC = (Lam + inv(pC))^-1,  qE = qC (Lam m + inv(pC) pE)
# rather than by planting the true value in qE. This matters: PEB inverts the
# first-level prior out of each posterior via Bayesian model reduction, so a
# naively planted qE implies a likelihood far from the value planted.
# =============================================================================

# Fresh temp directory, auto-removed when the calling test frame exits.
new_tempdir <- function(env = parent.frame()) {
  d <- tempfile("peb")
  dir.create(d)
  withr::defer(unlink(d, recursive = TRUE), envir = env)
  d
}

# Two-region priors, and helpers to synthesise subjects with a known effect.
peb_fixture <- function() {
  n   <- 2L
  pri <- dcm_fmri_priors(matrix(1, n, n), array(0, c(n, n, 1)),
                         matrix(c(1, 0), n, 1), array(0, c(n, n, 0)),
                         options = list())
  q       <- dcm_find_pC(list(M = list(pE = pri$pE, pC = pri$pC)), "A")$i
  pC_full <- as.matrix(pri$pC)
  list(pri = pri, q = q, pC_full = pC_full,
       iPC_A = solve(pC_full[q, q]), pE_A = dcm_vec(pri$pE)[q])
}

make_subject <- function(fx, theta, noise = 0.1, lam = 100) {
  Lam <- diag(lam, length(theta))
  qC  <- solve(Lam + fx$iPC_A)
  m   <- theta + stats::rnorm(length(theta), 0, noise)
  qE  <- as.vector(qC %*% (Lam %*% m + fx$iPC_A %*% fx$pE_A))
  Ep  <- fx$pri$pE; Ep$A <- matrix(qE, 2, 2)
  Cp  <- fx$pC_full * 0.5; Cp[fx$q, fx$q] <- qC
  list(M = list(pE = fx$pri$pE, pC = fx$pri$pC), Ep = Ep, Cp = Cp,
       F = -100 + stats::rnorm(1))
}

test_that("dcm_peb_design builds intercept and covariate designs", {
  X <- dcm_peb_design(4)
  expect_equal(dim(X), c(4L, 1L))
  expect_equal(colnames(X), "Intercept")
  expect_true(all(X == 1))

  Xc <- dcm_peb_design(4, covariates = data.frame(age = c(21, 34, 46, 58)))
  expect_equal(nrow(Xc), 4L)
  expect_equal(ncol(Xc), 2L)          # intercept + age
})

test_that("dcm_peb_files parses subject numbers and errors on an empty dir", {
  d <- new_tempdir()
  saveRDS(list(), file.path(d, "PEB_sub-01.rds"))
  saveRDS(list(), file.path(d, "PEB_sub-02.rds"))

  f <- dcm_peb_files(d)
  expect_equal(sort(f$subjects), c(1L, 2L))
  expect_length(f$files, 2L)

  empty <- new_tempdir()
  expect_error(dcm_peb_files(empty), "No PEB files found")
})

test_that("dcm_peb recovers a planted group effect", {
  set.seed(7)
  fx     <- peb_fixture()
  TRUE_A <- c(-0.20, 0.50, -0.30, 0.10)
  Ns     <- 16L
  P      <- lapply(seq_len(Ns), function(i) make_subject(fx, TRUE_A))

  prep <- dcm_peb_prepare(P, M = list(X = dcm_peb_design(Ns), Q = "single",
                                      maxit = 64), field = "A")
  expect_equal(prep$Ns, Ns)
  expect_equal(prep$Np, 4L)
  expect_equal(prep$Nx, 1L)

  res <- suppressMessages(dcm_peb_run(prep, verbose = FALSE))

  expect_true(all(is.finite(res$PEB$Ep)))
  expect_true(is.finite(res$PEB$F))
  expect_length(res$P, Ns)                      # updated first-level DCMs
  # Group estimate tracks the planted effect (PEB shrinks, so allow slack).
  expect_gt(stats::cor(as.numeric(res$PEB$Ep), TRUE_A), 0.95)
})

test_that("dcm_peb separates a between-subject group difference", {
  set.seed(11)
  fx   <- peb_fixture()
  Ns   <- 20L
  grp  <- rep(c(0, 1), each = Ns / 2)
  BASE <- c(-0.2, 0.5, -0.3, 0.1)
  DIFF <- c(0, 0.4, 0, 0)                       # effect on parameter 2 only
  P    <- lapply(seq_len(Ns), function(i) make_subject(fx, BASE + grp[i] * DIFF))

  X    <- dcm_peb_design(Ns, covariates = data.frame(group = grp))
  prep <- dcm_peb_prepare(P, M = list(X = X, Q = "single", maxit = 64),
                          field = "A")
  res  <- suppressMessages(dcm_peb_run(prep, verbose = FALSE))

  expect_equal(dim(res$PEB$Ep), c(4L, 2L))      # params x covariates
  # The group slope loads on parameter 2, where the difference was planted.
  expect_equal(which.max(abs(res$PEB$Ep[, 2])), 2L)
  # Design column names propagate into Xnames.
  expect_equal(res$PEB$Xnames, c("(Intercept)", "group"))
})

test_that("all M$Q options run, including 'none' (no hyperparameters)", {
  set.seed(5)
  fx <- peb_fixture()
  P  <- lapply(1:8, function(i) make_subject(fx, c(-0.2, 0.5, -0.3, 0.1)))
  X  <- dcm_peb_design(8)

  for (opt in c("single", "fields", "all", "none")) {
    prep <- dcm_peb_prepare(P, M = list(X = X, Q = opt, maxit = 16),
                            field = "A")
    res  <- suppressMessages(dcm_peb_run(prep, verbose = FALSE))
    expect_true(all(is.finite(res$PEB$Ep)), info = paste("Q =", opt))
    expect_true(is.finite(res$PEB$F),       info = paste("Q =", opt))
  }

  # 'none' means no covariance components, hence no hyperparameters at all.
  prep_none <- dcm_peb_prepare(P, M = list(X = X, Q = "none", maxit = 16),
                               field = "A")
  expect_equal(prep_none$Ng, 0L)
  expect_length(prep_none$g, 0L)
})

test_that("dcm_peb_of_pebs runs a third level and honours save_path", {
  set.seed(3)
  fx     <- peb_fixture()
  TRUE_A <- c(-0.2, 0.5, -0.3, 0.1)
  d      <- new_tempdir()

  # Level 2: one PEB per subject, written as PEB_sub-XX.rds
  for (s in 1:5) {
    Pi   <- lapply(1:6, function(i) make_subject(fx, TRUE_A, noise = 0.08))
    prep <- dcm_peb_prepare(Pi, M = list(X = dcm_peb_design(6), Q = "single",
                                         maxit = 32), field = "A")
    peb  <- suppressMessages(dcm_peb_run(prep, verbose = FALSE))$PEB
    saveRDS(peb, file.path(d, sprintf("PEB_sub-%02d.rds", s)))
  }

  # field = "all" on PEB inputs: M$pE is unnamed, so it resolves positionally
  res <- suppressMessages(dcm_peb_of_pebs(d, field = "all", verbose = FALSE))

  expect_true(is.finite(res$group_PEB$F))
  expect_true(all(is.finite(res$group_PEB$Ep)))
  expect_equal(sort(res$subjects), 1:5)
  expect_s3_class(res$summary, "data.frame")
  # Nothing written unless asked (same policy as dcm_estimate).
  expect_null(res$save_path)
  expect_length(list.files(d, pattern = "GROUP"), 0L)
  # Subject labels come from the filenames, not the inner PEB's first subject.
  expect_true(all(grepl("^sub", res$group_PEB$Snames)))

  # With save_path, the group PEB is written.
  sp <- file.path(new_tempdir(), "PEB_GROUP.rds")
  res2 <- suppressMessages(dcm_peb_of_pebs(d, field = "all", save_path = sp,
                                           verbose = FALSE))
  expect_true(file.exists(sp))
  expect_equal(res2$save_path, sp)
})
