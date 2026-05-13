# =============================================================================
# data-raw/make-toy-dcm.R
#
# Build the toy_dcm dataset shipped with the package by loading
# data-raw/DCM_model.rds (a real 2-region DCM contributed by the package
# authors), stripping any posterior-fit fields so the example is a clean
# input specification, and writing data/toy_dcm.rda.
#
# Usage (from package root, once):
#   Rscript data-raw/make-toy-dcm.R
#
# Requirements: an R session with no extra packages — uses only base R.
# Run after installing the package the first time, then commit the
# resulting data/toy_dcm.rda alongside it.
# =============================================================================

src_path <- "data-raw/DCM_model.rds"
if (!file.exists(src_path)) {
  stop(sprintf("Source file not found at %s. ", src_path),
       "Check you are running this script from the package root.")
}

DCM <- readRDS(src_path)
if (!is.list(DCM)) {
  stop("DCM_model.rds did not contain a list. Got: ",
       paste(class(DCM), collapse = "/"))
}

# -- Validate the input fields a DCM specification needs ---------------------

required <- c("a", "b", "c", "Y", "U")
missing  <- setdiff(required, names(DCM))
if (length(missing) > 0L) {
  stop("DCM is missing required fields: ", paste(missing, collapse = ", "))
}
if (is.null(DCM$Y$y) || !is.matrix(DCM$Y$y)) {
  stop("DCM$Y$y must be a numeric matrix.")
}
if (is.null(DCM$U$u)) {
  stop("DCM$U$u must be present (input design).")
}

# -- Strip posterior-fit fields so the example is a clean specification ------
# If the source DCM was the output of dcm_estimate(), it carries Ep/Cp/Vp/Pp
# and a populated M slot. We don't want those in the bundled toy dataset:
# users should be able to start from `data(toy_dcm); fit <- dcm_estimate(toy_dcm)`
# and see the inversion run from scratch.

posterior_fields <- c("Ep", "Cp", "Vp", "Pp", "Eh", "F", "L",
                      "H1", "K1", "R", "y", "T", "ID", "Ce",
                      "M",                                # rebuilt by estimate
                      "dFdp", "dFdpp")
toy_dcm <- DCM[setdiff(names(DCM), posterior_fields)]

# -- Make sure the options block exists and has sensible defaults ------------
if (is.null(toy_dcm$options)) toy_dcm$options <- list()
defaults <- list(
  nonlinear  = 0,
  two_state  = 0,
  stochastic = 0,
  centre     = 1,
  induced    = 0,
  maxit      = 32,
  maxnodes   = 8,
  hE         = 6,
  hC         = 1 / 128
)
for (nm in names(defaults)) {
  if (is.null(toy_dcm$options[[nm]])) toy_dcm$options[[nm]] <- defaults[[nm]]
}

# -- Make sure n and v are present -------------------------------------------
if (is.null(toy_dcm$n)) toy_dcm$n <- nrow(toy_dcm$a)
if (is.null(toy_dcm$v)) toy_dcm$v <- nrow(toy_dcm$Y$y)
if (is.null(toy_dcm$TE)) toy_dcm$TE <- 0.04

# -- Report ------------------------------------------------------------------

message(sprintf("Loaded DCM with %d regions, %d samples.",
                toy_dcm$n, toy_dcm$v))
message(sprintf("Fields kept: %s",
                paste(names(toy_dcm), collapse = ", ")))
stripped <- intersect(names(DCM), posterior_fields)
if (length(stripped) > 0L) {
  message(sprintf("Stripped posterior fields: %s",
                  paste(stripped, collapse = ", ")))
}

# -- Save --------------------------------------------------------------------

out_path <- file.path("data", "toy_dcm.rda")
if (!dir.exists("data")) dir.create("data")
save(toy_dcm, file = out_path, compress = "xz", compression_level = 9)

fsize <- file.info(out_path)$size
message(sprintf("Wrote %s (%.1f KB)", out_path, fsize / 1024))
if (fsize > 50 * 1024) {
  warning(sprintf("toy_dcm.rda is %.1f KB, above the 50 KB target. ",
                  fsize / 1024),
          "Consider whether all kept fields are necessary.")
}
