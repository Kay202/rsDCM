# rsDCM: build & release notes

These steps take the package from this scaffold to a working install,
then to a CRAN submission. Everything in this file is for you (the
maintainer); none of it ships to users.

## One-time setup

```r
install.packages(c("devtools", "roxygen2", "testthat",
                   "knitr", "rmarkdown", "Matrix", "expm", "MASS"))
```

### Pandoc: required, and the usual cause of a failed build

`R CMD build` knits `vignettes/introduction.Rmd`, which needs **pandoc**. If
pandoc is not on the `PATH`, the build does not skip the vignette; it
**aborts entirely** with:

```
Error: processing vignette 'introduction.Rmd' failed with diagnostics:
Pandoc is required to build R Markdown vignettes but not available.
```

RStudio bundles pandoc and puts it on the `PATH` itself, so this works from
the RStudio IDE and fails from a plain terminal. Check with:

```r
rmarkdown::pandoc_available()   # must be TRUE
rmarkdown::pandoc_version()
```

If it is `FALSE`, either build from RStudio, or point R at the copy RStudio
already ships before building:

```r
Sys.setenv(RSTUDIO_PANDOC =
  "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")
```

Do **not** work around a failed build by running `tar czf` on the source
directory by hand. That produces a tarball with an empty `man/` and an
empty `data/`, so whoever installs it gets no help pages and no `toy_dcm`.
See "Step 5: build the tarball to send" below for the correct command.

## Step 1: generate the toy dataset

The package ships `data/toy_dcm.rda` but the file isn't in this scaffold
because it has to be built from R. The build script reads a real DCM
specification from `data-raw/DCM_model.rds` (which **is** in the repo),
strips any posterior-fit fields, normalizes the options block, and
saves the result as `data/toy_dcm.rda`. Run once from the package root:

```bash
Rscript data-raw/make-toy-dcm.R
```

This writes `data/toy_dcm.rda` (target: < 50 KB, xz-compressed).

## Step 2: regenerate documentation and NAMESPACE

The hand-written `NAMESPACE` and the missing `man/*.Rd` files come from
the roxygen comments. Regenerate them:

```r
devtools::document()
```

This rewrites `NAMESPACE` to match every `@export` tag and creates one
`.Rd` per documented function in `man/`.

## Step 3: install locally and run the tests

```r
devtools::install()
devtools::test()
```

The fast tests (`test-vec.R`, `test-smoke.R`) should pass in a couple of
seconds. `test-estimate.R` runs the full Gauss-Newton inversion on
`toy_dcm` (a few seconds) and is skipped on CRAN.

## Step 4: local CRAN check

```r
devtools::check()              # equivalent to R CMD check --as-cran
```

You want **0 errors, 0 warnings, 0 notes**. Common things that fall out
on first check:

* "no visible binding for global variable": fix with `utils::globalVariables()`
  or a `@importFrom` tag.
* "Undocumented arguments": add `@param` for any missing arg in the
  roxygen block.
* Vignette build failure: usually a missing package; install it.
* `data/toy_dcm.rda` not found: run Step 1.

## Step 5: build the tarball to send

To hand the package to someone else (a reviewer, a supervisor, CRAN), build
it properly. **Always use `R CMD build`, never `tar`**:

```bash
cd ..                 # the directory ABOVE the package root
R CMD build rsDCM      # writes rsDCM_0.1.0.tar.gz
```

or from R:

```r
devtools::build()
```

`R CMD build` is what applies `.Rbuildignore` (stripping `.github/`,
`data-raw/`, `.gitignore`, and this file), knits the vignette into
`inst/doc/`, and includes the generated `man/*.Rd` and `data/toy_dcm.rda`.
A hand-rolled `tar czf` does none of that.

Sanity-check the result before sending it:

```bash
tar -tzf rsDCM_0.1.0.tar.gz | grep -c 'man/.*\.Rd'   # expect 45
tar -tzf rsDCM_0.1.0.tar.gz | grep 'data/'           # expect toy_dcm.rda
tar -tzf rsDCM_0.1.0.tar.gz | grep 'inst/doc/'       # expect the built vignette
```

The recipient installs it with:

```r
install.packages("rsDCM_0.1.0.tar.gz", repos = NULL, type = "source")
?dcm_estimate                  # help pages should now resolve
browseVignettes("rsDCM")
```

## Step 6: strengthen the regression test (optional but recommended)

The current `test-estimate.R` only checks that estimation completes and
returns finite values. Once you've installed and run the package once,
capture the reference posterior:

```r
data(toy_dcm)
fit <- suppressMessages(dcm_estimate(toy_dcm))
dput(fit$Ep, file = "tests/testthat/_ref_Ep.R")
```

then add a tolerance comparison to `test-estimate.R`:

```r
ref_Ep <- dget("_ref_Ep.R")
expect_equal(fit$Ep$A, ref_Ep$A, tolerance = 1e-4)
```

This catches any future numerical regression.

## Step 6: CRAN submission

1. Bump `Version:` in `DESCRIPTION` to `0.1.0` (already done).
2. Update `NEWS.md` with the release notes.
3. Final check on R-devel and R-oldrel:
   ```r
   devtools::check_win_devel()
   devtools::check_win_oldrelease()
   ```
4. Submit:
   ```r
   devtools::release()
   ```

## Continuous integration (already set up)

`.github/workflows/R-CMD-check.yaml` runs `R CMD check --as-cran` on every
push and pull request, across a 5-cell matrix:

| OS              | R version  |
| --------------- | ---------- |
| macOS-latest    | release    |
| windows-latest  | release    |
| ubuntu-latest   | devel      |
| ubuntu-latest   | release    |
| ubuntu-latest   | oldrel-1   |

Two pre-check steps regenerate `data/toy_dcm.rda` and the roxygen
documentation, so CI is green on the first push even before you've run
`devtools::document()` locally. Once you're routinely committing those
generated artifacts yourself, you can drop those steps from the workflow
to make CI strictly verify the committed state.

The build-status badge in `README.md` will auto-update once your first
push completes.

## File-by-file map

| Path | Purpose |
| --- | --- |
| `DESCRIPTION` | Package metadata. Update `URL`/`BugReports` once the GitHub repo exists. |
| `LICENSE`, `LICENSE.note` | GPL-2 declaration and SPM25 acknowledgment. |
| `NAMESPACE` | Generated by roxygen2; edit the tags in `R/rsDCM-package.R`, not this file. |
| `R/rsDCM-package.R` | Package-level docs + `.onLoad` + `rsdcm_options()`. |
| `R/utils-vec.R` | Flatten/unflatten primitives. |
| `R/utils-matrix.R` | `dcm_inv`, `dcm_logdet`, `dcm_svd`, `dcm_cat`, ... |
| `R/utils-misc.R` | `dcm_funcheck`, `dcm_data_id`, `dcm_Ncdf`, ... |
| `R/diff.R` | `dcm_diff` numerical Jacobian. |
| `R/expm.R` | `dcm_dx`, `expv`, `padm`, `dcm_expm`. |
| `R/cov-basis.R` | `dcm_Ce`, `dcm_Q` error covariance bases. |
| `R/fmri-model.R` | `dcm_fx_fmri`, `dcm_gx_fmri`, priors, `dcm_phi`. |
| `R/bireduce.R` | `dcm_bireduce`, `dcm_kernels`. |
| `R/integrate.R` | Optimized `dcm_int`. |
| `R/nlsi-gn.R` | `dcm_nlsi_GN` Gauss-Newton inversion. |
| `R/estimate.R` | `dcm_estimate`, model-evidence helpers. |
| `data-raw/make-toy-dcm.R` | Builds `data/toy_dcm.rda`. |
| `tests/testthat/` | Unit + smoke + estimation tests. |
| `vignettes/introduction.Rmd` | End-to-end worked example. |
| `inst/CITATION` | Citation entries. |
| `man/` | Generated by `devtools::document()`; 45 `.Rd` files. Do not edit by hand. |
