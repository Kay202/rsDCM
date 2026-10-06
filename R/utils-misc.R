# Miscellaneous utilities (internal).

#' Coerce a function-like object to a function
#'
#' Accepts a function, a function name (as character), or a one-line body and
#' returns a real \code{function}. Mirrors SPM25's \code{spm_funcheck}.
#'
#' @param f A function or character string.
#' @return A \code{function}.
#' @keywords internal
#' @export
dcm_funcheck <- function(f) {
  if (is.function(f) || is.null(f)) return(f)
  if (is.character(f)) {
    if (exists(f, mode = "function")) return(get(f, mode = "function"))
    tryCatch(eval(parse(text = paste0("function(...) { ", f, " }"))),
             error = function(e) stop(paste("Cannot convert", f, "to a function")))
  } else stop("Input must be a function or character string")
}

#' Compute a stable identifier for a numeric structure
#'
#' Used by SPM25 to fingerprint a dataset; returned in the DCM result.
#'
#' @param ... Numeric data (any shape, possibly nested).
#' @return Numeric scalar.
#' @keywords internal
#' @export
dcm_data_id <- function(...) {
  X <- list(...); if (length(X) == 1) X <- X[[1]]
  ID <- 0
  if (is.character(X)) { ID <- 0
  } else if (is.numeric(X)) {
    Y <- as.numeric(X)
    ID <- sum(abs(Y[!is.nan(Y) & !is.infinite(Y)]))
  } else if (is.list(X)) {
    for (elem in X) ID <- ID + dcm_data_id(elem)
  }
  if (ID > 0) ID <- 10^(-(floor(log10(ID)) - 2)) * ID
  ID
}

#' Univariate Normal CDF (SPM-style)
#'
#' Vectorized Normal CDF with \code{NaN} for non-positive variances.
#'
#' @param x Quantile.
#' @param u Mean.
#' @param v Variance.
#' @return Numeric vector of CDF values.
#' @examples
#' dcm_Ncdf(0)                     # 0.5
#' round(dcm_Ncdf(c(-1.96, 0, 1.96)), 4)
#'
#' # Note the third argument is a VARIANCE, not a standard deviation
#' dcm_Ncdf(2, u = 0, v = 4)       # == pnorm(2, mean = 0, sd = 2)
#' @keywords internal
#' @export
dcm_Ncdf <- function(x, u = 0, v = 1) {
  if (missing(x)) return(numeric(0))

  lx <- length(x); lu <- length(u); lv <- length(v)
  target_len <- max(lx, lu, lv)

  if (lx < target_len) x <- rep_len(x, target_len)
  if (lu < target_len) u <- rep_len(u, target_len)
  if (lv < target_len) v <- rep_len(v, target_len)

  F  <- numeric(target_len)
  md <- v > 0
  if (any(!md)) {
    F[!md] <- NaN
    warning("Returning NaN for variance <= 0.")
  }
  if (any(md)) F[md] <- stats::pnorm(x[md], mean = u[md], sd = sqrt(v[md]))
  F
}

#' Locate or describe entries in a structured parameter
#'
#' Given a structured parameter object \code{X} (with named fields), returns
#' the linear indices of one or more named fields, or the structural location
#' of a given linear index.
#'
#' @param X A structured parameter (list with numeric fields).
#' @param ... One or more field names or numeric indices.
#' @return Integer indices or a character description.
#' @keywords internal
#' @export
dcm_fieldindices <- function(X, ...) {
  varargin <- list(...)
  if (length(varargin) > 0 && is.numeric(varargin[[1]]) &&
      length(varargin[[1]]) > 1) {
    return(lapply(varargin[[1]], function(j) dcm_fieldindices(X, j)))
  }
  X0 <- dcm_zeros(X); ix <- dcm_vec(X0)
  for (i in seq_along(varargin)) {
    if (is.character(varargin[[i]])) {
      field <- varargin[[i]]; x <- X0
      if (field %in% names(x)) {
        f <- x[[field]]; f <- dcm_unvec(dcm_vec(f) + 1, f); x[[field]] <- f
        ix <- ix + dcm_vec(x)
      } else warning(paste("Field", field, "not found"))
    } else if (is.numeric(varargin[[i]])) {
      ind <- varargin[[i]]; x <- ix; x[ind] <- 1; x <- dcm_unvec(x, X)
      for (nm in names(x)) {
        if (any(dcm_vec(x[[nm]]) != 0)) {
          v <- x[[nm]]
          if (is.numeric(v) || is.matrix(v)) {
            idx <- which(v != 0, arr.ind = TRUE)
            if (is.null(dim(v)) || min(dim(v)) == 1) {
              return(sprintf("%s[%i]", nm, which(v != 0)[1]))
            } else if (length(dim(v)) < 3) {
              return(sprintf("%s[%i,%i]", nm, idx[1, 1], idx[1, 2]))
            } else {
              return(sprintf("%s[%i,%i,%i]", nm, idx[1, 1], idx[1, 2], idx[1, 3]))
            }
          }
        }
      }
    }
  }
  which(ix != 0)
}

#' Find indices of free parameters in a DCM
#'
#' @param ... Either \code{(DCM, fields)} or \code{(pC, pE, fields)}.
#' @return List with \code{i} (free indices), plus \code{pC}, \code{pE}, \code{Np}.
#' @keywords internal
#' @export
dcm_find_pC <- function(...) {
  varargin <- list(...)
  pC <- NULL; pE <- NULL; fields <- NULL
  if (length(varargin) > 2) {
    pC <- varargin[[1]]; pE <- varargin[[2]]; fields <- varargin[[3]]
  } else if (length(varargin) > 1) {
    DCM <- varargin[[1]]; fields <- varargin[[2]]
  } else { DCM <- varargin[[1]]; fields <- NULL }

  if (is.null(pC) || is.null(pE)) {
    if (is.character(DCM)) DCM <- readRDS(DCM)
    if (any(c("options", "M") %in% names(DCM))) {
      tryCatch({
        result <- dcm_find_rC(DCM); pC <- result$pC; pE <- result$pE
      }, error = function(e) {})
    }
  }

  q <- if (is.list(pC) && !is.matrix(pC)) dcm_vec(pC) else diag(pC)
  threshold <- mean(q[q < 1024]) / 1024
  i <- which(q > threshold); Np <- length(q)

  if (!is.null(fields)) {
    if (is.character(fields) && length(fields) == 1) fields <- list(fields)
    if (is.list(pE)) {
      j <- do.call(dcm_fieldindices, c(list(pE), fields))
      if (!length(j) && !(length(fields) > 0 && fields[[1]] == "none")) {
        warning(sprintf("%s not found. Returning all fields",
                        paste(unlist(fields), collapse = ",")))
      } else i <- j[j %in% i]
    }
  }
  list(i = i, pC = pC, pE = pE, Np = Np)
}

#' Recover priors from a DCM container
#'
#' @param DCM A DCM list.
#' @return List with \code{pC}, \code{pE}.
#' @keywords internal
#' @export
dcm_find_rC <- function(DCM) {
  if (!is.null(DCM$M$pC) && !is.null(DCM$M$pE)) {
    return(list(pC = DCM$M$pC, pE = DCM$M$pE))
  }
  if (is.null(DCM$d)) DCM$d <- array(0, dim = c(nrow(DCM$a), nrow(DCM$a), 0))
  result <- dcm_fmri_priors(DCM$a, DCM$b, DCM$c, DCM$d, DCM$options)
  list(pC = result$pC, pE = result$pE)
}
