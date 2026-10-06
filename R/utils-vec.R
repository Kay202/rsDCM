# Vector / list flattening utilities (internal). Performance-critical.

#' Flatten a nested numeric structure into a vector
#'
#' Recursively walks a numeric, logical, or nested list structure and returns
#' a single flat numeric vector. Internal helper; mirrors SPM25's \code{spm_vec}.
#'
#' @param X A numeric, logical, or list (possibly nested).
#' @param ... Additional structures to flatten and concatenate.
#' @return A numeric vector.
#' @examples
#' # A DCM parameter structure flattens in field order
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' dcm_vec(P)
#'
#' # Nested lists are walked recursively
#' dcm_vec(list(a = 1, b = list(c = 2:3, d = 4)))
#' @seealso \code{\link{dcm_unvec}} for the inverse operation.
#' @keywords internal
#' @export
dcm_vec <- function(X, ...) {
  args <- list(...)
  if (length(args) > 0) X <- c(list(X), args)

  if (is.numeric(X))  return(as.vector(as.numeric(X)))
  if (is.logical(X))  return(as.numeric(X))

  .fill <- function(X, buf, off) {
    if (is.numeric(X) || is.logical(X)) {
      v <- as.numeric(X)
      n <- length(v)
      if (n > 0) buf[off + seq_len(n)] <- v
      return(list(buf = buf, off = off + n))
    }
    if (is.list(X)) {
      items <- if (!is.null(names(X)) && all(names(X) != "")) X[names(X)] else X
      for (elem in items) {
        r   <- .fill(elem, buf, off)
        buf <- r$buf
        off <- r$off
      }
      return(list(buf = buf, off = off))
    }
    list(buf = buf, off = off)
  }

  total <- dcm_length(X)
  if (total == 0L) return(numeric(0))
  buf <- numeric(total)
  r   <- .fill(X, buf, 0L)
  r$buf
}

#' Total number of numeric entries in a nested structure
#'
#' Returns the length \code{dcm_vec(X)} would produce, without allocating it.
#'
#' @param X A numeric, logical, or list (possibly nested).
#' @return Integer length.
#' @examples
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' dcm_length(P)            # 6
#' length(dcm_vec(P))       # same, but allocates the vector
#' @keywords internal
#' @export
dcm_length <- function(X) {
  if (is.numeric(X) || is.logical(X)) return(length(X))
  if (!is.list(X)) return(0L)

  stack <- X
  total <- 0L
  while (length(stack) > 0) {
    item  <- stack[[length(stack)]]
    stack <- stack[-length(stack)]
    if (is.numeric(item) || is.logical(item)) {
      total <- total + length(item)
    } else if (is.list(item)) {
      stack <- c(stack, item)
    }
  }
  total
}

#' Reshape a flat vector back into a template structure
#'
#' Inverse of \code{\link{dcm_vec}}: splits a flat numeric vector across the
#' shape implied by one or more template arguments.
#'
#' @param vX Flat numeric vector.
#' @param ... One or more templates whose structure is copied; \code{vX} is
#'   distributed across them in order.
#' @return A structure (or list of structures) matching the templates.
#' @examples
#' # Round-trip: flatten a structure and rebuild it
#' P <- list(A = matrix(1:4, 2), C = c(5, 6))
#' identical(dcm_unvec(dcm_vec(P), P), P)
#'
#' # The template supplies the shape; the vector supplies the values
#' tmpl <- list(A = matrix(0, 2, 2), C = numeric(2))
#' dcm_unvec(1:6, tmpl)
#' @seealso \code{\link{dcm_vec}} for the forward operation.
#' @keywords internal
#' @export
dcm_unvec <- function(vX, ...) {
  vX <- as.numeric(vX)
  args <- list(...)

  if (length(args) > 1) {
    out <- vector("list", length(args))
    off <- 0L
    for (k in seq_along(args)) {
      nk <- dcm_length(args[[k]])
      out[[k]] <- dcm_unvec(if (nk > 0) vX[off + seq_len(nk)] else numeric(0),
                            args[[k]])
      off <- off + nk
    }
    return(out)
  }

  X <- if (length(args) == 1L) args[[1L]] else args
  if (!(is.numeric(vX) || is.logical(vX))) vX <- dcm_vec(vX)

  if (is.numeric(X) || is.logical(X)) {
    X[] <- vX
    return(X)
  }

  env <- new.env(parent = emptyenv())
  env$off <- 0L

  .fill_unvec <- function(X) {
    if (is.numeric(X) || is.logical(X)) {
      n <- length(X)
      if (n > 0) {
        X[] <- vX[env$off + seq_len(n)]
        env$off <- env$off + n
      }
      return(X)
    }
    if (is.list(X)) {
      for (i in seq_along(X)) X[[i]] <- .fill_unvec(X[[i]])
      return(X)
    }
    X
  }

  .fill_unvec(X)
}

#' Zero-fill a structure
#'
#' Returns a copy of \code{X} with the same shape but every numeric entry
#' replaced by 0.
#'
#' @param X A nested structure.
#' @return A structure with zeros, same shape as \code{X}.
#' @examples
#' dcm_zeros(list(A = matrix(1:4, 2), C = c(5, 6)))
#' @keywords internal
#' @export
dcm_zeros <- function(X) dcm_unvec(rep(0, dcm_length(X)), X)
