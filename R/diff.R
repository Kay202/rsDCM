# =============================================================================
# Numerical differentiation (internal). dcm_diff is a hot-path function:
# called once per parameter direction in every Gauss-Newton iteration.
# =============================================================================

#' High-order numerical Jacobian
#'
#' Forward-difference numerical Jacobian of a (possibly nested-output)
#' function with respect to one or more of its arguments. Mirrors SPM12's
#' \code{dcm_diff}.
#'
#' @param ... Function, its arguments, the index (or vector of indices) of
#'   the argument(s) to differentiate with respect to, and optionally a list
#'   of projection matrices.
#' @return A list containing the Jacobian \code{J} (or higher-order
#'   derivatives) and the function value \code{f0}.
#' @keywords internal
#' @export
dcm_diff <- function(...) {
  args <- list(...)
  if (length(args) < 2) stop("Improper call: need at least f and one x argument")

  dx <- tryCatch({
    if (exists("GLOBAL_DX", envir = .dcmR_env, inherits = FALSE) &&
        !is.null(get("GLOBAL_DX", envir = .dcmR_env)))
      get("GLOBAL_DX", envir = .dcmR_env)
    else exp(-8)
  }, error = function(e) exp(-8))

  f <- args[[1]]
  if (is.character(f)) {
    if (!exists(f, mode = "function")) stop("Function name not found: ", f)
    f <- get(f, mode = "function")
  } else if (!is.function(f)) {
    stop("First argument must be a function or name of a function")
  }
  if (!is.function(f)) f <- dcm_funcheck(f)

  q_flag <- TRUE; V_list <- NULL; n <- NULL
  last <- args[[length(args)]]

  if (is.list(last)) {
    V_list <- last; n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]; q_flag <- TRUE
  } else if (is.character(last)) {
    n <- as.integer(args[[length(args) - 1]])
    xs <- args[2:(length(args) - 2)]
    V_list <- vector("list", length(xs)); q_flag <- FALSE
  } else if (is.numeric(last)) {
    n <- as.integer(last); xs <- args[2:(length(args) - 1)]
    V_list <- vector("list", length(xs)); q_flag <- TRUE
  } else stop("Improper call.")

  if (any(n < 1L) || any(n > length(xs))) stop("n contains invalid argument indices")
  if (is.null(V_list)) V_list <- vector("list", length(xs))
  if (length(V_list) < length(xs))
    V_list <- c(V_list, vector("list", length(xs) - length(V_list)))

  for (i in seq_along(xs)) {
    if (is.null(V_list[[i]]) && any(n == i))
      V_list[[i]] <- diag(1, dcm_length(xs[[i]]))
  }

  m  <- n[length(n)]
  xm <- dcm_vec(xs[[m]])
  Vm_obj <- V_list[[m]]
  ncols_m <- if (!is.null(Vm_obj)) ncol(as.matrix(Vm_obj)) else 0L
  Jcols <- vector("list", ncols_m)

  V_is_identity <- FALSE
  if (ncols_m > 0L) {
    if (is.matrix(Vm_obj) && nrow(Vm_obj) == ncols_m && ncols_m == length(xm)) {
      V_is_identity <- isTRUE(all(diag(Vm_obj) == 1)) &&
                       isTRUE(sum(abs(Vm_obj)) == ncols_m)
    }
  }

  if (length(n) == 1L) {
    f0 <- do.call(f, xs)
    if (ncols_m > 0L) {
      if (V_is_identity) {
        for (i in seq_len(ncols_m)) {
          xm_p <- xm; xm_p[i] <- xm_p[i] + dx
          xi <- xs; xi[[m]] <- dcm_unvec(xm_p, xs[[m]])
          fi <- do.call(f, xi)
          Jcols[[i]] <- dcm_dfdx(fi, f0, dx)
        }
      } else {
        Vm <- as.matrix(Vm_obj)
        for (i in seq_len(ncols_m)) {
          xi <- xs
          xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
          fi <- do.call(f, xi)
          Jcols[[i]] <- dcm_dfdx(fi, f0, dx)
        }
      }
    }
    fv <- dcm_vec(f0)
    if (length(xm) == 0L) {
      J <- matrix(numeric(0), nrow = length(fv), ncol = 0)
    } else if (length(fv) == 0L) {
      J <- matrix(numeric(0), nrow = 0, ncol = length(xm))
    } else if (is.numeric(f0) && ncols_m > 0L && q_flag) {
      J <- dcm_dfdx_cat(Jcols)
    } else {
      J <- Jcols
    }
    return(list(J = J, f0 = f0))
  }

  if (length(n) == 2L) {
    d1 <- do.call(dcm_diff, c(list(f), xs, list(n[1], V_list)))
    if (ncols_m > 0L) {
      Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
      for (i in seq_len(ncols_m)) {
        xi <- xs
        xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
        fi <- do.call(dcm_diff, c(list(f), xi, list(n[1], V_list)))
        Jcols[[i]] <- dcm_dfdx(fi$J, d1$J, dx)
        p_numeric <- p_numeric && is.numeric(Jcols[[i]])
      }
      d2 <- if (p_numeric && q_flag) dcm_dfdx_cat(Jcols) else Jcols
    } else {
      d2 <- if (q_flag) matrix(0, nrow = length(dcm_vec(d1$f0)), ncol = 0)
            else list()
    }
    return(list(d2 = d2, d1 = d1$J, f0 = d1$f0))
  }

  lower <- do.call(dcm_diff, c(list(f), xs,
                               list(n[1:(length(n) - 1)], V_list)))
  if (ncols_m > 0L) {
    Vm <- as.matrix(V_list[[m]]); p_numeric <- TRUE
    for (i in seq_len(ncols_m)) {
      xi <- xs
      xi[[m]] <- dcm_unvec(xm + as.numeric(Vm[, i]) * dx, xs[[m]])
      fi <- do.call(dcm_diff, c(list(f), xi,
                                list(n[1:(length(n) - 1)], V_list)))
      Jcols[[i]] <- dcm_dfdx(fi$J, lower$J, dx)
      p_numeric <- p_numeric && is.numeric(Jcols[[i]])
    }
    J <- if (p_numeric && q_flag) dcm_dfdx_cat(Jcols) else Jcols
  } else {
    J <- if (q_flag) matrix(0, nrow = length(dcm_vec(lower$f0)), ncol = 0)
         else list()
  }
  list(J = J, f0 = lower$f0)
}

# Internal: forward-difference quotient that handles list-valued outputs.
dcm_dfdx <- function(f, f0, dx) {
  if (is.list(f)) {
    dfdx <- f
    for (i in seq_along(f)) dfdx[[i]] <- dcm_dfdx(f[[i]], f0[[i]], dx)
    return(dfdx)
  } else if (is.list(f0) && !is.data.frame(f0)) {
    return((dcm_vec(f) - dcm_vec(f0)) / dx)
  } else {
    return((f - f0) / dx)
  }
}

# Internal: concatenate per-column Jacobian pieces into a single matrix.
dcm_dfdx_cat <- function(J) {
  if (is.vector(J[[1]])) {
    if (is.null(dim(J[[1]])) || ncol(as.matrix(J[[1]])) == 1)
      return(dcm_cat(J))
    else
      return(t(dcm_cat(lapply(J, function(x) t(as.matrix(x))))))
  }
  J
}
