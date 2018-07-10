#' Inverse AI wavelet transform for curve of HPD matrices
#'
#' \code{InvWavTransf1D} computes an inverse intrinsic average-interpolation (AI) wavelet
#' transform mapping an array of coarsest-scale HPD midpoints combined with a pyramid of Hermitian
#' wavelet coefficients to a curve in the manifold of HPD matrices equipped with a metric specified by the user,
#' as described in \insertCite{CvS17}{pdSpecEst} and Chapter 3 of \insertCite{C18}{pdSpecEst}. This is
#' the inverse operation of the function \code{\link{WavTransf1D}}.
#'
#' The input list of arrays \code{D} and array \code{M0} correspond to a pyramid of wavelet coefficients and
#' the coarsest-scale HPD midpoints respectively, both are structured in the same way as in the output of
#' \code{WavTransf1D}. As in the forward AI wavelet transform, if the refinement order is an odd integer smaller or
#' equal to 9, the function computes the inverse wavelet transform using a fast wavelet refinement scheme based on
#' weighted intrinsic averages with pre-determined weights as explained in \insertCite{CvS17}{pdSpecEst} and Chapter 3 of
#' \insertCite{C18}{pdSpecEst}. If the refinement order is an odd integer larger than 9, the wavelet refinement
#' scheme uses intrinsic polynomial prediction based on Neville's algorithm in the Riemannian manifold (via \code{\link{pdNeville}}).\cr
#' The function computes the inverse intrinsic AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics.
#'
#' @param D a list of arrays containing the pyramid of wavelet coefficients, where each array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to the finest
#' wavelet scale \code{j = jmax}. This is the same format as the \code{$D} component given as output by
#'  \code{\link{WavTransf1D}}.
#' @param M0 a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the midpoint pyramid.
#' This is the same format as the \code{$M0} component given as output by \code{\link{WavTransf1D}}.
#' @param order an odd integer larger or equal to 1 corresponding to the order of the intrinsic AI refinement scheme,
#' defaults to \code{order = 5}. Note that if \code{order > 9}, the computational cost
#' significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param jmax the maximum scale (resolution) up to which the HPD midpoints (i.e. scaling coefficients) are reconstructed.
#' If \code{jmax} is not specified it is set equal to the resolution in the finest wavelet scale \code{jmax = length(D)}.
#' @param periodic a logical value determining whether the curve of HPD matrices can be reflected at the boundary for
#' improved wavelet refinement schemes near the boundaries of the domain. This is useful for spectral matrix estimation,
#' where the spectral matrix is a symmetric and periodic curve in the frequency domain. Defaults to \code{periodic = F}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. See also the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples1D(2^8, example = "bumps")
#' P.wt <- WavTransf1D(P$f) ## forward transform
#' P.f <- InvWavTransf1D(P.wt$D, P.wt$M0) ## backward transform
#' all.equal(P.f, P$f)
#'
#' @return Returns a (\eqn{d, d, m})-dimensional array corresponding to a length \eqn{m} curve of
#' (\eqn{d,d})-dimensional HPD matrices.
#'
#' @seealso \code{\link{WavTransf1D}}, \code{\link{pdSpecEst1D}}, \code{\link{pdNeville}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
InvWavTransf1D <- function(D, M0, order = 5, jmax, periodic = F, metric = "Riemannian", ...) {

  ## Initialize variables
  dots <- list(...)
  return_val <- (if(is.null(dots$return_val)) "f" else dots$return_val)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  chol_bias <- (if(is.null(dots$chol_bias)) F else dots$chol_bias)
  if (!(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order <- 5
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  L <- (order - 1) / 2
  L_round <- 2 * ceiling(L / 2)
  d <- nrow(D[[1]][, , 1])
  J <- (if(missing(jmax)) length(D) else jmax)
  if(!isTRUE((J > 0) & isTRUE(all.equal(as.integer(J), J)))) {
    stop("'jmax' should be an integer larger than zero")
  }

  ## Reconstruct midpoints
  m1 <- M0

  for (j in 0:(J - 1)) {
    n_M <- dim(m1)[3]
    L1 <- ifelse(order > n_M, floor((n_M - 1) / 2), L)
    ## Impute C++
    tm1 <- impute_C(m1, W_1D[[min(L1 + 1, 5)]], L1, T, metric, method)
    if(periodic){
      tm1 <- tm1[, , ifelse(j > 0, L_round, 2 * floor(L / 2)) + 1:(2^(j + 1) + 2 * L_round), drop = F]
    }
    L1 <- ifelse(periodic, L_round / 2 + ifelse((j > 0) | (L %% 2 == 0), 0, -1), 0)
    ## Reconstruct C++
    Dj <- (if(j < length(D)) D[[j + 1]] else array(dim = c(d, d, as.integer(dim(tm1)[3]/2))))
    m1 <- reconstr_C(tm1, m1, Dj, j, ifelse(j < length(D), dim(D[[j + 1]])[3], as.integer(dim(tm1)[3]/2)),
                     isTRUE(j < length(D)), L1, metric)
  }

  ## Transform back to manifold
  if(return_val == "f") {
    if(metric == "logEuclidean" | metric == "Cholesky" | metric == "rootEuclidean") {
      m1 <- Ptransf2D_C(m1, T, chol_bias, metric)
    }
  }
  return((if(periodic) m1[, , L_round + 1:2^J] else m1))
}

#' Inverse AI wavelet transform for surface of HPD matrices
#'
#' \code{InvWavTransf2D} computes the inverse intrinsic average-interpolation (AI) wavelet
#' transform mapping an array of coarsest-scale HPD midpoints combined with a 2D pyramid of Hermitian
#' wavelet coefficients to a surface in the manifold of HPD matrices equipped with a metric specified by the
#' user, as described in Chapter 5 of \insertCite{C18}{pdSpecEst}. This is the inverse operation of the
#' function \code{\link{WavTransf2D}}.
#'
#' The input list of arrays \code{D} and array \code{M0} correspond to a 2D pyramid of wavelet coefficients and
#' the coarsest-scale HPD midpoints respectively, both are structured in the same way as in the output of
#' \code{WavTransf2D}. As in the forward AI wavelet transform, the marginal refinement orders should be smaller
#' or equal to 9, and the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
#' intrinsic averages with pre-determined weights as explained in Chapter 5 of \insertCite{C18}{pdSpecEst}. By default
#' \code{WavTransf2D} computes the inverse intrinsic 2D AI wavelet transform equipping the space of HPD matrices with (i)
#' the affine-invariant Riemannian metric as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}.
#' Instead, the space of HPD matrices can also be equipped with one of the following metrics; (ii) the Log-Euclidean metric, the
#' Euclidean inner product between matrix logarithms; (iii) the Cholesky metric, the Euclidean inner product between Cholesky
#' decompositions; (iv) the Euclidean metric and (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian)
#' satisfies several useful properties not shared by the other metrics, see \insertCite{C18}{pdSpecEst} for more details. Note that this
#' comes at the cost of increased computation time in comparison to one of the other metrics.
#'
#' @param D a list of arrays containing the 2D pyramid of wavelet coefficients, where each array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to the finest
#' wavelet scale \code{j = jmax}. This is the same format as the \code{$D} component given as output by
#'  \code{\link{WavTransf2D}}.
#' @param M0 a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.
#' This is the same format as the \code{$M0} component given as output by \code{\link{WavTransf2D}}.
#' @param order a 2-dimensional numeric vector \eqn{(1,1) \le} \code{order} \eqn{\le (9,9)} corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}.
#' @param jmax the maximum scale (resolution) up to which the 2D surface of HPD midpoints (i.e. scaling coefficients) are
#' reconstructed. If \code{jmax} is not specified it is set equal to the resolution in the finest wavelet scale
#' \code{jmax = length(D)}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. See also the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples2D(c(2^4, 2^4), 2, example = "tvar")
#' P.wt <- WavTransf2D(P$f) ## forward transform
#' P.f <- InvWavTransf2D(P.wt$D, P.wt$M0) ## backward transform
#' all.equal(P.f, P$f)
#'
#' @return Returns a (\eqn{d, d, n_1, n_2})-dimensional array corresponding to a rectangular surface of size \eqn{n_1} by
#' \eqn{n_2} of (\eqn{d,d})-dimensional HPD matrices.
#'
#' @seealso \code{\link{WavTransf2D}}, \code{\link{pdSpecEst2D}}, \code{\link{pdNeville}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
InvWavTransf2D <- function(D, M0, order = c(3, 3), jmax, metric = "Riemannian", ...) {

  ## Set variables
  dots <- list(...)
  return_val <- (if(is.null(dots$return_val)) "f" else dots$return_val)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  chol_bias <- (if(is.null(dots$chol_bias)) F else dots$chol_bias)
  if (!isTRUE((order[1] %% 2 == 1) & (order[2] %% 2 == 1))) {
    warning("Refinement orders in both directions should be odd integers, by default set to c(5,5).")
    order <- c(3, 3)
  }
  if (!isTRUE((order[1] <= 9) & (order[2] <= 9))) {
    stop(paste0("Refinement orders in both directions should be smaller or equal to 9, please change ",
                order, " to be upper bounded by c(9, 9)." ))
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  L <- (order - 1) / 2
  d <- dim(D[[1]])[1]
  J <- (if(missing(jmax)) length(D) else jmax)
  J0_2D <- sum(sapply(1:length(D), function(j) any(dim(D[[j]]) == 1)))
  if(!isTRUE((J > 0) & isTRUE(all.equal(as.integer(J), J)))) {
    stop("'jmax' should be an integer larger than zero")
  }

  ## Reconstruct midpoints
  m1 <- M0
  for (j in 0:(J - 1)) {
    if(j < length(D)) {
      if(dim(D[[j + 1]])[3] == 1){
        ## Refine 1D
        L1 <- ifelse(order[2] > 2^j, floor((2^j - 1) / 2), L[2])
        tm1 <- impute_C(array(m1[, , 1, ], dim = c(d, d, 2^j)), W_1D[[min(L1 + 1, 5)]], L1, T, metric, method)
        ## Reconstruct midpoints 1D
        m1 <- reconstr2D_C(tm1, array(D[[j + 1]], dim = dim(tm1)), J0_2D + j, c(1, 2^(j + 1)), T, metric)
        m1 <- array(m1, dim = c(d, d, 1, 2^(j + 1)))
      } else if(dim(D[[j + 1]])[4] == 1){
        ## Refine 1D
        L1 <- ifelse(order[1] > 2^j, floor((2^j - 1) / 2), L[1])
        tm1 <- impute_C(array(m1[, , , 1], dim = c(d, d, 2^j)), W_1D[[min(L1 + 1, 5)]], L1, T, metric, method)
        ## Reconstruct midpoints 1D
        m1 <- reconstr2D_C(tm1, array(D[[j + 1]], dim = dim(tm1)), J0_2D + j, c(2^(j + 1), 1), T, metric)
        m1 <- array(m1, dim = c(d, d, 2^(j + 1), 1))
      } else {
        ## Refine 2D
        tm1 <- impute2D_R(m1, L, metric, method)
        ## Reconstruct midpoints 2D
        m1 <- reconstr2D_C(array(tm1, dim = c(d, d, dim(tm1)[3] * dim(tm1)[4])),
                           array(D[[j + 1]], dim = c(d, d, dim(D[[j + 1]])[3] * dim(D[[j + 1]])[4])), 2 * j,
                           c(dim(tm1)[3], dim(tm1)[4]), T, metric)
        m1 <- array(m1, dim = dim(tm1))
      }
    } else {
      ## Refine 2D outside of range of scales
      tm1 <- impute2D_R(m1, L, metric, method)
      ## Reconstruct midpoints 2D
      m1 <- reconstr2D_C(array(tm1, dim = c(d, d, dim(tm1)[3] * dim(tm1)[4])),
                         array(dim = c(d, d, dim(tm1)[3] * dim(tm1)[4])), 2 * j,
                         c(dim(tm1)[3], dim(tm1)[4]), F, metric)
      m1 <- array(m1, dim = dim(tm1))
    }
  }

  ## Transform back to manifold
  if(return_val == "f"){
    if(metric == "logEuclidean" | metric == "Cholesky" | metric == "rootEuclidean") {
      m1 <- array(Ptransf2D_C(array(m1, dim = c(d, d, dim(m1)[3] * dim(m1)[4])), T, chol_bias, metric), dim = dim(m1))
    }
  }

  return(m1)
}


