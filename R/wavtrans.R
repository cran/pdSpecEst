#' Forward AI wavelet transform for curve of HPD matrices
#'
#' \code{WavTransf1D} computes a forward intrinsic average-interpolating (AI) wavelet transform for a
#' curve in the manifold of HPD matrices equipped with a metric specified by the user, such as the
#' affine-invariant Riemannian metric, as described in \insertCite{CvS17}{pdSpecEst} and Chapter 3 of
#' \insertCite{C18}{pdSpecEst}.
#'
#' The input array \code{P} corresponds to a discretized curve of \eqn{(d,d)}-dimensional HPD matrices of
#' dyadic length. \code{WavTransf1D} then computes the intrinsic AI wavelet transform of \code{P} based on
#' the given refinement order and the chosen metric. If the refinement order is an odd integer smaller or
#' equal to 9, the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
#' intrinsic averages with pre-determined weights as explained in \insertCite{CvS17}{pdSpecEst} and Chapter 3 of
#' \insertCite{C18}{pdSpecEst}. If the refinement order is an odd integer larger than 9, the wavelet refinement
#' scheme uses intrinsic polynomial prediction based on Neville's algorithm in the Riemannian manifold (via \code{\link{pdNeville}}).\cr
#' The function computes the intrinsic AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of HPD matrices, corresponding to a sequence of \eqn{(d,d)}-dimensional HPD matrices
#' of length \eqn{m}, with \eqn{m = 2^J} for some \eqn{J > 0}.
#' @param order an odd integer larger or equal to 1 corresponding to the order of the intrinsic AI refinement scheme,
#' defaults to \code{order = 5}. Note that if \code{order > 9}, the computational cost
#' significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified, it is set equal to the maximum possible scale \code{jmax = J-1}, where \code{J = log2(m)}.
#' @param periodic a logical value determining whether the curve of HPD matrices can be reflected at the boundary for
#' improved wavelet refinement schemes near the boundaries of the domain. This is useful for spectral matrix estimation,
#' in which case the spectral matrix is a symmetric and periodic curve in the frequency domain. Defaults to \code{periodic = FALSE}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"},
#' \code{"Euclidean"} or \code{"Riemannian-Rahman"}. See also the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples1D(2^8, example = "bumps")
#' P.wt <- WavTransf1D(P$f, periodic = FALSE)
#'
#' @return The function returns a list with three components:
#' \item{D }{ the pyramid of wavelet coefficients. This is a list of arrays, where each array contains the
#' (\eqn{d,d})-dimensional Hermitian wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to
#' the finest wavelet scale \code{j = jmax}}.
#' \item{D.white }{ the pyramid of whitened wavelet coefficients. The structure of \code{D.white} is the same as
#' \code{D}, but with the wavelet coefficients replaced by their whitened counterparts as explained in
#' \insertCite{CvS17}{pdSpecEst}.}
#' \item{M0 }{ a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the midpoint pyramid.}
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and (depending on the
#' specified metric) may fail if matrices are close to being singular.
#'
#' @seealso \code{\link{InvWavTransf1D}}, \code{\link{pdSpecEst1D}}, \code{\link{pdNeville}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
WavTransf1D <- function(P, order = 5, jmax, periodic = FALSE, metric = "Riemannian", ...) {

  ## Initialize parameters
  n <- dim(P)[3]
  J <- log2(n)
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", n,
                " to dyadic number."))
  }
  if (!isTRUE(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order <- 5
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean", "Riemannian-Rahman"))
  dots <- list(...)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  d <- dim(P)[1]
  L <- (order - 1) / 2
  L_round <- 2 * ceiling(L / 2)
  N <- (2 * L + 1) * n
  Nj <- as.integer(ifelse(periodic & (order > 1), N, n) / (2^(0:J)))

  ## Compute midpoint pyramid
  Mper <- wavPyr_C(P, ifelse(periodic & (order > 1), L, 0), J, Nj, ifelse(metric == "Riemannian-Rahman", "Riemannian", metric))
  M <- list()
  M[[1]] <- Mper[, , sum(head(Nj, J)) + 1:Nj[J + 1], drop = FALSE]
  for(j in 1:J) {
    if(periodic & (order > 1)) {
      M[[j + 1]] <- Mper[, , ifelse(j < J, sum(Nj[1:(J - j)]), 0) + L * 2^j -
                           L_round + 1:(2^j + 2 * L_round), drop = FALSE]
    } else {
      M[[j + 1]] <- Mper[, , ifelse(j < J, sum(Nj[1:(J - j)]), 0) + 1:2^j, drop = FALSE]
    }
  }

  ## Create empty lists for predicted midpoints and wav. coeff's
  D <- Dw <- list()
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("'jmax' cannot exceed maximum scale j = ", J - 1))
    jmax <- J - 1
  }
  
  for (j in 0:jmax) {

    ## Compute predicted midpoints C++
    n_M <- dim(M[[j + 1]])[3]
    L1 <- ifelse(order > n_M, floor((n_M - 1) / 2), L)
    tm1 <- impute_C(M[[j + 1]], W_1D[[min(L1 + 1, 5)]], L1, FALSE, metric, method)[, , 2 * (1:n_M), drop = FALSE]
    tM <- (if(periodic){ tm1[, , L_round / 2 + ifelse(j > 0 | L %% 2 == 0, 0, -1) +
                               1:(2^j + L_round), drop = FALSE] } else tm1)
    ## Compute wavelet coefficients C++
    n_W <- dim(tM)[3]
    W <- wavCoeff_C(tM, M[[j + 2]][, , 2 * (1:n_W), drop = FALSE], j, ifelse(metric == "Riemannian-Rahman", "Riemannian", metric))
    Dw[[j + 1]] <- W[, , 1:n_W, drop = FALSE]
    D[[j + 1]] <- W[, , n_W + 1:n_W, drop = FALSE]
    names(D)[j + 1] <- names(Dw)[j + 1] <- paste0("D.scale", j)
  }
  return(list(D = D, D.white = Dw, M0 = M[[1]]))
}

#' Forward AI wavelet transform for surface of HPD matrices
#'
#' \code{WavTransf2D} computes a forward intrinsic average-interpolation (AI) wavelet transform for a
#' rectangular surface in the manifold of HPD matrices equipped with a metric specified by the user, such as the
#' affine-invariant Riemannian metric, as described in Chapter 5 of \insertCite{C18}{pdSpecEst}.
#'
#' The 4-dimensional array \code{P} corresponds to a discretized rectangular surface of \eqn{(d,d)}-dimensional
#' HPD matrices. The rectangular surface is of size \eqn{n_1} by \eqn{n_2}, where both \eqn{n_1} and
#' \eqn{n_2} are supposed to be dyadic numbers. \code{WavTransf2D} then computes the intrinsic AI wavelet transform
#' of \code{P} based on the given refinement orders and the chosen metric. The marginal refinement orders should be
#' smaller or equal to 9, and the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
#' intrinsic averages with pre-determined weights as explained in Chapter 5 of \insertCite{C18}{pdSpecEst}. By default \code{WavTransf2D}
#' computes the intrinsic 2D AI wavelet transform equipping the space of HPD matrices with (i) the affine-invariant Riemannian metric as
#' detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}. Instead, the space of HPD matrices
#' can also be equipped with one of the following metrics; (ii) the Log-Euclidean metric, the Euclidean inner product between matrix
#' logarithms; (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric and
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{C18}{pdSpecEst} for more details. Note that this comes at the cost of increased computation
#' time in comparison to one of the other metrics.
#'
#' @param P a (\eqn{d,d,n1,n2})-dimensional array of HPD matrices corresponding to a rectangular surface of \eqn{(d,d)}-dimensional HPD matrices
#' of size \eqn{n_1 \times n_2}, with \eqn{n_1 = 2^{J_1}} and \eqn{n_2 = 2^{J_2}} for some \eqn{J_1, J_2 > 0}.
#' @param order a 2-dimensional numeric vector \eqn{(1,1) \le} \code{order} \eqn{\le (9,9)} corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified, it is set equal to the maximum possible scale \code{jmax = max(J1, J2) - 1}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can be one of: \code{"Riemannian"}, \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. See also the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and (depending on the
#' specified metric) may fail if matrices are close to being singular.
#'
#' @examples
#' P <- rExamples2D(c(2^4, 2^4), 2, example = "tvar")
#' P.wt <- WavTransf2D(P$f)
#'
#' @return The function returns a list with three components:
#' \item{D }{ the 2D pyramid of wavelet coefficients. This is a list of arrays, where each 4-dimensional array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients in a 2D grid of locations from the coarsest wavelet scale \code{j = 0}
#' up to the finest wavelet scale \code{j = jmax}.}
#' \item{D.white }{ the 2D pyramid of whitened wavelet coefficients. The structure of \code{D.white} is the same as
#' \code{D}, but with the wavelet coefficients replaced by their whitened counterparts as explained in Chapter 5 of
#' \insertCite{C18}{pdSpecEst}.}
#' \item{M0 }{ a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.}
#'
#' @seealso \code{\link{InvWavTransf2D}}, \code{\link{pdSpecEst2D}}, \code{\link{pdNeville}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
WavTransf2D <- function(P, order = c(3, 3), jmax, metric = "Riemannian", ...) {

  ## Set variables
  J1 <- log2(dim(P)[3])
  J2 <- log2(dim(P)[4])
  J <- max(J1, J2)
  J0_2D <- abs(J1 - J2)
  if (!isTRUE(all.equal(as.integer(J1), J1) & all.equal(as.integer(J2), J2))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3], " or ",
                dim(P)[4], " to dyadic number."))
  }
  if (!isTRUE((order[1] %% 2 == 1) & (order[2] %% 2 == 1))) {
    warning("Refinement orders in both directions should be odd integers, by default set to c(3,3).")
    order <- c(3, 3)
  }
  if (!isTRUE((order[1] <= 9) & (order[2] <= 9))) {
    stop(paste0("Refinement orders in both directions should be smaller or equal to 9, please change ",
                order, " to be upper bounded by c(9, 9)." ))
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  dots <- list(...)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  L <- (order - 1) / 2
  d <- dim(P)[1]

  ## Transform surface
  if(metric == "logEuclidean" | metric == "Cholesky" | metric == "rootEuclidean") {
    P <- array(Ptransf2D_C(array(P, dim = c(d, d, dim(P)[3] * dim(P)[4])), FALSE, FALSE, metric), dim = dim(P))
  }

  ## Construct 2D midpoint pyramid C++
  grid_n <- cbind(2^((J1:0)[1:(J + 1)]), 2^((J2:0)[1:(J + 1)]))
  grid_n[which(is.na(grid_n))] <- 0
  M <- list()
  for (j in J:0) {
    if(j == J){
      Mper <- array(c(P), dim = c(d, d, grid_n[1, 1] * grid_n[1, 2]))
      M[[j + 1]] <- P
    } else {
      Mper <- wavPyr2D_C(Mper, max(grid_n[J - j, 1], 1), max(grid_n[J - j, 2], 1), metric)
      M[[j + 1]] <- array(c(Mper), dim = c(d, d, max(grid_n[J + 1 - j, 1], 1), max(grid_n[J + 1 - j, 2], 1)))
    }
  }

  ## 2D AI Wavelet transform
  D <- Dw <- list()
  W <- lapply(1:length(W_2D), function(i) array(c(aperm(W_2D[[i]], c(3, 4, 1, 2))),
                                  dim = c(dim(W_2D[[i]])[3] * dim(W_2D[[i]])[4], 4)))
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("'jmax' cannot exceed maximum scale j = ", J - 1))
    jmax <- J - 1
  }

  for (j in 0:jmax) {
    if(grid_n[J + 1 - j, 1] < 1) {
      ## Refine 1D C++
      n_M <- dim(M[[j + 1]])[4]
      L1 <- ifelse(order[2] > n_M, floor((n_M - 1) / 2), L[2])
      tm1 <- impute_C(array(M[[j + 1]][, , 1, ], dim = c(d, d, 2^j)), W_1D[[min(L1 + 1, 5)]], L1, TRUE, metric, method)
      ## Compute (correctly scaled) wavelet coefficients C++
      W <- wavCoeff_C(tm1, M[[j + 2]][, , 1, ], J0_2D + j, metric)
      Dw[[j + 1]] <- array(W[, , 1:(2 * n_M)], dim = dim(M[[j + 2]]))
      D[[j + 1]] <- array(W[, , (2 * n_M) + 1:(2 * n_M)], dim = dim(M[[j + 2]]))
    } else if(grid_n[J + 1 - j, 2] < 1) {
      ## Refine 1D C++
      n_M <- dim(M[[j + 1]])[3]
      L1 <- ifelse(order[1] > n_M, floor((n_M - 1) / 2), L[1])
      tm1 <- impute_C(array(M[[j + 1]][, , , 1], dim = c(d, d, 2^j)), W_1D[[min(L1 + 1, 5)]], L1, TRUE, metric, method)
      ## Compute (correctly scaled) wavelet coefficients C++
      W <- wavCoeff_C(tm1, M[[j + 2]][, , , 1], J0_2D + j, metric)
      Dw[[j + 1]] <- array(W[, , 1:(2 * n_M)], dim = dim(M[[j + 2]]))
      D[[j + 1]] <- array(W[, , (2 * n_M) + 1:(2 * n_M)], dim = dim(M[[j + 2]]))
    } else {
      ## Refine 2D
      tm1 <- impute2D_R(M[[j + 1]], L, metric, method)
      ## Compute (correctly scaled) wavelet coefficients C++
      n_tM <- dim(tm1)[3] * dim(tm1)[4]
      W <- wavCoeff_C(array(tm1, dim = c(d, d, n_tM)), array(M[[j + 2]], dim = c(d, d, n_tM)), 2 * j, metric)
      Dw[[j + 1]] <- array(W[, , 1:n_tM], dim = dim(M[[j + 2]]))
      D[[j + 1]] <- array(W[, , n_tM + 1:n_tM], dim = dim(M[[j + 2]]))
    }
    names(D)[j + 1] <- names(Dw)[j + 1] <- paste0("D.scale", j)
  }
  return(list(D = D, D.white = Dw, M0 = M[[1]]))
}
