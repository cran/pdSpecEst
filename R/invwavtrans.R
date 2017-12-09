#' Inverse average-interpolation 1D wavelet transform
#'
#' \code{InvWavTransf1D} computes the inverse intrinsic average-interpolation (AI) wavelet
#' transform of an array of coarsest-scale HPD midpoints combined with a pyramid of Hermitian
#' wavelet coefficients as explained in (Chau and von Sachs, 2017). This is the inverse operation
#' of the function \code{\link{WavTransf1D}}.
#'
#' The input list of arrays \code{D} and array \code{M0} correspond to a pyramid of wavelet coefficients and
#' the coarsest-scale HPD midpoints respectively, both are structured in the same way as in the output of
#' \code{WavTransf1D}. As in the forward AI wavelet transform, if the refinement order is an odd integer smaller or
#' equal to 9, the function computes the inverse wavelet transform using a fast wavelet refinement scheme based on
#' weighted geometric averages with pre-determined weights as explained in (Chau and von Sachs, 2017a). If the
#' refinement order is an odd integer larger than 9, the wavelet refinement scheme is based on intrinsic
#' polynomial prediction using Neville's algorithm on the Riemannian manifold. The function computes the inverse
#' intrinsic AI wavelet transform equipping the space of HPD matrices with one of the following metrics: (i) Riemannian
#' metric (default) as in (Bhatia, 2009, Chapter 6), (ii) log-Euclidean metric, the Euclidean inner product between
#' matrix logarithms, (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv)
#' Euclidean metric and (v) root-Euclidean metric. The default choice (Riemannian) has several appealing properties
#' not shared by the other metrics, see (Chau and von Sachs, 2017a) for more details.
#'
#' @param D a list of arrays containing the pyramid of wavelet coefficients, where each array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients from the finest wavelet scale \code{j = jmax} up to the coarsest
#' wavelet scale \code{j = 0}. This is the same format as the \code{$D} component given as output by
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
#' \code{"Euclidean"}. The inverse intrinsic AI wavelet transform fundamentally relies on the chosen metric.
#' @param progress should a console progress bar be displayed? Defaults to \code{progress = F}.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples(2^8, example = "bumps")
#' P.wt <- WavTransf1D(P$f) ## forward transform
#' P.f <- InvWavTransf1D(P.wt$D, P.wt$M0) ## backward transform
#' all.equal(P.f, P$f)
#'
#' @return Returns a (\eqn{d, d, m})-dimensional array corresponding to a curve of length \eqn{m} of
#' (\eqn{d,d})-dimensional HPD matrices.
#'
#' @seealso \code{\link{WavTransf1D}}, \code{\link{pdSpecEst1D}}, \code{\link{pdNeville}}
#'
#' @references Chau, J. and von Sachs, R. (2017) \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
InvWavTransf1D <- function(D, M0, order = 5, jmax, periodic = F, metric = "Riemannian", progress = F, ...) {

  ## Set variables
  dots = list(...)
  return_val = (if(is.null(dots$return_val)) "manifold" else dots$return_val)
  method = (if(is.null(dots$method)) ifelse(order <= 9, "weights",  "neville") else dots$method)
  chol.bias = (if(is.null(dots$chol.bias)) F else dots$chol.bias)

  if (!(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order = 5
  }
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  L = (order - 1) / 2
  L_round = 2 * ceiling(L / 2)
  d = nrow(D[[1]][, , 1])
  J = (if(missing(jmax)) length(D) else jmax)

  ## Reconstruct scaling coefficients
  m1 <- M0

  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:(J - 1)) {
    tm1 <- AIRefine1D(m1, L, method = method, inverse = T, metric = metric)
    if(periodic){
      tm1 <- tm1[, , ifelse(j > 0, L_round, 2 * floor(L / 2)) + 1:(2^(j + 1) + 2 * L_round), drop = F]
    }
    reconstr_even <- function(i) {
      if((j + 1) <= length(D)){
        if (any(c(D[[j + 1]][, , i]) != 0)) {
          if(metric == "Riemannian"){
            m1_i <- Expm(tm1[, , 2 * i], 2^(j/2) * D[[j + 1]][, , i])
          } else{
            m1_i <- 2^(j/2) * D[[j + 1]][, , i] + tm1[, , 2 * i]
          }
        } else {
          m1_i <- tm1[, , 2 * i]
        }
      } else {
        m1_i <- tm1[, , 2 * i]
      }
      return(m1_i)
    }

    grid_j <- (if((j + 1) <= length(D)) 1:dim(D[[j + 1]])[3] else 1:(dim(tm1)[3]/2))
    m2 <- array(dim = c(d, d, dim(tm1)[3]))
    m2[, , c(F, T)] <- sapply(grid_j, reconstr_even, simplify = "array")
    L_i <- ifelse(periodic, L_round / 2 + ifelse((j > 0) | (L %% 2 == 0), 0, -1), 0)

    if(metric == "Riemannian"){
      m2[, , c(T, F)] <- sapply(grid_j, function(i) (m1[, , i + L_i] %*%
                                                       solve(m2[, , 2 * i])) %*% m1[, , i + L_i], simplify = "array")
    } else{
      m2[, , c(T, F)] <- sapply(grid_j, function(i) 2 * m1[, , i + L_i] - m2[, , 2 * i],
                                simplify = "array")
    }
    m1 <- m2

    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / J))
    }
  }
  if(progress){
    close(pb)
  }

  if(return_val == "manifold"){
    m1 <- (if(metric == "logEuclidean"){
      sapply(1:dim(m1)[3], function(i) Expm(diag(d), m1[, , i]), simplify = "array")
    } else if(metric == "Cholesky"){
      sapply(1:dim(m1)[3], function(i) Chol(m1[, , i], inverse = T, bias.corr = chol.bias), simplify = "array")
    } else if(metric == "rootEuclidean") {
      sapply(1:dim(m1)[3], function(i) t(Conj(m1[, , i])) %*% m1[, , i], simplify = "array")
    } else m1)
  }

  return((if(periodic) m1[, , L_round + 1:2^J] else m1))
}

#' Inverse average-interpolation 2D wavelet transform
#'
#' \code{InvWavTransf2D} computes the inverse intrinsic average-interpolation (AI) wavelet
#' transform of an array of coarsest-scale HPD midpoints combined with a 2D pyramid of Hermitian
#' wavelet coefficients. This is the inverse operation of the function \code{\link{WavTransf2D}}.
#'
#' The input list of arrays \code{D} and array \code{M0} correspond to a 2D pyramid of wavelet coefficients and
#' the coarsest-scale HPD midpoints respectively, both are structured in the same way as in the output of
#' \code{WavTransf2D}. As in the forward AI wavelet transform, if both marginal refinement orders are
#' smaller or equal to 9, the function computes the inverse wavelet transform using a fast wavelet refinement scheme based
#' on weighted geometric averages with pre-determined weights. If one of the marginal refinement order is an odd integer
#' larger than 9, the wavelet refinement scheme is based on intrinsic polynomial surface prediction using Neville's algorithm on the
#' Riemannian manifold (\code{\link{pdNeville}}). By default \code{InvWavTransf2D} computes the inverse intrinsic 2D AI wavelet transform
#' equipping the space of HPD matrices with (i) the Riemannian metric. Instead, the space of HPD matrices can also be
#' equipped with one of the following metrics; (ii) log-Euclidean metric, the Euclidean inner product between matrix logarithms,
#' (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv) Euclidean metric and
#' (v) root-Euclidean metric. The default choice (Riemannian) has several appealing properties not shared by the
#' other metrics, see (Chau and von Sachs, 2017a) for more details.
#'
#' @param D a list of arrays containing the 2D pyramid of wavelet coefficients, where each array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients from the finest wavelet scale \code{j = jmax} up to the coarsest
#' wavelet scale \code{j = 0}. This is the same format as the \code{$D} component given as output by
#'  \code{\link{WavTransf2D}}.
#' @param M0 a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.
#' This is the same format as the \code{$M0} component given as output by \code{\link{WavTransf2D}}.
#' @param order a 2-dimensional numeric vector of odd integers larger or equal to 1 corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(5, 5)}. Note that if \code{max(order) > 9},
#' the computational cost significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param jmax the maximum scale (resolution) up to which the 2D surface of HPD midpoints (i.e. scaling coefficients) are
#' reconstructed. If \code{jmax} is not specified it is set equal to the resolution in the finest wavelet scale
#' \code{jmax = length(D)}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The inverse intrinsic 2D AI wavelet transform fundamentally relies on the chosen metric.
#' @param progress should a console progress bar be displayed? Defaults to \code{progress = T}.
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
#' @references Chau, J. and von Sachs, R. (2017) \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
InvWavTransf2D <- function(D, M0, order = c(3, 3), jmax, metric = "Riemannian", progress = T, ...) {

  ## Set variables
  dots = list(...)
  return_val = (if(is.null(dots$return_val)) "manifold" else dots$return_val)
  method = (if(is.null(dots$method)) "weights" else dots$method)
  chol.bias = (if(is.null(dots$chol.bias)) F else dots$chol.bias)
  if (!isTRUE((order[1] %% 2 == 1) & (order[2] %% 2 == 1))) {
    warning("Refinement orders in both directions should be odd integers, by default set to c(5,5).")
    order = c(3, 3)
  }
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  L = (order - 1) / 2
  d = dim(D[[1]])[1]
  J = (if(missing(jmax)) length(D) else jmax)

  ## Reconstruct scaling coefficients
  m1 <- M0
  J0_2D <- sum(sapply(1:length(D), function(j) any(dim(D[[j]]) == 1)))

  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:(J - 1)) {

    if((j + 1) <= length(D)){
      if(dim(D[[j + 1]])[3] == 1){
        ## Refine 1D
        tm1 <- array(AIRefine1D(array(m1[, , 1, ], dim = c(d, d, 2^j)), L[2], method = method,
                                inverse = T, metric = metric), dim = c(d, d, 1, 2^(j + 1)))
      } else if (dim(D[[j + 1]])[4] == 1){
        ## Refine 1D
        tm1 <- array(AIRefine1D(array(m1[, , , 1], dim = c(d, d, 2^j)), L[1], method = method,
                                inverse = T, metric = metric), dim = c(d, d, 2^(j + 1), 1))
      }else{
        ## Refine 2D
        tm1 <- AIRefine2D(m1, L, method = method, metric = metric)
      }
    } else{
      ## Refine 2D
      tm1 <- AIRefine2D(m1, L, method = method, metric = metric)
    }

    reconstr <- function(i1, i2) {
      if((j + 1) <= length(D)) {
        if(any(c(D[[j + 1]][, , i1, i2]) != 0)){
          if(metric == "Riemannian"){
            m1_i <- Expm(tm1[, , i1, i2], ifelse(any(dim(D[[j + 1]]) == 1),
                                        2^((J0_2D + j)/2), 2^j)  * D[[j + 1]][, , i1, i2])
          } else{
            m1_i <- ifelse(any(dim(D[[j + 1]]) == 1), 2^((J0_2D + j)/2), 2^j) *
              D[[j + 1]][, , i1, i2] + tm1[, , i1, i2]
          }
        } else{
          m1_i <- tm1[, , i1, i2]
        }
      } else {
        m1_i <- tm1[, , i1, i2]
      }
      return(m1_i)
    }

    if((j + 1) <= length(D)){
      grid_j <- expand.grid(1:dim(D[[j + 1]])[3], 1:dim(D[[j + 1]])[4])
    } else{
      grid_j <- expand.grid(1:(dim(D[[length(D)]])[3] * 2^((j + 1) - length(D))),
                            1:(dim(D[[length(D)]])[4] * 2^((j + 1) - length(D))))
    }
    m1 <- array(c(mapply(function(i1, i2) reconstr(i1, i2), grid_j$Var1, grid_j$Var2, SIMPLIFY = "array")),
                dim = c(d, d, attributes(grid_j)$out.attrs$dim[1], attributes(grid_j)$out.attrs$dim[2]))

    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / J))
    }
  }
  if(progress){
    close(pb)
  }

  ## Back transform to manifold
  if(return_val == "manifold"){
    m1 <- (if(metric == "logEuclidean"){
      array(apply(m1, c(3, 4), function(M) Expm(diag(d), M)),
            dim = c(d, d, dim(m1)[3], dim(m1)[4]))
    } else if(metric == "Cholesky"){
      array(apply(m1, c(3, 4), function(M) Chol(M, inverse = T, bias.corr = chol.bias)),
            dim = c(d, d, dim(m1)[3], dim(m1)[4]))
    } else if(metric == "rootEuclidean") {
      array(apply(m1, c(3, 4), function(M) t(Conj(M)) %*% M),
            dim = c(d, d, dim(m1)[3], dim(m1)[4]))
    } else m1)
  }

  return(m1)
}

