#' Forward average-interpolation 1D wavelet transform
#'
#' \code{WavTransf1D} computes the forward intrinsic average-interpolation (AI) wavelet transform of a
#' curve in the manifold of HPD matrices equipped with a metric specified by the user (e.g. the Riemannian
#' metric) as described in (Chau and von Sachs, 2017).
#'
#' The input array \code{P} corresponds to a discretized curve of \eqn{(d,d)}-dimensional HPD matrices of
#' dyadic length. \code{WavTransf1D} then computes the intrinsic AI wavelet transform of \code{P} based on
#' the given refinement order and the chosen metric. If the refinement order is an odd integer smaller or
#' equal to 9, the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
#' geometric averages with pre-determined weights as explained in (Chau and von Sachs, 2017a). If the
#' refinement order is an odd integer larger than 9, the wavelet refinement scheme is based on intrinsic
#' polynomial prediction using Neville's algorithm on the Riemannian manifold.
#' The function computes the intrinsic AI wavelet transform equipping the space of HPD matrices with
#' one of the following metrics: (i) Riemannian metric (default) as in (Bhatia, 2009, Chapter 6),
#' (ii) log-Euclidean metric, the Euclidean inner product between matrix logarithms,
#' (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv) Euclidean metric and
#' (v) root-Euclidean metric. The default choice (Riemannian) has several appealing properties not shared by the
#' other metrics, see (Chau and von Sachs, 2017a) for more details.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of HPD matrices, with \eqn{m = 2^J} for some \eqn{J > 0}.
#' @param order an odd integer larger or equal to 1 corresponding to the order of the intrinsic AI refinement scheme,
#' defaults to \code{order = 5}. Note that if \code{order > 9}, the computational cost
#' significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified it is set equal to the maximum possible scale \code{jmax = J-1}.
#' @param periodic a logical value determining whether the curve of HPD matrices can be reflected at the boundary for
#' improved wavelet refinement schemes near the boundaries of the domain. This is useful for spectral matrix estimation,
#' where the spectral matrix is a symmetric and periodic curve in the frequency domain. Defaults to \code{periodic = F}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic AI wavelet transform fundamentally relies on the chosen metric.
#' @param progress should a console progress bar be displayed? Defaults to \code{progress = F}.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples(2^8, example = "bumps")
#' P.wt <- WavTransf1D(P$f, periodic = FALSE)
#'
#' @return The function returns a list with three components:
#' \item{D }{the pyramid of wavelet coefficients. This is a list of arrays, where each array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients from the finest wavelet scale \code{j = jmax} up to the coarsest
#' wavelet scale \code{j = 0}.}
#' \item{D.white }{the pyramid of whitened wavelet coefficients. The structure of \code{D.white} is the same as
#' \code{D}, but with the wavelet coefficients replaced by their whitened counterparts as explained in
#' (Chau and von Sachs, 2017).}
#' \item{M0 }{a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the midpoint pyramid.}
#'
#' @seealso \code{\link{InvWavTransf1D}}, \code{\link{pdSpecEst1D}}, \code{\link{pdNeville}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
WavTransf1D <- function(P, order = 5, jmax, periodic = F, metric = "Riemannian", progress = F, ...) {

  ## Set variables
  J = log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                " to dyadic number."))
  }
  if (!isTRUE(order %% 2 == 1)) {
    warning("Refinement order should be an odd integer, by default set to 5")
    order = 5
  }
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  dots = list(...)
  method = (if(is.null(dots$method)) ifelse(order <= 9, "weights", "neville") else dots$method)
  d = dim(P)[1]
  L = (order - 1) / 2
  L_round = 2 * ceiling(L / 2)

  P <- (if(metric == "logEuclidean"){
    sapply(1:2^J, function(i) Logm(diag(d), P[, , i]), simplify = "array")
  } else if(metric == "Cholesky"){
    sapply(1:2^J, function(i) Chol(P[, , i]), simplify = "array")
  } else if(metric == "rootEuclidean"){
    sapply(1:2^J, function(i) Sqrt(P[, , i]), simplify = "array")
  } else P)

  if(periodic & (order > 1)){
    P_per <- array(if(L %% 2 == 0) {
      c(rep(c(P, P[, , c(2^J, 2^J:2)]), times = L), P)
    } else {
      c(P[, , c(2^J, 2^J:2)], rep(c(P, P[, , c(2^J, 2^J:2)]), times = L))
    }, dim = c(d, d, (2 * L + 1) * 2^J))
  } else{
    P_per <- P
  }

  M <- list()
  for (j in J:0) {
    if (j == J) {
      Mper <- P_per
    } else {
      if(!(metric == "Riemannian")){
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) 0.5 * (Mper[, , 2 * i - 1] +  Mper[, , 2 * i]),
                       simplify = "array")
      } else {
        Mper <- sapply(1:(dim(Mper)[3]/2), function(i) Mid(Mper[, , 2 * i - 1], Mper[, , 2 * i]),
                       simplify = "array")
      }
    }
    M[[j + 1]] <- if(periodic & (j > 0)){
      Mper[, , L * (2^j) - L_round + 1:(2^j + 2 * L_round)]
    } else{
      Mper
    }
  }

  D <- tM <- D.white <- list()
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("'jmax' cannot exceed maximum scale j = ", J - 1))
    jmax <- J - 1
  }

  ## Compute wavelet transform
  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:jmax) {
    tm1 <- AIRefine1D(M[[j + 1]], L, method, inverse = F, metric = metric)
    tM[[j + 1]] <- (if(periodic) tm1[, , L_round / 2 + ifelse(j > 0 | L %% 2 == 0, 0, -1) + 1:(2^j + L_round), drop = F] else tm1)
    if(!(metric == "Riemannian")){
      D[[j + 1]] <- D.white[[j + 1]] <- sapply(1:dim(tM[[j + 1]])[3], function(l) 2^(-j/2) * (M[[j + 2]][, , 2 * l] - tM[[j + 1]][, , l]),
                                               simplify = "array")
    } else{
      iSqrt_tm1 <- sapply(1:dim(tM[[j + 1]])[3], function(l) iSqrt(tM[[j + 1]][, , l]), simplify = "array")
      D[[j + 1]] <- sapply(1:dim(tM[[j + 1]])[3], function(l) 2^(-j/2) * Logm(tM[[j + 1]][, , l], M[[j + 2]][, , 2 * l]),
                           simplify = "array")
      D.white[[j + 1]] <- sapply(1:dim(D[[j + 1]])[3], function(l) (iSqrt_tm1[, , l] %*% D[[j + 1]][, , l]) %*%
                                   iSqrt_tm1[, , l], simplify = "array")
    }
    names(D)[j + 1] <- names(D.white)[j + 1] <- paste0("D.scale", j)
    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / (jmax + 1)))
    }
  }
  if(progress){
    close(pb)
  }
  return(list(D = D, D.white = D.white, M0 = M[[1]]))
}

#' Forward average-interpolation 2D wavelet transform
#'
#' \code{WavTransf2D} computes the forward intrinsic average-interpolation (AI) wavelet transform of a
#' rectangular surface in the manifold of HPD matrices equipped with a metric specified by the user
#' (e.g. the Riemannian metric).
#'
#' The 4-dimensional array \code{P} corresponds to a discretized rectangular surface of \eqn{(d,d)}-dimensional
#' HPD matrices. The rectangular surface is of size \eqn{n_1} by \eqn{n_2}, where both \eqn{n_1} and
#' \eqn{n_2} are supposed to be dyadic numbers. \code{WavTransf2D} then computes the intrinsic AI wavelet transform
#' of \code{P} based on the given refinement orders and the chosen metric. If both marginal refinement orders are
#' smaller or equal to 9, the function computes the wavelet transform using a fast wavelet refinement scheme based on weighted
#' geometric averages with pre-determined weights. If one of the marginal refinement order is an odd integer larger than 9,
#' the wavelet refinement scheme is based on intrinsic polynomial surface prediction using Neville's algorithm on the
#' Riemannian manifold (\code{\link{pdNeville}}). By default \code{WavTransf2D} computes the intrinsic 2D AI wavelet transform
#' equipping the space of HPD matrices with (i) the Riemannian metric. Instead, the space of HPD matrices can also be
#' equipped with one of the following metrics; (ii) log-Euclidean metric, the Euclidean inner product between matrix logarithms,
#' (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv) Euclidean metric and
#' (v) root-Euclidean metric. The default choice (Riemannian) has several appealing properties not shared by the
#' other metrics, see (Chau and von Sachs, 2017a) for more details.
#'
#' @param P a (\eqn{d,d,n1,n2})-dimensional array of Hermitian PD matrices, with \eqn{n_1 = 2^{J_1}} and \eqn{n_2 = 2^{J_2}}
#' for some \eqn{J_1, J_2 > 0}.
#' @param order a 2-dimensional numeric vector of odd integers larger or equal to 1 corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}. Note that the computational cost
#' significantly increases if \code{max(order) > 9} as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified it is set equal to the maximum possible scale \code{jmax = max(J1, J2) - 1}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can be one of: \code{"Riemannian"}, \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic AI wavelet transform fundamentally relies on the chosen metric.
#' @param progress should a console progress bar be displayed? Defaults to \code{progress = T}.
#' @param ... additional arguments for internal use.
#'
#' @examples
#' P <- rExamples2D(c(2^4, 2^4), 2, example = "tvar")
#' P.wt <- WavTransf2D(P$f)
#'
#' @return The function returns a list with three components:
#' \item{D }{the 2D pyramid of wavelet coefficients. This is a list of arrays, where each 4-dimensional array contains the
#' (\eqn{d,d})-dimensional wavelet coefficients in a 2D grid of locations from the finest wavelet scale \code{j = jmax}
#' up to the coarsest wavelet scale \code{j = 0}.}
#' \item{D.white }{the 2D pyramid of whitened wavelet coefficients. The structure of \code{D.white} is the same as
#' \code{D}, but with the wavelet coefficients replaced by their whitened counterparts as explained in
#' (Chau and von Sachs, 2017).}
#' \item{M0 }{a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.}
#'
#' @seealso \code{\link{InvWavTransf2D}}, \code{\link{pdSpecEst2D}}, \code{\link{pdNeville}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
WavTransf2D <- function(P, order = c(3, 3), jmax, metric = "Riemannian", progress = T, ...) {

  ## Set variables
  J1 = log2(dim(P)[3])
  J2 = log2(dim(P)[4])
  J = max(J1, J2)
  J0_2D = abs(J1 - J2)
  if (!isTRUE(all.equal(as.integer(J1), J1) & all.equal(as.integer(J2), J2))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3], " or ",
                dim(P)[4], " to dyadic number."))
  }
  if (!isTRUE((order[1] %% 2 == 1) & (order[2] %% 2 == 1))) {
    warning("Refinement orders in both directions should be odd integers, by default set to c(5,5).")
    order = c(3, 3)
  }
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  dots = list(...)
  method = (if(is.null(dots$method)) "weights" else dots$method)
  L = (order - 1) / 2
  d = dim(P)[1]

  P <- (if(metric == "logEuclidean"){
    array(apply(P, c(3, 4), function(Pi) Logm(diag(d), Pi)), dim = c(d, d, 2^J1, 2^J2))
  } else if(metric == "Cholesky"){
    array(apply(P, c(3, 4), function(Pi) Chol(Pi)), dim = c(d, d, 2^J1, 2^J2))
  } else if(metric == "rootEuclidean"){
    array(apply(P, c(3, 4), function(Pi) Sqrt(Pi)), dim = c(d, d, 2^J1, 2^J2))
  } else P)

  Mper_2D <- function(Mper, j1, j2) {
    l1 <- t(sapply(1:2^j1, function(i) c(1, 2) + 2 * (i - 1)))
    l2 <- t(sapply(1:2^j2, function(i) c(1, 2) + 2 * (i - 1)))
    grid <- expand.grid(1:2^j1, 1:2^j2)
    Mper <- array(c(mapply(function(i1, i2) pdMean(array(c(Mper[,,l1[i1,],l2[i2,]]), dim = c(d, d, 4)),
                                                   metric = ifelse(metric == "Riemannian", "Riemannian", "Euclidean")),
                           grid$Var1, grid$Var2, SIMPLIFY = "array")), dim = c(d, d, 2^j1, 2^j2))
    return(Mper)
  }

  grid_j <- cbind((J1:0)[1:(J + 1)], (J2:0)[1:(J + 1)])

  M <- list()
  for (j in J:0) {
    if(j == J){
      Mper <- unname(P)
    } else if(j >= J0_2D) {
      Mper <- Mper_2D(Mper, grid_j[J + 1 - j, 1], grid_j[J + 1 - j, 2])
    } else {
      if(is.na(grid_j[J + 1 - j, 1])) {
        j2 <- grid_j[J + 1 - j, 2]
        Mper <- array(c(sapply(1:2^j2, function(i){
          if(metric == "Riemannian") Mid(Mper[, , , 2 * i - 1], Mper[, , , 2 * i]) else{
            0.5 * (Mper[, , , 2 * i - 1] + Mper[, , , 2 * i])
          }})), dim = c(d, d, 1, 2^j2))
      } else if(is.na(grid_j[J + 1 - j, 2])) {
        j1 <- grid_j[J + 1 - j, 1]
        Mper <- array(c(sapply(1:2^j1, function(i){
          if(metric == "Riemannian") Mid(Mper[, , 2 * i - 1, ], Mper[, , 2 * i, ]) else{
            0.5 * (Mper[, , 2 * i - 1, ] + Mper[, , 2 * i, ])
          }})), dim = c(d, d, 2^j1, 1))
      }
    }
    M[[j + 1]] <- Mper
  }

  D <- tM <- D.white <- list()
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("'jmax' cannot exceed maximum scale j = ", J - 1))
    jmax <- J - 1
  }

  ## Compute 2D wavelet transform
  if(progress){
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for (j in 0:jmax) {

    if(is.na(grid_j[J + 1 - j, 1])) {

      ## Refine 1D
      tm1 <- array(AIRefine1D(array(M[[j + 1]][, , 1, ], dim = c(d, d, 2^j)), L[2], method = method, inverse = T, metric = metric),
                              dim = c(d, d, 1, 2^(j + 1)))

      ## Construct (correctly scaled) wavelet coefficients
      if(!(metric == "Riemannian")){
        D[[j + 1]] <- D.white[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) 2^((-J0_2D - j)/2) *
                                                  (M[[j + 2]][, , 1, i] - tm1[, , 1, i]))), dim = c(d, d, 1, 2^(j + 1)))
      } else{
        iSqrt_tm1 <- sapply(1:2^(j + 1), function(i) iSqrt(tm1[, , 1, i]), simplify = "array")
        D[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) 2^((-J0_2D - j)/2) *
                                       Logm(tm1[, , 1, i], M[[j + 2]][, , 1, i]))), dim = c(d, d, 1, 2^(j + 1)))
        D.white[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) (iSqrt_tm1[, , i] %*% D[[j + 1]][, , 1, i]) %*%
                                             iSqrt_tm1[, , i])), dim = c(d, d, 1, 2^(j + 1)))
      }
    } else if(is.na(grid_j[J + 1 - j, 2])) {

      ## Refine 1D
      tm1 <- array(AIRefine1D(array(M[[j + 1]][, , , 1], dim = c(d, d, 2^j)), L[1], method = method, inverse = T, metric = metric),
                              dim = c(d, d, 2^(j + 1), 1))

      ## Construct (correctly scaled) wavelet coefficients
      if(!(metric == "Riemannian")){
        D[[j + 1]] <- D.white[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) 2^((-J0_2D - j)/2) *
                                                           (M[[j + 2]][, , i, 1] - tm1[, , i, 1]))), dim = c(d, d, 2^(j + 1), 1))
      } else{
        iSqrt_tm1 <- sapply(1:2^(j + 1), function(i) iSqrt(tm1[, , i, 1]), simplify = "array")
        D[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) 2^((-J0_2D - j)/2) *
                                       Logm(tm1[, , i, 1], M[[j + 2]][, , i, 1]))), dim = c(d, d, 2^(j + 1), 1))
        D.white[[j + 1]] <- array(c(sapply(1:2^(j + 1), function(i) (iSqrt_tm1[, , i] %*% D[[j + 1]][, , i, 1]) %*%
                                             iSqrt_tm1[, , i])), dim = c(d, d, 2^(j + 1), 1))
      }
    } else {

      ## Refine 2D
      tm1 <- AIRefine2D(M[[j + 1]], L, method = method, metric = metric)
      grid <- expand.grid(1:dim(tm1)[3], 1:dim(tm1)[4])

      if(!(metric == "Riemannian")){
        ## Euclidean wavelet coeff's
        D[[j + 1]] <- D.white[[j + 1]] <- array(c(mapply(function(i1, i2) 2^(-j) * (M[[j + 2]][, , i1, i2] -
                                                tm1[, , i1, i2]), grid$Var1, grid$Var2)), dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
      } else{
        ## Riemannian wavelet coeff's
        iSqrt_tm1 <- array(apply(tm1, c(3, 4), function(tM) iSqrt(tM)), dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
        D[[j + 1]] <- array(c(mapply(function(i1, i2) 2^(-j) * Logm(tm1[, , i1, i2], M[[j + 2]][, , i1, i2]),
                                     grid$Var1, grid$Var2)), dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
        D.white[[j + 1]] <- array(c(mapply(function(i1, i2) (iSqrt_tm1[, , i1, i2] %*% D[[j + 1]][, , i1, i2]) %*%
                                             iSqrt_tm1[, , i1, i2], grid$Var1, grid$Var2)), dim = c(d, d, dim(tm1)[3], dim(tm1)[4]))
      }
    }
    tM[[j + 1]] <- tm1
    names(D)[j + 1] <- names(D.white)[j + 1] <- paste0("D.scale", j)
    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (j + 1) / (jmax + 1)))
    }
  }
  if(progress){
    close(pb)
  }
  return(list(D = D, D.white = D.white, M0 = M[[1]]))
}
