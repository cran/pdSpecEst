#' Intrinsic 1D wavelet-based spectral matrix estimation
#'
#' \code{pdSpecEst1D} calculates a \eqn{(d,d)}-dimensional HPD wavelet-denoised spectral matrix estimator
#' by: (i) applying an intrinsic 1D AI wavelet transform (\code{\link{WavTransf1D}}) to an initial noisy
#' HPD spectral estimate, (ii) (tree-structured) thresholding of the wavelet coefficients (\code{\link{pdCART}})
#' and (iii) applying an intrinsic inverse 1D AI wavelet transform (\code{\link{InvWavTransf1D}}). The complete
#' estimation procedure is described in detail in (Chau and von Sachs, 2017a).
#'
#' The input array \code{P} corresponds to an initial noisy HPD spectral estimate of the (\eqn{d, d})-dimensional
#' spectral matrix at \code{m} different frequencies, with \eqn{m = 2^J} for some \eqn{J > 0}. This can be e.g.
#' a multitaper HPD periodogram given as output by the function \code{\link{pdPgram}}.\cr
#' \code{P} is transformed to the wavelet domain by the function \code{\link{WavTransf1D}}, which applies an intrinsic
#' 1D AI wavelet transform based on e.g. the Riemannian metric. The noise is removed by tree-structured thresholding
#' of the wavelet coefficients based on the trace of the whitened coefficients as in \code{\link{pdCART}} by
#' minimization of a \emph{complexity penalized residual sum of squares} (CPRESS) criterion in (Donoho, 1997),
#' via a fast tree-pruning algorithm. As in \code{\link{pdCART}}, the sparsity parameter is set equal to \code{alpha}
#' times the universal threshold where the noise variance of the traces of the whitened wavelet
#' coefficients determined from the finest wavelet scale. If the thresholding policy is set to \code{policy = "universal"},
#' the sparsity parameter is set equal to the universal threshold. If the thresholding policy is set to \code{policy = "cv"},
#' a data-adaptive sparsity parameter is computed via two-fold cross-validation as in (Nason, 1996) based on the chosen metric.\cr
#' If \code{return == 'f'} the thresholded wavelet coefficients are transformed back to the frequency domain by
#' the inverse intrinsic 1D AI wavelet transform via \code{\link{InvWavTransf1D}} giving the wavelet-denoised
#' HPD spectral estimate.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of HPD matrices, with \eqn{m} a dyadic number.
#' @param order an odd integer larger or equal to 1 corresponding to the order of the intrinsic AI refinement scheme,
#' defaults to \code{order = 5}. Note that if \code{order > 9}, the computational cost
#' significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param policy a character, one of \code{"universal"} or \code{"cv"}, defaults to \code{policy = "universal"}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic AI wavelet transform fundamentally relies on the chosen metric.
#' @param alpha an optional tuning parameter in the wavelet the thresholding procedure. If \code{policy = "universal"},
#' the sparsity parameter in the tree-structured wavelet thresholding procedure is set to \code{alpha} times the
#' universal threshold, defaults to \code{alpha = 1}.
#' @param return an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with four components:
#' \item{f }{a (\eqn{d,d,m})-dimensional array corresponding to the wavelet-denoised HPD (\eqn{d,d})-dimensional
#' spectral estimate at the \code{m} different frequencies. If \code{!(return == 'f')}, the inverse wavelet transform
#' of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{the pyramid of threshold wavelet coefficients. This is a list of arrays, where each array contains the
#' (\eqn{d,d})-dimensional thresholded wavelet coefficients from the finest wavelet scale \code{j = jmax} up to the coarsest
#' wavelet scale \code{j = 0}.}
#' \item{M0 }{a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the midpoint pyramid.}
#' \item{tree.weights }{a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale.}
#' \item{alpha.opt }{the wavelet thresholding tuning parameter equal to the input argument \code{alpha} if
#' \code{policy = "universal"}; or determined data-adaptively via two-fold cross-validation if \code{policy = "cv"}.}
#'
#' @examples
#' P <- rExamples(2^8, example = "bumps")$per
#' f <- pdSpecEst1D(P)
#'
#' @seealso \code{\link{pdPgram}}, \code{\link{WavTransf1D}}, \code{\link{InvWavTransf1D}}, \code{\link{pdCART}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Nason, G.P. (1996). \emph{Wavelet shrinkage using cross-validation}. Journal of the
#' Royal Statistical Society: Series B, 58, 463-479.
#'
#' @export
pdSpecEst1D <- function(P, order = 5, policy = "universal", metric = "Riemannian", alpha = 1, return = "f", ...) {

  ## Set variables
  dots = list(...)
  tol = (if(is.null(dots$tol)) 0.01 else dots$tol)
  alpha.range = (if(is.null(dots$alpha.range)) c(0.5, 2) else dots$alpha.range)
  tree = (if(is.null(dots$tree)) T else dots$tree)
  periodic = (if(is.null(dots$periodic)) T else dots$periodic)

  policy = match.arg(policy, c("universal", "cv"))
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  J = log2(dim(P)[3])
  d = dim(P)[1]
  B = (if(is.null(dots$B)) d else dots$B)
  progress = (if(is.null(dots$progress)) F else dots$progress)
  J.out = (if(is.null(dots$J.out)) J else dots$J.out)
  jmax = min((if(is.null(dots$jmax)) J - 3 else dots$jmax), J.out - 1)
  jmax.cv = min((if(is.null(dots$jmax.cv)) J - 3 else dots$jmax.cv), J.out - 1)
  bias.corr = (if(is.null(dots$bias.corr)) T else dots$bias.corr)
  return.D = (if(is.null(dots$return.D)) NA else dots$return.D)

  # Manifold bias-correction
  P <- (if((metric == "Riemannian" | metric == "logEuclidean") & bias.corr) {
    B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * P } else P)

  if(policy == "cv"){

    ## Two-fold cross-validation (Nason, 1996)
    Pmid <- (if(metric == "logEuclidean") {
      sapply(1:2^J, function(i) Logm(diag(d), P[, , i]), simplify = "array")
    } else if(metric == "Cholesky") {
      sapply(1:2^J, function(i) Chol(P[, , i]), simplify = "array")
    } else if(metric == "rootEuclidean"){
      sapply(1:2^J, function(i) Sqrt(P[, , i]), simplify = "array")
    } else P)

    for (j in J:(jmax.cv + 1)) {
      Pmid <- sapply(1:(dim(Pmid)[3]/2), function(i) (if(metric == "Riemannian"){
        Mid(Pmid[, , 2 * i - 1], Pmid[, , 2 * i]) } else { 0.5 * (Pmid[, , 2 * i - 1] +
                                                                    Pmid[, , 2 * i]) }), simplify = "array")
    }

    P1 <- list(odd = Pmid[, , c(T, F)], even = Pmid[, , c(F, T)])

    if(metric == "Riemannian"){
      for(m in c(1,2)){
        P1[[m]] <- array(apply(P1[[m]], 3, function(P) Logm(diag(d), P)), dim = dim(P1[[m]]))
      }
    }

    coeff.odd <- WavTransf1D(P1$odd, order, periodic = periodic, metric = "Euclidean")
    coeff.even <- WavTransf1D(P1$even, order, periodic = periodic, metric = "Euclidean")

    cv <- function(alpha){

      D <- list(odd = pdCART(coeff.odd$D, coeff.odd$D.white, order = order, B = B, tree = F,
                             alpha = alpha, periodic = periodic)$D_w,
                even = pdCART(coeff.even$D, coeff.even$D.white, order = order, B = B, tree = F,
                              alpha = alpha, periodic = periodic)$D_w)

      f.hat <- list(odd = InvWavTransf1D(D$odd, coeff.odd$M0, order,
                                         periodic = periodic, metric = "Euclidean", return_val = "tangent"),
                    even = InvWavTransf1D(D$even, coeff.even$M0, order,
                                          periodic = periodic, metric = "Euclidean", return_val = "tangent"))
      ## Predicted points
      f.pred <- list(even = array(c(sapply(1:(dim(f.hat$odd)[3] - 1), function(k) 0.5 * (f.hat$odd[, , k] + f.hat$odd[, , k + 1]),
                                           simplify = "array"), f.hat$odd[, , dim(f.hat$odd)[3]]), dim = dim(f.hat$odd)),
                     odd = array(c(f.hat$even[, , 1], sapply(1:(dim(f.hat$even)[3] - 1), function(k) 0.5 * (f.hat$even[, , k] +
                                                                                                              f.hat$even[, , k + 1]), simplify = "array")), dim = dim(f.hat$even)))

      return(mean(sapply(1:dim(f.hat$odd)[3], function(k) NormF(f.pred$even[, , k] - P1$even[, , k])^2 +
                           NormF(f.pred$odd[, , k] - P1$odd[, , k])^2)))
    }

    ## Find minimum and rescale twice number of data points
    alpha.opt <- gss(alpha.range, cv, tol)

  } else {
    alpha.opt <- alpha
  }

  ## Threshold full data using 'alpha.opt'
  coeff <- (if(policy == "cv"){
    WavTransf1D(P, order, jmax = jmax.cv, periodic = periodic, metric = metric, progress = progress)
  } else WavTransf1D(P, order, jmax = jmax, periodic = periodic, metric = metric, progress = progress))

  coeff.opt <- pdCART(coeff$D, coeff$D.white, alpha = alpha.opt, tree = tree, periodic = periodic,
                      order = order, B = B, return.D = return.D)

  ## Return 'f' or not
  f <- (if(return == "f"){
    InvWavTransf1D(coeff.opt$D_w, coeff$M0, order, jmax = J.out, periodic = periodic, metric = metric, progress = progress)
  } else NULL)

  ## Return whitened coeff's or not
  if(!isTRUE(return.D == "D.white")){
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, alpha.opt = alpha.opt)
  } else{
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, alpha.opt = alpha.opt,
                D.white = coeff.opt$D.white_w)
  }

  return(res)
}

#' Intrinsic 2D wavelet-based time-varying spectral matrix estimation
#'
#' \code{pdSpecEst2D} calculates a \eqn{(d,d)}-dimensional HPD wavelet-denoised time-varying spectral matrix estimator
#' by: (i) applying an intrinsic 2D AI wavelet transform (\code{\link{WavTransf2D}}) to an initial noisy
#' HPD spectral estimate, (ii) (tree-structured) thresholding of the wavelet coefficients (\code{\link{pdCART}})
#' and (iii) applying an intrinsic inverse 2D AI wavelet transform (\code{\link{InvWavTransf2D}}).
#'
#' The input array \code{P} corresponds to an initial noisy HPD time-varying spectral estimate of the (\eqn{d, d})-dimensional
#' spectral matrix at \eqn{m_1 \times m_2} different time-frequency points, with \eqn{m_1, m_2} dyadic numbers. This can be e.g.
#' a multitaper HPD time-varying periodogram given as output by the function \code{\link{pdPgram2D}}.\cr
#' \code{P} is transformed to the wavelet domain by the function \code{\link{WavTransf2D}}, which applies an intrinsic
#' 2D AI wavelet transform based on e.g. the Riemannian metric. The noise is removed by tree-structured thresholding
#' of the wavelet coefficients based on the trace of the whitened coefficients as in \code{\link{pdCART}} by
#' minimization of a \emph{complexity penalized residual sum of squares} (CPRESS) criterion in (Donoho, 1997),
#' via a fast tree-pruning algorithm. As in \code{\link{pdCART}}, the sparsity parameter is set equal to \code{alpha}
#' times the universal threshold where the noise variance of the traces of the whitened wavelet
#' coefficients determined from the finest wavelet scale. \cr
#' If \code{return == 'f'} the thresholded wavelet coefficients are transformed back to the frequency domain by
#' the inverse intrinsic 2D AI wavelet transform via \code{\link{InvWavTransf2D}} giving the wavelet-denoised
#' HPD time-varying spectral estimate.
#'
#' @param P a (\eqn{d,d,m_1,m_2})-dimensional array of HPD matrices, with \eqn{m_1, m_2} both dyadic numbers.
#' @param order a 2-dimensional numeric vector of odd integers larger or equal to 1 corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}. Note that the computational cost
#' significantly increases if \code{max(order) > 9} as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. The intrinsic AI wavelet transform fundamentally relies on the chosen metric.
#' @param alpha an optional tuning parameter in the wavelet the thresholding procedure. If \code{policy = "universal"},
#' the sparsity parameter in the tree-structured wavelet thresholding procedure is set to \code{alpha} times the
#' universal threshold, defaults to \code{alpha = 1}.
#' @param return an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with four components:
#' \item{f }{a (\eqn{d,d,m_1,m_2})-dimensional array corresponding to the wavelet-denoised HPD (\eqn{d,d})-dimensional
#' nonstationary spectral estimate at the \eqn{m1 \times m2} different time-frequency points. If \code{!(return == 'f')},
#' the inverse wavelet transform of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{the 2D pyramid of threshold wavelet coefficients. This is a list of arrays, where each array contains the 2D grid of
#' (\eqn{d,d})-dimensional thresholded wavelet coefficients from the finest wavelet scale \code{j = jmax} up to the coarsest
#' wavelet scale \code{j = 0}.}
#' \item{M0 }{a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.}
#' \item{tree.weights }{a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale.}
#' \item{D.raw }{the 2D pyramid of non-thresholded wavelet coefficients in the same format as the component \code{$D}.}
#'
#' @examples
#' \dontrun{
#'  P <- rExamples2D(c(2^7, 2^7), 3, example = "tvar")$per
#'  f <- pdSpecEst2D(P)
#' }
#'
#' @seealso \code{\link{pdPgram2D}}, \code{\link{WavTransf2D}}, \code{\link{InvWavTransf2D}}, \code{\link{pdCART}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
pdSpecEst2D <- function(P, order = c(3, 3), metric = "Riemannian", alpha = 1, return = "f", ...) {

  ## Set variables
  dots = list(...)
  tree = (if(is.null(dots$tree)) T else dots$tree)
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  J1 = log2(dim(P)[3])
  J2 = log2(dim(P)[4])
  J = max(J1, J2)
  d = dim(P)[1]
  B = (if(is.null(dots$B)) d else dots$B)
  progress = (if(is.null(dots$progress)) T else dots$progress)
  J.out = (if(is.null(dots$J.out)) J else dots$J.out)
  jmax = min((if(is.null(dots$jmax)) J - 2 else dots$jmax), J.out - 1)
  bias.corr = (if(is.null(dots$bias.corr)) T else dots$bias.corr)
  return.D = (if(is.null(dots$return.D)) NA else dots$return.D)

  # Manifold bias-correction
  P <- (if((metric == "Riemannian" | metric == "logEuclidean") & bias.corr) {
    B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * P } else P)

  ## Threshold full data using 'alpha'
  coeff <- WavTransf2D(P, order = order, jmax = jmax, metric = metric, progress = progress)
  coeff.opt <- pdCART(coeff$D, coeff$D.white, alpha = alpha, tree = tree, order = order, B = B,
                      return.D = return.D)

  ## Return 'f' or not
  f <- (if(return == "f"){
    InvWavTransf2D(coeff.opt$D_w, coeff$M0, order = order, jmax = J.out, metric = metric,
                   progress = progress, chol.bias = T)
  } else NULL)

  ## Return whitened coeff's or not
  if(!isTRUE(return.D == "D.white")){
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, D.raw = coeff$D)
  } else{
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, D.raw = coeff$D,
                D.white = coeff.opt$D.white_w)
  }

  return(res)
}

#' Tree-structured trace thresholding of wavelet coefficients
#'
#' \code{pdCart()} performs hard tree-structured thresholding of the wavelet coefficients obtained with \code{\link{WavTransf1D}}
#' or \code{\link{WavTransf2D}} based on the trace of the whitened wavelet coefficients, see e.g. (Chau and von Sachs, 2017).
#'
#' Depending on the structure of the input list of arrays \code{D} the function performs 1D or 2D tree-structured thresholding of wavelet coefficients.
#' The optimal tree of wavelet coefficients is found by minimization of the \emph{complexity penalized residual sum of squares} (CPRESS) criterion
#' in (Donoho, 1997), via a fast tree-pruning algorithm. By default, the sparsity parameter is set equal to \code{alpha} times
#' the universal threshold \eqn{\sigma_w\sqrt(2\log(n))}, where \eqn{\sigma_w^2} is the noise variance of the traces of the whitened wavelet
#' coefficients determined from the finest wavelet scale and \eqn{n} is the total number of coefficients. By default, \code{alpha = 1},
#' with \code{alpha = 0}, the sparsity parameter is zero and we do not threshold any coefficients.
#'
#' @note For thresholding of 1D wavelet coefficients, the noise
#' variance of the traces of the whitened wavelet coefficients is constant across scales as shown in (Chau and von Sachs, 2017a). For thresholding of 2D
#' wavelet coefficients, there is a discrepancy between the constant noise variance of the traces of the whitened wavelet coefficients of the first
#' \code{abs(J1 - J2)} scales and the remaining scales, where \eqn{J_1 = \log_2(n_1)} and \eqn{J_2 = \log_2(n_2)} with \eqn{n_1} and \eqn{n_2}
#' the dyadic number of observations in each marginal direction of the 2D rectangular tensor grid.  The reason is that the variances of the traces of
#' the whitened coefficients are not homogeneous between: (i) scales at which the 1D wavelet refinement scheme is applied and (ii) scales at which the
#' 2D wavelet refinement scheme is applied. To correct for this discrepancy, the variances of the coefficients at the 2D wavelet scales are normalized
#' by the noise variance determined from the finest wavelet scale. The variances of the coefficients at the 1D wavelet scales are normalized using the
#' analytic noise variance of the traces of the whitened coefficients for a grid of complex random Wishart matrices, which corresponds to the distributional
#' behavior of the pseudo HPD periodogram matrices, or the asymptotic distributional behavior of the actual HPD periodogram matrices. Note that if the
#' 2D time-frequency grid of is a square grid, i.e. \eqn{n_1 = n_2}, the variances of the traces of the whitened coefficients are again homogeneous across
#' all wavelet scales.
#'
#' @param D a list of wavelet coefficients as obtained from \code{\link{WavTransf1D}} or \code{\link{WavTransf2D}}.
#' @param D.white a list of whitened wavelet coefficients as obtained from \code{\link{WavTransf1D}} or \code{\link{WavTransf2D}}.
#' @param alpha tuning parameter specifying the sparsity parameter as \code{alpha} times the universal threshold.
#' @param tree logical value, if \code{tree = T} performs tree-structured thresholding, otherwise performs
#'  non-tree-structured hard thresholding of the coefficients.
#' @param order the order(s) of the intrinsic 1D or 2D AI refinement scheme as in \code{\link{WavTransf1D}} and \code{\link{WavTransf2D}}.
#' @param ... additional parameters for internal use.
#'
#' @return Returns a list with two components:
#'    \item{\code{w} }{ a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale.}
#'    \item{\code{D_w} }{ the list of thresholded wavelet coefficients, with each list component corresponding
#'    to an individual wavelet scale.}
#'
#' @examples
#' ## 1D
#' X <- rExamples(256, example = "bumps")
#' Coeffs <- WavTransf1D(X$per)
#' pdCART(Coeffs$D, Coeffs$D.white, order = 5)$w ## logical tree of non-zero coefficients
#'
#' \dontrun{
#' ## 2D
#' P <- rExamples2D(c(2^7, 2^7), 3, example = "tvar")$per
#' Coeffs <- WavTransf2D(P, jmax = 5)
#' pdCART(Coeffs$D, Coeffs$D.white, order = c(3, 3))$w
#' }
#'
#' @seealso \code{\link{WavTransf1D}}, \code{\link{InvWavTransf1D}}, \code{\link{WavTransf2D}}, \code{\link{InvWavTransf2D}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Donoho, D.L. (1997). \emph{CART and best-ortho-basis: a connection}. Annals of Statistics,
#' 25(5), 1870-1911.
#'
#' @export
pdCART <- function(D, D.white, order, alpha = 1, tree = T, ...) {

  ## Set variables
  J = length(D)
  d = dim(D[[1]])[1]
  is_2D = ifelse(length(dim(D[[1]])) == 4, T, F)
  dots = list(...)
  B = (if(is.null(dots$B)) d else dots$B)
  periodic = (if(is.null(dots$periodic) | is_2D) F else dots$periodic)
  return.D = (if(is.null(dots$return.D)) NA else dots$return.D)
  if(periodic){
    L = (order - 1) / 2
    L_b = ceiling(L / 2)
  } else{
    L_b = 0
  }
  lam = (if(is.null(dots$lam)) NA else dots$lam)

  if(is_2D){
    D_trace <- lapply(2:J, function(j) apply(D.white[[j]], c(3, 4), function(A) Re(sum(diag(A)))))
    J_tr <- length(D_trace)
    s_e2D <- stats::mad(c(D_trace[[J_tr]]))

    J0_2D <- sum(sapply(1:J, function(j) any(dim(D[[j]]) == 1)))
    if(J0_2D > 1){
      n_2D <- which.max(dim(D[[J0_2D]])[c(3,4)])
      if((max(order > 9))){
        warning(paste0('The first ', J0_2D[1], ' are not thresholded, since the maximum marginal
                     refinement order is larger than 9. To include all wavelet scales in the
                     thresholding procedure choose marginal refinement order smaller or equal to 9.'))
      }
      for(j in 1:length(D_trace)){
        if(j < J0_2D){
          n <- dim(D_trace[[j]])[n_2D]
          N <- ifelse(order[n_2D] > n, 2 * floor((n - 1) / 2) + 1, order[n_2D])
          L <- (N - 1) / 2
          l <- sapply(1:n, function(k) if((k - L) < 1) 2 * k else if((k + L) >= n) 2 * (N - (n - k)) else N + 1)
          s_e1D <- sqrt(sapply(1:n, function(k) 4^(-J) * sum(W_1D[[L + 1]][l[k], ]^2) * sum(trigamma(B - (d - 1:d)))))
          D_trace[[j]] <- D_trace[[j]] / s_e1D
        } else{
          D_trace[[j]] <- D_trace[[j]] / s_e2D
        }
      }
    } else {
      D_trace <- lapply(1:J_tr, function(j) D_trace[[j]] / s_e2D)
    }
  } else {
    D_trace_full <- lapply(2:J, function(j) apply(D.white[[j]], 3, function(A) Re(sum(diag(A)))))
    J_tr <- length(D_trace_full)
    s_e <- max(stats::mad(c(D_trace_full[[J_tr]])),
               sqrt(2^(-J) * sum(W_1D[[(order + 1)/2]][order + 1, ]^2) * sum(trigamma(B - (d - 1:d)))))
    D_trace <- lapply(1:J_tr, function(j) D_trace_full[[j]][L_b + 1:2^j] / s_e)
  }
  if(is.na(lam)){
    lam <- alpha * sqrt(2 * log(length(unlist(D_trace))))
  }

  if(tree){

    ## Dyadic CART
    w <- D_trace

    for(j in J_tr:1){
      if(j == J_tr){
        w[[j]] <- ifelse(abs(D_trace[[j]]) > lam, T, F)
        R <- pmin(D_trace[[j]]^2, lam^2)
        V <- D_trace[[j]]^2
      } else{
        if(is_2D){
          dims <- dim(D_trace[[j]])
          if(all(dim(D_trace[[min(j + 1, J_tr)]]) > 1)){
            l1 <- t(sapply(1:dims[1], function(i) c(1, 2) + 2 * (i - 1)))
            l2 <- t(sapply(1:dims[2], function(i) c(1, 2) + 2 * (i - 1)))
            grid <- expand.grid(1:dims[1], 1:dims[2])
            V <- array(c(mapply(function(i1, i2) sum(V[l1[i1, ], l2[i2, ]]), grid$Var1, grid$Var2)),
                     dim = dims) + D_trace[[j]]^2
            R <- array(c(mapply(function(i1, i2) sum(R[l1[i1, ], l2[i2, ]]), grid$Var1, grid$Var2)),
                     dim = dims) + lam^2
          } else if(dim(D_trace[[j + 1]])[1] == 1){
            V <- sapply(1:dims[2], function(i) V[1, 2 * i - 1] + V[1, 2 * i]) + D_trace[[j]]^2
            R <- sapply(1:dims[2], function(i) R[1, 2 * i - 1] + R[1, 2 * i]) + lam^2
          } else if(dim(D_trace[[j + 1]])[2] == 1){
            V <- sapply(1:dims[1], function(i) V[2 * i - 1, 1] + V[2 * i, 1]) + D_trace[[j]]^2
            R <- sapply(1:dims[1], function(i) R[2 * i - 1, 1] + R[2 * i, 1]) + lam^2
          }
        } else {
          dims <- length(D_trace[[j]])
          V <- sapply(1:dims, function(i) V[2 * i - 1] + V[2 * i]) + D_trace[[j]]^2
          R <- sapply(1:dims, function(i) R[2 * i - 1] + R[2 * i]) + lam^2
        }
        w[[j]] <- ifelse(R < V, T, F)
        R <- pmin(V, R)
      }
    }
  } else{
    w <- lapply(1:J_tr, function(j) abs(D_trace[[j]]) > lam)
  }

  ## Threshold wavelet coefficients with weights 'w'
  D_w <- D
  if(isTRUE(return.D == "D.white")){
    D.white_w <- D.white
  }
  if(is_2D){
    ## 2D
    for(j in 2:J){
      if(tree){
        if(j > 2){
          dims <- dim(w[[j - 1]])
          roots <- matrix(rep(matrix(rep(t(w[[j - 2]]), each = ifelse(dims[2] > 1, 2, 1)),
                              byrow = T, ncol = dims[2]), each = ifelse(dims[1] > 1, 2, 1)), nrow = dims[1])
          w[[j - 1]] <- w[[j - 1]] & roots
        }
      }
      D0 <- array(D_w[[j]], dim = c(d, d, dim(D_w[[j]])[3] * dim(D_w[[j]])[4]))
      D0[, , !(w[[j - 1]])] <- 0
      D_w[[j]] <- array(D0, dim = c(d, d, dim(D_w[[j]])[3], dim(D_w[[j]])[4]))

      if(isTRUE(return.D == "D.white")){
        D0 <- array(D.white[[j]], dim = c(d, d, dim(D.white_w[[j]])[3] * dim(D.white_w[[j]])[4]))
        D0[, , !(w[[j - 1]])] <- 0
        D.white_w[[j]] <- array(D0, dim = c(d, d, dim(D.white_w[[j]])[3], dim(D.white_w[[j]])[4]))
      }
    }
  } else {
    ## 1D
    for(j in 2:J){
      w[[j - 1]] <- (if(tree) w[[j - 1]] & rep((if(j == 2) T else w[[j - 2]]), each = 2) else w[[j - 1]])
      if(periodic & (L_b > 0)){
        zeros <- !(c(abs(D_trace_full[[j - 1]][1:L_b]) > lam, w[[j - 1]], abs(D_trace_full[[j - 1]][2^(j - 1) + L_b + 1:L_b]) > lam))
      } else{
        zeros <- !(w[[j - 1]])
      }
      D_w[[j]][, , zeros] <- 0
      if(isTRUE(return.D == "D.white")){
        D.white_w[[j]][, , zeros] <- 0
      }
    }
  }

  res <- (if(!isTRUE(return.D == "D.white")) list(w = w, D_w = D_w) else list(w = w, D_w = D_w, D.white_w = D.white_w))

  return(res)
}


