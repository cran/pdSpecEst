#' Intrinsic wavelet HPD spectral estimation
#'
#' \code{pdSpecEst1D} calculates a \eqn{(d,d)}-dimensional HPD wavelet-denoised spectral matrix estimator
#' by applying the following steps to an initial noisy HPD spectral estimate (obtained with e.g., \code{\link{pdPgram}}):
#' \enumerate{
#'     \item a forward intrinsic AI wavelet transform, with \code{\link{WavTransf1D}},
#'     \item (tree-structured) thresholding of the wavelet coefficients, with \code{\link{pdCART}},
#'     \item an inverse intrinsic AI wavelet transform, with \code{\link{InvWavTransf1D}}.
#' }
#' The complete estimation procedure is described in more detail in \insertCite{CvS17}{pdSpecEst} or Chapter 3 of
#' \insertCite{C18}{pdSpecEst}.
#'
#' The input array \code{P} corresponds to an initial noisy HPD spectral estimate of the (\eqn{d,d})-dimensional
#' spectral matrix at \code{m} different frequencies, with \eqn{m = 2^J} for some \eqn{J > 0}. This can be e.g.,
#' a multitaper HPD periodogram given as output by the function \code{\link{pdPgram}}.\cr
#' \code{P} is transformed to the wavelet domain by the function \code{\link{WavTransf1D}}, which applies an intrinsic
#' 1D AI wavelet transform based on a metric specified by the user. The noise is removed by tree-structured
#' thresholding of the wavelet coefficients based on the trace of the whitened coefficients with \code{\link{pdCART}} by
#' minimization of a \emph{complexity penalized residual sum of squares} (CPRESS) criterion via the fast tree-pruning algorithm
#' in \insertCite{D97}{pdSpecEst}. The penalty or sparsity parameter in the optimization procedure is set equal to \code{alpha}
#' times the universal threshold, where the noise variance of the traces of the whitened wavelet
#' coefficients are determined from the finest wavelet scale. See \insertCite{CvS17}{pdSpecEst} and Chapter 3 of \insertCite{C18}{pdSpecEst}
#' for further details. \cr
#' The function computes the forward and inverse intrinsic AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics. \cr
#' If \code{return_val = 'f'} the thresholded wavelet coefficients are transformed back to the frequency domain by
#' the inverse intrinsic 1D AI wavelet transform via \code{\link{InvWavTransf1D}}, returning the wavelet-denoised
#' HPD spectral estimate.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of HPD matrices, corresponding to a sequence of \eqn{(d,d)}-dimensional HPD matrices
#' of length \eqn{m}, with \eqn{m = 2^J} for some \eqn{J > 0}.
#' @param order an odd integer larger or equal to 1 corresponding to the order of the intrinsic AI refinement scheme,
#' defaults to \code{order = 5}. Note that if \code{order > 9}, the computational cost
#' significantly increases as the wavelet transform no longer uses a fast wavelet refinement scheme based
#' on pre-determined weights.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"},
#' \code{"Euclidean"} or \code{"Riemannian-Rahman"}. See also the Details section below.
#' @param alpha an optional tuning parameter in the wavelet thresholding procedure. The penalty (or sparsity)
#' parameter in the tree-structured wavelet thresholding procedure in \code{\link{pdCART}} is set to \code{alpha}
#' times the estimated universal threshold, defaults to \code{alpha = 1}.
#' @param return_val an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not. See the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with the following five components:
#' \item{f }{ a (\eqn{d,d,m})-dimensional array of HPD matrices, corresponding to the HPD wavelet-denoised estimate
#' of the same resolution as the input array \code{P}. If \code{return_val != 'f'}, the inverse wavelet transform
#' of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{ the pyramid of threshold wavelet coefficients. This is a list of arrays, where each array contains the
#' (\eqn{d,d})-dimensional thresholded wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to the finest
#' wavelet scale \code{j = jmax}.}
#' \item{M0 }{ a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the midpoint pyramid.}
#' \item{tree.weights }{a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale starting from the coarsest wavelet scale \code{j = 0}.}
#' \item{D.raw }{ the pyramid of non-thresholded wavelet coefficients in the same format as the component \code{$D}.}
#'
#' @examples
#' P <- rExamples1D(2^8, example = "bumps")$P
#' f <- pdSpecEst1D(P)
#'
#' @seealso \code{\link{pdPgram}}, \code{\link{WavTransf1D}}, \code{\link{InvWavTransf1D}}, \code{\link{pdCART}}
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and (depending on the
#' specified metric) may fail if matrices are close to being singular.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdSpecEst1D <- function(P, order = 5, metric = "Riemannian", alpha = 1, return_val = "f", ...) {

  ## Set variables
  dots <- list(...)
  tree <- (if(is.null(dots$tree)) TRUE else dots$tree)
  w.tree <- (if(is.null(dots$w.tree)) NULL else dots$w.tree)
  periodic <- (if(is.null(dots$periodic)) TRUE else dots$periodic)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean", "Riemannian-Rahman"))
  J <- log2(dim(P)[3])
  d <- dim(P)[1]
  B <- (if(is.null(dots$B)) d else dots$B)
  J.out <- (if(is.null(dots$J.out)) J else dots$J.out)
  jmax <- min((if(is.null(dots$jmax)) J - 2 else dots$jmax), J.out - 1)
  bias.corr <- (if(is.null(dots$bias.corr)) TRUE else dots$bias.corr)
  return.D <- (if(is.null(dots$return.D)) NA else dots$return.D)

  ## Wishart bias-correction
  P <- (if((grepl("Riemannian", metric) | metric == "logEuclidean") & bias.corr) {
    B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * P } else P)

  ## (1) Transform data to wavelet domain
  coeff <- WavTransf1D(P, order, jmax = jmax, periodic = periodic, metric = metric, method = method)

  ## (2) Threshold coefficients in wavelet domain
  coeff.thresh <- pdCART(coeff$D, coeff$D.white, alpha = alpha, tree = tree, periodic = periodic,
                      w.tree = w.tree, order = order, B = B, return.D = return.D)

  ## (3) Transform back to HPD space
  f <- (if(return_val == "f"){
    InvWavTransf1D(coeff.thresh$D_w, coeff$M0, order = order, jmax = J.out, periodic = periodic,
                   metric = metric, method = method, chol_bias = bias.corr)
  } else NULL)

  ## Return whitened coeff's or not
  if(!isTRUE(return.D == "D.white")){
    res <- list(f = f, D = coeff.thresh$D_w, M0 = coeff$M0, tree.weights = coeff.thresh$w, D.raw = coeff$D)
  } else{
    res <- list(f = f, D = coeff.thresh$D_w, M0 = coeff$M0, tree.weights = coeff.thresh$w, D.raw = coeff$D,
                D.white = coeff.thresh$D.white_w)
  }
  return(res)
}

#' Intrinsic wavelet HPD time-varying spectral estimation
#'
#' \code{pdSpecEst2D} calculates a \eqn{(d,d)}-dimensional HPD wavelet-denoised time-varying spectral matrix estimator
#' by applying the following steps to an initial noisy HPD time-varying spectral estimate (obtained with e.g., \code{\link{pdPgram2D}}):
#' \enumerate{
#'     \item a forward intrinsic AI wavelet transform, with \code{\link{WavTransf2D}},
#'     \item (tree-structured) thresholding of the wavelet coefficients, with \code{\link{pdCART}},
#'     \item an inverse intrinsic AI wavelet transform, with \code{\link{InvWavTransf2D}}.
#' }
#' The complete estimation procedure is described in more detail in Chapter 5 of \insertCite{C18}{pdSpecEst}.
#'
#' The input array \code{P} corresponds to an initial noisy HPD time-varying spectral estimate of the (\eqn{d, d})-dimensional
#' spectral matrix at a time-frequency grid of size \eqn{m_1 \times m_2}, with \eqn{m_1, m_2} dyadic numbers. This can be e.g.,
#' a multitaper HPD time-varying periodogram given as output by the function \code{\link{pdPgram2D}}.\cr
#' \code{P} is transformed to the wavelet domain by the function \code{\link{WavTransf2D}}, which applies an intrinsic
#' 2D AI wavelet transform based on a metric specified by the user. The noise is removed by tree-structured
#' thresholding of the wavelet coefficients based on the trace of the whitened coefficients with \code{\link{pdCART}} by
#' minimization of a \emph{complexity penalized residual sum of squares} (CPRESS) criterion via the fast tree-pruning algorithm
#' in \insertCite{D97}{pdSpecEst}. The penalty (i.e., sparsity) parameter in the optimization procedure is set equal to \code{alpha}
#' times the universal threshold, where the noise variance of the traces of the whitened wavelet
#' coefficients are determined from the finest wavelet scale. See Chapter 5 of \insertCite{C18}{pdSpecEst}
#' for further details. \cr
#' The function computes the forward and inverse intrinsic 2D AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics. \cr
#' If \code{return_val = 'f'} the thresholded wavelet coefficients are transformed back to the time-frequency domain by
#' the inverse intrinsic 2D AI wavelet transform via \code{\link{InvWavTransf2D}}, returning the wavelet-denoised
#' HPD time-varying spectral estimate.
#'
#' @param P a (\eqn{d,d,n1,n2})-dimensional array of HPD matrices corresponding to a rectangular surface of \eqn{(d,d)}-dimensional HPD matrices
#' of size \eqn{n_1 \times n_2}, with \eqn{n_1 = 2^{J_1}} and \eqn{n_2 = 2^{J_2}} for some \eqn{J_1, J_2 > 0}.
#' @param order a 2-dimensional numeric vector \eqn{(1,1) \le} \code{order} \eqn{\le (9,9)} corresponding to the marginal
#' orders of the intrinsic 2D AI refinement scheme, defaults to \code{order = c(3, 3)}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. See also the Details section below.
#' @param alpha an optional tuning parameter in the wavelet thresholding procedure. The penalty (or sparsity)
#' parameter in the tree-structured wavelet thresholding procedure in \code{\link{pdCART}} is set to \code{alpha}
#' times the estimated universal threshold, defaults to \code{alpha = 1}.
#' @param return_val an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not. See the Details section below.
#' @param ... additional arguments for internal use.
#'
#' @return The function returns a list with the following five components:
#' \item{f }{ a (\eqn{d,d,m1,m2})-dimensional array of HPD matrices, corresponding to the HPD wavelet-denoised estimate
#' on the same resolution grid of size \eqn{m_1 \times m_2} as specified by the input array \code{P}. If \code{return_val != 'f'}, the
#' inverse wavelet transform of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{ the 2D pyramid of threshold wavelet coefficients. This is a list of arrays, where each array contains the rectangular grid
#' (\eqn{d,d})-dimensional thresholded wavelet coefficients from the coarsest wavelet scale \code{j = 0} up to the finest
#' wavelet scale \code{j = jmax}.}
#' \item{M0 }{ a numeric array containing the midpoint(s) at the coarsest scale \code{j = 0} in the 2D midpoint pyramid.}
#' \item{tree.weights }{ a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale starting from the coarsest wavelet scale \code{j = 0}.}
#' \item{D.raw }{ the 2D pyramid of non-thresholded wavelet coefficients in the same format as the component \code{$D}.}
#'
#' @examples
#' \dontrun{
#' P <- rExamples2D(c(2^6, 2^6), 2, example = "tvar")$P
#' f <- pdSpecEst2D(P)
#' }
#'
#' @seealso \code{\link{pdPgram2D}}, \code{\link{WavTransf2D}}, \code{\link{InvWavTransf2D}}, \code{\link{pdCART}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdSpecEst2D <- function(P, order = c(3, 3), metric = "Riemannian", alpha = 1, return_val = "f", ...) {

  ## Set variables
  dots <- list(...)
  tree <- (if(is.null(dots$tree)) TRUE else dots$tree)
  w.tree <- (if(is.null(dots$w.tree)) NULL else dots$w.tree)
  method <- (if(is.null(dots$method)) "fast" else dots$method)
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  J1 <- log2(dim(P)[3])
  J2 <- log2(dim(P)[4])
  J <- max(J1, J2)
  d <- dim(P)[1]
  B <- (if(is.null(dots$B)) d else dots$B)
  J.out <- (if(is.null(dots$J.out)) J else dots$J.out)
  jmax <- min((if(is.null(dots$jmax)) J - 2 else dots$jmax), J.out - 1)
  bias.corr <- (if(is.null(dots$bias.corr)) TRUE else dots$bias.corr)
  return.D <- (if(is.null(dots$return.D)) NA else dots$return.D)

  # Wishart bias-correction
  P <- (if((metric == "Riemannian" | metric == "logEuclidean") & bias.corr) {
    B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * P } else P)

  ## (1) Transform data to wavelet domain
  coeff <- WavTransf2D(P, order = order, jmax = jmax, metric = metric, method = method)

  ## (2) Threshold coefficients in wavelet domain
  coeff.opt <- pdCART(coeff$D, coeff$D.white, alpha = alpha, tree = tree, w.tree = w.tree,
                      order = order, B = B, return.D = return.D)

  ## (3) Transform back to HPD space
  f <- (if(return_val == "f"){
    InvWavTransf2D(coeff.opt$D_w, coeff$M0, order = order, jmax = J.out, metric = metric,
                   method = method, return_val = return_val, chol_bias = bias.corr)
  } else NULL)

  ## Return whitened coeff's or not
  if(!isTRUE(return.D == "D.white")){
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, D.raw = coeff$D)
  } else{
    res <- list(f = f, D = coeff.opt$D_w, M0 = coeff$M0, tree.weights = coeff.opt$w, D.raw = coeff$D,
                D.white = coeff.opt$D.white_w, D.raw_white = coeff$D.white)
  }
  return(res)
}

#' Tree-structured trace thresholding of wavelet coefficients
#'
#' \code{pdCART} performs hard tree-structured thresholding of the Hermitian matrix-valued wavelet coefficients obtained with
#' \code{\link{WavTransf1D}} or \code{\link{WavTransf2D}} based on the trace of the whitened wavelet coefficients, as explained in
#' \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst}. This function is primarily written for internal use in other functions and
#' is typically not used as a stand-alone function.
#'
#' Depending on the structure of the input list of arrays \code{D} the function performs 1D or 2D tree-structured thresholding of wavelet coefficients.
#' The optimal tree of wavelet coefficients is found by minimization of the \emph{complexity penalized residual sum of squares} (CPRESS) criterion
#' in \insertCite{D97}{pdSpecEst}, via a fast tree-pruning algorithm. By default, the penalty parameter in the optimization procedure is set equal to
#' \code{alpha} times the universal threshold \eqn{\sigma_w\sqrt(2\log(n))}, where \eqn{\sigma_w^2} is the noise variance of the traces of the whitened
#' wavelet coefficients determined from the finest wavelet scale and \eqn{n} is the total number of coefficients. By default, \code{alpha = 1},
#' if \code{alpha = 0}, the penalty parameter is zero and the coefficients remain untouched.
#'
#' @note For thresholding of 1D wavelet coefficients, the noise
#' variance of the traces of the whitened wavelet coefficients is constant across scales as seen in \insertCite{CvS17}{pdSpecEst}. For thresholding of 2D
#' wavelet coefficients, there is a discrepancy between the constant noise variance of the traces of the whitened wavelet coefficients at the first
#' \code{abs(J1 - J2)} scales and the remaining scales, as discussed in Chapter 5 of \insertCite{C18}{pdSpecEst}, where \eqn{J_1 = \log_2(n_1)} and
#' \eqn{J_2 = \log_2(n_2)} with \eqn{n_1} and \eqn{n_2} the dyadic number of observations in each marginal direction of the 2D rectangular tensor grid.
#' The reason is that the variances of the traces of the whitened coefficients are not homogeneous between: (i) scales at which the 1D wavelet refinement
#' scheme is applied and (ii) scales at which the 2D wavelet refinement scheme is applied. To correct for this discrepancy, the variances of the coefficients
#' at the 2D wavelet scales are normalized by the noise variance determined from the finest wavelet scale. The variances of the coefficients at the 1D wavelet
#' scales are normalized using the analytic noise variance of the traces of the whitened coefficients for a grid of complex random Wishart matrices, which
#' corresponds to the asymptotic distributional behavior of the HPD periodogram matrices obtained with e.g., \code{\link{pdPgram2D}}. Note that if the
#' time-frequency grid is square, i.e., \eqn{n_1 = n_2}, the variances of the traces of the whitened coefficients are again homogeneous across all wavelet scales.
#'
#' @param D a list of wavelet coefficients as obtained from the \code{$D} component of \code{\link{WavTransf1D}} or \code{\link{WavTransf2D}} .
#' @param D.white a list of whitened wavelet coefficients as obtained from the \code{$D.white} component of \code{\link{WavTransf1D}} or \code{\link{WavTransf2D}}.
#' @param alpha tuning parameter specifying the penalty/sparsity parameter as \code{alpha} times the universal threshold.
#' @param tree logical value, if \code{tree = TRUE} performs tree-structured thresholding, otherwise performs
#'  non-tree-structured hard thresholding of the coefficients.
#' @param order the order(s) of the intrinsic 1D or 2D AI refinement scheme as in \code{\link{WavTransf1D}} and \code{\link{WavTransf2D}}.
#' @param ... additional arguments for internal use.
#'
#' @return Returns a list with two components:
#'    \item{\code{w} }{ a list of logical values specifying which coefficients to keep, with each list component
#'    corresponding to an individual wavelet scale starting from the coarsest wavelet scale \code{j = 0}.}
#'    \item{\code{D_w} }{ the list of thresholded wavelet coefficients, with each list component corresponding
#'    to an individual wavelet scale.}
#'
#' @examples
#' ## 1D tree-structured trace thresholding
#' P <- rExamples1D(2^8, example = "bumps")$P
#' Coeffs <- WavTransf1D(P)
#' pdCART(Coeffs$D, Coeffs$D.white, order = 5)$w ## logical tree of non-zero coefficients
#'
#' \dontrun{
#' ## 2D tree-structured trace thresholding
#' P <- rExamples2D(c(2^6, 2^6), 2, example = "tvar")$P
#' Coeffs <- WavTransf2D(P)
#' pdCART(Coeffs$D, Coeffs$D.white, order = c(3, 3))$w
#' }
#'
#' @seealso \code{\link{WavTransf1D}}, \code{\link{InvWavTransf1D}}, \code{\link{WavTransf2D}}, \code{\link{InvWavTransf2D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdCART <- function(D, D.white, order, alpha = 1, tree = TRUE, ...) {

  ## Set variables
  J <- length(D)
  d <- dim(D[[1]])[1]
  is_2D <- ifelse(length(dim(D[[1]])) == 4, TRUE, FALSE)
  dots <- list(...)
  B <- (if(is.null(dots$B)) d else dots$B)
  periodic <- (if(is.null(dots$periodic) | is_2D) FALSE else dots$periodic)
  return.D <- (if(is.null(dots$return.D)) NA else dots$return.D)
  w.tree <- (if(is.null(dots$w.tree)) NULL else dots$w.tree)
  if(periodic){
    L <- (order - 1) / 2
    L_b <- ceiling(L / 2)
  } else{
    L <- L_b <- 0
  }
  lam <- (if(is.null(dots$lam)) NA else dots$lam)

  ## Compute traces and standardize variance across wavelet scales
  if(is_2D){
    D_trace <- lapply(2:J, function(j) apply(D.white[[j]], c(3, 4), function(A) Re(sum(diag(A)))))
    J_tr <- length(D_trace)
    s_e2D <- stats::mad(c(D_trace[[J_tr]]))

    J0_2D <- sum(sapply(1:J, function(j) any(dim(D[[j]]) == 1)))
    if(J0_2D > 1){
      n_2D <- which.max(dim(D[[J0_2D]])[c(3,4)])
      if((max(order > 9))){
        stop("The marginal refinement orders in 'order' should all be smaller or equal to 9")
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

  if(tree && is.null(w.tree)){

    ## Dyadic CART
    w <- D_trace

    for(j in J_tr:1){
      if(j == J_tr){
        w[[j]] <- ifelse(abs(D_trace[[j]]) > lam, TRUE, FALSE)
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
        w[[j]] <- ifelse(R < V, TRUE, FALSE)
        R <- pmin(V, R)
      }
    }
  } else if(!isTRUE(tree) && is.null(w.tree))  {
    w <- lapply(1:J_tr, function(j) abs(D_trace[[j]]) > lam)
  } else if(!is.null(w.tree)){
    w <- w.tree
  }

  ## Threshold wavelet coefficients with weights 'w'
  D_w <- D
  if(isTRUE(return.D == "D.white")){
    D.white_w <- D.white
  }
  if(is_2D){
    ## 2D
    for(j in 2:J){
      if(isTRUE(tree)){
        if(j > 2){
          dims <- dim(w[[j - 1]])
          roots <- matrix(rep(matrix(rep(t(w[[j - 2]]), each = ifelse(dims[2] > 1, 2, 1)),
                              byrow = TRUE, ncol = dims[2]), each = ifelse(dims[1] > 1, 2, 1)), nrow = dims[1])
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
      w[[j - 1]] <- (if(isTRUE(tree)) w[[j - 1]] & rep((if(j == 2) TRUE else w[[j - 2]]), each = 2) else w[[j - 1]])
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

