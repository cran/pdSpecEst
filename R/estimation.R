#' Wavelet-thresholded multivariate spectral estimator.
#'
#' \code{pdSpecEst} calculates a \eqn{(d \times d)}-dimensional Hermitian PD wavelet-denoised multivariate
#' spectral estimator by thresholding wavelet coefficients in the manifold wavelet domain. The
#' estimation procedure is described in detail (Chau and von Sachs, 2017).
#'
#' The input array \code{P} corresponds to an initial noisy Hermitian PD spectral estimate of the
#' (\eqn{d \times d})-dimensional spectral matrix at \code{m} different frequencies, with \eqn{m = 2^J} for some
#' \eqn{J > 0}. This can be e.g. the curve of averaged periodograms given as output in \code{\link{pdPgram}}.\cr
#' \code{P} is transformed to the manifold wavelet domain by the function \code{\link{WavTransf}}, and the
#' wavelet coefficients are decomposed in terms of an orthonormal basis of the space of Hermitian matrices. \cr
#' The variances of the components of the wavelet coefficients are standardized across scales via weighted
#' log-linear regression, where the weights increase exponentially across scales at a rate depending on the
#' tuning parameter \code{alpha}. Without prior knowledge of the sparsity of the signal, a choice \code{alpha}
#' in \eqn{[0.5, 1)} is reasonable in most settings. If \code{alpha} is not specified, by default \code{alpha = 0.75} \cr
#' The components of the wavelet coefficients are thresholded based on the hard (keep-or-kill)
#' threshold \code{lam}. If \code{lam} is unspecified, the threshold is determined in a data-adaptive manner
#' by a twofold cross-validation procedure, which is described in detail in (Chau and von Sachs, 2017). \cr
#' If \code{return == 'f'} the thresholded wavelet coefficients are transformed back to the frequency domain by
#' the function \code{\link{InvWavTransf}} giving the wavelet-denoised Hermitian PD spectral estimate.
#'
#' @param P a (\eqn{d,d,m})-dimensional array of Hermitian PD matrices, with \eqn{m} a dyadic number.
#' @param lam an optional argument specifying the wavelet threshold, if \code{lam}
#'  is not specified the threshold is calculated by a twofold cross-validation procedure.
#' @param order an odd integer between 1 and 9 corresponding to the refinement order of the MI wavelet transform.
#' @param return an optional argument that specifies whether the denoised spectral estimator
#'  is returned or not.
#' @param alpha an optional argument tuning the weights used to normalize the variances
#' of the wavelet coefficients across scales. By default, \code{alpha} is set to \code{0.75}.
#'
#' @return The function returns a list with four components:
#' \item{f }{a (\eqn{d,d,m})-dimensional array corresponding to the wavelet-denoised Hermitian PD (\eqn{d \times d})-dimensional
#' spectral estimate at the \code{m} different frequencies. If \code{!(return == 'f')}, the inverse wavelet transform
#' of the thresholded wavelet coefficients is not computed and \code{f} is set equal to \code{NULL}.}
#' \item{D }{a list of arrays, each (\eqn{d, d, 2^j})-dimensional array contains the thresholded
#' (\eqn{d \times d})-dimensional wavelet coefficents at the \eqn{2^j} different locations in the given wavelet scale
#' \eqn{j}. The first list element contains the midpoints at the coarsest scale in the
#' midpoint pyramid \eqn{j=1}, see (Chau and von Sachs, 2017) for more details.}
#' \item{lam }{the hard threshold used to threshold the components of the wavelet coefficients.}
#' \item{components }{a list of arrays, each (\eqn{d^2, 2^j})-dimensional array contains the components
#' of the thresholded (\eqn{d \times d})-dimensional wavelet coefficients in terms of an orthonormal basis of the
#' space of (\eqn{d \times d})-dimensional Hermitian matrices. The columns correspond to the \eqn{d^2} basis components
#' at each of the \eqn{2^j} different locations at wavelet scale \eqn{j}.}
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(2^10, 2, Phi, Theta, Sigma)
#' ts.plot(ts.sim$X) # plot generated time series traces.
#'
#' pgram <- pdPgram(ts.sim$X)
#' f <- pdSpecEst(pgram$P)
#'
#' @seealso \code{\link{pdPgram}}, \code{\link{WavTransf}}, \code{\link{InvWavTransf}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @importFrom stats mad
#' @export
pdSpecEst <- function(P, lam = NULL, order = 5, return = "f", alpha = 0.75) {

  ## Define variables
  J <- log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    warning(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                 " to dyadic number."))
  }
  stopifnot(isTRUE(all.equal(as.integer(J), J)))
  if (!(order %in% c(1, 3, 5, 7, 9))) {
    warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  dim <- dim(P)[1]
  E <- E_basis(dim)

  ## Find optimal threshold
  if (is.null(lam)) {

    ## Wavelet transforms
    P.half <- list(odd = P[, , c(T, F)], even = P[, , c(F, T)])
    D.half <- list(odd = WavTransf(P.half$odd, order)$D, even = WavTransf(P.half$even, order)$D)
    d <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d$odd[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$odd[[j + 1]][, , k], E))
      d$even[[j]] <- sapply(1:2^j, function(k) E_coeff(D.half$even[[j + 1]][, , k], E))
    }

    ## Normalize variances wavelet coefficients
    W <- diag(exp(alpha * (1:(J - 2))))
    X <- matrix(c(rep(1, J - 2), 1:(J - 2)), nrow = J - 2)
    mads <- lapply(1:2, function(l) sapply(1:(J - 2), function(j) apply(d[[l]][[j]], 1, stats::mad)))
    Y <- colMeans(rbind(log(mads[[1]]), log(mads[[2]])))
    beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
    sd.j <- function(j) exp(sum(beta * c(1, j)))

    d.vec <- list(odd = NULL, even = NULL)
    d.new <- list(odd = list(), even = list())
    for (j in 1:(J - 2)) {
      d.new$odd[[j]] <- d$odd[[j]]/sd.j(j)
      d.vec$odd <- cbind(d.vec$odd, d.new$odd[[j]])
      d.new$even[[j]] <- d$even[[j]]/sd.j(j)
      d.vec$even <- cbind(d.vec$even, d.new$even[[j]])
    }

    ## Cross-validation scores
    cv <- function(lam) {
      D.lam <- D.half
      d.lam <- d
      for (j in 3:(J - 2)) {
        zero.odd <- sapply(1:2^j, function(k) (abs(d.new$odd[[j]][, k]) <
                                                 lam) | (d.lam$odd[[j - 1]][, ceiling(k/2)] == 0))
        zero.even <- sapply(1:2^j, function(k) (abs(d.new$even[[j]][, k]) <
                                                  lam) | (d.lam$even[[j - 1]][, ceiling(k/2)] == 0))
        d.lam$odd[[j]][zero.odd] <- 0
        d.lam$even[[j]][zero.even] <- 0
      }
      for (j in 1:(J - 2)) {
        D.lam$odd[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$odd[[j]][, k], E),
                                      simplify = "array")
        D.lam$even[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d.lam$even[[j]][, k], E),
                                       simplify = "array")
      }
      f.hat <- list(odd = InvWavTransf(D.lam$odd, order), even = InvWavTransf(D.lam$even, order))

      ## Predicted points
      f.t.even <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$odd[, , k], f.hat$odd[, , k + 1]),
                          simplify = "array")
      f.t.even <- array(c(f.t.even, f.hat$odd[, , 2^(J - 1)]), dim = c(dim, dim, 2^(J - 1)))
      f.t.odd <- sapply(1:(2^(J - 1) - 1), function(k) Mid(f.hat$even[, , k], f.hat$even[, , k + 1]),
                         simplify = "array")
      f.t.odd <- array(c(f.hat$even[, , 1], f.t.odd), dim = c(dim, dim, 2^(J - 1)))

      return(sum(sapply(1:2^(J - 1), function(k) RiemmDist(f.t.even[, , k], P.half$even[, , k])^2 +
                                                  RiemmDist(f.t.odd[, , k], P.half$odd[, , k])^2)))
    }

    ## Golden section search
    lam.range <- sort(abs(c(d.vec$even)), decreasing = T)[c(10, round(0.25 * length(c(d.vec$even))))]
    lam.cv <- gss(lam.range, cv)
  }

  ## Rescale threshold to twice number of observations
  lam.cv <- ifelse(is.null(lam), 1/sqrt((1 - log(2)/log(2^J * dim^2))) * lam.cv, lam)

  ## Transform original noisy data
  D <- WavTransf(P, order)$D
  d <- list()
  for (j in 1:(J - 1)) {
    d[[j]] <- sapply(1:2^j, function(k) E_coeff(D[[j + 1]][, , k], E))
  }

  ## Normalize variances wavelet coefficients
  W <- diag(exp(alpha * (1:(J - 1))))
  X <- matrix(c(rep(1, J - 1), 1:(J - 1)), nrow = J - 1)
  Y <- colMeans(log(sapply(1:(J - 1), function(j) apply(d[[j]], 1, stats::mad))))
  beta <- ((solve((t(X) %*% W) %*% X) %*% t(X)) %*% W) %*% Y
  sd.j <- function(j) exp(sum(beta * c(1, j)))

  d.vec <- NULL
  d.new <- list()
  for (j in 1:(J - 1)) {
    d.new[[j]] <- d[[j]]/sd.j(j)
    d.vec <- cbind(d.vec, d.new[[j]])
  }

  ## Threshold coefficients
  for (j in 3:(J - 1)) {
    zero <- sapply(1:2^j, function(k) (abs(d.new[[j]][, k]) < lam.cv) |
                                       (d[[j - 1]][, ceiling(k/2)] == 0))
    d[[j]][zero] <- 0
  }

  ## Inverse transform denoised data
  for (j in 1:(J - 1)) {
    D[[j + 1]] <- sapply(1:2^j, function(k) E_coeff_inv(d[[j]][, k], E), simplify = "array")
  }

  if (return == "f") {
    f <- InvWavTransf(D, order)
  } else {
    f <- NULL
  }
  return(list(f = f, D = D, lam = lam.cv, components = d))
}
