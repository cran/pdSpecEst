#' Multitaper HPD periodogram matrix
#'
#' Given a multivariate time series, \code{pdPgram} computes a multitapered HPD periodogram matrix based on
#' averaging raw Hermitian PSD periodogram matrices of tapered multivariate time series segments.
#'
#' If \code{method = "multitaper"}, \code{pdPgram} calculates a \eqn{(d,d)}-dimensional multitaper
#' periodogram matrix based on \eqn{B} DPSS (Discrete Prolate Spheroidal Sequence or Slepian) orthogonal tapering functions
#' as in \code{\link[multitaper]{dpss}} applied to the \eqn{d}-dimensional time series \code{X}. If \code{method = "bartlett"},
#' \code{pdPgram} computes a Bartlett spectral estimator by averaging the periodogram matrices of \code{B} non-overlapping
#' segments of the \eqn{d}-dimensional time series \code{X}. Note that Bartlett's spectral estimator is a
#' specific (trivial) case of a multitaper spectral estimator with uniform orthogonal tapering windows.\cr
#' In the case of subsequent periodogram matrix denoising in the space of HPD matrices equipped with the
#' affine-invariant Riemannian metric, one should set \code{bias.corr = T}, thereby correcting for the asymptotic
#' bias of the periodogram matrix in the manifold of HPD matrices equipped with the affine-invariant metric as explained in
#' \insertCite{CvS17}{pdSpecEst} and Chapter 3 of \insertCite{C18}{pdSpecEst}. The pre-smoothed HPD periodogram matrix
#' (i.e., an initial noisy HPD spectral estimator) can be given as input to the function \code{\link{pdSpecEst1D}} to perform
#' intrinsic wavelet-based spectral matrix estimation. In this case, set \code{bias.corr = F} (the default) as the appropriate
#' bias-corrections are applied internally by the function \code{\link{pdSpecEst1D}}.
#'
#' @param X  an (\eqn{n,d})-dimensional matrix corresponding to a multivariate time series,
#' with the \code{d} columns corresponding to the components of the time series.
#' @param B  depending on the argument \code{method}, either the number of orthogonal DPSS tapers, or the number of
#' non-overlapping segments to compute Bartlett's averaged periodogram. By default,
#' \code{B = d}, such that the averaged periodogram is guaranteed to be positive definite.
#' @param method the tapering method, either \code{"multitaper"} or \code{"bartlett"} explained in the Details
#' section below. Defaults to \code{"multitaper"}.
#' @param bias.corr should an asymptotic bias-correction under the affine-invariant Riemannian metric be applied to
#' the HPD periodogram matrix? Defaults to \code{FALSE}.
#' @param nw a positive numeric value corresponding to the time-bandwidth parameter of the DPSS tapering functions,
#' see also \code{\link[multitaper]{dpss}}, defaults to \code{nw = 3}.
#'
#' @return A list containing two components:
#'    \item{\code{freq} }{ vector of \eqn{n/2} frequencies in the range \eqn{[0,0.5)} at which the periodogram is evaluated.}
#'    \item{\code{P} }{ a \eqn{(d, d, n/2)}-dimensional array containing the
#'      (\eqn{d,d})-dimensional multitaper periodogram matrices at frequencies corresponding
#'      to \code{freq}.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{pdPgram2D}}, \code{\link[multitaper]{dpss}}
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(200, 2, Phi, Theta, Sigma)
#' ts.plot(ts.sim$X) # plot generated time series traces
#' pgram <- pdPgram(ts.sim$X)
#'
#' @importFrom multitaper dpss
#'
#' @export
pdPgram <- function(X, B, method = c("multitaper", "bartlett"), bias.corr = F, nw = 3) {

  ## Initialize parameters
  if(missing(method)){
    method <- "multitaper"
  }
  method <- match.arg(method, c("multitaper", "bartlett"))
  d <- ncol(X)
  n <- nrow(X)
  if (missing(B)) {
    B <- d
  }
  if(B < d){
    warning("The number of tapers 'B' is smaller than the dimension of the time series 'ncol(X)';
            the periodogram matrix is not positive definite.")
  }
  ## Taper weights
  h <- switch(method,
              bartlett = matrix(0, 1, 1),
              multitaper = multitaper::dpss(n, B, nw = nw, returnEigenvalues = F)$v * sqrt(n))
  ## Compute periodogram C++
  P <- pgram_C(X, B, h, method, F)
  ## Bias-correction
  bias <- ifelse(bias.corr, B * exp(-1/d * sum(digamma(B - (d - 1:d)))), 1)

  return(list(freq = pi * (0:(dim(P)[3]-1))/dim(P)[3], P = bias * P))
}

#' Multitaper HPD time-varying periodogram matrix
#'
#' Given a multivariate time series, \code{pdPgram2D} computes a multitapered HPD time-varying periodogram matrix based on
#' averaging raw Hermitian PSD time-varying periodogram matrices of tapered multivariate time series segments.
#'
#' If \code{method = "dpss"}, \code{pdPgram2D} calculates a \eqn{(d,d)}-dimensional multitaper time-varying
#' periodogram matrix based on sliding \eqn{B} DPSS (Discrete Prolate Spheroidal Sequence or Slepian) orthogonal tapering functions
#' as in \code{\link[multitaper]{dpss}} applied to the \eqn{d}-dimensional time series \code{X}. If \eqn{B \ge d}, the
#' multitaper time-varying periodogram matrix is guaranteed to be positive definite at each time-frequency point in the
#' grid \code{expand.grid(tf.grid$time, tf.grid$frequency)}. In short, the function \code{pdPgram2D} computes a multitaper
#' periodogram matrix (as in \code{\link{pdPgram}}) in each of a number of non-overlapping time series
#' segments of \code{X}, with the time series segments centered around the (rescaled) time points in \code{tf.grid$time}.
#' If \code{method = "hermite"}, the function calculates a multitaper time-varying periodogram matrix replacing the DPSS
#' tapers by orthogonal Hermite tapering functions as in e.g., \insertCite{BB96}{pdSpecEst}. \cr
#' In the case of subsequent periodogram matrix denoising in the space of HPD matrices equipped with the
#' affine-invariant Riemannian metric, one should set \code{bias.corr = T}, thereby correcting for the asymptotic
#' bias of the periodogram matrix in the manifold of HPD matrices equipped with the affine-invariant metric as explained in
#' \insertCite{CvS17}{pdSpecEst} and Chapter 3 and 5 of \insertCite{C18}{pdSpecEst}. The pre-smoothed HPD periodogram matrix
#' (i.e., an initial noisy HPD spectral estimator) can be given as input to the function \code{\link{pdSpecEst2D}} to perform
#' intrinsic wavelet-based time-varying spectral matrix estimation. In this case, set \code{bias.corr = F} (the default) as the
#' appropriate bias-corrections are applied internally by the function \code{\link{pdSpecEst2D}}.
#'
#' @param X  an (\eqn{n,d})-dimensional matrix corresponding to a multivariate time series,
#' with the \code{d} columns corresponding to the components of the time series.
#' @param B  depending on the argument \code{method}, either the number of orthogonal DPSS or Hermite tapering functions.
#' By default, \code{B = d}, such that the multitaper periodogram is guaranteed to be positive definite.
#' @param tf.grid a list with two components \code{tf.grid$time} and \code{tf.grid$frequency} specifying the
#' rectangular grid of time-frequency points at which the multitaper periodogram is evaluated. \code{tf.grid$time}
#' should be a numeric vector of rescaled time points in the range \code{(0,1)}. \code{tf.grid$frequency} should be a numeric
#' vector of frequency points in the range \code{(0,0.5)}, with 0.5 corresponding to the Nyquist frequency.
#' @param method the tapering method, either \code{"dpss"} or \code{"hermite"} explained in the Details
#' section below. Defaults to \code{method = "dpss"}.
#' @param nw a positive numeric value corresponding to the time-bandwidth parameter of the tapering functions,
#' see also \code{\link[multitaper]{dpss}}, defaults to \code{nw = 3}. Both the DPSS and Hermite tapers are
#' rescaled with the same time-bandwidth parameter.
#' @param bias.corr should an asymptotic bias-correction under the affine-invariant Riemannian metric be applied to
#' the HPD periodogram matrix? Defaults to \code{FALSE}.
#'
#' @return A list containing two components:
#'    \item{\code{tf.grid} }{ a list with two components corresponding to the rectangular grid of time-frequency points
#'    at which the multitaper periodogram is evaluated.}
#'    \item{\code{P} }{ a \eqn{(d,d,m_1,m_2)}-dimensional array with \code{m_1 = length(tf.grid$time)} and
#'    \code{m_2 = length(tf.grid$frequency)} corresponding to the (\eqn{d,d})-dimensional tapered periodogram matrices
#'    evaluated at the time-frequency points in \code{tf.grid}.}
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{pdPgram}}, \code{\link[multitaper]{dpss}}
#'
#' @examples
#' ## Coefficient matrices
#' Phi1 <- array(c(0.4, 0, 0, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Phi2 <- array(c(0.8, 0, 0, 0.4, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#'
#' ## Generate piecewise stationary time series
#' ts.Phi <- function(Phi) rARMA(2^9, 2, Phi, Theta, Sigma)$X
#' ts <- rbind(ts.Phi(Phi1), ts.Phi(Phi2))
#'
#' pgram <- pdPgram2D(ts)
#'
#' @export
pdPgram2D <- function(X, B, tf.grid, method = c("dpss", "hermite"), nw = 3, bias.corr = F){

  ## Set variables
  d <- ncol(X)
  n <- nrow(X)
  if(missing(method)) {
    method <- "dpss"
  }
  method <- match.arg(method, c("dpss", "hermite"))
  if (missing(B)) {
    B <- d
  }
  if(B < d){
    warning("The number of tapers 'B' is smaller than the dimension of the time series 'ncol(X)';
             the periodogram matrix is not positive definite!")
  }
  L <- ifelse(missing(tf.grid), 2^round(log2(sqrt(n))), length(tf.grid$time))
  m <- 2 * floor(n / (2 * L))
  if(missing(tf.grid)){
    tf.grid <- list(time = seq(L + 1, n - L - 1, length = L)/n, freq = pi * 0:(m - 1) / m)
  }
  bias.corr <-ifelse(isTRUE(bias.corr & B >= d), B * exp(-1/d * sum(digamma(B - (d - 1:d)))), 1)

  ## Taper weights
  if(method == "dpss"){
    ## DPSS tapers
    h <- multitaper::dpss(m, B, nw, returnEigenvalues = F)$v * sqrt(m)
  } else if(method == "hermite"){
    ## Hermite tapers
    Hermite <- function(k, t){
      res <- numeric(length(t))
      for(ti in t){
        ## Hermite polynomials
        H <- numeric(k + 1)
        for(i in 0:k){
          H[i + 1] <- (if(i == 0) 1 else if(i == 1) 2 * ti else 2 * ti * H[i] - 2 * (i - 1) * H[i - 1])
        }
        res[which(t == ti)] <- exp(-ti^2 / 2) * H[k + 1] / sqrt(sqrt(pi) * 2^k * factorial(k))
      }
      return(res)
    }
    h0 <- sapply(0:(B-1), function(b) Hermite(b, seq(-nw, nw, length = m)))
    norm.h0 <- sqrt(mean(h0[, 1]^2))
    h <- apply(h0, 2, function(h.col) h.col / norm.h0)
  }

  ## Sliding dpss or hermite multitaper
  Per <- sapply(tf.grid$time, function(ti) {
                pgram_C(X[round(ti * n) + (-floor(m/2 - 1)):(m / 2), ], B, h, "multitaper", T)
                }, simplify = "array")

  return(list(tf.grid = tf.grid, P = bias.corr * aperm(Per, c(1, 2, 4, 3))))
}
