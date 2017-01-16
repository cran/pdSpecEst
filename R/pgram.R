#' Bias-corrected averaged periodogram
#'
#' \code{pdPgram} calculates an averaged (\eqn{d \times d})-dimensional Hermitian PD periodogram matrix
#' based on Bartlett's method by averaging the periodograms of non-overlapping segments
#' of the \eqn{d}-dimensional time series \code{X}. The averaged periodogram is rescaled by
#' applying the bias-correction as described in (Chau and von Sachs, 2017)
#'
#' @param X  a multivariate time series, with the \code{d} columns corresponding to the components
#'  of the time series.
#' @param B  the number of segments over which the averaged periodogram is computed. By default, \code{B = d},
#'  such that the averaged periodogram is guaranteed to be positive-definite.
#'
#' @return A list containing two components:
#'    \item{\code{freq} }{ vector of frequencies at which the periodogram is computed.}
#'    \item{\code{P} }{ a \code{(d, d, length(freq))}-dimensional array containing the
#'      (\eqn{d \times d})-dimensional averaged periodogram matrices at frequencies corresponding
#'      to \code{freq}.}
#'
#' @references Bartlett, M.S. (1950). \emph{Periodogram analysis and continuous spectra}.
#' Biometrika 37 (1-2): 1-16.
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @note The curve of HPD periodogram matrices obtained from \code{pdPgram(X)$P} can be
#' used as input in the functions \code{\link{WavTransf}} or \code{\link{pdSpecEst}}.
#'
#' @seealso \code{\link{rARMA}}
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' ts.sim <- rARMA(200, 2, Phi, Theta, Sigma)
#' ts.plot(ts.sim$X) # plot generated time series traces
#' pgram <- pdPgram(ts.sim$X)
#'
#' @importFrom astsa mvspec
#' @export
pdPgram <- function(X, B) {

  d <- ncol(X)
  n <- nrow(X)
  if (missing(B)) {
    B <- d
  }
  freq <- pi * (1:floor(n/B))/floor(n/B)
  Per <- sapply(1:B, function(b) 1/(2 * pi) * astsa::mvspec(X[floor(n/B) * (b - 1) + 1:floor(n/B), ],
                                                             plot = F)$fxx, simplify = "array")
  P <- B * exp(-1/d * sum(digamma(B - (d - 1:d)))) * apply(Per, c(1, 2, 3), mean)
  freq <- pi * (1:dim(P)[3])/dim(P)[3]

  return(list(freq = freq, P = P))
}

