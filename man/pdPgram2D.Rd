% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pgram.R
\name{pdPgram2D}
\alias{pdPgram2D}
\title{Multitaper HPD time-varying periodogram matrix}
\usage{
pdPgram2D(X, B, tf.grid, method = c("dpss", "hermite"), nw = 3,
  bias.corr = F)
}
\arguments{
\item{X}{an (\eqn{n,d})-dimensional matrix corresponding to a multivariate time series,
with the \code{d} columns corresponding to the components of the time series.}

\item{B}{depending on the argument \code{method}, either the number of orthogonal DPSS or Hermite tapering functions.
By default, \code{B = d}, such that the multitaper periodogram is guaranteed to be positive definite.}

\item{tf.grid}{a list with two components \code{tf.grid$time} and \code{tf.grid$frequency} specifying the
rectangular grid of time-frequency points at which the multitaper periodogram is evaluated. \code{tf.grid$time}
should be a numeric vector of rescaled time points in the range \code{(0,1)}. \code{tf.grid$frequency} should be a numeric
vector of frequency points in the range \code{(0,0.5)}, with 0.5 corresponding to the Nyquist frequency.}

\item{method}{the tapering method, either \code{"dpss"} or \code{"hermite"} explained in the Details
section below. Defaults to \code{method = "dpss"}.}

\item{nw}{a positive numeric value corresponding to the time-bandwidth parameter of the tapering functions,
see also \code{\link[multitaper]{dpss}}, defaults to \code{nw = 3}. Both the DPSS and Hermite tapers are
rescaled with the same time-bandwidth parameter.}

\item{bias.corr}{should an asymptotic bias-correction under the affine-invariant Riemannian metric be applied to
the HPD periodogram matrix? Defaults to \code{FALSE}.}
}
\value{
A list containing two components:
   \item{\code{tf.grid} }{ a list with two components corresponding to the rectangular grid of time-frequency points
   at which the multitaper periodogram is evaluated.}
   \item{\code{P} }{ a \eqn{(d,d,m_1,m_2)}-dimensional array with \code{m_1 = length(tf.grid$time)} and
   \code{m_2 = length(tf.grid$frequency)} corresponding to the (\eqn{d,d})-dimensional tapered periodogram matrices
   evaluated at the time-frequency points in \code{tf.grid}.}
}
\description{
Given a multivariate time series, \code{pdPgram2D} computes a multitapered HPD time-varying periodogram matrix based on
averaging raw Hermitian PSD time-varying periodogram matrices of tapered multivariate time series segments.
}
\details{
If \code{method = "dpss"}, \code{pdPgram2D} calculates a \eqn{(d,d)}-dimensional multitaper time-varying
periodogram matrix based on sliding \eqn{B} DPSS (Discrete Prolate Spheroidal Sequence or Slepian) orthogonal tapering functions
as in \code{\link[multitaper]{dpss}} applied to the \eqn{d}-dimensional time series \code{X}. If \eqn{B \ge d}, the
multitaper time-varying periodogram matrix is guaranteed to be positive definite at each time-frequency point in the
grid \code{expand.grid(tf.grid$time, tf.grid$frequency)}. In short, the function \code{pdPgram2D} computes a multitaper
periodogram matrix (as in \code{\link{pdPgram}}) in each of a number of non-overlapping time series
segments of \code{X}, with the time series segments centered around the (rescaled) time points in \code{tf.grid$time}.
If \code{method = "hermite"}, the function calculates a multitaper time-varying periodogram matrix replacing the DPSS
tapers by orthogonal Hermite tapering functions as in e.g., \insertCite{BB96}{pdSpecEst}. \cr
In the case of subsequent periodogram matrix denoising in the space of HPD matrices equipped with the
affine-invariant Riemannian metric, one should set \code{bias.corr = T}, thereby correcting for the asymptotic
bias of the periodogram matrix in the manifold of HPD matrices equipped with the affine-invariant metric as explained in
\insertCite{CvS17}{pdSpecEst} and Chapter 3 and 5 of \insertCite{C18}{pdSpecEst}. The pre-smoothed HPD periodogram matrix
(i.e., an initial noisy HPD spectral estimator) can be given as input to the function \code{\link{pdSpecEst2D}} to perform
intrinsic wavelet-based time-varying spectral matrix estimation. In this case, set \code{bias.corr = F} (the default) as the
appropriate bias-corrections are applied internally by the function \code{\link{pdSpecEst2D}}.
}
\examples{
## Coefficient matrices
Phi1 <- array(c(0.4, 0, 0, 0.8, rep(0, 4)), dim = c(2, 2, 2))
Phi2 <- array(c(0.8, 0, 0, 0.4, rep(0, 4)), dim = c(2, 2, 2))
Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)

## Generate piecewise stationary time series
ts.Phi <- function(Phi) rARMA(2^9, 2, Phi, Theta, Sigma)$X
ts <- rbind(ts.Phi(Phi1), ts.Phi(Phi2))

pgram <- pdPgram2D(ts)

}
\references{
\insertAllCited{}
}
\seealso{
\code{\link{pdPgram}}, \code{\link[multitaper]{dpss}}
}
