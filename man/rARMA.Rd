% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/examples.R
\name{rARMA}
\alias{rARMA}
\title{Simulate vARMA(2,2) time series}
\usage{
rARMA(n, d, Phi, Theta, Sigma, burn = 100, freq = NULL)
}
\arguments{
\item{n}{number of time series observations to be generated.}

\item{d}{dimension of the multivariate time series.}

\item{Phi}{a (\eqn{d, d, 2})-dimensional array, with \code{Phi[, , 1]} and \code{Phi[, , 2]}
the autoregressive (AR) coefficient matrices.}

\item{Theta}{a (\eqn{d, d, 2})-dimensional array, with \code{Theta[, , 1]} and \code{Theta[, , 2]}
the moving-average (MA) coefficient matrices.}

\item{Sigma}{the covariance matrix of the Gaussian white noise component.}

\item{burn}{a burn-in period when generating the time series observations, by default \code{burn = 100}.}

\item{freq}{an optional vector of frequencies, if \code{!is.null(freq)} the function also returns the
underlying Fourier spectral matrix of the stationary generating process evaluated at the frequencies in \code{freq}.}
}
\value{
The function returns a list with two components:
   \item{\code{X} }{ generated time series observations, the \code{d} columns correspond to the components of
    the multivariate time series.}
   \item{\code{f} }{ if \code{!is.null(freq)}, \code{f} is a \code{(d, d, length(freq))}-dimensional array corresponding
   to the underlying Fourier spectral matrix curve of \eqn{(d,d)}-dimensional HPD matrices evaluated at the frequencies
   in \code{freq}. If \code{is.null(freq)}, \code{f} is set to \code{NULL}.}
}
\description{
\code{rARMA} generates \code{d}-dimensional time series observations from a vARMA(2,2)
(vector-autoregressive-moving-average) process based on Gaussian white noise for testing and simulation
purposes.
}
\examples{
## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
freq <- seq(from = pi / 100, to = pi, length = 100)
Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
ts.sim <- rARMA(200, 2, Phi, Theta, Sigma, freq = freq)
ts.plot(ts.sim$X) # plot generated time series traces.

}
\references{
\insertRef{BD06}{pdSpecEst}
}
