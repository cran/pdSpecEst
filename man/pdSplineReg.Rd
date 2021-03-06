% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/benchmark.R
\name{pdSplineReg}
\alias{pdSplineReg}
\title{Cubic smoothing spline regression for HPD matrices}
\usage{
pdSplineReg(P, f0, lam = 1, Nd, ini_step = 1, max_iter = 100,
  eps = 0.001, ...)
}
\arguments{
\item{P}{a \eqn{(d,d,n)}-dimensional array corresponding to a length \eqn{n} sequence of (\eqn{d, d})-dimensional
noisy HPD matrices.}

\item{f0}{a \eqn{(d,d,n)}-dimensional array corresponding to an initial estimate of the smooth
target curve of (\eqn{d, d})-dimensional HPD matrices.}

\item{lam}{a smoothness penalty, defaults to \code{lam = 1}. If \code{lam = 0}, the penalized curve estimate
coincides with geodesic interpolation of the data points with respect to the Riemannian metric.
If \code{lam} increases to \eqn{\infty}, the penalized regression estimator is approximately a fitted geodesic curve.}

\item{Nd}{a numeric value (\code{Nd <= n}) determining a lower resolution of the cubic spline regression estimator to speed up
computation time, defaults to \eqn{n}.}

\item{ini_step}{initial candidate step size in a backtracking line search based on the Armijo-Goldstein
condition, defaults to \code{ini_step = 1}.}

\item{max_iter}{maximum number of gradient descent iterations, defaults to \code{max_iter = 100}.}

\item{eps}{optional tolerance parameter in gradient descent algorithm. The gradient descent procedure exits if the
absolute difference between consecutive evaluations of the objective function is smaller than \code{eps},
defaults to \code{eps = 1E-3}.}

\item{...}{additional arguments for internal use.}
}
\value{
A list with three components:
\item{f }{a \eqn{(d, d, N_d)}-dimensional array corresponding to a length \code{Nd} estimated cubic smoothing spline
curve of (\eqn{d, d})-dimensional HPD matrices.}
\item{cost }{a numeric vector containing the costs of the objective function at each gradient descent iteration.}
\item{total_iter }{total number of gradient descent iterations.}
}
\description{
\code{pdSplineReg()} performs cubic smoothing spline regression in the space of HPD matrices equipped with the
affine-invariant Riemannian metric through minimization of a penalized regression objective function using a
geometric conjugate gradient descent method as outlined in \insertCite{BA11}{pdSpecEst} and \insertCite{BA11b}{pdSpecEst}.
This is a specific implementation of the more general algorithm in \insertCite{BA11}{pdSpecEst} and \insertCite{BA11b}{pdSpecEst},
setting the part in the objective function based on the first-order finite geometric differences to zero, such that the solutions
of the regression problem are approximating cubic splines.
}
\note{
This function does not check for positive definiteness of the matrices given as input, and may fail
if matrices are close to being singular.
}
\examples{
\dontrun{
set.seed(2)
P <- rExamples1D(50, example = 'gaussian', noise.level = 0.1)
P.spline <- pdSplineReg(P$P, P$P, lam = 0.5, Nd = 25)

## Examine matrix-component (1,1)
plot((1:50)/50, Re(P$P[1, 1, ]), type = "l", lty = 2) ## noisy observations
lines((1:25)/25, Re(P.spline$f[1, 1, ])) ## estimate
lines((1:50)/50, Re(P$f[1, 1, ]), col = 2, lty = 2) ## smooth target
}
}
\references{
\insertAllCited{}
}
