#' Inverse MI wavelet transform
#'
#' \code{InvWavTransf} computes the inverse \emph{midpoint-interpolation} (MI) wavelet
#' transform of an array of coarse-scale Hermitian PD midpoints combined with a pyramid of Hermitian
#' matrix-valued wavelet coefficients as described in (Chau and von Sachs, 2017).
#'
#' @param D a list of arrays containing coarse-scale midpoints and Hermitian matrix-valued wavelet
#'  coefficients in the same format as the \code{$D} component given as output by the function
#'  \code{\link{WavTransf}}.
#' @param order an odd integer between 1 and 9 corresponding to the order of the MI refinement scheme.
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
#' D <- WavTransf(pgram$P)$D
#' P <- InvWavTransf(D)
#' all.equal(pgram$P, P)
#'
#' @return Returns a (\eqn{d, d, m})-dimensional array corresponding to a curve of length \eqn{m} of
#' (\eqn{d \times d})-dimensional Hermitian PD matrices.
#'
#' @seealso \code{\link{WavTransf}}, \code{\link{pdSpecEst}}
#'
#' @references Chau, J. and von Sachs, R. (2017) \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
InvWavTransf <- function(D, order = 5) {

  if (!(order %in% c(1, 3, 5, 7, 9))) {
    warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  d <- nrow(D[[1]][, , 1])
  J <- length(D)
  Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
             N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
             N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)
  m1 <- D[[1]]
  for (j in 1:(J - 1)) {
    tm1 <- Impute_man(m1, (order - 1)/2, Nw)
    m2 <- array(dim = c(d, d, 2^(j + 1)))
    reconstr <- function(i) {
      if (any(c(D[[j + 1]][, , i]) != 0)) {
        Sqrt_tm1 <- Sqrt(tm1[, , i])
        m2_even <- (Sqrt_tm1 %*% Expm(diag(d), D[[j + 1]][, , i])) %*% Sqrt_tm1
      } else {
        m2_even <- tm1[, , i]
      }
      return(m2_even)
    }
    m2[, , c(F, T)] <- sapply(1:2^j, reconstr, simplify = "array")
    m2[, , c(T, F)] <- sapply(1:2^j, function(i) solveMid(m2[, , 2 * i], m1[, , i]), simplify = "array")
    m1 <- m2
  }
  return(m1)
}
