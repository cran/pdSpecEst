#' Forward MI wavelet transform
#'
#' \code{WavTransf} computes the forward \emph{midpoint-interpolation} (MI) wavelet transform of a
#' curve of length \eqn{m} of (\eqn{d \times d})-dimensional Hermitian PD matrices as described in
#' (Chau and von Sachs, 2017).
#'
#' @param P a (\eqn{d,d,m})-dimensional array of Hermitian PD matrices, with \eqn{m = 2^J} for some \eqn{J > 0}.
#' @param order an odd integer between 1 and 9 corresponding to the order of the MI refinement scheme.
#' @param jmax the maximum scale up to which the wavelet coefficients are computed. If \code{jmax} is not
#' specified it is set equal to the maximum possible scale \code{jmax = J-1}.
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
#' D <- WavTransf(pgram$P)
#'
#' @return The function returns a list with two components:
#' \item{D }{a list of arrays, where each (\eqn{d, d, 2^j})-dimensional array contains the (\eqn{d \times d})-
#' dimensional wavelet coefficients at the \eqn{2^j} different locations in the given wavelet scale. The first
#' list element is a (\eqn{d, d, 2})-dimensional array containing the (\eqn{d \times d})-dimensional midpoints
#' at the coarsest scale (\eqn{j = 1}) in the midpoint pyramid.}
#' \item{M }{a list of arrays, where each (\eqn{d, d, 2^j})-dimensional array contains the (\eqn{d \times d})-
#' dimensional midpoints at the \eqn{2^j} different locations in the given midpoint scale. The first list
#' element is equivalent to the first list element in the \code{$D} component.}
#'
#' @seealso \code{\link{InvWavTransf}}, \code{\link{pdSpecEst}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
WavTransf <- function(P, order = 5, jmax) {

  ## Set variables
  J <- log2(dim(P)[3])
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                " to dyadic number."))
  }
  if (!(order %in% c(1, 3, 5, 7, 9))) {
    warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
    order <- 5
  }
  d <- dim(P)[1]
  M <- list()
  for (j in J:1) {
    if (j == J) {
      Mper <- unname(P)
    } else {
      Mper <- sapply(1:(dim(Mper)[3]/2), function(i) Mid(Mper[, , 2 * i - 1], Mper[, , 2 * i]),
                     simplify = "array")
    }
    M[[j]] <- array(Mper[, , 1:2^j], dim = c(d, d, 2^j))
  }
  names(M) <- paste0("M.scale", 1:J)
  D <- list(M.scale1 = M[[1]])
  Nw <- list(N1 = 1, N3 = c(-1, 8, 1)/8, N5 = c(3, -22, 128, 22, -3)/128,
             N7 = c(-5, 44, -201, 1024, 201, -44, 5)/1024,
             N9 = c(35, -370, 1898, -6922, 32768, 6922, -1898, 370, -35)/32768)
  if (missing(jmax)) {
    jmax <- J - 1
  }
  if (jmax > J - 1) {
    warning(paste0("jmax cannot exceed maximum scale j=", J - 1))
    jmax <- J - 1
  }

  ## Compute wavelet transform
  for (j in 1:jmax) {
    tm1 <- Impute_man(M[[j]], (order - 1)/2, Nw)
    iSqrt_tm1 <- sapply(1:2^j, function(l) iSqrt(tm1[, , l]), simplify = "array")
    D[[j + 1]] <- sapply(1:2^j, function(l) Logm(diag(d), (iSqrt_tm1[, , l] %*%
                                                             M[[j + 1]][, , 2 * l]) %*% iSqrt_tm1[, , l]), simplify = "array")
    names(D)[j + 1] <- paste0("D.scale", j)
  }

  return(list(D = D, M = M))
}
