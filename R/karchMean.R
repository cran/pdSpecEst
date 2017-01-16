#' Approximate weighted Karcher mean
#'
#' \code{KarchMean} calculates the approximate weighted Karcher mean of \eqn{S} different
#' \eqn{(d \times d)}-dimensional Hermitian PD matrices by the recursive algorithm described
#' in (Chau and von Sachs, 2017). By default, the unweighted Karcher mean is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' Ave <- KarchMean(M, w)
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive-definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @seealso \code{\link{Mid}}
#'
#' @export
KarchMean <- function(M, w) {
  if (missing(w)) {
    w <- rep(1/dim(M)[3], dim(M)[3])
  }
  return(kMean(do.call(rbind, lapply(1:dim(M)[3], function(s) M[, , s])), w))
}







