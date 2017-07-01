#' Compute distance between two HPD matrices
#'
#' \code{pdDist} calculates a distance between two Hermitian PD matrices.
#'
#' Available distance measures between two Hermitian PD matrices are (i) Riemannian distance (default) as in
#' (Bhatia, 2009, Chapter 6), (ii) log-Euclidean distance, the Euclidean distance between matrix logarithms,
#' (iii) Cholesky distance, the Euclidean distance between Cholesky decompositions, (iv) Euclidean distance,
#' and (v) Procrustes distance as in (Dryden et al., 2009). In particular, \code{pdDist} generalizes the function
#' \code{\link[shapes]{distcov}} to compute the distance between two symmetric positive definite matrices to the
#' distance between two Hermitian positive definite matrices.
#'
#' @param A,B Hermitian positive definite matrices (of equal dimension).
#' @param method the distance measure, one of \code{'Riemannian'},
#' \code{'logEuclidean'}, \code{'Cholesky'}, \code{'Euclidean'} or \code{'Procrustes'}. Defaults to \code{'Riemannian'}.
#'
#' @examples
#'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  A <- t(Conj(a)) %*% a
#'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  B <- t(Conj(b)) %*% b
#'  pdDist(A, B) ## Riemannian distance
#'
#' @references Bhatia, R. (2009). \emph{Positive Definite Matrices}. New Jersey: Princeton University Press.
#' @references Dryden, I.L., Koloydenko, A., Zhou, D. (2009). Non-Euclidean statistics for covariance matrices,
#' with applications to diffusion tensor imaging. \emph{Annals of Applied Statistics}, 3(3), 1102-1123.
#'
#' @seealso \code{\link[shapes]{distcov}}
#'
#' @export
pdDist <- function(A, B, method = "Riemannian") {

  if (!(isTRUE(all.equal(dim(A), dim(B)) & (dim(A)[1] == dim(A)[2]) & (length(dim(A)) == 2)))) {
    stop("Incorrect input dimensions for arguments: 'A' and/or 'B',
         consult the function documentation for the requested inputs.")
  }
  method <- match.arg(method, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "Procrustes"))
  d <- nrow(A)

  if (method == "Riemannian") {
    dd <- RiemmDist(A, B)
  }
  if (method == "logEuclidean") {
    dd <- NormF(Logm(diag(d), A) - Logm(diag(d), B))
  }
  if (method == "Cholesky") {
    dd <- NormF(Chol(A) - Chol(B))
  }
  if (method == "Euclidean") {
    dd <- NormF(A - B)
  }
  if (method == "Procrustes") {
    l1 <- Sqrt(A)
    l2 <- Sqrt(B)
    dd <- sqrt(NormF(l1)^2 + NormF(l2)^2 - 2 * sum(svd(t(Conj(l2)) %*% l1)$d))
  }
  return(dd)
}
