#' Compute distance between two HPD matrices
#'
#' \code{pdDist} calculates a distance between two Hermitian PD matrices.
#'
#' Available distance measures between two HPD matrices are: (i) the affine-invariant Riemannian distance (default) as in
#' e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}; (ii) the Log-Euclidean distance,
#' the Euclidean distance between matrix logarithms; (iii) the Cholesky distance, the Euclidean distance between Cholesky decompositions;
#' (iv) the Euclidean distance; (v) the root-Euclidean distance; and (vi) the Procrustes distance as in \insertCite{D09}{pdSpecEst}.
#' In particular, \code{pdDist} generalizes the function \code{shapes::distcov}, to compute the distance between two symmetric positive
#' definite matrices, in order to compute the distance between two Hermitian positive definite matrices.
#'
#' @param A,B Hermitian positive definite matrices (of equal dimension).
#' @param metric the distance measure, one of \code{'Riemannian'},
#' \code{'logEuclidean'}, \code{'Cholesky'}, \code{'Euclidean'}, \code{'rootEuclidean'} or \code{'Procrustes'}.
#' Defaults to \code{'Riemannian'}.
#'
#' @examples
#'  a <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  A <- t(Conj(a)) %*% a
#'  b <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#'  B <- t(Conj(b)) %*% b
#'  pdDist(A, B) ## Riemannian distance
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdDist <- function(A, B, metric = "Riemannian") {

  if (!(isTRUE(all.equal(dim(A), dim(B)) & (dim(A)[1] == dim(A)[2]) & (length(dim(A)) == 2)))) {
    stop("Incorrect input dimensions for arguments: 'A' and/or 'B',
         consult the function documentation for the requested inputs.")
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean", "Procrustes"))

  if (metric != "Procrustes") {
    ## Compute distance C++
    dd <- pdDist_C(A, B, metric)
  } else {
    l1 <- Sqrt(A)
    l2 <- Sqrt(B)
    dd <- sqrt(NormF(l1)^2 + NormF(l2)^2 - 2 * sum(svd(t(Conj(l2)) %*% l1)$d))
  }
  return(dd)
}

#' Weighted Karcher mean of HPD matrices
#'
#' \code{pdMean} calculates an (approximate) weighted Karcher or Frechet mean of a sample of
#' \eqn{(d,d)}-dimensional HPD matrices intrinsic to a user-specified metric. In the case of the
#' affine-invariant Riemannian metric as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or
#' \insertCite{PFA05}{pdSpecEst}, the weighted Karcher mean is either approximated via
#' the fast recursive algorithm in \insertCite{H13}{pdSpecEst} or computed via the slower, but more accurate,
#' gradient descent algorithm in \insertCite{P06}{pdSpecEst}. By default, the unweighted Karcher mean is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array corresponding to a sample of \eqn{(d,d)}-dimensional HPD matrices of
#' size \eqn{S}.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param metric the distance measure, one of \code{'Riemannian'}, \code{'logEuclidean'},
#' \code{'Cholesky'}, \code{'Euclidean'} or \code{'rootEuclidean'}. Defaults to \code{'Riemannian'}.
#' @param grad_desc if \code{metric = "Riemannian"}, a logical value indicating if the
#' gradient descent algorithm in \insertCite{P06}{pdSpecEst} should be used, defaults to \code{FALSE}.
#' @param maxit maximum number of iterations in gradient descent algorithm, only used if
#' \code{grad_desc = TRUE} and \code{metric = "Riemannian"}. Defaults to \code{1000}
#' @param reltol optional tolerance parameter in gradient descent algorithm, only used if
#' \code{grad_desc = TRUE} and \code{metric = "Riemannian"}. Defaults to \code{1E-10}.
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and (depending on the
#' specified metric) may fail if matrices are close to being singular.
#'
#' @examples
#' ## Generate random sample of HPD matrices
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' ## Generate random weight vector
#' w <- abs(z)/sum(abs(z))
#' ## Compute weighted (Riemannian) Karcher mean
#' pdMean(M, w)
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{Mid}}, \code{\link{pdMedian}}
#'
#' @export
pdMean <- function(M, w, metric = "Riemannian", grad_desc = FALSE, maxit = 1000, reltol) {

  if (!(isTRUE(is.array(M) & (dim(M)[1] == dim(M)[2]) & (length(dim(M)) == 3)))) {
    stop("Incorrect input dimensions for arguments: 'M',
         consult the function documentation for the requested inputs.")
  }
  ## Set parameters
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  w <- (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d <- dim(M)[1]

  if(dim(M)[3] == 1){
    Mean <- M[, , 1]
  } else{

    if(metric == "Riemannian"){
      ## Recursive algorithm
      Mean <- pdMean_C_approx(M, w)
      ## Gradient descent algorithm
      if(grad_desc){
        reltol <- (if(missing(reltol)) 1E-10 else reltol)
        if(!isTRUE(is.numeric(reltol) & is.numeric(maxit) & maxit > 0)){
          stop("Incorrect input for arguments: 'reltol' or 'maxit'.")
        }
        Mean <- pdMean_C(Mean, M, w, round(maxit), reltol)
      }
    } else {
      ## Transform
      if(metric %in% c("logEuclidean", "Cholesky", "rootEuclidean")) {
        M <- Ptransf2D_C(M, FALSE, FALSE, metric)
      }
      ## Euclidean mean
      Mean <- apply(sweep(M, 3, w, "*"), c(1, 2), sum)
      ## Transform back
      Mean <- switch(metric,
                       logEuclidean = Expm(diag(d), Mean),
                       Cholesky = Chol_C(Mean, FALSE, TRUE),
                       rootEuclidean = t(Conj(Mean)) %*% Mean,
                       Euclidean = Mean)
    }
  }
  return(Mean)
}

#' Weighted intrinsic median of HPD matrices
#'
#' \code{pdMedian} calculates a weighted intrinsic median of a sample of \eqn{(d,d)}-dimensional
#' HPD matrices based on a Weiszfeld algorithm intrinsic to the chosen metric.In the case of the
#' affine-invariant Riemannian metric as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or
#' \insertCite{PFA05}{pdSpecEst}, the intrinsic Weiszfeld algorithm in \insertCite{F09}{pdSpecEst} is used.
#' By default, the unweighted intrinsic median is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array corresponding to a sample of \eqn{(d,d)}-dimensional HPD matrices of
#' size \eqn{S}.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param metric the distance measure, one of \code{'Riemannian'}, \code{'logEuclidean'},
#' \code{'Cholesky'}, \code{'Euclidean'} or \code{'rootEuclidean'}. Defaults to \code{'Riemannian'}.
#' @param maxit maximum number of iterations in gradient descent algorithm. Defaults to \code{1000}
#' @param reltol optional tolerance parameter in gradient descent algorithm. Defaults to \code{1E-10}.
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and (depending on the
#' specified metric) may fail if matrices are close to being singular.
#'
#' @examples
#' ## Generate random sample of HPD matrices
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' ## Generate random weight vector
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' ## Compute weighted intrinsic (Riemannian) median
#' pdMedian(M, w)
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{pdMean}}
#'
#' @export
pdMedian <- function(M, w, metric = "Riemannian", maxit = 1000, reltol) {

  if (!(isTRUE(is.array(M) & (dim(M)[1] == dim(M)[2]) & (length(dim(M)) == 3)))) {
    stop("Incorrect input dimensions for arguments: 'M',
         consult the function documentation for the requested inputs.")
  }
  ## Set parameters
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  w <- (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d <- dim(M)[1]
  reltol <- (if(missing(reltol)) 1E-10 else reltol)
  if(!isTRUE(is.numeric(reltol) & is.numeric(maxit) & maxit > 0)){
    stop("Incorrect input for arguments: 'reltol' or 'maxit'.")
  }

  if(dim(M)[3] == 1){
    Med <- M[, , 1]
  } else{

    if(metric == "Riemannian"){
      ## Initial estimate
      Med0 <- pdMean_C_approx(M, w)
      ## Weiszfeld algorithm
      Med <- pdMedian_C(Med0, M, w, round(maxit), reltol)
    } else {
      ## Transform
      if(metric %in% c("logEuclidean", "Cholesky", "rootEuclidean")) {
        M <- Ptransf2D_C(M, FALSE, FALSE, metric)
      }
      ## Initial estimate
      Med0 <- apply(sweep(M, 3, w, "*"), c(1, 2), sum)
      ## Euclidean Weiszfeld algorithm
      Med <- Euclid_Median_C(Med0, M, w, round(maxit), reltol)
      ## Transform back
      Med <- switch(metric,
                     logEuclidean = Expm(diag(d), Med),
                     Cholesky = Chol_C(Med, FALSE, TRUE),
                     rootEuclidean = t(Conj(Med)) %*% Med,
                     Euclidean = Med)
    }
  }
  return(Med)
}


