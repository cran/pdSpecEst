#' K-means clustering for HPD matrices
#'
#' \code{pdkMeans} performs (fuzzy) k-means clustering for collections of HPD matrices, such as covariance or
#' spectral density matrices, based on a number of different metrics in the space of HPD matrices.
#'
#' The input array \code{X} corresponds to a collection of \eqn{(d,d)}-dimensional HPD matrices
#' for \eqn{S} different subjects. If the fuzziness parameter satisfies \code{m > 1}, the \eqn{S} subjects are assigned to
#' \eqn{K} different clusters in a probabilistic fashion according to a fuzzy k-means algorithm as detailed in classical texts,
#' such as \insertCite{BE81}{pdSpecEst}. If \code{m = 1}, the \eqn{S} subjects are assigned to the \eqn{K} clusters in a non-probabilistic
#' fashion according to a standard (hard) k-means algorithm. If not specified by the user, the \eqn{K} cluster
#' centers are initialized by random sampling without replacement from the input array of HPD matrices \code{X}.
#' The distance measure in the (fuzzy) k-means algorithm is induced by the metric on the space of HPD matrices specified by the user.
#' By default, the space of HPD matrices is equipped with (i) the affine-invariant Riemannian metric (\code{metric = 'Riemannian'})
#' as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}. Instead, this can also be one of:
#' (ii) the log-Euclidean metric (\code{metric = 'logEuclidean'}), the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric (\code{metric = 'Cholesky'}), the Euclidean inner product between Cholesky decompositions; (iv) the
#' Euclidean metric (\code{metric = 'Euclidean'}); or (v) the root-Euclidean metric (\code{metric = 'rootEuclidean'}). The default
#' choice of metric (affine-invariant Riemannian) satisfies several useful properties not shared by the other metrics, see e.g.,
#' \insertCite{C18}{pdSpecEst} for more details. Note that this comes at the cost of increased computation time in comparison to one
#' of the other metrics.
#'
#' @param X a (\eqn{d,d,S})-dimensional array of (\eqn{d,d})-dimensional HPD matrices for \eqn{S}
#' different subjects. Also accepts a (\eqn{d,d,n,S})-dimensional array, which is understood to be an array of
#' \eqn{n}-dimensional sequences of (\eqn{d,d})-dimensional HPD matrices for \eqn{S} different subjects.
#' @param K the number of clusters, a positive integer larger than 1.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. Additional details are given below.
#' @param m a fuzziness parameter larger or equal to \eqn{1}. If \eqn{m = 1} the cluster assignments are no longer fuzzy,
#' i.e., the procedure performs hard clustering. Defaults to \code{m = 1}.
#' @param eps an optional tolerance parameter determining the stopping criterion. The k-means algorithm
#' terminates if the intrinsic distance between cluster centers is smaller than \code{eps}, defaults to \code{eps = 1e-05}.
#' @param max_iter an optional parameter tuning the maximum number of iterations in the
#' k-means algorithm, defaults to \code{max_iter = 100}.
#' @param centroids an optional (\eqn{d,d,K})- or (\eqn{d,d,n,K})-dimensional array depending on the input array \code{X}
#' specifying the initial cluster centroids. If not specified, \code{K} initial cluster centroids are randomly sampled without
#' replacement from the input array \code{X}.
#'
#' @return Returns a list with two components:
#' \describe{
#'   \item{cl.assignments }{ an (\eqn{S,K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#'   matrix corresponds to the (probabilistic or binary) cluster membership assignment of subject \eqn{s} with respect
#'   to cluster \eqn{k}.}
#'   \item{cl.centroids }{ either a (\eqn{d,d,K})- or (\eqn{d,d,n,K})-dimensional array depending on the input array \code{X}
#'   corresponding respectively to the \code{K} \eqn{(d,d)}- or (\eqn{d,d,n})-dimensional final cluster centroids.
#'   }
#' }
#'
#' @examples
#' ## Generate 20 random HPD matrices in 2 groups
#' m <- function(rescale){
#'  x <- matrix(complex(real = rescale * rnorm(9), imaginary = rescale * rnorm(9)), nrow = 3)
#'  t(Conj(x)) %*% x
#' }
#' X <- array(c(replicate(10, m(0.25)), replicate(10, m(1))), dim = c(3, 3, 20))
#'
#' ## Compute fuzzy k-means cluster assignments
#' cl <- pdkMeans(X, K = 2, m = 2)$cl.assignments
#'
#' @seealso \code{\link{pdDist}}, \code{\link{pdSpecClust1D}}, \code{\link{pdSpecClust2D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdkMeans <- function(X, K, metric = "Riemannian", m = 1, eps = 1e-05, max_iter = 100, centroids) {

  ## Initialize parameters
  d <- dim(X)[1]
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  if(!isTRUE(length(dim(X)) %in% c(3, 4) & dim(X)[1] == dim(X)[2])){
    stop("Incorrect input dimensions for argument: 'X',
                  consult the function documentation for the requested inputs.")
  }
  S <- tail(dim(X), 1)
  n <- ifelse(length(dim(X)) == 3, 1, dim(X)[3])
  ## Initialize centers and reshape X if necessary
  if(length(dim(X)) == 3) {
    init_cent <- (if(missing(centroids)) X[, , sample(1:S, K), drop = FALSE] else centroids)
  } else {
    init_cent <- (if(missing(centroids)) X[, , , sample(1:S, K), drop = FALSE] else centroids)
    init_cent <- array(init_cent, dim = c(d, d, n * K))
    X <- array(aperm(X, c(1, 2, 4, 3)), dim = c(d, d, S * n))
  }

  ## C++ fuzzy c-means algorithm
  clust <- cMeans_C(X, init_cent, S, K, m, eps, max_iter, metric, matrix(1, S, K))[, , 1]
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  ## Recompute final centers
  centroids <- array(dim = c(d, d, n, K))
  for(k in 1:K) {
    w <- clust[, k]^m / sum(clust[, k]^m)
    if(length(dim(X)) == 3) {
      centroids[, , , k] <- pdMean(X, w, metric, grad_desc = TRUE)
    } else {
      centroids[, , , k] <- sapply(1:n, function(i) pdMean(X[, , i, ], w, metric, grad_desc = TRUE),
                                   simplify = "array")
    }
  }
  ## Return result
  return(list(cl.assignments = clust, cl.centroids = centroids[, , 1:n, ]))
}

#' Intrinsic wavelet HPD spectral matrix clustering
#'
#' \code{pdSpecClust1D} performs clustering of HPD spectral matrices corrupted by noise (e.g. HPD periodograms)
#' by combining wavelet thresholding and fuzzy clustering in the intrinsic wavelet coefficient domain according to
#' the following steps:
#' \enumerate{
#'   \item Transform a collection of noisy HPD spectral matrices to the intrinsic wavelet domain and denoise the
#'   HPD matrix curves by (tree-structured) thresholding of wavelet coefficients with \code{\link{pdSpecEst1D}}.
#'   \item Apply an intrinsic fuzzy c-means algorithm to the coarsest midpoints at scale \code{j = 0} across subjects.
#'   \item Taking into account the fuzzy cluster assignments in the previous step, apply a weighted fuzzy c-means
#'   algorithm to the nonzero thresholded wavelet coefficients across subjects from scale \code{j = 1} up to \code{j = jmax}.
#' }
#' More details can be found in Chapter 3 of \insertCite{C18}{pdSpecEst} and the accompanying vignettes.
#'
#' The input array \code{P} corresponds to a collection of initial noisy HPD spectral estimates of the \eqn{(d,d)}-dimensional
#' spectral matrix at \code{n} different frequencies, with \eqn{n = 2^J} for some \eqn{J > 0}, for \eqn{S} different subjects.
#' These can be e.g., multitaper HPD periodograms given as output by the function \code{\link{pdPgram}}.\cr
#' First, for each subject \eqn{s = 1,\ldots,S}, thresholded wavelet coefficients in the intrinsic wavelet domain are
#' calculated by \code{\link{pdSpecEst1D}}, see the function documentation for additional details on the wavelet thresholding
#' procedure.\cr
#' The maximum wavelet scale taken into account in the clustering procedure is determined by the arguments
#' \code{jmax} and \code{d.jmax}. The maximum scale is set to the minimum of \code{jmax} and the wavelet
#' scale \eqn{j} for which the proportion of nonzero thresholded wavelet coefficients (averaged
#' across subjects) is smaller than \code{d.jmax}.\cr
#' The \eqn{S} subjects are assigned to \eqn{K} different clusters in a probabilistic fashion according to a
#' two-step procedure:
#' \enumerate{
#' \item In the first step, an intrinsic fuzzy c-means algorithm, with fuzziness parameter \eqn{m} is applied to the
#' \eqn{S} coarsest midpoints at scale \code{j = 0} in the subject-specific midpoint pyramids. Note that the distance
#' function in the intrinsic c-means algorithm relies on the chosen metric on the space of HPD matrices.
#' \item In the second step, a weighted fuzzy c-means algorithm based on the Euclidean
#' distance function, also with fuzziness parameter \eqn{m}, is applied to the nonzero thresholded wavelet
#' coefficients of the \eqn{S} different subjects. The tuning parameter \code{tau} controls the weight given
#' to the cluster assignments obtained in the first step of the clustering algorithm.
#' }
#' The function computes the forward and inverse intrinsic AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{CvS17}{pdSpecEst} or \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics. \cr
#' If \code{return.centers = TRUE}, the function also returns the \code{K} HPD spectral matrix curves corresponding to
#' the cluster centers based on the given metric by applying the intrinsic inverse AI wavelet transform (
#' \code{\link{InvWavTransf1D}}) to the cluster centers in the wavelet domain.
#'
#' @param P a (\eqn{d,d,n,S})-dimensional array of HPD matrices, corresponding to a collection of sequences of
#' \eqn{(d,d)}-dimensional HPD matrices of length \eqn{n}, with \eqn{n = 2^J} for some \eqn{J > 0}, for \eqn{S} different
#' subjects.
#' @param K the number of clusters, a positive integer larger than 1.
#' @param jmax an upper bound on the maximum wavelet scale to be considered in the clustering procedure. If
#' \code{jmax} is not specified, it is set equal to the maximum (i.e., finest) wavelet scale minus 2.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. Additional details are given below.
#' @param m the fuzziness parameter for both fuzzy c-means algorithms. \code{m} should be larger or equal to \eqn{1}.
#' If \eqn{m = 1} the cluster assignments are no longer fuzzy, i.e., the procedure
#' performs hard clustering.
#' @param d.jmax a proportion that is used to determine the maximum wavelet scale to be considered in the clustering
#' procedure. A larger value \code{d.jmax} leads to less wavelet coefficients being taken into account, and therefore
#' lower computational effort in the procedure. If \code{d.jmax} is not specified, by default \code{d.jmax = 0.1}.
#' @param eps an optional vector with two components determining the stopping criterion. The first step in the cluster procedure
#' terminates if the (integrated) intrinsic distance between cluster centers is smaller than \code{eps[1]}.
#' The second step in the cluster procedure terminates if the (integrated) Euclidean distance between cluster centers is smaller
#' than \code{eps[2]}. By default \code{eps = c(1e-04, 1e-04)}.
#' @param tau an optional argument tuning the weight given to the cluster assignments obtained in the first step of
#' the clustering algorithm. If \code{tau} is not specified, by default \code{tau = 0.5}.
#' @param max_iter an optional argument tuning the maximum number of iterations in both the first and second step of the
#' clustering algorithm, defaults to \code{max_iter = 50}.
#' @param return.centers should the cluster centers transformed back the space of HPD matrices also be returned?
#' Defaults to \code{return.centers = FALSE}.
#' @param ... additional arguments passed on to \code{\link{pdSpecEst1D}}.
#'
#' @return Depending on the input the function returns a list with five or six components:
#' \describe{
#'   \item{cl.prob }{ an (\eqn{S,K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#'   matrix corresponds to the probabilistic cluster membership assignment of subject \eqn{s} with respect
#'   to cluster \eqn{k}.}
#'   \item{cl.centers.D }{ a list of \code{K} wavelet coefficient pyramids, where each pyramid of wavelet
#'   coefficients is associated to a cluster center.}
#'   \item{cl.centers.M0 }{ a list of \code{K} arrays of coarse-scale midpoints at scale \code{j = 0}, where each
#'   array is associated to a cluster center.}
#'   \item{cl.centers.f }{ only available if \code{return.centers = TRUE}, returning a list of \code{K} \eqn{(d,d,n)}-dimensional arrays,
#'   where each array corresponds to a length \eqn{n} curve of \eqn{(d,d)}-dimensional HPD matrices associated to a cluster center.}
#'   \item{cl.jmax }{ the maximum wavelet scale taken into account in the clustering procedure determined by
#'   the input arguments \code{jmax} and \code{d.jmax}.}
#' }
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#' Phi1 <- array(c(0.5, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Phi2 <- array(c(0.7, 0, 0, 0.4, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#'
#' ## Generate periodogram data for 10 subjects in 2 groups
#' pgram <- function(Phi) pdPgram(rARMA(2^9, 2, Phi, Theta, Sigma)$X)$P
#' P <- array(c(replicate(5, pgram(Phi1)), replicate(5, pgram(Phi2))), dim=c(2,2,2^8,10))
#'
#' cl <- pdSpecClust1D(P, K = 2, metric = "logEuclidean")
#'
#' @seealso \code{\link{pdSpecEst1D}}, \code{\link{pdSpecClust2D}}, \code{\link{pdkMeans}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdSpecClust1D <- function(P, K, jmax, metric = "Riemannian", m = 2, d.jmax = 0.1, eps = c(1e-04, 1e-04),
                          tau = 0.5, max_iter = 50, return.centers = FALSE, ...) {

  ## Set variables
  d <- dim(P)[1]
  S <- dim(P)[4]
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  dots <- list(...)
  order <- (if(is.null(dots$order)) 5 else dots$order)
  periodic <- (if(is.null(dots$periodic)) TRUE else dots$periodic)
  alpha <- (if(is.null(dots$alpha)) 1 else dots$alpha)
  bias.corr <- (if(is.null(dots$bias.corr)) TRUE else dots$bias.corr)
  jmax <- (if(missing(jmax)) log2(dim(P)[3]) - 2 else jmax)
  if(periodic){
    L <- (order - 1) / 2
    L_b <- ceiling(L / 2)
  } else{
    L <- L_b <- 0
  }

  ## Compute denoised wavelet coefficients
  D <- D.est <- list()
  M0 <- array(dim = c(d, d, S))
  M0.est <- array(dim = c(d, d, order, S))
  D.nzero <- matrix(NA, S, jmax + 1)
  for(s in 1:S) {
    D.s <- pdSpecEst1D(P[, , , s], order, metric, alpha, return_val = "D", jmax = jmax,
                       periodic = periodic, return.D = "D.white", bias.corr = bias.corr)
    D[[s]] <- lapply(0:jmax, function(j) D.s$D.white[[j + 1]][, , L_b + 1:2^j, drop = FALSE])
    M0[, , s] <- D.s$M0[, , L + 1]
    M0.est[, , , s] <- D.s$M0
    D.nzero[s, ] <- sapply(0:jmax, function(j) ifelse(j == 0, TRUE, sum(D.s$tree.weights[[j]])/2^j))
    D.est[[s]] <- D.s$D
  }

  ## C++ c-medoids algorithm coarsest midpoints
  clust <- cMeans_C(M0, M0[, , sample(1:S, K)], S, K, m, eps[1], max_iter,
                    ifelse(metric == "Riemannian", "Riemannian", "Euclidean"), matrix(1, S, K))

  ## Set up variables for weighted c-means in wavelet domain
  jmax <- min(jmax, sum(colMeans(D.nzero) > d.jmax) - 1)
  n_jmax <- 2^(jmax + 1) - 1
  DD <- array(aperm(array(sapply(D, function(D) unlist(D[1:(jmax + 1)])),
                          dim = c(d, d, n_jmax, S)), c(1, 2, 4, 3)), dim = c(d, d, n_jmax * S))
  cent <- array(dim = c(d, d, n_jmax * K))
  for(k in 1:K) {
    w <- clust[, k, 1]^m / sum(clust[, k, 1]^m)
    for(i in 1:n_jmax) {
      cent[, , (k - 1) * n_jmax + i] <- apply(sweep(DD[, , (i - 1) * S + 1:S], 3, w, "*"), c(1, 2), sum)
    }
  }
  dist_weights <- (1 - exp(-tau * clust[, , 2])) / (1 + exp(-tau * clust[, , 2]))

  ## C++ weighted c-means algorithm wavelet domain
  clust <- cMeans_C(DD, cent, S, K, m, eps[2], max_iter, "Euclidean", dist_weights)[, , 1]
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  ## Compute wavelet domain centroids
  cent.M0 <- cent.D <- list()
  for(k in 1:K) {
    w <- clust[, k]^m / sum(clust[, k]^m)
    cent.M0[[k]] <- sapply(1:dim(M0.est)[3], function(i) pdMean(M0.est[, , i,], w, ifelse(metric == "Riemannian",
                            "Riemannian", "Euclidean"), grad_desc = TRUE), simplify = "array")
    cent.D[[k]] <- lapply(1:length(D.est[[1]]), function(j) {
      apply(sweep(sapply(D.est, "[[", j, simplify = "array"), 4, w, "*"), c(1, 2, 3), sum)
    })
  }
  names(cent.D) <- names(cent.M0) <- paste0("cluster.center.", 1:K)

  ## Inverse wavelet transforms and return output
  if(isTRUE(return.centers)){
    cent.f <- lapply(1:K, function(k) InvWavTransf1D(cent.D[[k]], cent.M0[[k]], order, periodic = periodic, metric = metric))
    names(cent.f) <- paste0("cluster.center.", 1:K)
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.centers.f = cent.f, cl.jmax = jmax)
  } else {
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.jmax = jmax)
  }
  return(res)
}

#' Intrinsic wavelet HPD time-varying spectral clustering
#'
#' \code{pdSpecClust2D} performs clustering of HPD time-varying spectral matrices corrupted by noise (e.g. HPD time-varying
#' periodograms) by combining wavelet thresholding and fuzzy clustering in the intrinsic wavelet coefficient domain according to
#' the following steps:
#' \enumerate{
#'   \item Transform a collection of noisy HPD time-varying spectral matrices to the intrinsic wavelet domain and denoise the
#'   HPD matrix surfaces by (tree-structured) thresholding of wavelet coefficients with \code{\link{pdSpecEst2D}}.
#'   \item Apply an intrinsic fuzzy c-means algorithm to the coarsest midpoints at scale \code{j = 0} across subjects.
#'   \item Taking into account the fuzzy cluster assignments in the previous step, apply a weighted fuzzy c-means
#'   algorithm to the nonzero thresholded wavelet coefficients across subjects from scale \code{j = 1} up to \code{j = jmax}.
#' }
#' More details can be found in Chapter 3 of \insertCite{C18}{pdSpecEst} and the accompanying vignettes.
#'
#' The input array \code{P} corresponds to a collection of initial noisy HPD time-varying spectral estimates of the
#' \eqn{(d,d)}-dimensional time-varying spectral matrix at \eqn{n_1 \times n_2} time-frequency points, with \eqn{n_1, n_2}
#' dyadic numbers, for \eqn{S} different subjects. These can be e.g., multitaper HPD time-varying periodograms given as output
#' by the function \code{\link{pdPgram2D}}.\cr
#' First, for each subject \eqn{s = 1,\ldots,S}, thresholded wavelet coefficients in the intrinsic wavelet domain are
#' calculated by \code{\link{pdSpecEst2D}}, see the function documentation for additional details on the wavelet thresholding
#' procedure.\cr
#' The maximum wavelet scale taken into account in the clustering procedure is determined by the arguments
#' \code{jmax} and \code{d.jmax}. The maximum scale is set to the minimum of \code{jmax} and the wavelet
#' scale \eqn{j} for which the proportion of nonzero thresholded wavelet coefficients (averaged
#' across subjects) is smaller than \code{d.jmax}.\cr
#' The \eqn{S} subjects are assigned to \eqn{K} different clusters in a probabilistic fashion according to a
#' two-step procedure:
#' \enumerate{
#' \item In the first step, an intrinsic fuzzy c-means algorithm, with fuzziness parameter \eqn{m} is applied to the
#' \eqn{S} coarsest midpoints at scale \code{j = 0} in the subject-specific 2D midpoint pyramids. Note that the distance
#' function in the intrinsic c-means algorithm relies on the chosen metric on the space of HPD matrices.
#' \item In the second step, a weighted fuzzy c-means algorithm based on the Euclidean
#' distance function, also with fuzziness parameter \eqn{m}, is applied to the nonzero thresholded wavelet
#' coefficients of the \eqn{S} different subjects. The tuning parameter \code{tau} controls the weight given
#' to the cluster assignments obtained in the first step of the clustering algorithm.
#' }
#' The function computes the forward and inverse intrinsic 2D AI wavelet transform in the space of HPD matrices equipped with
#' one of the following metrics: (i) the affine-invariant Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}; (ii) the log-Euclidean metric, the Euclidean inner product between matrix logarithms;
#' (iii) the Cholesky metric, the Euclidean inner product between Cholesky decompositions; (iv) the Euclidean metric; or
#' (v) the root-Euclidean metric. The default choice of metric (affine-invariant Riemannian) satisfies several useful properties
#' not shared by the other metrics, see \insertCite{C18}{pdSpecEst} for more details. Note that this comes
#' at the cost of increased computation time in comparison to one of the other metrics. \cr
#' If \code{return.centers = TRUE}, the function also returns the \code{K} HPD time-varying spectral matrices corresponding to
#' the cluster centers based on the given metric by applying the intrinsic inverse 2D AI wavelet transform (
#' \code{\link{InvWavTransf2D}}) to the cluster centers in the wavelet domain.
#'
#' @param P a (\code{d,d,n[1],n[2],S})-dimensional array of HPD matrices, corresponding to a collection of surfaces of
#' \eqn{(d,d)}-dimensional HPD matrices of size \eqn{n_1 \times n_2}, with \eqn{n_1 = 2^{J_1}} and \eqn{n_2 = 2^{J_2}}
#' for some \eqn{J_1,J_2 > 0}, for \eqn{S} different subjects.
#' @param K the number of clusters, a positive integer larger than 1.
#' @param jmax an upper bound on the maximum wavelet scale to be considered in the clustering procedure. If
#' \code{jmax} is not specified, it is set equal to the maximum (i.e., finest) wavelet scale minus 2.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. Additional details are given below.
#' @param m the fuzziness parameter for both fuzzy c-means algorithms. \code{m} should be larger or equal to \eqn{1}.
#' If \eqn{m = 1} the cluster assignments are no longer fuzzy, i.e., the procedure
#' performs hard clustering.
#' @param d.jmax a proportion that is used to determine the maximum wavelet scale to be considered in the clustering
#' procedure. A larger value \code{d.jmax} leads to less wavelet coefficients being taken into account, and therefore
#' lower computational effort in the procedure. If \code{d.jmax} is not specified, by default \code{d.jmax = 0.1}.
#' @param eps an optional vector with two components determining the stopping criterion. The first step in the cluster procedure
#' terminates if the (integrated) intrinsic distance between cluster centers is smaller than \code{eps[1]}.
#' The second step in the cluster procedure terminates if the (integrated) Euclidean distance between cluster centers is smaller
#' than \code{eps[2]}. By default \code{eps = c(1e-04, 1e-04)}.
#' @param tau an optional argument tuning the weight given to the cluster assignments obtained in the first step of
#' the clustering algorithm. If \code{tau} is not specified, by default \code{tau = 0.5}.
#' @param max_iter an optional argument tuning the maximum number of iterations in both the first and second step of the
#' clustering algorithm, defaults to \code{max_iter = 50}.
#' @param return.centers should the cluster centers transformed back the space of HPD matrices also be returned?
#' Defaults to \code{return.centers = FALSE}.
#' @param ... additional arguments passed on to \code{\link{pdSpecEst2D}}.
#'
#' @return Depending on the input the function returns a list with five or six components:
#' \describe{
#'   \item{cl.prob }{ an (\eqn{S,K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#'   matrix corresponds to the probabilistic cluster membership assignment of subject \eqn{s} with respect
#'   to cluster \eqn{k}.}
#'   \item{cl.centers.D }{ a list of \code{K} wavelet coefficient pyramids, where each 2D pyramid of wavelet
#'   coefficients is associated to a cluster center.}
#'   \item{cl.centers.M0 }{ a list of \code{K} arrays of coarse-scale midpoints at scale \code{j = 0}, where each
#'   array is associated to a cluster center.}
#'   \item{cl.centers.f }{ only available if \code{return.centers = TRUE}, returning a list of \code{K} \code{(d,d,n[1],n[2])}-dimensional
#'   arrays, where each array corresponds to an\eqn{n_1 \times n_2}-sized surface of \eqn{(d,d)}-dimensional HPD matrices associated
#'   to a cluster center.}
#'   \item{cl.jmax }{ the maximum wavelet scale taken into account in the clustering procedure determined by
#'   the input arguments \code{jmax} and \code{d.jmax}.}
#' }
#'
#' @examples
#' \dontrun{
#' ## Generate noisy HPD surfaces for 6 subjects in 2 groups
#' n <- c(2^5, 2^5)
#' P <- array(c(rExamples2D(n, example = "tvar", replicates = 3)$P,
#'              rExamples2D(n, example = "tvar", replicates = 3)$P), dim = c(2, 2, n, 6))
#' cl <- pdSpecClust2D(P, K = 2, metric = "logEuclidean")
#' }
#'
#' @seealso \code{\link{pdSpecEst2D}}, \code{\link{WavTransf2D}}, \code{\link{pdDist}}, \code{\link{pdPgram2D}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdSpecClust2D <- function(P, K, jmax, metric = "Riemannian", m = 2, d.jmax = 0.1,
                          eps = c(1e-04, 1e-04), tau = 0.5, max_iter = 50, return.centers = FALSE, ...) {

  ## Set variables
  d <- dim(P)[1]
  S <- dim(P)[5]
  jmax <- (if(missing(jmax)) max(log2(dim(P)[3]), log2(dim(P)[4])) - 2 else jmax)
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))
  dots <- list(...)
  order <- (if(is.null(dots$order)) c(3, 3) else dots$order)
  alpha <- (if(is.null(dots$alpha)) 1 else dots$alpha)
  bias.corr <- (if(is.null(dots$bias.corr)) TRUE else dots$bias.corr)

  ## Compute denoised wavelet coefficients
  D <- D.est <- list()
  M0 <- array(dim = c(d, d, S))
  D.nzero <- D.nzero1 <- matrix(NA, S, jmax + 1)
  for(s in 1:S) {
    D.s <- pdSpecEst2D(P[, , , , s], order = order, metric = metric, alpha = alpha,
                       return_val = "D", jmax = jmax, return.D = "D.white", bias.corr = bias.corr)

    D[[s]] <- D.s$D.white
    M0[, , s] <- D.s$M0[, , 1, 1]
    D.nzero[s, ] <- sapply(0:jmax, function(j) ifelse(j == 0, 1, mean(D.s$tree.weights[[j]])))
    D.est[[s]] <- D.s$D
  }

  ## C++ c-medoids algorithm coarsest midpoints
  clust <- cMeans_C(M0, M0[, , sample(1:S, K)], S, K, m, eps[1], max_iter,
                    ifelse(metric == "Riemannian", "Riemannian", "Euclidean"), matrix(1, S, K))

  ## Set up variables for weighted c-means in wavelet domain
  jmax <- min(jmax, sum(colMeans(D.nzero) > d.jmax) - 1)
  n_jmax <- sum(apply(sapply(D[[1]], dim)[c(3,4), 1:(jmax + 1)], 2, prod))
  DD <- array(aperm(array(sapply(D, function(D) unlist(D[1:(jmax + 1)])),
                          dim = c(d, d, n_jmax, S)), c(1, 2, 4, 3)), dim = c(d, d, n_jmax * S))
  cent <- array(dim = c(d, d, n_jmax * K))
  for(k in 1:K) {
    w <- clust[, k, 1]^m / sum(clust[, k, 1]^m)
    for(i in 1:n_jmax) {
      cent[, , (k - 1) * n_jmax + i] <- apply(sweep(DD[, , (i - 1) * S + 1:S], 3, w, "*"), c(1, 2), sum)
    }
  }
  dist_weights <- (1 - exp(-tau * clust[, , 2])) / (1 + exp(-tau * clust[, , 2]))

  ## C++ weighted c-means algorithm wavelet domain
  clust <- cMeans_C(DD, cent, S, K, m, eps[2], max_iter, "Euclidean", dist_weights)[, , 1]
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  ## Compute wavelet domain centroids
  cent.M0 <- cent.D <- list()
  for(k in 1:K) {
    w <- clust[, k]^m / sum(clust[, k]^m)
    cent.M0[[k]] <- array(pdMean(M0, w, ifelse(metric == "Riemannian", "Riemannian", "Euclidean"),
                                 grad_desc = TRUE), dim = c(d, d, 1, 1))
    cent.D[[k]] <- lapply(1:length(D.est[[1]]), function(j) {
      apply(sweep(sapply(D.est, "[[", j, simplify = "array"), 5, w, "*"), 1:4, sum)
    })
  }
  names(cent.D) <- names(cent.M0) <- paste0("cluster.center.", 1:K)

  ## Inverse wavelet transforms and return output
  if(isTRUE(return.centers)){
    cent.f <- lapply(1:K, function(k) InvWavTransf2D(cent.D[[k]], cent.M0[[k]], order = order, metric = metric))
    names(cent.f) <- paste0("cluster.center.", 1:K)
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.centers.f = cent.f, cl.jmax = jmax)
  } else {
    res <- list(cl.prob = clust, cl.centers.D = cent.D, cl.centers.M0 = cent.M0, cl.jmax = jmax)
  }
  return(res)
}



