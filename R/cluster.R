#' Wavelet-based clustering of multivariate spectra.
#'
#' \code{pdSpecClust} performs clustering of multivariate spectral matrices via a two-step fuzzy
#' clustering algorithm in the manifold wavelet domain. The clustering procedure is described in
#' detail in (Chau and von Sachs, 2017).
#'
#' The input array \code{P} contains initial noisy Hermitian PD spectral estimates of the
#' (\eqn{d \times d})-dimensional spectral matrices at \eqn{m} different frequencies for \eqn{S}
#' different subjects. Here \eqn{m = 2^J} for some \eqn{J > 0}. The initial spectral estimates can
#' be e.g. the averaged periodograms given as output in \code{\link{pdPgram}}. \cr
#' For each subject \eqn{s}, the thresholded wavelet coefficients in the manifold wavelet domain are
#' calculated by the function \code{\link{pdSpecEst}}. Instead of providing the noisy input data \code{P},
#' it is possible to directly provide the thresholded wavelet coefficients via the argument \code{D.hat}.
#' \code{D.hat} should be a list of \eqn{S} elements, where each element corresponds to an object obtained
#' from the function \code{\link{pdSpecEst}}.\cr
#' The maximum wavelet scale taken into account in the clustering procedure is determined by the arguments
#' \code{jmax} and \code{d.jmax}. The maximum scale is set to the minimum of \code{jmax} and the wavelet
#' scale \eqn{j} for which the proportion of nonzero components of thresholded wavelet coefficients (averaged
#' across subjects) is smaller than \code{d.jmax}.\cr
#' The \eqn{S} subjects are assigned to \eqn{K} different fuzzy clusters according to a two-step procedure.
#' In the first step, a fuzzy c-medoids algorithm based on the Riemannian distance function as in
#' \code{\link{pdDist}}, with fuzziness parameter \eqn{m}, is applied to the \eqn{S} coarsest-scale midpoint
#' vectors, i.e. scale \code{j = 1}. \cr
#' In the second step, a weighted fuzzy c-means algorithm based on the Euclidean
#' distance function, also with fuzziness parameter \eqn{m}, is applied to the nonzero thresholded wavelet
#' coefficients for the \eqn{S} different subjects. The tuning parameter \code{tau} controls the weight given
#' to the cluster assignments obtained in the first step of the clustering algorithm. For more details on the
#' procedure we refer to (Chau and von Sachs, 2017).
#'
#' @param P a (\eqn{d,d,m,S})-dimensional array, with \eqn{m} a dyadic number.
#' @param D.hat a list with \eqn{S} elements, where each element is an object given as output by the function
#' \code{\link{pdSpecEst}}. If both arguments \code{P} and \code{D.hat} are provided, only the values in \code{D.hat}
#' are taken into account.
#' @param K the number of clusters to be considered.
#' @param m the fuzziness parameter for both the fuzzy c-medoids and the weighted fuzzy c-means algorithm. \code{m}
#' should be larger or equal to \eqn{1}. If \eqn{m = 1} the cluster assignments are no longer fuzzy (i.e. the procedure
#' performs hard clustering).
#' @param jmax an upper bound on the maximum wavelet scale to be considered in the clustering procedure. If
#' \code{jmax} is not specified, it is set equal to the largest (i.e. finest) possible wavelet scale.
#' @param d.jmax a proportion that is used to determine the maximum wavelet scale to be considered in the clustering
#' procedure. A larger value \code{d.jmax} leads to less wavelet coefficients being taken into account, and therefore
#' lower computational effort in the procedure. If \code{d.jmax} is not specified, by default \code{d.jmax = 0.1}.
#' @param eps an optional vector with two elements, \code{eps[1]} determines the termination criterion in the fuzzy c-medoids
#' algorithm (i.e. first clustering step), and \code{eps[2]} determines the termination criterion in the weighted fuzzy c-means algorithm
#' (i.e. second clustering step). If \code{eps} is not specified, by default \code{eps = c(1E-4, 1E-4)}.
#' @param tau an optional argument tuning the weight given to the cluster assignments obtained in the first step of
#' the clustering algorithm. If \code{tau} is not specified, by default \code{tau = 0.5}.
#' @param return.D an optional argument specifying whether to return also the list of coarsest midpoints and wavelet coefficients
#' (i.e. the features) for each individual subject, by default \code{return.D = FALSE}.
#' @param ... additional arguments passed on to \code{\link{pdSpecEst}}. These arguments are only used if
#' \code{is.null(D.hat)}, otherwise the function \code{\link{pdSpecEst}} is not called.
#'
#' @return The function returns an (\eqn{S, K})-dimensional matrix, where the value at position (\eqn{s,k}) in the
#' matrix corresponds to the probabilistic cluster membership assignment of subject \eqn{s} with respect to
#' cluster \eqn{k}. If \code{isTRUE(return.D)} also returns the list of coarsest midpoints and wavelet coefficients
#' of each subject used in the clustering procedure.
#'
#' @examples
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#'
#' Phi1 <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Phi2 <- array(c(0.7, 0, 0, 0.5, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#'
#' ## Generate periodogram data for 10 subjects
#' pgram <- function(Phi) pdPgram(rARMA(2^9, 2, Phi, Theta, Sigma)$X)$P
#' P <- array(c(replicate(5, pgram(Phi1)), replicate(5, pgram(Phi2))), dim=c(2,2,2^7,10))
#'
#' cl <- pdSpecClust(P, K=2)
#'
#' @seealso \code{\link{pdPgram}}, \code{\link{pdSpecEst}}, \code{\link{pdDist}}
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
pdSpecClust <- function(P, D.hat = NULL, K, m = 2, jmax, d.jmax = 0.1, eps = c(1e-04, 1e-04), tau = 0.5, return.D = FALSE, ...) {

  ## missing arguments
  if (missing(P)) {
    P <- NULL
  }
  dots <- list(...)
  order <- dots$order
  if (is.null(order)) {
    order <- 5
  }
  alpha <- dots$alpha
  if (is.null(alpha)) {
    alpha <- 0.75
  }
  lam <- dots$lam

  if (is.null(D.hat)) {

    ## Define variables
    J <- log2(dim(P)[3])
    if (!isTRUE(all.equal(as.integer(J), J))) {
      stop(paste0("Input length is non-dyadic, please change length ", dim(P)[3],
                   " to dyadic number."))
    }
    if (!(order %in% c(1, 3, 5, 7, 9))) {
      warning("Refinement order should be an odd integer between 1 and 9, by default set to 5")
      order <- 5
    }
    if (missing(jmax)) {
      jmax <- J - 1
    }
    dim <- dim(P)[1]
    S <- dim(P)[4]
    d.nzero <- matrix(1, ncol = jmax, nrow = S)
    D <- list()
    for (s in 1:S) {
      D.hat <- pdSpecEst(P[, , , s], lam, order, return = "D", alpha)
      D[[s]] <- D.hat$D
      d.nzero[s, ] <- sapply(1:jmax, function(j) sum(as.logical(D.hat$components$thresholded[[j]]))/(dim^2 *  2^j))
    }
    jmax <- min(jmax, sum(colMeans(d.nzero) > d.jmax))
  } else {

    ## Define variables
    J <- length(D.hat[[1]]$D)
    dim <- dim(D.hat[[1]]$D[[1]])[1]
    S <- length(D.hat)
    if (missing(jmax))
      jmax <- J - 1
    D <- list()
    d.nzero <- matrix(1, ncol = jmax, nrow = S)
    for (s in 1:S) {
      D[[s]] <- D.hat[[s]]$D
      d.nzero[s, ] <- sapply(1:jmax, function(j) sum(as.logical(D.hat[[s]]$components$thresholded[[j]]))/(dim^2 * 2^j))
    }
    jmax <- min(jmax, sum(colMeans(d.nzero) > d.jmax))
  }

  D0 <- D

  ## c-medoids algorithm
  M <- sapply(1:S, function(s) D[[s]][[1]], simplify = "array")
  cent <- M[, , , sample(1:S, K)]
  stopit <- F
  i <- 0
  distM <- function(M1, M2) {
    RiemmDist(M1[, , 1], M2[, , 1])^2 + RiemmDist(M1[, , 2], M2[, , 2])^2
  }
  while ((!stopit) & (i < 50)) {
    dist <- t(sapply(1:S, function(s) t(sapply(1:K, function(k) distM(M[, , , s], cent[, , , k])))))
    mu <- function(s) {
      if (!any(dist[s, ] < .Machine$double.eps)) {
        sapply(1:K, function(k) 1/sum((dist[s, k]/dist[s, ])^(1/(m - 1))))
      } else {
        as.numeric(dist[s, ] < .Machine$double.eps)/sum(dist[s, ] < .Machine$double.eps)
      }
    }
    clust <- t(sapply(1:S, mu))
    cent1 <- array(dim = c(dim, dim, 2, K))
    for (k in 1:K) {
      w <- clust[, k]^m/sum(clust[, k]^m)
      cent1[, , , k] <- sapply(1:2, function(i) KarchMean(M[, , i, ], w), simplify = "array")
    }
    stopit <- ifelse(isTRUE(sum(sapply(1:K, function(k) distM(cent[, , , k],
                                        cent1[, , , k]))) > eps[1]), F, T)
    cent <- cent1
    i <- i + 1
  }

  ## Weighted c-means algorithm
  D <- lapply(1:jmax, function(j) sapply(1:S, function(s) D[[s]][[j + 1]], simplify = "array"))
  cent <- lapply(1:jmax, function(j) array(dim = c(dim, dim, 2^j, K)))
  for (k in 1:K) {
    w <- clust[, k]^m/sum(clust[, k]^m)
    for (j in 1:jmax) {
      cent[[j]][, , , k] <- sapply(1:2^j, function(i) apply(array(rep(w, each = dim^2),
                                    dim = c(dim, dim, S)) * D[[j]][, , i, ], c(1, 2), sum), simplify = "array")
    }
  }
  stopit <- F
  i <- 0
  dist0 <- dist
  clust0 <- clust
  distj <- function(D1, D2, j) {
    sum(sapply(1:2^j, function(i) 2^(jmax - j) * NormF(D1[, , i] - D2[, , i])^2))
  }
  while ((!stopit) & (i < 50)) {
    dist <- t(sapply(1:S, function(s) t(sapply(1:K, function(k) sum(sapply(1:jmax,
                                                function(j) (1 - exp(-tau * dist0[s, k]))/(1 + exp(-tau * dist0[s, k])) *
                                                  distj(D[[j]][, , , s], cent[[j]][, , , k], j)))))))
    mu <- function(s) {
      if (!any(dist[s, ] < .Machine$double.eps)) {
        sapply(1:K, function(k) 1/sum((dist[s, k]/dist[s, ])^(1/(m - 1))))
      } else if (all(dist[s, ] < .Machine$double.eps)) {
        clust0[s, ]
      } else {
        as.numeric(dist[s, ] < .Machine$double.eps)/sum(dist[s, ] < .Machine$double.eps)
      }
    }
    clust <- t(sapply(1:S, mu))
    cent1 <- lapply(1:jmax, function(j) array(dim = c(dim, dim, 2^j, K)))
    for (k in 1:K) {
      w <- clust[, k]^m/sum(clust[, k]^m)
      for (j in 1:jmax) {
        cent1[[j]][, , , k] <- sapply(1:2^j, function(i) apply(array(rep(w,  each = dim^2),
                                       dim = c(dim, dim, S)) * D[[j]][, , i, ], c(1, 2), sum),
                                         simplify = "array")
      }
    }
    stopit <- ifelse(isTRUE(sum(sapply(1:K, function(k) sum(sapply(1:jmax, function(j)
                                        distj(cent[[j]][, , , k], cent1[[j]][, , , k], j))))) > eps[2]), F, T)
    cent <- cent1
    i <- i + 1
  }
  rownames(clust) <- paste0("Subject", 1:S)
  colnames(clust) <- paste0("Cluster", 1:K)

  if(!return.D){
    cl <- clust
  } else{
    cl <- list(clust = clust, D = D0)
  }
  return(cl)
}
