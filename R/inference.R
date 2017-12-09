#' Intrinsic depth-based bootstrap confidence regions
#'
#' \code{pdConfInt1D} constructs depth-based parametric bootstrap confidence regions as explained in
#' (Chau and von Sachs, 2017a) for a wavelet-based spectral matrix estimator obtained with \code{\link{pdSpecEst1D}}
#' based on intrinsic manifold data depths  as in \code{\link{pdDepth}}.
#'
#' The parametric bootstrap procedure exploits the data generating process of a stationary time series
#' via its Cramer representation, and is equivalent to e.g. (Fiecas and Ombao, 2016) among others.
#' Given a consistent wavelet-based spectral estimator \code{f} obtained with \code{\link{pdSpecEst1D}} corresponding
#' to a sequence of \eqn{m} HPD matrices along the frequency range \eqn{(0, \pi]}, where w.l.o.g. \eqn{m} is assumed
#' to be even, the bootstrap confidence regions are constructed as:
#' \enumerate{
#'    \item Generate a bootstrapped time series trace \code{X_b} of length \eqn{n} via its Cramer representation using
#'    complex normal random variates and transfer functions given by the Hermitian square root of the estimate \code{f}.
#'    \item Compute a bootstrapped wavelet-based spectral estimate \code{f_b} based on \code{X_b}.
#'    \item Repeat 1. and 2. many times. Construct \eqn{100(1-\alpha)\%}-pointwise or -simultaneous confidence sets
#'    by taking the \eqn{100(1-\alpha)\%}-most central depth-based region in the cloud of bootstrapped spectral estimates
#'    based on one of the intrinsic manifold data depths in \code{\link{pdDepth}}.
#' }
#' The depth-based confidence balls (maximum/minimum depths and radii) are given by the component \code{depth.CI}. In particular,
#' a given curve of HPD matrices is covered by the confidence ball if its manifold data depth with respect to the cloud of
#' bootstrapped spectral estimates is above the minimum depth in the given confidence ball. \cr
#' If \code{return.f = TRUE}, the function also returns the \eqn{100(1-\alpha)\%}-most central bootstrapped spectral estimates,
#' by default \code{return.f = FALSE} to save memory use of the returned object. \cr
#' If, in addition, we supply a target spectrum \code{f.0} (i.e. a numerical
#' array with the same dimensions as \code{f}), the function checks whether the target \code{f.0} is covered by the
#' constructed confidence regions or not by computing the depth of \code{f.0} with respect to the cloud of bootstrapped spectral
#' estimates.
#'
#' @param f a (\eqn{d,d,m})-dimensional array corresponding to a wavelet-denoised HPD (\eqn{d,d})-dimensional
#' spectral estimate at \code{m} different frequencies given as output with \code{\link{pdSpecEst1D}}.
#' @param alpha a numerical vector of quantiles (between 0 and 1) determining the \eqn{100(1-\alpha)\%}-confidence levels,
#' defaults to \code{alpha = 0.05}.
#' @param ci.region a 2-dimensional numeric vector \code{c(min.ci, max.ci)} specifying the domain of the simultaneous confidence region.
#' The frequency domain \eqn{[0, \pi]} is normalized to a unit interval, e.g. \code{ci.region = c(0.5, 1)} constructs a simultaneous
#' confidence region over the second half of the frequency domain. \code{ci.region} can also be a \eqn{(2,L)}-dimensional matrix, where each
#' column of the matrix specificies the domain of an individual simultaneous confidence region. If \code{ci.region} is not specified, it defaults
#' to the unit interval, i.e. the entire frequency domain.
#' @param boot.samples number of bootstrap spectral estimates, defaults to \code{boot.samples = 1e3}.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"logEuclidean"},
#' but this can also be one of: \code{"Riemannian"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. This argument is passed on to the spectral estimator in \code{\link{pdSpecEst1D}} and the depth calculation in
#' \code{\link{pdDepth}}. The default choice is the Log-Euclidean metric as the computational effort is significantly reduced
#' in comparison to e.g. the Riemannian metric.
#' @param depth the data depth measure, one of \code{'zonoid'}, \code{'gdd'}, or \code{'spatial'} corresponding to
#' the manifold zonoid depth, geodesic distance depth, and manifold spatial depth respectively.
#' @param return.f a logical value, if \code{return.f = TRUE} the function also returns the \eqn{100(1-\alpha)\%}-most central bootstrapped
#' spectral estimates. By default \code{return.f = FALSE} to save memory use of the returned object.
#' @param f.0 optional target spectrum in the same format as \code{f}. The function computes the depth of the target spectrum with respect to
#' the cloud of bootstrap spectral estimates and checks whether the target \code{f.0} is covered by the constructed depth-based confidence region.
#' @param ... additional arguments passed on to the functions \code{\link{pdPgram}} and \code{\link{pdSpecEst1D}} used internally.
#'
#' @return The function returns a list with the following components:
#' \item{depth.CI }{ a matrix or a list of matrices, with each matrix giving the depth-based confidence ball (maximum depth, minimum depth and radius)
#' at an individual region in the frequency domain as specified by the argument \code{ci.region} and at the various \eqn{100(1-\alpha)\%}-confidence
#'  levels specified by the argument \code{alpha}. }
#' \item{f.CI }{ If \code{return.f = TRUE} returns the \eqn{100(1-\alpha)\%}-most central bootstrapped spectral estimates at each of the individual regions
#' in the frequency domain specified by the argument \code{ci.region} at the level corresponding to the first argument in the input vector \eqn{\alpha}.
#' If \code{return.f = FALSE} returns \code{NULL}.}
#' If a target spectrum \code{f.0} is supplied, i.e. \code{!is.null(f.0)}, the returned list includes the additional components:
#' \item{cover.f }{ a list of logical vectors checking whether the constructed confidence balls cover the target spectrum \code{f.0}. Each vector corresponds
#' to an individual region in the frequency domain as specified by the argument \code{ci.region} at the various \eqn{100(1-\alpha)\%}-confidence levels
#' specified by the argument \code{alpha}.}
#' \item{depth.f }{ a numeric vector of depth values of the target spectrum \code{f.0} with respect to the cloud of bootstrapped spectral estimates at each of
#' the individual regions in the frequency domain as specified by the argument \code{ci.region}.}
#'
#' @examples
#' \dontrun{
#' set.seed(123)
#' example <- rExamples(2^8, example = "gaussian")
#' f.hat <- pdSpecEst1D(example$per, metric = "logEuclidean")
#' boot.ci <- pdConfInt1D(f.hat$f, alpha = c(0.1, 0.05, 0.01), ci.region = c(0.45, 0.55),
#'                        boot.samples = 1E3, f.0 = example$f)
#' boot.ci
#' }
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Chau, J., Ombao, H., and von Sachs, R. (2017b). \emph{Data depth and rank-based
#' tests for covariance and spectral density matrices}. Available at \url{http://arxiv.org/abs/1706.08289}.
#' @references Fiecas, M., Ombao, H. (2016). \emph{Modeling the evolution of dynamic brain processes during
#' an associative learning experiment}. Journal of the American Statistical Association, 111(516), 1440-1453.
#'
#' @export
pdConfInt1D <- function(f, alpha = 0.05, ci.region, boot.samples = 1000, metric = "logEuclidean",
                      depth = "spatial", return.f = F, f.0 = NULL, ...) {

  ## Set variables
  J = log2(dim(f))[3]
  n = 2 * dim(f)[3]
  d = dim(f)[1]
  depth_vec = match.arg(depth, c("gdd", "zonoid", "spatial"), several.ok = T)
  metric = match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "rootEuclidean", "Euclidean"))

  if(!isTRUE(all(alpha > 0 & alpha < 1))) {
    stop(paste0("alpha = ", alpha, " should be a numeric vector (of quantiles) between 0 and 1"))
  }
  if (!isTRUE(all.equal(as.integer(J), J))) {
    stop(paste0("Input length is non-dyadic, please change length ", dim(f)[3],
                " to dyadic number."))
  }

  dots = list(...)
  B = (if(is.null(dots$B)) d else dots$B)
  method = (if(is.null(dots$method)) "multitaper" else dots$method)
  bias.corr = (if(is.null(dots$bias.corr)) F else dots$bias.corr)
  nw = (if(is.null(dots$nw)) pi else dots$nw)
  policy = (if(is.null(dots$policy)) "universal" else dots$policy)
  order = (if(is.null(dots$order)) 5 else dots$order)
  lam = (if(is.null(dots$lam)) 1 else dots$lam)
  periodic = (if(is.null(dots$periodic)) T else dots$periodic)
  J.out = (if(is.null(dots$J.out)) J else dots$J.out)
  return.f = (if(is.null(dots$return.f)) F else dots$return.f)
  jmax = min((if(is.null(dots$jmax)) J - 3 else dots$jmax), J.out - 1)
  jmax.cv = min(if(is.null(dots$jmax.cv)) J - 3 else dots$jmax.cv, J.out - 1)
  progress = (if(is.null(dots$progress)) T else dots$progress)

  if(missing(ci.region)) {
    ci.region <- list(1:2^J.out)
    warning("'ci.region' is missing, by default it is set to the entire domain")
  } else if(is.vector(ci.region)){
    ci.region <- list(max(round(ci.region[1] * 2^J.out), 1):min(round(ci.region[2] * 2^J.out), 2^J.out))
  } else if(isTRUE(is.matrix(ci.region) & nrow(ci.region) == 2)){
    ci.region <- lapply(1:ncol(ci.region), function(L) max(round(ci.region[1, L] * 2^J.out),
                          1):min(round(ci.region[2, L] * 2^J.out), 2^J.out))
  } else {
    stop("'ci.region' does not have the correct format, see the documentation for more details.")
  }

  if(metric == "logEuclidean"){
    f.hat.J <- sapply(1:2^J, function(i) Logm(diag(d), f[, , i]), simplify = "array")
    if(!is.null(f.0)){
      f.0.J <- sapply(1:2^J, function(i) Logm(diag(d), f.0[, , i]), simplify = "array")
    }
  } else if(metric == "Cholesky"){
    f.hat.J <- sapply(1:2^J, function(i) Chol(f[, , i]), simplify = "array")
    if(!is.null(f.0)){
      f.0.J <- sapply(1:2^J, function(i) Chol(f.0[, , i]), simplify = "array")
    }
  } else if(metric == "rootEuclidean"){
    f.hat.J <- sapply(1:2^J, function(i) Sqrt(f[, , i]), simplify = "array")
    if(!is.null(f.0)){
      f.0.J <- sapply(1:2^J, function(i) Sqrt(f.0[, , i]), simplify = "array")
    }
  }

  if(J.out < J){
    for (j in J:(J.out + 1)) {
      if(!(metric == "Riemannian")){
        f.hat.J <- sapply(1:(dim(f.hat.J)[3]/2), function(i) 0.5 * (f.hat.J[, , 2 * i - 1] +  f.hat.J[, , 2 * i]),
                       simplify = "array")
      } else {
        f.hat.J <- sapply(1:(dim(f.hat.J)[3]/2), function(i) Mid(f.hat.J[, , 2 * i - 1], f.hat.J[, , 2 * i]),
                       simplify = "array")
      }
      if(!is.null(f.0)){
        if(!(metric == "Riemannian")){
          f.0.J <- sapply(1:(dim(f.0.J)[3]/2), function(i) 0.5 * (f.0.J[, , 2 * i - 1] +  f.0.J[, , 2 * i]),
                            simplify = "array")
        } else {
          f.0.J <- sapply(1:(dim(f.0.J)[3]/2), function(i) Mid(f.0.J[, , 2 * i - 1], f.0.J[, , 2 * i]),
                            simplify = "array")
        }
      }
    }
  }

  if(metric == "logEuclidean"){
    f.hat.J <- array(apply(f.hat.J, 3, function(f) Expm(diag(d), f)), dim = dim(f.hat.J))
    if(!is.null(f.0)){
      f.0.J <- array(apply(f.0.J, 3, function(f) Expm(diag(d), f)), dim = dim(f.0.J))
    }
  } else if(metric == "Cholesky"){
    f.hat.J <- array(apply(f.hat.J, 3, function(f) Chol(f, inverse = T, bias.corr = T)), dim = dim(f.hat.J))
    if(!is.null(f.0)){
      f.0.J <- array(apply(f.0.J, 3, function(f) Chol(f, inverse = T, bias.corr = T)), dim = dim(f.0.J))
    }
  } else if(metric == "rootEuclidean"){
    f.hat.J <- array(apply(f.hat.J, 3, function(f) t(Conj(f)) %*% f), dim = dim(f.hat.J))
    if(!is.null(f.0)){
      f.0.J <- array(apply(f.0.J, 3, function(f) t(Conj(f)) %*% f), dim = dim(f.0.J))
    }
  }

  ## Generate bootstrap samples
  f.sqrt <- sapply(1:(n/2), function(i) Sqrt(f[, , i]), simplify = "array")
  f.B <- array(dim = c(d, d, 2^J.out, boot.samples))
  if(progress){
    cat("1. Generating bootstrap samples...")
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for(b in 1:boot.samples){
    ## Generate time series via Cramer
    chi <- matrix(nrow = d, ncol = n)
    chi[, 1:(n/2 - 1)] <- replicate(n/2 - 1, complex(d, rnorm(d, sd = sqrt(1/2)), rnorm(d, sd = sqrt(1/2))))
    chi[, c(n/2, n)] <- replicate(2, rnorm(d))
    chi[, (n/2 + 1):(n - 1)] <- Conj(chi[, 1:(n/2 - 1)])

    f.sqrt1 <- array(c(f.sqrt, Conj(f.sqrt[, , (n/2):1])), dim = c(d, d, n))
    ts <- sqrt(2 * pi) / sqrt(n) * mvfft(t(sapply(1:n, function(i) f.sqrt1[, , i] %*% chi[, i])), inverse = T)

    ## Compute b-th pre-smoothed periodogram
    per.b <- pdPgram(ts, B = B, method = method, bias.corr = bias.corr, nw = nw)$P

    ## Compute b-th spectral estimator
    coeffs <- pdSpecEst1D(per.b, order = order, policy = policy, alpha = lam, metric = metric,
                          periodic = periodic, jmax = jmax, jmax.cv = jmax.cv, return = "coeff")
    f.B[, , , b] <- InvWavTransf1D(coeffs$D, coeffs$M0, order = order, periodic = periodic,
                                   jmax = J.out, metric = metric)

    if(progress) utils::setTxtProgressBar(pb, round(100 * b / boot.samples))
  }
  if(progress) close(pb)

  ## Compute data depths and confidence regions
  depth.CI <- lapply(1:length(depth_vec), function(j) array(dim = c(length(alpha), 3, length(ci.region))))
  f.CI <- lapply(1:length(depth_vec), function(j) lapply(1:length(ci.region), function(I) NA))
  cover.f0 <- lapply(1:length(depth_vec), function(j) matrix(nrow = length(alpha), ncol = length(ci.region)))
  depth.f0 <- lapply(1:length(depth_vec), function(j) numeric(length(ci.region)))
  names(depth.CI) <- names(f.CI) <- names(cover.f0) <- names(depth.f0) <- depth_vec

  if(progress){
    cat("2. Computing integrated depth values and constructing confidence regions...")
    pb <- utils::txtProgressBar(1, 100, style = 3)
  }
  for(d_i in 1:length(depth_vec)){
    dimnames(depth.CI[[d_i]]) <- list(sapply(1:length(alpha), function(i) paste0(100 * (1 - alpha[i]), "%-CI")),
                                                            c('max-depth', 'min-depth', 'radius'),
                                                            sapply(1:length(ci.region), function(i) paste0("ci.region.", i)))
    dimnames(cover.f0[[d_i]]) <- dimnames(depth.CI[[d_i]])[c(1, 3)]
    names(depth.f0[[d_i]]) <- dimnames(depth.CI[[d_i]])[[3]]
    names(f.CI[[d_i]]) <- dimnames(depth.CI[[d_i]])[[3]]
    dd <- matrix(nrow = length(ci.region), ncol = boot.samples)

    for(I in 1:length(ci.region)){

      ## Compute data depths
      dd[I, ] <- pdDepth(X = f.B[, , ci.region[[I]], ], method = depth_vec[d_i], metric = metric)

      ## Construct confidence regions
      depth.CI[[d_i]][, 1:2, I] <- t(sapply(1:length(alpha), function(i) c(1, unname(quantile(dd[I, ], alpha[i])))))
      if(return.f){
        f.CI[[d_i]][[I]] <- f.B[, , , which(dd[I, ] > depth.CI[[d_i]][1, 2, I])]
      }
      depth.CI[[d_i]][, 3, I] <- sapply(1:length(alpha), function(j) mean(sapply(ci.region[[I]], function(i)
              pdDist(f.hat.J[, , i], f.B[, , i, which.min(abs(dd[I, ] - depth.CI[[d_i]][j, 2, I]))], method = metric))))

      if(!is.null(f.0)){
        depth.f0[[d_i]][I] <- pdDepth(f.0.J[, , ci.region[[I]], drop = F], f.B[, , ci.region[[I]], , drop = F],
                                      method = depth_vec[d_i], metric = metric)
        cover.f0[[d_i]][, I] <- (depth.f0[[d_i]][I] > depth.CI[[d_i]][, 2, I])
      }
    }
    if(progress){
      utils::setTxtProgressBar(pb, round(100 * (length(ci.region) * (d_i - 1) + I) /
                                           (length(depth_vec) * length(ci.region))))
    }
  }

  if(!is.null(f.0)){
    res <- list(depth.CI = depth.CI, f.CI = (if(return.f) f.CI else NULL),
                cover.f = cover.f0, depth.f = depth.f0)
  } else{
    res <- list(depth.CI = depth.CI, f.CI = (if(return.f) f.CI else NULL))
  }

  return(res)
}


