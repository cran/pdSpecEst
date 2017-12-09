#' Rank-based hypothesis tests for HPD matrices
#'
#' \code{pdRankTests} performs generalized rank-based hypothesis testing in the metric space of HPD matrices equipped
#' with the Riemannian or Log-Euclidean metric for samples of Hermitian PD matrices or samples of sequences (curves)
#' of Hermitian PD matrices as described in (Chau, Ombao, and von Sachs, 2017b).
#'
#' For samples of \eqn{(d \times d)}-dimensional Hermitian PD matrices with pooled sample size \eqn{S}, the argument
#' \code{data} is a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices, where the individual samples are
#' combined along the third array dimension. For samples of sequences of \eqn{(d \times d)}-dimensional Hermitian PD
#' matrices with pooled sample size \eqn{S}, the argument \code{data} is a \eqn{(d,d,n,S)}-dimensional array of sequences
#' of Hermitian PD matrices, where the individual samples are combined along the fourth array dimension. The argument
#' \code{sample.sizes} specifies the sizes of the individual samples so that \code{sum(sample.sizes)} is equal to \code{S}. \cr
#' The available generalized rank-based testing procedures (specified by the argument \code{test}) are:
#' \describe{
#' \item{\code{"rank.sum"}}{Manifold Wilcoxon rank-sum test to test for homogeneity of distributions of two independent
#' samples of Hermitian PD matrices or samples of sequences of Hermitian PD matrices. The usual univariate ranks are replaced by data depth
#' induced ranks via \code{\link{pdDepth}}.}
#' \item{\code{"krusk.wall"}}{Manifold Kruskal-Wallis test to test for homogeneity of distributions of more than two independent
#' samples of Hermitian PD matrices or samples of sequences of Hermitian PD matrices. The usual univariate ranks are replaced by data depth
#' induced ranks via \code{\link{pdDepth}}.}
#' \item{\code{"signed.rank"}}{Manifold signed-rank test to test for homogeneity of distributions of independent paired or matched samples
#' of Hermitian PD matrices. The manifold signed-rank test is \emph{not} based on data depth induced ranks, but on a specific difference score on the Riemannian
#' manifold of Hermitian PD matrices.}
#' \item{\code{"bartels"}}{Manifold Bartels-von Neumann test to test for randomness (i.e. exchangeability) within a single independent sample of
#' Hermitian PD matrices or a sample of sequences of Hermitian PD matrices. The usual univariate ranks are replaced by data depth induced
#' ranks via \code{\link{pdDepth}}.}
#' }
#' The function computes the generalized rank-based test statistics in the \emph{complete} metric space of HPD matrices equipped with one of the following metrics:
#' (i) Riemannian metric (default) as in (Bhatia, 2009, Chapter 6), or (ii) Log-Euclidean metric, the Euclidean inner product between matrix logarithms.
#' The default Riemannian metric is invariant under congruence transformation by any invertible matrix, whereas the Log-Euclidean metric is only
#' invariant under congruence transformation by unitary matrices, see (Chau, Ombao and von Sachs 2017b) for more details.
#'
#' @note The manifold signed-rank test also provides a valid test for equivalence of spectral matrices of two multivariate stationary time
#' series based on the Hermitian PD periodogram matrices obtained via \code{\link{pdPgram}}, see (Chau, Ombao, and von Sachs, 2017b) for the details.
#'
#' @param data either a \eqn{(d,d,S)}-dimensional array corresponding to an array of pooled individual samples of Hermitian PD matrices, or a
#' \eqn{(d,d,n,S)}-dimensional array corresponding to an array of pooled individual samples of sequences of Hermitian PD matrices.
#' @param sample.sizes a numeric vector corresponding to the individual sample sizes in the argument \code{data}, such that \code{sum(sample.sizes)} is
#' equal to \code{S}. Not required for tests \code{"signed-rank"} and \code{"bartels"}, as the sample sizes are automatically determined from \code{data}.
#' @param test rank-based hypothesis testing procedure, one of \code{"rank.sum"}, \code{"krusk.wall"}, \code{"signed.rank"}, \code{"bartels"} explained
#' in the Details section below.
#' @param depth data depth measure used in the rank-based tests, one of \code{"gdd"}, \code{"zonoid"}, or \code{"spatial"} corresponding to the
#' geodesic distance depth, manifold zonoid depth and manifold spatial depth respectively. Defaults to \code{"gdd"}. Not required for test
#' \code{"signed.rank"}.
#' @param metric the metric that the space of HPD matrices is equipped with, either \code{"Riemannian"} or \code{"logEuclidean"}. Defaults to
#' \code{"Riemannian"}.
#'
#' @return The function returns a list with five components:
#' \item{test }{name of the rank-based test}
#' \item{p.value }{p-value of the test}
#' \item{statistic }{computed test statistic}
#' \item{null.distr }{distribution of the test statistic under the null hypothesis}
#' \item{depth.values }{computed data depth values (if available)}
#'
#' @examples
#' ## null hypothesis is true
#' data <- replicate(100, Expm(diag(2), H.coeff(rnorm(4), inverse = TRUE)))
#' pdRankTests(data, sample.sizes = c(50, 50), test = "rank.sum") ## homogeneity 2 samples
#' pdRankTests(data, sample.sizes = rep(25, 4), test = "krusk.wall") ## homogeneity 4 samples
#' pdRankTests(data, test = "bartels") ## randomness
#'
#' ## null hypothesis is false
#' data1 <- array(c(data, replicate(50, Expm(diag(2), H.coeff(0.5 * rnorm(4), inverse = TRUE)))),
#'                  dim = c(2,2,150))
#' pdRankTests(data1, sample.sizes = c(100, 50), test = "rank.sum")
#' pdRankTests(data1, sample.sizes = rep(50, 3), test = "krusk.wall")
#' pdRankTests(data1, test = "bartels")
#'
#' ## signed-rank test for equivalence of spectra
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' pgram <- function(Sigma) pdPgram(rARMA(2^9, 2, Phi, Theta, Sigma)$X)$P
#'
#' ## null is true
#' pdRankTests(array(c(pgram(Sigma), pgram(Sigma)), dim = c(2,2,2^8)), test = "signed.rank")
#' ## null is false
#' pdRankTests(array(c(pgram(Sigma), pgram(0.5 * Sigma)), dim = c(2,2,2^8)), test = "signed.rank")
#'
#' @seealso \code{\link{pdDepth}}, \code{\link{pdPgram}}
#'
#' @references Chau, J., Ombao, H., and von Sachs, R. (2017b). \emph{Data depth and rank-based
#' tests for covariance and spectral density matrices}. Available at \url{http://arxiv.org/abs/1706.08289}.
#' @references Brockwell, P.J. and Davis, R.A. (1991). \emph{Time series: Theory and Methods}. New York: Springer.
#'
#' @export
pdRankTests <- function(data, sample.sizes, test = c("rank.sum", "krusk.wall", "signed.rank", "bartels"),
                        depth = c("gdd", "zonoid", "spatial"), metric = c("Riemannian", "logEuclidean")) {
  if (missing(depth)) {
    depth <- "gdd"
  }
  ddim <- dim(data)
  if (missing(sample.sizes)) {
    sample.sizes <- NA
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean"))
  test <- match.arg(test, c("rank.sum", "krusk.wall", "signed.rank", "bartels"))
  depth <- match.arg(depth, c("gdd", "zonoid", "spatial"))
  err.message <- "Incorrect input lenghts for arguments: 'samples' and/or 'sample.sizes',
                                    consult the function documentation for the requested inputs."
  n <- sample.sizes
  if ((test == "krusk.wall") & (length(n) == 2)) {
    warning("Argument 'test' changed to 'rank.sum' to test for homogeneity of
                        distributions of two independent samples of HPD matrices.")
    test <- "rank.sum"
  }

  ## Manifold rank-sum test
  if (test == "rank.sum") {
    if (!isTRUE((((length(ddim) == 3) & (ddim[3] == sum(n))) | ((length(ddim) == 4) &
                      (ddim[4] == sum(n)))) & (ddim[1] == ddim[2]) & (length(n) == 2))) {
      stop(err.message)
    }

    dd <- pdDepth(X = data, method = depth, metric = metric)
    T1 <- (sum(rank(dd, ties.method = "random")[1:n[1]]) - n[1] * (sum(n) + 1)/2) /
                                                    sqrt(n[1] * n[2] * (sum(n) + 1)/12)

    output <- list(test = "Manifold Wilcoxon rank-sum", p.value = 2 * stats::pnorm(abs(T1), lower.tail = F), statistic = T1,
                                                  null.distr = "Standard normal distribution", depth.values = dd)
  }

  ## Manifold Kruskal-Wallis test
  if (test == "krusk.wall") {
    N <- sum(n)
    if (!isTRUE((((length(ddim) == 3) & (ddim[3] == N)) | ((length(ddim) == 4) &
                            (ddim[4] == N))) & (ddim[1] == ddim[2]) & (length(n) > 2))) {
      stop(err.message)
    }
    dd <- pdDepth(X = data, method = depth, metric = metric)
    R_bar <- unname(unlist(lapply(split(rank(dd, ties.method = "random"),
                                      f = rep(1:length(n), times = n)), mean)))
    T2 <- 12/(N * (N + 1)) * sum(n * (R_bar - (N + 1)/2)^2)

    output <- list(test = "Manifold Kruskal-Wallis", p.value = min(stats::pchisq(T2, df = 2, lower.tail = T),
                            pchisq(T2, df = 2, lower.tail = F)), statistic = T2,
                                    null.distr = "Chi-squared distribution (df = 2)", depth.values = dd)
  }

  ## Manifold signed-rank test
  if (test == "signed.rank") {
    if (!isTRUE((length(ddim) == 3) & (ddim[1] == ddim[2]) & (ddim[3]%%2 == 0))) {
      stop(err.message)
    }
    n <- ddim[3]/2
    d <- ddim[1]
    if(metric == "Riemannian"){
      ast <- function(A, B) t(Conj(A)) %*% B %*% A
      diff <- sapply(1:n, function(i) Re(sum(diag(Logm(diag(d), ast(iSqrt(data[, , n + i]), data[, , i]))))))
    } else{
      diff <- sapply(1:n, function(i) Re(sum(diag(Logm(diag(d), data[, , n + i]) - Logm(diag(d), data[, , i])))))
    }
    T3 <- stats::wilcox.test(x = diff, y = rep(0, n), paired = T, correct = T)
    output <- list(test = "Manifold Wilcoxon signed-rank", p.value = T3$p.value, statistic = T3$statistic, null.distr = T3$method)
  }

  ## Manifold Bartels-von Neumann test
  if (test == "bartels") {
    if (!isTRUE(((length(ddim) == 3) | ((length(ddim) == 4))) & (ddim[1] == ddim[2]))) {
      stop(err.message)
    }
    n <- utils::tail(ddim, 1)
    dd <- pdDepth(X = data, method = depth, metric = metric)
    T4 <- sum(diff(rank(dd, ties.method = "random"))^2)/(n * (n^2 - 1)/12)
    sigma <- sqrt(4 * (n - 2) * (5 * n^2 - 2 * n - 9)/(5 * n * (n + 1) * (n - 1)^2))

    output <- list(test = "Manifold Bartels-von Neumann", p.value = 2 * pnorm(abs((T4 - 2)/sigma), lower.tail = F),
                      statistic = (T4 - 2)/sigma, null.distr = "Standard normal distribution",
                   depth.values = dd)
  }

  return(output)
}
