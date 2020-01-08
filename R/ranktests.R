#' Rank-based hypothesis tests for HPD matrices
#'
#' \code{pdRankTests} performs a number of generalized rank-based hypothesis tests in the metric space of HPD matrices equipped
#' with the affine-invariant Riemannian metric or Log-Euclidean metric for samples of HPD matrices or samples of sequences
#' (curves) of HPD matrices as described in Chapter 4 of \insertCite{C18}{pdSpecEst}.
#'
#' For samples of \eqn{(d,d)}-dimensional HPD matrices with pooled sample size \eqn{S}, the argument
#' \code{data} is a \eqn{(d,d,S)}-dimensional array of \eqn{(d,d)}-dimensional HPD matrices, where the individual samples are
#' combined along the third array dimension. For samples of sequences of \eqn{(d,d)}-dimensional HPD matrices with pooled sample
#' size \eqn{S}, the argument \code{data} is a \eqn{(d,d,n,S)}-dimensional array of length \eqn{n} sequences
#' of \eqn{(d,d)}-dimensional HPD matrices, where the individual samples are combined along the fourth array dimension. The argument
#' \code{sample_sizes} specifies the sizes of the individual samples so that \code{sum(sample_sizes)} is equal to \code{S}. \cr
#' The available generalized rank-based testing procedures (specified by the argument \code{test}) are:
#' \describe{
#' \item{\code{"rank.sum"}}{Intrinsic Wilcoxon rank-sum test to test for homogeneity of distributions of two independent
#' samples of HPD matrices or samples of sequences of HPD matrices. The usual univariate ranks are replaced by data depth
#' induced ranks obtained with \code{\link{pdDepth}}.}
#' \item{\code{"krusk.wall"}}{Intrinsic Kruskal-Wallis test to test for homogeneity of distributions of more than two independent
#' samples of HPD matrices or samples of sequences of HPD matrices. The usual univariate ranks are replaced by data depth
#' induced ranks obtained with \code{\link{pdDepth}}.}
#' \item{\code{"signed.rank"}}{Intrinsic signed-rank test to test for homogeneity of distributions of independent paired or matched samples
#' of HPD matrices. The intrinsic signed-rank test is \emph{not} based on data depth induced ranks, but on a specific difference score in the Riemannian
#' manifold of HPD matrices equipped with either the affine-invariant Riemannian or Log-Euclidean metric.}
#' \item{\code{"bartels"}}{Intrinsic Bartels-von Neumann test to test for randomness (i.e., exchangeability) within a single independent sample of
#' HPD matrices or a sample of sequences of HPD matrices. The usual univariate ranks are replaced by data depth induced
#' ranks obtained with \code{\link{pdDepth}}.}
#' }
#' The function computes the generalized rank-based test statistics in the \emph{complete} metric space of HPD matrices equipped with one of the following metrics:
#' (i) the Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or \insertCite{PFA05}{pdSpecEst}; or (ii) the Log-Euclidean metric,
#' the Euclidean inner product between matrix logarithms. The default Riemannian metric is invariant under congruence transformation by any invertible matrix,
#' whereas the Log-Euclidean metric is only invariant under congruence transformation by unitary matrices, see \insertCite{C18}{pdSpecEst}[Chapter 4] for more details.
#'
#' @note The intrinsic signed-rank test also provides a valid test for equivalence of spectral matrices of two multivariate stationary time
#' series based on the HPD periodogram matrices obtained via \code{\link{pdPgram}}, see \insertCite{C18}{pdSpecEst}[Chapter 4] for the details.
#'
#' @note The function does not check for positive definiteness of the input matrices, and may fail
#' if matrices are close to being singular.
#'
#' @note The data depth computations under the Riemannian metric are more involved than under the Log-Euclidean
#' metric, and may therefore result in (significantly) higher computation times.
#'
#' @param data either a \eqn{(d,d,S)}-dimensional array corresponding to an array of pooled individual samples of \eqn{(d,d)}-dimensional
#' HPD matrices, or a \eqn{(d,d,n,S)}-dimensional array corresponding to an array of pooled individual samples of length \eqn{n} sequences
#' of \eqn{(d,d)}-dimensional HPD matrices.
#' @param sample_sizes a numeric vector specifying the individual sample sizes in the pooled sample \code{data}, such that \code{sum(sample_sizes)} is
#' equal to \code{S}. Not required for tests \code{"signed-rank"} and \code{"bartels"}, as the sample sizes are automatically determined from the input array
#' \code{data}.
#' @param test rank-based hypothesis testing procedure, one of \code{"rank.sum"}, \code{"krusk.wall"}, \code{"signed.rank"}, \code{"bartels"} explained
#' in the Details section below.
#' @param depth data depth measure used in the rank-based tests, one of \code{"gdd"}, \code{"zonoid"}, or \code{"spatial"} corresponding to the
#' geodesic distance depth, intrinsic zonoid depth and intrinsic spatial depth respectively. Defaults to \code{"gdd"}. Not required for test
#' \code{"signed.rank"}. See the documentation of the function \code{\link{pdDepth}} for additional details about the different depth measures.
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
#' pdRankTests(data, sample_sizes = c(50, 50), test = "rank.sum") ## homogeneity 2 samples
#' pdRankTests(data, sample_sizes = rep(25, 4), test = "krusk.wall") ## homogeneity 4 samples
#' pdRankTests(data, test = "bartels") ## randomness
#'
#' ## null hypothesis is false
#' data1 <- array(c(data, replicate(50, Expm(diag(2), H.coeff(0.5 * rnorm(4), inverse = TRUE)))),
#'                  dim = c(2,2,150))
#' pdRankTests(data1, sample_sizes = c(100, 50), test = "rank.sum")
#' pdRankTests(data1, sample_sizes = rep(50, 3), test = "krusk.wall")
#' pdRankTests(data1, test = "bartels")
#'
#' \dontrun{
#' ## signed-rank test for equivalence of spectra of multivariate time series
#' ## ARMA(1,1) process: Example 11.4.1 in (Brockwell and Davis, 1991)
#' Phi <- array(c(0.7, 0, 0, 0.6, rep(0, 4)), dim = c(2, 2, 2))
#' Theta <- array(c(0.5, -0.7, 0.6, 0.8, rep(0, 4)), dim = c(2, 2, 2))
#' Sigma <- matrix(c(1, 0.71, 0.71, 2), nrow = 2)
#' pgram <- function(Sigma) pdPgram(rARMA(2^8, 2, Phi, Theta, Sigma)$X)$P
#'
#' ## null is true
#' pdRankTests(array(c(pgram(Sigma), pgram(Sigma)), dim = c(2,2,2^8)), test = "signed.rank")
#' ## null is false
#' pdRankTests(array(c(pgram(Sigma), pgram(0.5 * Sigma)), dim = c(2,2,2^8)), test = "signed.rank")
#' }
#' @seealso \code{\link{pdDepth}}, \code{\link{pdPgram}}
#'
#' @references
#' \insertAllCited{}
#'
#' @export
pdRankTests <- function(data, sample_sizes, test = c("rank.sum", "krusk.wall", "signed.rank", "bartels"),
                        depth = c("gdd", "zonoid", "spatial"), metric = c("Riemannian", "logEuclidean")) {
  if (missing(depth)) {
    depth <- "gdd"
  }
  ddim <- dim(data)
  if (missing(sample_sizes)) {
    sample_sizes <- NA
  }
  metric <- match.arg(metric, c("Riemannian", "logEuclidean"))
  test <- match.arg(test, c("rank.sum", "krusk.wall", "signed.rank", "bartels"))
  depth <- match.arg(depth, c("gdd", "zonoid", "spatial"))
  err.message <- "Incorrect input lenghts for arguments: 'samples' and/or 'sample_sizes',
                                    consult the function documentation for the requested inputs."
  n <- sample_sizes
  if ((test == "krusk.wall") & (length(n) == 2)) {
    warning("Argument 'test' changed to 'rank.sum' to test for homogeneity of
                        distributions of two independent samples of HPD matrices.")
    test <- "rank.sum"
  }

  ## Intrinsic rank-sum test
  if (test == "rank.sum") {
    if (!isTRUE((((length(ddim) == 3) & (ddim[3] == sum(n))) | ((length(ddim) == 4) &
                      (ddim[4] == sum(n)))) & (ddim[1] == ddim[2]) & (length(n) == 2))) {
      stop(err.message)
    }

    dd <- pdDepth(X = data, method = depth, metric = metric)
    T1 <- (sum(rank(dd, ties.method = "random")[1:n[1]]) - n[1] * (sum(n) + 1)/2) /
                                                    sqrt(n[1] * n[2] * (sum(n) + 1)/12)

    output <- list(test = "Intrinsic Wilcoxon rank-sum", p.value = 2 * stats::pnorm(abs(T1), lower.tail = FALSE), statistic = T1,
                                                  null.distr = "Standard normal distribution", depth.values = dd)
  }

  ## Intrinsic Kruskal-Wallis test
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

    output <- list(test = "Intrinsic Kruskal-Wallis", p.value = min(stats::pchisq(T2, df = 2, lower.tail = TRUE),
                            pchisq(T2, df = 2, lower.tail = FALSE)), statistic = T2,
                                    null.distr = "Chi-squared distribution (df = 2)", depth.values = dd)
  }

  ## Intrinsic signed-rank test
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
    T3 <- stats::wilcox.test(x = diff, y = rep(0, n), paired = TRUE, correct = TRUE)
    output <- list(test = "Intrinsic Wilcoxon signed-rank", p.value = T3$p.value, statistic = T3$statistic, null.distr = T3$method)
  }

  ## Intrinsic Bartels-von Neumann test
  if (test == "bartels") {
    if (!isTRUE(((length(ddim) == 3) | ((length(ddim) == 4))) & (ddim[1] == ddim[2]))) {
      stop(err.message)
    }
    n <- utils::tail(ddim, 1)
    dd <- pdDepth(X = data, method = depth, metric = metric)
    T4 <- sum(diff(rank(dd, ties.method = "random"))^2)/(n * (n^2 - 1)/12)
    sigma <- sqrt(4 * (n - 2) * (5 * n^2 - 2 * n - 9)/(5 * n * (n + 1) * (n - 1)^2))

    output <- list(test = "Intrinsic Bartels-von Neumann", p.value = 2 * pnorm(abs((T4 - 2)/sigma), lower.tail = FALSE),
                      statistic = (T4 - 2)/sigma, null.distr = "Standard normal distribution",
                   depth.values = dd)
  }

  return(output)
}
