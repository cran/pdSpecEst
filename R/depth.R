#' Data depth for HPD matrices
#'
#' \code{pdDepth} calculates the data depth of a Hermitian PD matrix with respect
#' to a given data cloud (i.e. a sample or collection) of Hermitian PD matrices, or the integrated
#' data depth of a sequence (curve) of Hermitian PD matrices with respect to a given data cloud of
#' sequences (curves) of Hermitian PD matrices.
#'
#' Available pointwise or integrated manifold data depth functions for samples of Hermitian PD matrices are (i)
#' manifold zonoid depth, (ii) geodesic distance depth and (iii) manifold spatial depth.
#' The various data depth measures and their theoretical properties are described in
#' (Chau, Ombao, and von Sachs, 2017b). If \code{y} is a \eqn{(d \times d)}-dimensional Hermitian PD matrix (i.e.
#' \eqn{(d,d)}-dimensional array), \code{X} should be a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices and the pointwise
#' data depth values are computed. If \code{y} is a sequence of \eqn{(d \times d)}-dimensional Hermitian PD matrices of length \code{n}
#' (i.e. \eqn{(d,d,n)}-dimensional array), \code{X} should be a \eqn{(d,d,n,S)}-dimensional array of sequences of Hermitian PD matrices
#' and the integrated data depth values according to (Chau, Ombao, and von Sachs, 2017b) are computed. If \code{is.null(y)}, the data depth
#' of each individual object (i.e. a Hermitian PD matrix or a sequence of Hermitian PD matrices) in \code{X} is computed with
#' respect to the data cloud \code{X}. \cr
#' The function computes the intrinsic data depths based on the metric space of HPD matrices equipped with
#' one of the following metrics: (i) Riemannian metric (default) as in (Bhatia, 2009, Chapter 6),
#' (ii) log-Euclidean metric, the Euclidean inner product between matrix logarithms,
#' (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv) Euclidean metric and
#' (v) root-Euclidean metric. The default choice (Riemannian) has several appealing properties not shared by the
#' other metrics, see (Chau, Ombao and von Sachs, 2017b) for more details.
#'
#' @param y either a \eqn{(d,d)}-dimensional array corresponding to a \eqn{(d \times d)}-dimensional Hermitian PD matrix,
#'  or a \eqn{(d, d, n)}-dimensional array corresponding to a sequence or curve of Hermitian PD matrices.
#'  Defaults to \code{NULL}, in which case the data depth of each individual object in \code{X} with respect to the
#'  data cloud \code{X} itself is calculated.
#' @param X depending on the input \code{y} either a \eqn{(d,d,S)}-dimensional array corresponding to a data cloud of
#'  \eqn{S} individual Hermitian PD matrices, or a \eqn{(d,d,n,S)}-dimensional array corresponding to a data cloud of \eqn{S}
#'  sequences or curves of \eqn{n} individual Hermitian PD matrices.
#' @param method the data depth measure, one of \code{'zonoid'}, \code{'gdd'}, or \code{'spatial'} corresponding to
#' the manifold zonoid depth, geodesic distance depth, and manifold spatial depth respectively.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}.
#'
#' @return If \code{!is.null(y)}, \code{pdDepth} returns the numeric depth value of \code{y} with
#' respect to \code{X}. If \code{is.null(y)}, \code{pdDepth} returns a numeric vector of length \code{S} corresponding to the depth values
#' of each individual object in \code{X} with respect to \code{X} itself.
#'
#' @examples
#' ## Pointwise depth
#' X1 <- replicate(50, Expm(diag(2), H.coeff(rnorm(4), inverse = TRUE)))
#' pdDepth(y = diag(2), X1, method = "gdd") ## depth of one point
#' pdDepth(X = X1, method = "gdd") ## depth of each point in the data cloud
#'
#' ## Integrated depth
#' X2 <- replicate(50, replicate(5, Expm(diag(2), H.coeff(rnorm(4), inverse = TRUE))))
#' pdDepth(y = replicate(5, diag(2)), X2, method = "zonoid", metric = "logEuclidean")
#' pdDepth(X = X2, method = "zonoid", metric = "logEuclidean")
#'
#' @seealso \code{\link{pdDist}}, \code{\link{pdRankTests}}
#'
#' @note Note that the data depth computations in the metric space equipped with the Riemannian
#' metric may be significantly slower than the depth computations based on one of the alternative
#' metrics.
#'
#' @references Chau, J., Ombao, H., and von Sachs, R. (2017b). \emph{Data depth and rank-based
#' tests for covariance and spectral density matrices}. Available at \url{http://arxiv.org/abs/1706.08289}.
#'
#' @importFrom ddalpha depth.zonoid
#' @importFrom ddalpha depth.spatial
#'
#' @export
pdDepth <- function(y = NULL, X, method = c("zonoid", "gdd", "spatial"), metric = "Riemannian") {

  method <- match.arg(method, c("zonoid", "gdd", "spatial"))
  metric <- match.arg(metric, c("Riemannian", "logEuclidean", "Cholesky", "Euclidean", "rootEuclidean"))
  d <- dim(X)[1]
  err.message <- "Incorrect input dimensions for arguments: 'y' and/or 'X',
                      consult the function documentation for the requested inputs."
  if (length(dim(X)) == 3) {
    if (!is.null(y)) {
      if (!isTRUE(all.equal(dim(y), dim(X)[c(1, 2)]))) {
        stop(err.message)
      }
    }
    S <- dim(X)[3]
  } else if (length(dim(X)) == 4) {
    if (!is.null(y)) {
      if (!isTRUE(all.equal(dim(y), dim(X)[c(1, 2, 3)]))) {
        stop(err.message)
      }
    }
    S <- dim(X)[4]
    n <- dim(X)[3]
  } else {
    stop(err.message)
  }

  ## Manifold zonoid depth
  if (method == "zonoid") {
    ZD <- function(y, X) {
      point <- ifelse(length(dim(y)) == 2, T, F)
      if(metric == "Riemannian" & point){
        y.coeff <- as.matrix(rep(0, d^2))
        X.coeff <- apply(X, 3, function(X) E_coeff(Logm(y, X)))
      } else if(metric == "logEuclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(Logm(diag(d), X)))
        y.coeff <- (if(point) E_coeff(Logm(diag(d), y)) else X.coeff)
      } else if(metric == "Cholesky"){
        X.coeff <- apply(X, 3, function(X) E_chol(Chol(X)))
        y.coeff <- (if(point) E_chol(Chol(y)) else X.coeff)
      } else if(metric == "Euclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(X))
        y.coeff <- (if(point) E_coeff(y) else X.coeff)
      } else if(metric == "rootEuclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(Sqrt(X)))
        y.coeff <- (if(point) E_coeff(Sqrt(y)) else X.coeff)
      }
      ddalpha::depth.zonoid(t(y.coeff), t(X.coeff))
    }

    ## Pointwise zonoid depth
    if (length(dim(X)) == 3) {
      if (!is.null(y)) {
        depth <- ZD(y, X)
      } else {
        depth <- (if(metric == "Riemannian") sapply(1:S, function(s) ZD(X[, , s], X)) else{
            ZD(X, X) })
      }
    }

    ## Integrated zonoid depth
    if (length(dim(X)) == 4) {
      iZD <- function(y, X) {
        point <- ifelse(length(dim(y)) == 3, T, F)
        if(point){
          dd <- mean(sapply(1:n, function(t) ZD(y[, , t], X[, , t, ])))
        } else{
          if(metric == "Riemannian"){
            dd <- sapply(1:S, function(s) mean(sapply(1:n, function(t) ZD(y[, , t, s], X[, , t, ]))))
          } else{
            dd <- rowMeans(sapply(1:n, function(t) ZD(y[, , t, ], X[, , t, ])))
          }
        }
        dd
      }
      depth <- iZD(y = (if(!is.null(y)) y else X), X = X)
    }
  }

  ## Geodesic distance depth
  if (method == "gdd") {
    if (!is.null(y)) {
      if (length(dim(X)) == 3) {
        depth <- exp(-mean(sapply(1:S, function(s) pdDist(y, X[, , s], method = metric))))
      } else if (length(dim(X)) == 4) {
        depth <- exp(-mean(sapply(1:S, function(s) mean(sapply(1:n,
                                            function(t) pdDist(y[, , t], X[, , t, s], method = metric))))))
      }
    } else {
      dist <- matrix(0, nrow = S, ncol = S)
      grid <- expand.grid(1:S, 1:S)
      grid <- grid[grid$Var1 > grid$Var2, ]
      if (length(dim(X)) == 3) {
        dist[lower.tri(dist)] <- mapply(function(i, j) pdDist(X[, , i], X[, , j], method = metric),
                                        grid$Var1, grid$Var2)
      } else if (length(dim(X)) == 4) {
        dist[lower.tri(dist)] <- mapply(function(i, j) mean(sapply(1:n,
                                            function(t) pdDist(X[, , t, i], X[, , t, j], method = metric))),
                                        grid$Var1, grid$Var2)
      }
      dist[upper.tri(dist)] <- t(dist)[upper.tri(dist)]
      depth <- exp(-colMeans(dist))
    }
  }

  ## Manifold spatial depth
  if (method == "spatial") {
    SD <- function(y, X) {
      point <- ifelse(length(dim(y)) == 2, T, F)
      if(metric == "Riemannian"){
        y.isqrt <- iSqrt(y)
        log.yx <- sapply(1:S, function(s) Logm(diag(d), (y.isqrt %*% X[, , s]) %*% y.isqrt), simplify = "array")
        dd <- 1 - NormF(apply(sapply(1:S,
                                                 function(s) log.yx[, , s]/max(NormF(log.yx[, , s]), .Machine$double.eps),
                                                 simplify = "array"), c(1, 2), mean))
      } else {
        if(metric == "logEuclidean"){
          X.coeff <- apply(X, 3, function(X) E_coeff(Logm(diag(d), X)))
          y.coeff <- (if(point) E_coeff(Logm(diag(d), y)) else X.coeff)
        } else if(metric == "Cholesky"){
          X.coeff <- apply(X, 3, function(X) E_chol(Chol(X)))
          y.coeff <- (if(point) E_chol(Chol(y)) else X.coeff)
        } else if(metric == "Euclidean"){
          X.coeff <- apply(X, 3, function(X) E_coeff(X))
          y.coeff <- (if(point) E_coeff(y) else X.coeff)
        } else if(metric == "rootEuclidean"){
          X.coeff <- apply(X, 3, function(X) E_coeff(Sqrt(X)))
          y.coeff <- (if(point) E_coeff(Sqrt(y)) else X.coeff)
        }
        dd <- ddalpha::depth.spatial(t(y.coeff), t(X.coeff), mah.estimate = "none")
      }
      dd
    }

    ## Pointwise spatial depth
    if (length(dim(X)) == 3) {
      if (!is.null(y)) {
        depth <- SD(y, X)
      } else {
        depth <- (if(metric == "Riemannian") sapply(1:S, function(s) SD(X[, , s], X)) else{
          SD(X, X) })
      }
    }

    ## Integrated spatial depth
    if (length(dim(X)) == 4) {
      iSD <- function(y, X) {
        point <- ifelse(length(dim(y)) == 3, T, F)
        if(point){
          dd <- mean(sapply(1:n, function(t) SD(y[, , t], X[, , t, ])))
        } else{
          if(metric == "Riemannian"){
            dd <- sapply(1:S, function(s) mean(sapply(1:n, function(t) SD(y[, , t, s], X[, , t, ]))))
          } else{
            dd <- rowMeans(sapply(1:n, function(t) SD(y[, , t, ], X[, , t, ])))
          }
        }
        dd
      }
      depth <- iSD(y = (if(!is.null(y)) y else X), X = X)
    }
  }
  return(depth)
}

