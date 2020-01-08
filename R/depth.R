#' Data depth for HPD matrices
#'
#' \code{pdDepth} calculates the data depth of a HPD matrix with respect
#' to a given data cloud (i.e., a sample or collection) of HPD matrices, or the integrated
#' data depth of a sequence (curve) of HPD matrices with respect to a given data cloud of
#' sequences (curves) of HPD matrices as detailed in \insertCite{COvS17}{pdSpecEst}.
#'
#' Available pointwise or integrated intrinsic data depth functions for samples of HPD matrices are: (i)
#' geodesic distance depth, (ii) intrinsic zonoid depth and (iii) intrinsic spatial depth.
#' The various data depth measures and their theoretical properties are described in
#' \insertCite{COvS17}{pdSpecEst}. If \code{y} is a \eqn{(d,d)}-dimensional HPD matrix, \code{X} should be a \eqn{(d,d,S)}-dimensional array
#' corresponding to a length \code{S} sequence of \eqn{(d,d)}-dimensional HPD matrices and the pointwise
#' data depth values are computed. If \code{y} is a sequence of \eqn{(d,d)}-dimensional HPD matrices of length \code{n}
#' (i.e., \eqn{(d,d,n)}-dimensional array), \code{X} should be a \eqn{(d,d,n,S)}-dimensional array of replicated sequences of HPD matrices
#' and the integrated data depth values according to \insertCite{COvS17}{pdSpecEst} are computed. If \code{is.null(y)}, the data depth
#' of each individual object (i.e., a HPD matrix or a sequence of HPD matrices) in \code{X} is computed with
#' respect to the data cloud \code{X}. \cr
#' The function computes the intrinsic data depth values based on the metric space of HPD matrices equipped with
#' one of the following metrics: (i) Riemannian metric (default) as detailed in e.g., \insertCite{B09}{pdSpecEst}[Chapter 6] or
#' \insertCite{PFA05}{pdSpecEst}, (ii) log-Euclidean metric, the Euclidean inner product between matrix logarithms,
#' (iii) Cholesky metric, the Euclidean inner product between Cholesky decompositions, (iv) Euclidean metric and
#' (v) root-Euclidean metric. The default choice (Riemannian) has several properties not shared by the
#' other metrics, see \insertCite{COvS17}{pdSpecEst} for more details.
#'
#' @param y either a \eqn{(d,d)}-dimensional HPD matrix, or a \eqn{(d, d, n)}-dimensional array corresponding to a sequence
#'  or curve of HPD matrices. Defaults to \code{NULL}, in which case the data depth of each individual object in \code{X}
#'  with respect to the data cloud \code{X} itself is calculated.
#' @param X depending on the input \code{y}, \code{X} is either a \eqn{(d,d,S)}-dimensional array corresponding to a data cloud of
#'  \eqn{S} individual HPD matrices, or a \eqn{(d,d,n,S)}-dimensional array corresponding to a data cloud of \eqn{S}
#'  sequences or curves of \eqn{n} individual Hermitian PD matrices.
#' @param method the data depth measure, one of \code{'gdd'}, \code{'zonoid'} or \code{'spatial'} corresponding to
#' the geodesic distance depth, intrinsic zonoid depth, and intrinsic spatial depth respectively.
#' @param metric the metric that the space of HPD matrices is equipped with. The default choice is \code{"Riemannian"},
#' but this can also be one of: \code{"logEuclidean"}, \code{"Cholesky"}, \code{"rootEuclidean"} or
#' \code{"Euclidean"}. See also the Details section below.
#'
#' @return If \code{!is.null(y)}, \code{pdDepth} returns the numeric depth value of \code{y} with
#' respect to \code{X}. If \code{is.null(y)}, \code{pdDepth} returns a numeric vector of length \code{S} corresponding to
#' the vector of depth values for each individual object in \code{X} with respect to \code{X} itself.
#'
#' @examples
#' ## Pointwise depth
#' X1 <- replicate(50, Expm(diag(2), H.coeff(rnorm(4), inverse = TRUE)))
#' pdDepth(y = diag(2), X = X1) ## depth of one point
#' pdDepth(X = X1) ## depth of each point in the data cloud
#'
#' ## Integrated depth
#' X2 <- replicate(50, replicate(5, Expm(diag(2), H.coeff(rnorm(4), inverse = TRUE))))
#' pdDepth(y = replicate(5, diag(2)), X2, method = "zonoid", metric = "logEuclidean")
#' pdDepth(X = X2, method = "zonoid", metric = "logEuclidean")
#'
#' @seealso \code{\link{pdDist}}, \code{\link{pdRankTests}}
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and may fail
#' if matrices are close to being singular.
#'
#' @note
#' The data depth computations under the Riemannian metric are more involved than under the other
#' metrics, and may therefore result in (significantly) higher computation times.
#'
#' @references
#' \insertAllCited{}
#'
#' @importFrom ddalpha depth.zonoid
#' @importFrom ddalpha depth.spatial
#'
#' @export
pdDepth <- function(y = NULL, X, method = "gdd", metric = "Riemannian") {

  method <- match.arg(method, c("gdd", "zonoid", "spatial"))
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
    n <- 1
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

  if (method == "gdd") {

    ## Geodesic distance depth C++
    X <- array(X, dim = c(d, d, n * S))
    if(!is.null(y)){
      depth <- c(gdd_C(X, array(y, dim = c(d, d, n)), metric, S, "depth.y"))
    } else {
      depth <- c(gdd_C(X, array(0, dim = c(1, 1, 1)), metric, S, "depth.X"))
    }

  } else if (method == "zonoid") {

    ## intrinsic zonoid depth
    ZD <- function(y, X) {
      point <- ifelse(isTRUE(length(dim(y)) == 2), TRUE, FALSE)
      if(metric == "Riemannian" & point){
        X.coeff <- apply(X, 3, function(X) E_coeff(Logm(y, X)))
        X.coeff <- matrix(X.coeff[apply(X.coeff, 1, function(X){isTRUE(diff(range(X)) > .Machine$double.eps)}),],
                          ncol = dim(X)[3])
        y.coeff <- as.matrix(rep(0, nrow(X.coeff)))
      } else if(metric == "logEuclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(Logm(diag(d), X)))
        X.rows <- apply(X.coeff, 1, function(X){isTRUE(diff(range(X)) > .Machine$double.eps)})
        X.coeff <- X.coeff[X.rows, , drop = FALSE]
        y.coeff <- (if(point) E_coeff(Logm(diag(d), y))[X.rows] else X.coeff)
      } else if(metric == "Cholesky"){
        X.coeff <- apply(X, 3, function(X) E_chol(Chol_C(X, FALSE, FALSE)))
        X.rows <- apply(X.coeff, 1, function(X){isTRUE(diff(range(X)) > .Machine$double.eps)})
        X.coeff <- X.coeff[X.rows, , drop = FALSE]
        y.coeff <- (if(point) E_chol(Chol_C(y, FALSE, FALSE))[X.rows] else X.coeff)
      } else if(metric == "Euclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(X))
        X.rows <- apply(X.coeff, 1, function(X){isTRUE(diff(range(X)) > .Machine$double.eps)})
        X.coeff <- X.coeff[X.rows, , drop = FALSE]
        y.coeff <- (if(point) E_coeff(y)[X.rows] else X.coeff)
      } else if(metric == "rootEuclidean"){
        X.coeff <- apply(X, 3, function(X) E_coeff(Sqrt(X)))
        X.rows <- apply(X.coeff, 1, function(X){isTRUE(diff(range(X)) > .Machine$double.eps)})
        X.coeff <- X.coeff[X.rows, , drop = FALSE]
        y.coeff <- (if(point) E_coeff(Sqrt(y))[X.rows] else X.coeff)
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
        point <- ifelse(length(dim(y)) == 3, TRUE, FALSE)
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

  } else if (method == "spatial") {

    ## intrinsic spatial depth
    SD <- function(y, X) {
      point <- ifelse(length(dim(y)) == 2, TRUE, FALSE)
      if(metric == "Riemannian"){
        y.isqrt <- iSqrt(y)
        log.yx <- sapply(1:S, function(s) Logm(diag(d), (y.isqrt %*% X[, , s]) %*% y.isqrt), simplify = "array")
        dd <- 1 - NormF(apply(sapply(1:S, function(s) log.yx[, , s]/max(NormF(log.yx[, , s]), .Machine$double.eps),
                                     simplify = "array"), c(1, 2), mean))
      } else {
        if(metric == "logEuclidean"){
          X.coeff <- apply(X, 3, function(X) E_coeff(Logm(diag(d), X)))
          y.coeff <- (if(point) E_coeff(Logm(diag(d), y)) else X.coeff)
        } else if(metric == "Cholesky"){
          X.coeff <- apply(X, 3, function(X) E_chol(Chol_C(X, FALSE, FALSE)))
          y.coeff <- (if(point) E_chol(Chol_C(y, FALSE, FALSE)) else X.coeff)
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
        point <- ifelse(length(dim(y)) == 3, TRUE, FALSE)
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
