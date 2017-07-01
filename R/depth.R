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
#' respect to the data cloud \code{X}.
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
#'
#' @return If \code{!is.null(y)}, \code{pdDepth} returns the numeric depth value of \code{y} with
#' respect to \code{X}. If \code{is.null(y)}, \code{pdDepth} returns a numeric vector of length \code{S} corresponding to the depth values
#' of each individual object in \code{X} with respect to \code{X} itself.
#'
#' @examples
#' ## Pointwise depth
#' X1 <- replicate(50, Expm(diag(2), pdSpecEst:::E_coeff_inv(rnorm(4))))
#' pdDepth(y = diag(2), X1, method = "gdd") ## depth of one point
#' pdDepth(X = X1, method = "gdd") ## depth of each point in the data cloud
#'
#' ## Integrated depth
#' X2 <- replicate(50, replicate(5, Expm(diag(2), pdSpecEst:::E_coeff_inv(rnorm(4)))))
#' pdDepth(y = replicate(5, diag(2)), X2, method = "gdd") ## depth of one curve
#' pdDepth(X = X2, method = "gdd") ## depth of each curve in the data cloud
#'
#' @seealso \code{\link{pdDist}}, \code{\link{pdRankTests}}
#'
#' @references Chau, J., Ombao, H., and von Sachs, R. (2017b). \emph{Data depth and rank-based
#' tests for covariance and spectral density matrices}. Available at \url{http://arxiv.org/abs/1706.08289}.
#'
#' @importFrom ddalpha depth.zonoid
#' @export
pdDepth <- function(y = NULL, X, method = c("zonoid", "gdd", "spatial")) {

  method <- match.arg(method, c("zonoid", "gdd", "spatial"))
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
    if (length(dim(X)) == 3) {
      ZD <- function(y, X) {
        # E_y <- T_basis(E, y)
        return(ddalpha::depth.zonoid(t(as.matrix(rep(0, d^2))), t(sapply(1:S,
                                        function(s) E_coeff(Logm(y, X[, , s]))))))
      }
      if (!is.null(y)) {
        depth <- ZD(y, X)
      } else {
        depth <- sapply(1:S, function(s) ZD(X[, , s], X))
      }
    } else if (length(dim(X)) == 4) {
      iZD <- function(y, X) {
        depth.t <- numeric(n)
        for (t in 1:n) {
          # E_y <- T_basis(E, y[, , t])
          depth.t[t] <- ddalpha::depth.zonoid(t(as.matrix(rep(0, d^2))), t(sapply(1:S,
                                        function(s) E_coeff(Logm(y[, , t], X[, , t, s])))))
        }
        return(mean(depth.t))
      }
      if (!is.null(y)) {
        depth <- iZD(y, X)
      } else {
        depth <- sapply(1:S, function(s) iZD(X[, , , s], X))
      }
    }
  }

  ## Geodesic distance depth
  if (method == "gdd") {
    if (!is.null(y)) {
      if (length(dim(X)) == 3) {
        depth <- exp(-mean(sapply(1:S, function(s) pdDist(y, X[, , s]))))
      } else if (length(dim(X)) == 4) {
        depth <- exp(-mean(sapply(1:S, function(s) mean(sapply(1:n,
                                            function(t) pdDist(y[, , t], X[, , t, s]))))))
      }
    } else {
      dist <- matrix(0, nrow = S, ncol = S)
      grid <- expand.grid(1:S, 1:S)
      grid <- grid[grid$Var1 > grid$Var2, ]
      if (length(dim(X)) == 3) {
        dist[lower.tri(dist)] <- mapply(function(i, j) pdDist(X[, , i], X[, , j]), grid$Var1, grid$Var2)
      } else if (length(dim(X)) == 4) {
        dist[lower.tri(dist)] <- mapply(function(i, j) mean(sapply(1:n,
                                            function(t) pdDist(X[, , t, i], X[, , t, j]))), grid$Var1, grid$Var2)
      }
      dist[upper.tri(dist)] <- t(dist)[upper.tri(dist)]
      depth <- exp(-colMeans(dist))
    }
  }

  ## Manifold spatial depth
  if (method == "spatial") {
    if (length(dim(X)) == 3) {
      SD <- function(y, X) {
        y.isqrt <- iSqrt(y)
        log.yx <- sapply(1:S, function(s) Logm(diag(d), (y.isqrt %*% X[, , s]) %*% y.isqrt), simplify = "array")
        return(1 - NormF(apply(sapply(1:S,
            function(s) log.yx[, , s]/NormF(log.yx[, , s]), simplify = "array"), c(1, 2), mean)))
      }
      if (!is.null(y)) {
        depth <- SD(y, X)
      } else {
        depth <- sapply(1:S, function(s) SD(X[, , s], X))
      }
    } else if (length(dim(X)) == 4) {
      iSD <- function(y, X) {
        depth.t <- numeric(n)
        for (t in 1:n) {
          y.isqrt <- iSqrt(y[, , t])
          log.yx <- sapply(1:S, function(s) Logm(diag(d), (y.isqrt %*% X[, , t, s]) %*% y.isqrt), simplify = "array")
          depth.t[t] <- 1 - NormF(apply(sapply(1:S,
                                function(s) log.yx[, , s]/NormF(log.yx[, , s]), simplify = "array"), c(1, 2), mean))
        }
        return(mean(depth.t))
      }
      if (!is.null(y)) {
        depth <- iSD(y, X)
      } else {
        depth <- sapply(1:S, function(s) iSD(X[, , , s], X))
      }
    }
  }
  return(depth)
}

