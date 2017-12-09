#' Polynomial interpolation of curves (1D) or surfaces (2D) of HPD matrices
#'
#' \code{pdNeville()} performs intrinsic polynomial interpolation of curves or surfaces on the manifold of
#' HPD matrices equipped with the Riemannian metric via iterative geodesic interpolation, see e.g.
#' (Chau and von Sachs, 2017).
#'
#' For polynomial curve interpolation, given \eqn{N} control points (i.e. HPD
#' matrices), the degree of the interpolated polynomial is \eqn{N - 1}. For polynomial surface interpolation,
#' given \eqn{N_1 \times N_2} control points (i.e. HPD matrices) on a tensor product grid, the interpolated
#' polynomial has bi-degree \eqn{(N_1 - 1, N_2 - 1)}. The function \code{pdNeville()} determines whether
#' polynomial curve or polynomial surface interpolation has to be performed based on the function input.
#'
#' @param P for polynomial curve interpolation, a \eqn{(d, d, N)}-dimensional array corresponding to a sequence
#' of \eqn{(d, d)}-dimensional HPD matrices, i.e. control points, through which the interpolated polynomial
#' curve passes. For polynomial surface interpolation, a \eqn{(d, d, N_1, N_2)}-dimensional array corresponding
#' to a tensor product grid of \eqn{(d, d)}-dimensional matrices, i.e. control points, through which the interpolated
#' polynomial surface passes.
#' @param X for polynomial curve interpolation, a numeric vector of length \eqn{N} specifying the time points
#' the interpolated polynomial passes through the control points \code{P}. For polynomial surface interpolation, a list
#' with as elements two numeric vectors \code{x} and \code{y} of length \eqn{N_1} and \eqn{N_2} respectively. The numeric
#' vectors specify the time points on the tensor product grid \code{expand.grid(X$x, X$y)} the interpolated polynomial passes
#' trough the control points \code{P}.
#' @param x for polynomial curve interpolation, a numeric vector specifying the time grid (resolution) at which the
#' interpolated polynomial is estimated. For polynomial surface interpolation, a list with as elements two numeric vectors
#' \code{x} and \code{y} specifying the time tensor product grid (resolution) \code{expand.grid(x$x, x$y)} at which the
#' interpolated polynomial surface is estimated.
#' @param metric the metric the space of HPD gets equipped with, by default \code{metric = "Riemannian"}, but instead
#' this can also be set to \code{metric = "Euclidean"} to perform (standard) Euclidean polynomial interpolation of curves or
#' surfaces in the space of HPD matrices.
#'
#' @return For polynomial curve interpolation, a \code{(d, d, length(x))}-dimensional array containing the interpolated
#' polynomial of degree \eqn{N-1} at the time grid \code{x} passing through the control points \code{P} at times \code{X}.
#' For polynomial surface interpolation, a \code{(d, d, length(x$x), length(x$y))}-dimensional array containing the interpolated
#' polynomial of bi-degree \eqn{N_1 - 1, N_2 - 1} at the time grid \code{expand.grid(x$x, x$y)} passing through the control points
#' \code{P} at times \code{expand.grid(X$x, X$y)}.
#'
#' @examples
#' ### Polynomial curve interpolation
#' P <- rExamples(100, example = 'gaussian')$f[, , 10*(1:5)]
#' P.poly <- pdNeville(P, (1:5)/5, (1:50)/50)
#' ## Examine matrix-component (1,1)
#' plot((1:50)/50, Re(P.poly[1, 1, ]), type = "l") ## interpolated polynomial
#' lines((1:5)/5, Re(P[1, 1, ]), col = 2) ## control points
#'
#' ### Polynomial surface interpolation
#' P.surf <- array(P[, , 1:4], dim = c(2,2,2,2)) ## control points
#' P.poly <- pdNeville(P.surf, list(x = c(0, 1), y = c(0, 1)), list(x = (0:10)/10, y = (0:10)/10))
#'
#' @seealso \code{\link{pdPolynomial}}
#'
#' @references Chau, J. and von Sachs, R. (2017a). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#'
#' @export
pdNeville <- function(P, X, x, metric = "Riemannian"){

  ## Set variables
  if(!isTRUE((length(dim(P)) == 3 | length(dim(P)) == 4) & (is.complex(P) | is.numeric(P)))){
    stop("'P' should be a numeric or complex 3- or 4-dimensional array")
  }
  is_2D = (length(dim(P)) == 4)
  if(!is_2D){
    if(!isTRUE(dim(P)[3] == length(X))){
      stop("The number of control points 'dim(P)[3]' should be equal to 'length(X)'.")
    }
    n = dim(P)[3] - 1
  } else if(is_2D){
    if(!isTRUE(dim(P)[3] == length(X$x) & dim(P)[4] == length(X$y))){
      stop("The dimensions of the grid of control points in 'P' and 'X' do not match.")
    }
    n = c(dim(P)[3] - 1, dim(P)[4] - 1)
  }
  d = dim(P)[1]
  metric = match.arg(metric, match.arg(metric, c("Riemannian", "Euclidean")))

  ## 1D Neville's algorithm via geodesic curve interpolation
  pdNeville1D <- function(P, X, x){

    nn <- dim(P)[3] - 1
    geo <- function(Pi, t.range, t){
      if(metric == "Riemannian"){
        res <- Expm(Pi[, , 1], (t - t.range[1]) / (t.range[2] - t.range[1]) * Logm(Pi[, , 1], Pi[, , 2]))
      } else if(metric == "Euclidean"){
        res <- (t.range[2] - t)/(t.range[2] - t.range[1]) * Pi[, , 1] +
                  (t - t.range[1])/(t.range[2] - t.range[1]) * Pi[, , 2]
      }
      res
    }
    ind <- expand.grid(1:nn, x)
    Li <- array(c(mapply(function(i, xj) geo(P[, , c(i, i + 1)], X[c(i, i + 1)], xj),
                         ind$Var1, ind$Var2)), dim = c(d, d, nn, length(x)))
    if(nn > 1){
      for(m in 2:nn){
        ind_m <- expand.grid(1:(nn - m  + 1), 1:length(x))
        Li <- array(c(mapply(function(i, j) geo(Li[, , c(i, i + 1), j], X[c(i, i + m)], x[j]),
                             ind_m$Var1, ind_m$Var2)), dim = c(d, d, nn - m + 1, length(x)))
      }
    }
    Li[, , 1, ]
  }

  if(!is_2D) {

    PP <- pdNeville1D(P, X, x)

  } else if(is_2D) {

    ## 2D Neville's algorithm via geodesic surface interpolation
    if(n[1] < 1 & n[2] < 1){
      PP <- array(P[, , 1, 1], dim = c(d, d, length(x$x), length(x$y)))
    } else if(n[1] < 1 | n[2] < 1){
      if(length(X$y) > 1){
        PP_i <- pdNeville1D(P[, , 1, ], X$y, x$y)
        PP <- aperm(sapply(1:length(x$x), function(i) PP_i, simplify = "array"), c(1, 2, 4, 3))
      } else if(length(X$x) > 1){
        PP_j <- pdNeville1D(P[, , , 1], X$x, x$x)
        PP <- sapply(1:length(x$y), function(j) PP_j, simplify = "array")
      }
    } else {
      PP_i <- sapply(1:(n[1] + 1), function(i) pdNeville1D(P[, , i, ], X$y, x$y), simplify = "array")
      PP <- sapply(1:length(x$y), function(j) pdNeville1D(PP_i[, , j, ], X$x, x$x), simplify = "array")
    }

  }
  return(PP)
}

#' Generate intrinsic HPD polynomial curves
#'
#' \code{pdPolynomial()} generates intrinsic polynomial curves on the manifold of HPD matrices
#' equipped with the Riemannian metric according to the numerical integration procedure described in (Hinkle et al., 2014).
#' Given an initial starting point \code{p0} (i.e. HPD matrix) on the Riemannian manifold and the covariant
#' derivatives up to order \eqn{k - 1} at \code{p0}, \code{pdPolynomial()} approximates the uniquely existing
#' intrinsic polynomial curve of degree \eqn{k} passing through \code{p0} with the given covariant derivatives up
#' to order \eqn{k - 1} and vanishing higher order covariant derivatives.
#'
#' @param p0 a \eqn{(d, d)}-dimensional HPD matrix specifying the starting point of the polynomial curve.
#' @param v0 a \eqn{(d, d, k)}-dimensional array corresponding to a sequence of covariant derivatives of
#' order zero up to order \eqn{k - 1} at the starting point \code{p0}.
#' @param delta.t a numeric value determining the incrementing step size in the numerical integration procedure.
#' A smaller step size results in a higher resolution and therefore a more accurate approximation of the polynomial curve,
#' defaults to \code{delta.t = 0.01}.
#' @param steps number of incrementing steps in the numerical integration procedure, defaults to \code{steps = 100}.
#'
#' @examples
#' ## First-order polynomial
#' p0 <- diag(3) ## HPD starting point
#' v0 <- array(H.coeff(rnorm(9), inverse = TRUE), dim = c(3, 3, 1)) ## zero-th order cov. derivative
#' P.poly <- pdPolynomial(p0, v0)
#'
#' ## First-order polynomials coincide with geodesic curves
#' geo <- function(A, B, t) Expm(A, t * Logm(A, B))
#' P.geo <- sapply(seq(0, 1, length = 100), function(t) geo(p0, P.poly[, , 100], t),
#'                 simplify = "array")
#' all.equal(P.poly, P.geo)
#'
#' @return A \code{(d, d, length(steps))}-dimensional array containing the approximated intrinsic polynomial
#' curve of degree \eqn{k} passing through \code{p0} with the given covariant derivatives up to order
#' \eqn{k - 1} and vanishing higher order covariant derivatives.
#'
#' @references
#' Hinkle, J., Fletcher, P. and Joshi, S. (2014). Intrinsic polynomials for regression on Riemannian manifolds.
#' \emph{Journal of Mathematical Imaging and Vision} 50, 32-52.
#'
#' @seealso \code{\link{pdNeville}}
#'
#' @export
pdPolynomial <- function(p0, v0, delta.t = 0.01, steps = 100) {

  if(!isTRUE(length(dim(v0)) == 3)){
    stop("'v0' should be a 3-dimensional array of covariant derivatives.")
  }
  d <- dim(p0)[1]
  k <- max(dim(v0)[3], 2)
  p <- array(dim = c(d, d, steps))
  p[, , 1] <- p0
  vi <- array(0, dim = c(d, d, k))
  vi[, , 1:dim(v0)[3]] <- v0

  for(ti in 1:(steps - 1)){
    w <- vi[, , 1]
    vi[, , 1:(k - 1)] <- sapply(1:(k - 1), function(i) ParTrans(p[, , ti], delta.t * w, vi[, , i] +
                                                                  delta.t * vi[, , i + 1]), simplify = "array")
    vi[, , k] <- ParTrans(p[, , ti], delta.t * w, vi[, , k])
    p[, , ti + 1] <- Expm(p[, , ti], delta.t * w)
  }

  return(p)
}

