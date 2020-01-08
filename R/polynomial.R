#' Polynomial interpolation of curves (1D) or surfaces (2D) of HPD matrices
#'
#' \code{pdNeville} performs intrinsic polynomial interpolation of curves or surfaces of HPD matrices
#' in the metric space of HPD matrices equipped with the affine-invariant Riemannian metric (see \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}) via Neville's algorithm based on iterative geodesic interpolation detailed
#' in \insertCite{CvS17}{pdSpecEst} and in Chapter 3 and 5 of \insertCite{C18}{pdSpecEst}.
#'
#' For polynomial curve interpolation, given \eqn{N} control points (i.e., HPD matrices), the degree of the
#' interpolated polynomial is \eqn{N - 1}. For polynomial surface interpolation, given \eqn{N_1 \times N_2} control points
#' (i.e., HPD matrices) on a tensor product grid, the interpolated polynomial surface is of bi-degree \eqn{(N_1 - 1, N_2 - 1)}.
#' Depending on the input array \code{P}, the function decides whether polynomial curve or polynomial surface interpolation
#' is performed.
#'
#' @param P for polynomial curve interpolation, a \eqn{(d, d, N)}-dimensional array corresponding to a length \eqn{N} sequence
#' of \eqn{(d, d)}-dimensional HPD matrices (control points) through which the interpolating polynomial
#' curve passes. For polynomial surface interpolation, a \eqn{(d, d, N_1, N_2)}-dimensional array corresponding
#' to a tensor product grid of \eqn{(d, d)}-dimensional HPD matrices (control points) through which the interpolating
#' polynomial surface passes.
#' @param X for polynomial curve interpolation, a numeric vector of length \eqn{N} specifying the time points at which
#' the interpolating polynomial passes through the control points \code{P}. For polynomial surface interpolation, a list
#' with as elements two numeric vectors \code{x} and \code{y} of length \eqn{N_1} and \eqn{N_2} respectively. The numeric
#' vectors specify the time points on the tensor product grid \code{expand.grid(X$x, X$y)} at which the interpolating
#' polynomial passes trough the control points \code{P}.
#' @param x for polynomial curve interpolation, a numeric vector specifying the time points (locations) at which the
#' interpolating polynomial is evaluated. For polynomial surface interpolation, a list with as elements two numeric vectors
#' \code{x} and \code{y} specifying the time points (locations) on the tensor product grid \code{expand.grid(x$x, x$y)} at which the
#' interpolating polynomial surface is evaluated.
#' @param metric the metric on the space of HPD matrices, by default \code{metric = "Riemannian"}, but instead of the Riemannian metric
#' this can also be set to \code{metric = "Euclidean"} to perform (standard) Euclidean polynomial interpolation of curves or
#' surfaces in the space of HPD matrices.
#'
#' @return For polynomial curve interpolation, a \code{(d, d, length(x))}-dimensional array corresponding to the interpolating polynomial
#' curve of \eqn{(d,d)}-dimensional matrices of degree \eqn{N-1} evaluated at times \code{x} and passing through the control points \code{P}
#' at times \code{X}. For polynomial surface interpolation, a \code{(d, d, length(x$x), length(x$y))}-dimensional array corresponding to the
#' interpolating polynomial surface of \eqn{(d,d)}-dimensional matrices of bi-degree \eqn{N_1 - 1, N_2 - 1} evaluated at times \code{expand.grid(x$x, x$y)}
#' and passing through the control points \code{P} at times \code{expand.grid(X$x, X$y)}.
#'
#' @note
#' If \code{metric = "Euclidean"}, the interpolating curve or surface may not be positive definite everywhere as the space of HPD
#' matrices equipped with the Euclidean metric has its boundary at a finite distance.
#'
#' @note
#' The function does not check for positive definiteness of the input matrices, and may fail if \code{metric = "Riemannian"} and
#' the input matrices are close to being singular.
#'
#' @examples
#' ### Polynomial curve interpolation
#' P <- rExamples1D(50, example = 'gaussian')$f[, , 10*(1:5)]
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
#' @references
#' \insertAllCited{}
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

  if(!is_2D) {
    ## 1D Neville's algorithm C++
    PP <- pdNeville_C(P, X, x, metric)
  } else if(is_2D) {
    ## 2D Neville's algorithm via geodesic surface interpolation
    if(n[1] < 1 & n[2] < 1){
      PP <- array(P[, , 1, 1], dim = c(d, d, length(x$x), length(x$y)))
    } else if(n[1] < 1 | n[2] < 1){
      if(length(X$y) > 1){
        PP_i <- pdNeville_C(P[, , 1, ], X$y, x$y, metric)
        PP <- aperm(replicate(length(x$x), PP_i), c(1, 2, 4, 3))
      } else if(length(X$x) > 1){
        PP_j <- pdNeville_C(P[, , , 1], X$x, x$x, metric)
        PP <- replicate(length(x$y), PP_j)
      }
    } else {
      PP_i <- sapply(1:(n[1] + 1), function(i) pdNeville_C(P[, , i, ], X$y, x$y, metric), simplify = "array")
      PP <- sapply(1:length(x$y), function(j) pdNeville_C(PP_i[, , j, ], X$x, x$x, metric), simplify = "array")
    }
  }
  return(PP)
}

#' Generate intrinsic HPD polynomial curves
#'
#' \code{pdPolynomial} generates intrinsic polynomial curves in the manifold of HPD matrices
#' equipped with the affine-invariant Riemannian metric (see \insertCite{B09}{pdSpecEst}[Chapter 6]
#' or \insertCite{PFA05}{pdSpecEst}) according to the numerical integration procedure in \insertCite{HFJ14}{pdSpecEst}.
#' Given an initial starting point \code{p0} (i.e., a HPD matrix) in the Riemannian manifold and covariant
#' derivatives up to order \eqn{k - 1} at \code{p0}, \code{pdPolynomial} approximates the uniquely existing
#' intrinsic polynomial curve of degree \eqn{k} passing through \code{p0} with the given covariant derivatives up
#' to order \eqn{k - 1} and vanishing higher order covariant derivatives.
#'
#' @param p0 a \eqn{(d, d)}-dimensional HPD matrix specifying the starting point of the polynomial curve.
#' @param v0 a \eqn{(d, d, k)}-dimensional array corresponding to a sequence of \eqn{(d,d)}-dimensional Hermitian matrix-valued
#' covariant derivatives from order zero up to order \eqn{k - 1} at the starting point \code{p0}.
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
#' P.geo <- sapply(seq(0, 1, length = 100), function(t) Expm(p0, t * Logm(p0, P.poly[, , 100])),
#'                simplify = "array")
#' all.equal(P.poly, P.geo)
#'
#' @return A \code{(d, d, length(steps))}-dimensional array corresponding to a generated (approximate)
#' intrinsic polynomial curve in the space of \eqn{(d,d)}-dimensional HPD matrices of degree \eqn{k}
#' passing through \code{p0} with the given covariant derivatives \code{v0} up to order \eqn{k - 1}
#' and vanishing higher order covariant derivatives.
#'
#' @references
#' \insertAllCited{}
#'
#' @seealso \code{\link{pdNeville}}, \code{\link{pdParTrans}}
#'
#' @export
pdPolynomial <- function(p0, v0, delta.t = 0.01, steps = 100) {

  if(!isTRUE(length(dim(v0)) == 3 & identical(dim(v0)[1:2], dim(p0)))) {
    stop("Incorrect input formats for arguments: 'p0' or 'v0', consult
         the function documentation for the correct inputs.")
  }
  if(!isTRUE(all.equal(apply(v0, 3, function(x) t(Conj(x))), apply(v0, 3, c)))) {
    stop("The matrix-valued covariant derivatives in 'v0' are not Hermitian, consult
          the function documentation for the correct inputs.")
  }
  if(!isTRUE(all(eigen(p0, symmetric = TRUE, only.values = TRUE)$values > 0))){
    stop("'p0' is not a positive definite matrix, consult the function documentation
         for the correct input.")
  }
  d <- dim(p0)[1]
  k <- max(dim(v0)[3], 2)
  p <- array(dim = c(d, d, steps))
  p[, , 1] <- p0
  vi <- array(0, dim = c(d, d, k))
  vi[, , 1:dim(v0)[3]] <- v0

  for(ti in 1:(steps - 1)){
    w <- vi[, , 1]
    vi[, , 1:(k - 1)] <- sapply(1:(k - 1), function(i) pdParTrans(p[, , ti], delta.t * w, vi[, , i] +
                                                                  delta.t * vi[, , i + 1]), simplify = "array")
    vi[, , k] <- pdParTrans(p[, , ti], delta.t * w, vi[, , k])
    p[, , ti + 1] <- Expm(p[, , ti], delta.t * w)
  }

  return(p)
}

