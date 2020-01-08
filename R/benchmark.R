#' Cubic smoothing spline regression for HPD matrices
#'
#' \code{pdSplineReg()} performs cubic smoothing spline regression in the space of HPD matrices equipped with the
#' affine-invariant Riemannian metric through minimization of a penalized regression objective function using a
#' geometric conjugate gradient descent method as outlined in \insertCite{BA11}{pdSpecEst} and \insertCite{BA11b}{pdSpecEst}.
#' This is a specific implementation of the more general algorithm in \insertCite{BA11}{pdSpecEst} and \insertCite{BA11b}{pdSpecEst},
#' setting the part in the objective function based on the first-order finite geometric differences to zero, such that the solutions
#' of the regression problem are approximating cubic splines.
#'
#' @param P a \eqn{(d,d,n)}-dimensional array corresponding to a length \eqn{n} sequence of (\eqn{d, d})-dimensional
#' noisy HPD matrices.
#' @param f0 a \eqn{(d,d,n)}-dimensional array corresponding to an initial estimate of the smooth
#' target curve of (\eqn{d, d})-dimensional HPD matrices.
#' @param lam a smoothness penalty, defaults to \code{lam = 1}. If \code{lam = 0}, the penalized curve estimate
#' coincides with geodesic interpolation of the data points with respect to the Riemannian metric.
#' If \code{lam} increases to \eqn{\infty}, the penalized regression estimator is approximately a fitted geodesic curve.
#' @param Nd a numeric value (\code{Nd <= n}) determining a lower resolution of the cubic spline regression estimator to speed up
#' computation time, defaults to \eqn{n}.
#' @param ini_step initial candidate step size in a backtracking line search based on the Armijo-Goldstein
#' condition, defaults to \code{ini_step = 1}.
#' @param max_iter maximum number of gradient descent iterations, defaults to \code{max_iter = 100}.
#' @param eps optional tolerance parameter in gradient descent algorithm. The gradient descent procedure exits if the
#' absolute difference between consecutive evaluations of the objective function is smaller than \code{eps},
#' defaults to \code{eps = 1E-3}.
#' @param ... additional arguments for internal use.
#'
#' @note This function does not check for positive definiteness of the matrices given as input, and may fail
#' if matrices are close to being singular.
#'
#' @return A list with three components:
#' \item{f }{a \eqn{(d, d, N_d)}-dimensional array corresponding to a length \code{Nd} estimated cubic smoothing spline
#' curve of (\eqn{d, d})-dimensional HPD matrices.}
#' \item{cost }{a numeric vector containing the costs of the objective function at each gradient descent iteration.}
#' \item{total_iter }{total number of gradient descent iterations.}
#'
#' @examples
#' \dontrun{
#' set.seed(2)
#' P <- rExamples1D(50, example = 'gaussian', noise.level = 0.1)
#' P.spline <- pdSplineReg(P$P, P$P, lam = 0.5, Nd = 25)
#'
#' ## Examine matrix-component (1,1)
#' plot((1:50)/50, Re(P$P[1, 1, ]), type = "l", lty = 2) ## noisy observations
#' lines((1:25)/25, Re(P.spline$f[1, 1, ])) ## estimate
#' lines((1:50)/50, Re(P$f[1, 1, ]), col = 2, lty = 2) ## smooth target
#' }
#' @references
#' \insertAllCited{}
#'
#' @export
pdSplineReg <- function(P, f0, lam = 1, Nd, ini_step = 1, max_iter = 100, eps = 1E-3, ...) {

  ## Set variables
  dots = list(...)
  tau = (if(is.null(dots$tau)) 0.5 else dots$tau)
  sigma = (if(is.null(dots$sigma)) 0.5 else dots$sigma)
  max_iter_a = (if(is.null(dots$max_iter_a)) 100 else dots$max_iter_a)

  d <- dim(P)[1]
  n <- dim(P)[3]
  if(missing(Nd)){
    Nd <- n
  }

  Nd.seq <- round(seq(from = 1, to = n, length = Nd))
  Nd.ind <- sapply(1:n, function(j) which.min(abs(Nd.seq - j)))
  delta.t <- mean(diff(Nd.seq))

  ast <- function(A, B) (A %*% B) %*% t(Conj(A))

  ## Compute directional derivatives and gradients as in (Boumal and Absil, 2011b)
  DLog <- function(A, H){
    e <- eigen(A, symmetric = TRUE)
    e_vec <- e$vectors
    e_val <- e$values
    H1 <- (t(Conj(e_vec)) %*% H) %*% e_vec
    grid <- expand.grid(1:d, 1:d)
    Z1 <- array(mapply(function(i1, i2) if(i1 == i2) 1/e_val[i1] else (log(e_val[i1]) - log(e_val[i2])) /
                         (e_val[i1] - e_val[i2]), grid$Var1, grid$Var2), dim = c(d, d))
    e_vec %*% (H1 * Z1) %*% t(Conj(e_vec))
  }

  grad_fA <- function(A, B) -2 * Logm(A, B)

  grad_gA <- function(A, B, C){
    B1 <- iSqrt(B)
    C1 <- iSqrt(C)
    ast(A, ast(B1, DLog(ast(B1, A), Logm(diag(d), ast(Sqrt(solve(ast(B1, C))), ast(B1, A)))))) +
      ast(A, ast(C1, DLog(ast(C1, A), Logm(diag(d), ast(Sqrt(ast(C1, A)), solve(ast(C1, B)))))))
  }

  grad_gB <- function(A, B, C){
    A1 <- iSqrt(A)
    ast(B, ast(A1, DLog(ast(A1, B), Logm(diag(d), ast(A1, C)))))
  }

  grad_E <- function(gamma){
    sapply(1:Nd, function(k) 0.5 * apply(sapply(which(Nd.ind == k),
                                                function(ki) grad_fA(gamma[, , k], P[, , ki]), simplify = "array"),
                                         c(1, 2), sum) + delta.t^3  * lam / 2 * (if(k %in% 2:(Nd-1)){
                                           2 * grad_fA(gamma[, , k], gamma[, , k + 1]) +
                                             2 * grad_fA(gamma[, , k], gamma[, , k - 1])  +
                                             2 * grad_gA(gamma[, , k], gamma[, , k + 1], gamma[, , k - 1]) +
                                             (if(k > 2) 2 * grad_gB(gamma[, , k - 1], gamma[, , k], gamma[, , k - 2]) else 0) +
                                             (if(k < Nd-1) 2 * grad_gB(gamma[, , k + 1], gamma[, , k], gamma[, , k + 2]) else 0)
                                         } else 0), simplify = "array")
  }

  ## Cost objective function
  E <- function(gamma){
    gamma.isqrt <- sapply(1:Nd, function(k) iSqrt(gamma[, , k]), simplify = "array")
    return(0.5 * sum(sapply(1:n, function(i) pdDist(P[, , i], gamma[, , which.min(abs(Nd.seq - i))])^2)) +
             delta.t^3 * lam / 2 * sum(sapply(2:(Nd - 1), function(k) NormF(ast(gamma.isqrt[, , k],
                                                                                Logm(gamma[, , k], gamma[, , k + 1]) + Logm(gamma[, , k], gamma[, , k - 1])))^2)))
  }

  ## Backtrack line search to determine optimal step size
  backtrack <- function(p, gamma, E0){
    alpha <- ini_step
    t <- mean(sapply(1:Nd, function(k) sigma * NormF(ast(iSqrt(gamma[, , k]), p[, , k]))^2))
    E1 <- NULL
    iter_a <- 0
    while(isTRUE(iter_a < max_iter_a) & isTRUE(is.null(E1) | isTRUE((E0 - E1) < (alpha * t)))){
      E1 <- tryCatch({ E(sapply(1:dim(p)[3], function(i) Expm(gamma[, , i], alpha * p[, , i]),
                                simplify = "array")) }, error = function(e) return(NULL))
      alpha <- tau * alpha
      iter_a <- iter_a + 1
    }
    return(c(alpha, iter_a))
  }

  ## Geometric conjugate gradient descent
  p <- -grad_E(f0[, , Nd.seq])
  cost <- E(f0[, , Nd.seq])
  gamma_0 <- f0[, , Nd.seq]
  iter <- 0
  cost_diff <- -1

  while(isTRUE(abs(cost_diff) > eps) & isTRUE(cost_diff < 0) & isTRUE(iter < max_iter)){
    alpha <- backtrack(p, gamma_0, tail(cost, 1))
    if(alpha[2] == max_iter_a){
      message("Backtracking line search not converging, increase 'max_iter_a' or choose smaller 'ini_step'")
      break
    }
    gamma_1 <- sapply(1:Nd, function(k) Expm(gamma_0[, , k], alpha[1] * p[, , k]), simplify = "array")
    grad_1 <- grad_E(gamma_1)
    beta_1 <- sapply(1:Nd, function(k) NormF(grad_1[, , k])^2/NormF(p[, , k])^2)
    p <- -grad_1 + array(rep(beta_1, each = d^2), dim = c(d, d, Nd)) * p
    cost <- c(cost, E(gamma_1))
    cost_diff <- tail(cost, 1) - tail(cost, 2)[1]
    gamma_0 <- gamma_1
    iter <- iter + 1
    if(iter == max_iter){
      message("Reached maximum number of iterations in gradient descent")
    }
  }

  return(list(f = gamma_0, cost = cost, total_iter = iter))

}

#' Orthonormal basis expansion of a Hermitian matrix
#'
#' \code{H.coeff} expands a \eqn{(d,d)}-dimensional Hermitian matrix \code{H}  with respect to
#' an orthonormal (in terms of the Frobenius inner product) basis of the space of Hermitian matrices.
#' That is, \code{H.coeff} transforms \code{H} into a numeric vector of \eqn{d^2} real-valued basis coefficients,
#' which is possible as the space of Hermitian matrices is a real vector space. Let \eqn{E_{nm}} be a
#' \eqn{(d,d)}-dimensional zero matrix with a 1 at location \eqn{(1, 1) \leq (n,m) \leq (d,d)}.
#' The orthonormal basis contains the following matrix elements; let  \eqn{1 \le n \le d} and
#' \eqn{1 \le m \le d},
#' \describe{
#'   \item{If \code{n == m}}{ the real matrix element \eqn{E_{nn}}}
#'   \item{If \code{n < m}}{ the complex matrix element \eqn{2i/\sqrt 2 E_{nm}}}
#'   \item{If \code{n > m}}{ the real matrix element \eqn{2/\sqrt 2 E_{nm}}}
#' }
#' The orthonormal basis coefficients are ordered by scanning through the matrix \code{H} in a row-by-row
#' fashion.
#'
#' @param H if \code{inverse = FALSE}, a \eqn{(d,d)}-dimensional Hermitian matrix; if \code{inverse = TRUE}, a numeric
#' vector of length \eqn{d^2} with \eqn{d} an integer.
#' @param inverse a logical value that determines whether the forward basis transform (\code{inverse = FALSE}) or the inverse
#' basis transform (\code{inverse = TRUE}) should be applied.
#'
#' @return If \code{inverse = FALSE} takes as input a \eqn{(d,d)}-dimensional Hermitian matrix and outputs a numeric
#' vector of length \eqn{d^2} containing the real-valued basis coefficients. If \code{inverse = TRUE} takes as input a
#' \eqn{d^2}-dimensional numeric vector of basis coefficients and outputs the corresponding \eqn{(d,d)}-dimensional
#' Hermitian matrix.
#'
#' @examples
#' ## random Hermitian matrix
#' H <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
#' diag(H) <- rnorm(3)
#' H[lower.tri(H)] <- t(Conj(H))[lower.tri(H)]
#'
#' ## orthonormal basis expansion
#' h <- H.coeff(H)
#' H1 <- H.coeff(h, inverse = TRUE) ## reconstructed Hermitian matrix
#' all.equal(H, H1)
#'
#' @export
H.coeff <- function(H, inverse = FALSE){
  if(!isTRUE(inverse)){
    if(!isTRUE(all.equal(t(Conj(H)), H))){
      stop("'H' should be an Hermitian matrix.")
    }
    HH <- E_coeff(H)
  } else{
    if(!isTRUE(is.numeric(H) & all.equal(sqrt(length(H)), round(sqrt(length(H)))))){
      stop("'H' should be a numeric vector of length 'd^2' with 'd' integer.")
    }
    HH <- E_coeff_inv(H)
  }
  return(HH)
}

## Orthonormal basis expansion Cholesky matrix
E_chol <- function(R, inverse = FALSE){
  d <- (if(!inverse) dim(R)[1] else round(sqrt(length(R))))

  if(!inverse){
    P <- c(Re(R)[lower.tri(R, diag = TRUE)], Im(R)[lower.tri(R)])
  } else {
    P1 <- array(0, dim = c(d, d))
    P2 <- array(0i, dim = c(d, d))
    P1[lower.tri(P1, diag = TRUE)] <- R[1:(d * (d + 1)/2)]
    P2[lower.tri(P2)] <- complex(imaginary = tail(R, -d * (d + 1)/2))
    P <- P1 + P2
  }
  return(P)
}

## Basis Hermitian matrices
E_basis <- function(d) {
  E <- function(i, j) {
    E_ij <- matrix(0, nrow = d, ncol = d)
    if (i == j) {
      E_ij[i, j] <- 1
    } else if (i < j) {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- 1/sqrt(2)
    } else {
      E_ij[c((j - 1) * d + i, (i - 1) * d + j)] <- complex(imaginary = c(-1, 1))/sqrt(2)
    }
    E_ij
  }

  indices <- expand.grid(1:d, 1:d)
  return(mapply(E, indices$Var1, indices$Var2, SIMPLIFY = "array"))
}

## Basis Tangent space
T_basis <- function(E, y) {
  d <- nrow(E)
  y.sqrt <- Sqrt(y)
  return(array(c(apply(E, 3, function(E) (y.sqrt %*% E) %*% y.sqrt)), dim = c(d, d, d^2)))
}

## Golden section search
gss <- function(range, f, tol = 0.01) {
  x1 <- range[1]
  x4 <- range[2]
  x2 <- x4 + (x1 - x4)/((1 + sqrt(5))/2)
  x3 <- x1 + (x4 - x1)/((1 + sqrt(5))/2)
  fx2 <- f(x2)
  fx3 <- f(x3)
  i <- 0
  while (!isTRUE(all.equal(fx2, fx3)) & (abs(x2 - x3) > tol) & (i <= 20)) {
    i <- i + 1
    if (fx2 < fx3) {
      x4 <- x3
      x3 <- x2
      fx3 <- fx2
      x2 <- x4 + (x1 - x4)/((1 + sqrt(5))/2)
      fx2 <- f(x2)
    } else {
      x1 <- x2
      x2 <- x3
      fx2 <- fx3
      x3 <- x1 + (x4 - x1)/((1 + sqrt(5))/2)
      fx3 <- f(x3)
    }
  }
  return(mean(c(x2, x3)))
}
