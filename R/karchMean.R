#' Weighted geometric mean of HPD matrices
#'
#' \code{pdMean} calculates an (approximate) weighted geometric mean of \eqn{S} different
#' \eqn{(d \times d)}-dimensional Hermitian PD matrices based on the Riemannian metric by
#' the fast recursive algorithm in (Chau and von Sachs, 2017) or the slower but more accurate
#' gradient descent algorithm in (Pennec, 2006). By default, the unweighted geometric mean is computed.
#'
#' @param M a \eqn{(d,d,S)}-dimensional array of Hermitian PD matrices.
#' @param w an \eqn{S}-dimensional nonnegative weight vector, such that \code{sum(w) = 1}.
#' @param grad_desc a logical value deciding if the gradient descent algorithm be used, defaults to
#' \code{FALSE}.
#' @param max_iter maximum number of iterations in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc)}.
#' @param tol optional tolerance parameter in gradient descent algorithm, only used if
#' \code{isTRUE(grad_desc)}, defaults to \code{.Machine$double.eps}.
#' @param ... additional arguments for internal usage.
#'
#' @examples
#' m <- function(){
#'  X <- matrix(complex(real=rnorm(9), imaginary=rnorm(9)), nrow=3)
#'  t(Conj(X)) %*% X
#' }
#' M <- replicate(100, m())
#' z <- rnorm(100)
#' w <- abs(z)/sum(abs(z))
#' Ave <- pdMean(M, w)
#'
#' @references Chau, J. and von Sachs, R. (2017). \emph{Positive definite multivariate spectral
#' estimation: a geometric wavelet approach}. Available at \url{http://arxiv.org/abs/1701.03314}.
#' @references Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
#' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
#'
#' @seealso \code{\link{Mid}}
#'
#' @export
pdMean <- function(M, w, grad_desc = F, max_iter = 1000, tol, ...) {

  ## Set parameters
  dots = list(...)
  metric = (if(is.null(dots$metric)) "Riemannian" else dots$metric)
  tol = (if(missing(tol)) NA else tol)
  w = (if(missing(w)) rep(1/dim(M)[3], dim(M)[3]) else w)
  d = dim(M)[1]

    if(dim(M)[3] == 1){
      Mean <- M[, , 1]
    } else{

      if(metric == "Riemannian"){

      ## Recursive algorithm
      Mean <- kMean(do.call(rbind, lapply(1:dim(M)[3], function(s) M[, , s])), w)

      ## Gradient descent algorithm
      if(grad_desc){
        tol <- (if(is.na(tol)) .Machine$double.eps else tol)
        if(!isTRUE(is.numeric(tol))){
          stop("'tol' should be NA (default) or a numeric value.")
        }
        Mean_new <- Mean
        i <- 0
        while((pdDist(Mean_new, Mean) > tol) & (i < max_iter)){
          Mean <- Mean_new
          Mean_new <- Expm(Mean, apply(sapply(1:dim(M)[3], function(i) w[i] *
                                                Logm(Mean, M[, , i]), simplify = "array"), c(1, 2), sum))
          i <- i+1
        }
        Mean <- Mean_new
      }
    } else if(metric == "Euclidean"){

      ## Euclidean weighted mean
      Mean <- apply(array(rep(w, each = d^2), dim = c(d, d, dim(M)[3])) * M, c(1, 2), sum)
    }
  }
  return(Mean)
}







