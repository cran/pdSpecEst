# include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]

//' Logarithmic map
//'
//' \code{Logm(P, Q)} computes the projection of a Hermitian PD matrix \code{Q} on the Riemannian manifold
//' to the tangent space attached at the Hermitian PD matrix \code{P} via the logarithmic map as in (Pennec, 2006).
//' This is the unique inverse of the exponential map \code{\link{Expm}}.
//'
//' @param P a Hermitian positive definite matrix.
//' @param Q a Hermitian positive definite matrix (of equal dimension as \code{P}).
//'
//' @examples
//'  q <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  Q <- t(Conj(q)) %*% q
//'  p <- matrix(complex(real = rnorm(9), imaginary = rnorm(9)), nrow = 3)
//'  P <- t(Conj(p)) %*% p
//'  Logm(P, Q)
//'
//' @references
//' Pennec, X. (2006). Intrinsic statistics on Riemannian manifolds: Basic tools for geometric
//' measurements. \emph{Journal of Mathematical Imaging and Vision} 25(1), 127-154.
//'
//' @seealso \code{\link{Expm}}
//'
//' @export
// [[Rcpp::export()]]

arma::cx_mat Logm(arma::cx_mat P, arma::cx_mat Q) {

  arma::cx_mat P1 = arma::sqrtmat_sympd(P);

  arma::cx_mat P2 = arma::inv_sympd(P1);

  arma::cx_mat P3 = arma::logmat_sympd(P2 * Q * P2);

  return P1 * P3 * P1;

}





